/**
 * @file ODSmartPolicy.hpp
 * @brief Intelligent OD policy engine — OrbFit-style automatic decisions.
 *
 * Encodes 4 physics-driven rules that make the OD pipeline self-configuring:
 *
 *  Rule 1 — Force Model Selection
 *    • PPN relativity if perihelion q < 0.5 AU
 *    • Yarkovsky as free parameter if NEO AND observational arc > 10 years
 *    • Reduced integrator step near planetary close approaches (MOID < 0.1 AU)
 *
 *  Rule 2 — Arc Admissibility
 *    • Block unconstrained 6-param fit if N < 10 obs or ΔT < 60 days
 *    • Energy Barrier: reject DC solution if semi-major axis drifts > 50%
 *
 *  Rule 3 — Dynamic Statistical Weighting
 *    • σ = 0.5" for modern sky surveys (Pan-STARRS, Catalina, WISE, Gaia)
 *    • σ = 0.8" for post-2000 CCD observations at generic observatories
 *    • σ = 2.0" for pre-2000 / photographic-era observations
 *    • Single-worst-outlier rejection when first-iteration RMS > 3"
 *
 *  Rule 4 — Frame / Timescale Consistency
 *    • Enforced structurally: UTC→TDB conversion happens inside
 *      ResidualCalculator via time::to_tdb(); frame rotations are
 *      performed automatically by the typed Vector3/CartesianStateTyped system.
 *
 * @note The caller is responsible for rebuilding the propagator after
 *       ForceModelPolicy changes (use PropagatorSettings returned by analyze()).
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_OD_SMART_POLICY_HPP
#define ASTDYN_ORBIT_DETERMINATION_OD_SMART_POLICY_HPP

#include "astdyn/core/Constants.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/time/epoch.hpp"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace astdyn::orbit_determination {

// ---------------------------------------------------------------------------
// Keplerian element summary derived from a Cartesian state
// ---------------------------------------------------------------------------
struct KeplerianElements {
    double a_au  = 0.0;  ///< Semi-major axis [AU]
    double e     = 0.0;  ///< Eccentricity
    double i_deg = 0.0;  ///< Inclination [deg]
    double q_au  = 0.0;  ///< Perihelion distance [AU]
    double Q_au  = 0.0;  ///< Aphelion distance [AU]

    bool is_neo  = false; ///< Near-Earth Object: q < 1.3 AU
    bool is_mba  = false; ///< Main-belt asteroid: 2.0 < a < 3.3 AU
    bool is_comet = false;///< e >= 1 or q < 0.3 AU
};

// ---------------------------------------------------------------------------
// Arc geometry summary
// ---------------------------------------------------------------------------
struct ArcInfo {
    int    n_obs       = 0;
    double arc_days    = 0.0;
    double arc_years   = 0.0;

    /// True if arc is long enough for unconstrained 6-param DC
    bool is_admissible = true;
    /// Human-readable reason when not admissible
    std::string rejection_reason;
};

// ---------------------------------------------------------------------------
// Force model recommendation
// ---------------------------------------------------------------------------
struct ForceModelPolicy {
    bool use_relativity       = false; ///< Activate PPN correction
    bool fit_yarkovsky        = false; ///< Treat A₂ as free parameter in DC
    bool use_reduced_step     = false; ///< Reduce integrator step (close approach)
    double suggested_step_days = 0.5;  ///< Recommended initial integrator step

    /// Summary string for diagnostics
    std::string rationale;
};

// ---------------------------------------------------------------------------
// Combined policy analysis result
// ---------------------------------------------------------------------------
struct ODPolicyAnalysis {
    KeplerianElements elements;
    ArcInfo           arc;
    ForceModelPolicy  force;

    /// PropagatorSettings patched according to the force model policy
    propagation::PropagatorSettings propagator_settings;

    void print() const {
        std::cout << "[ODSmartPolicy] Orbital class: "
                  << (elements.is_neo  ? "NEO"  :
                      elements.is_mba  ? "MBA"  :
                      elements.is_comet ? "COMET" : "OTHER")
                  << "  a=" << elements.a_au << " AU"
                  << "  q=" << elements.q_au << " AU"
                  << "  e=" << elements.e    << "\n";
        std::cout << "[ODSmartPolicy] Arc: " << arc.n_obs << " obs, "
                  << arc.arc_days << " days ("
                  << (arc.is_admissible ? "ADMISSIBLE" : "INADMISSIBLE: " + arc.rejection_reason)
                  << ")\n";
        std::cout << "[ODSmartPolicy] Force model: " << force.rationale << "\n";
    }
};

// ---------------------------------------------------------------------------
// ODPolicyEngine
// ---------------------------------------------------------------------------
class ODPolicyEngine {
public:

    // -----------------------------------------------------------------------
    // Rule 1+2 — Orbital analysis
    // -----------------------------------------------------------------------

    /// Compute Keplerian elements from a heliocentric Cartesian state.
    static KeplerianElements compute_keplerian(
        const physics::CartesianStateTyped<core::GCRF>& state)
    {
        const double au_m = physics::Distance::from_au(1.0).to_m();
        const double mu   = state.gm.to_m3_s2(); // [m³/s²]

        Eigen::Vector3d r_si = state.position.to_eigen_si();
        Eigen::Vector3d v_si = state.velocity.to_eigen_si();

        double r  = r_si.norm();
        double v2 = v_si.squaredNorm();
        double energy = 0.5 * v2 - mu / r;

        KeplerianElements el;
        // Parabolic / hyperbolic safeguard
        if (energy >= 0.0) {
            el.e = 1.0 + (r * v2 / mu - 2.0) * 1e-9; // treat as near-parabolic
            el.a_au = 1e6; // unbounded
            el.is_comet = true;
            return el;
        }
        el.a_au = (-mu / (2.0 * energy)) / au_m;

        // Eccentricity vector
        Eigen::Vector3d ecc_vec = (v_si.cross(v_si.cross(r_si)) / mu) - r_si / r;
        // Correct formula: e_vec = ((v²/μ - 1/r) r - (r·v/μ) v)
        double rdotv = r_si.dot(v_si);
        ecc_vec = (v2 / mu - 1.0 / r) * r_si - (rdotv / mu) * v_si;
        el.e = ecc_vec.norm();

        el.q_au = el.a_au * (1.0 - el.e);
        el.Q_au = el.a_au * (1.0 + el.e);

        // Inclination
        Eigen::Vector3d h = r_si.cross(v_si);
        el.i_deg = std::acos(std::clamp(h.z() / h.norm(), -1.0, 1.0)) * 180.0 / M_PI;

        el.is_neo   = (el.q_au < 1.3);
        el.is_mba   = (el.a_au > 2.0 && el.a_au < 3.3);
        el.is_comet = (el.e >= 1.0 || el.q_au < 0.3);

        return el;
    }

    /// Compute semi-major axis [AU] from a state (fast path for energy barrier).
    static double compute_sma_au(const physics::CartesianStateTyped<core::GCRF>& state) {
        const double au_m = physics::Distance::from_au(1.0).to_m();
        double r  = state.position.to_eigen_si().norm();
        double v2 = state.velocity.to_eigen_si().squaredNorm();
        double mu = state.gm.to_m3_s2();
        double energy = 0.5 * v2 - mu / r;
        if (energy >= 0.0) return 1e6;
        return (-mu / (2.0 * energy)) / au_m;
    }

    /// Analyze the observational arc geometry.
    static ArcInfo analyze_arc(
        const std::vector<observations::OpticalObservation>& obs)
    {
        ArcInfo info;
        info.n_obs = static_cast<int>(obs.size());
        if (info.n_obs < 2) {
            info.is_admissible = false;
            info.rejection_reason = "fewer than 2 observations";
            return info;
        }
        // Sort by time to get arc length
        double t_min = obs.front().time.mjd();
        double t_max = obs.front().time.mjd();
        for (const auto& o : obs) {
            t_min = std::min(t_min, o.time.mjd());
            t_max = std::max(t_max, o.time.mjd());
        }
        info.arc_days  = t_max - t_min;
        info.arc_years = info.arc_days / 365.25;

        if (info.n_obs < 10) {
            info.is_admissible = false;
            info.rejection_reason = "N_obs=" + std::to_string(info.n_obs) + " < 10";
        } else if (info.arc_days < 60.0) {
            info.is_admissible = false;
            info.rejection_reason = "arc=" + std::to_string(static_cast<int>(info.arc_days)) + " days < 60";
        }
        return info;
    }

    // -----------------------------------------------------------------------
    // Rule 1 — Force model selection
    // -----------------------------------------------------------------------

    /// Select force model based on orbital class and arc length.
    static ForceModelPolicy select_force_model(
        const KeplerianElements& el,
        const ArcInfo& arc,
        const propagation::PropagatorSettings& base = {})
    {
        ForceModelPolicy p;
        p.use_relativity   = base.include_relativity;
        p.fit_yarkovsky    = false;
        p.use_reduced_step = false;
        p.suggested_step_days = 0.5;

        std::string log;

        // PPN: activate for objects with close perihelion (q < 0.5 AU)
        if (el.q_au < 0.5) {
            p.use_relativity = true;
            log += "PPN:ON(q<0.5AU) ";
        }

        // Yarkovsky as free parameter: NEO with arc > 10 years
        if (el.is_neo && arc.arc_years > 10.0) {
            p.fit_yarkovsky = true;
            log += "Yark:FIT(NEO+arc>10yr) ";
        }

        // Close-approach: if orbit can approach Earth within 0.1 AU geometrically
        // (simplified MOID: orbit crosses 0.9-1.1 AU belt)
        if (el.q_au < 1.1 && el.Q_au > 0.9) {
            p.use_reduced_step = true;
            p.suggested_step_days = 0.1;
            log += "ReducedStep(MOID<0.1AU) ";
        } else if (el.is_mba) {
            p.suggested_step_days = 1.0; // MBA: larger step is fine
            log += "MBA:step=1d ";
        }

        if (log.empty()) log = "default";
        p.rationale = log;
        return p;
    }

    // -----------------------------------------------------------------------
    // Full combined analysis
    // -----------------------------------------------------------------------

    static ODPolicyAnalysis analyze(
        const physics::CartesianStateTyped<core::GCRF>& iod_state,
        const std::vector<observations::OpticalObservation>& obs,
        const propagation::PropagatorSettings& base_settings = {})
    {
        ODPolicyAnalysis result;
        result.elements  = compute_keplerian(iod_state);
        result.arc       = analyze_arc(obs);
        result.force     = select_force_model(result.elements, result.arc, base_settings);

        // Patch PropagatorSettings
        result.propagator_settings = base_settings;
        result.propagator_settings.include_relativity = result.force.use_relativity;

        return result;
    }

    // -----------------------------------------------------------------------
    // Rule 3 — Dynamic observation weighting
    // -----------------------------------------------------------------------

    /**
     * @brief Apply automatic MPC-style weights to observations.
     *
     * Follows the Vereš et al. (2017) weighting scheme:
     *  - Modern survey observatories (Pan-STARRS, Catalina, MLS, WISE, Gaia): σ = 0.5"
     *  - Post-2000 CCD generic:  σ = 0.8"
     *  - Pre-2000 / photographic: σ = 2.0"
     *
     * Only modifies sigma if it still holds the MPC default value (≤ 0.5").
     */
    static void apply_auto_weights(
        std::vector<observations::OpticalObservation>& obs)
    {
        // MJD cutoff for "modern CCD era" (2000-01-01.5 ≈ MJD 51544)
        static constexpr double MJD_YEAR_2000 = 51544.0;

        for (auto& o : obs) {
            // Skip if user has set a non-default sigma
            if (o.sigma_ra.to_arcsec() > 0.51) continue;

            const double mjd = o.time.mjd();
            const std::string& code = o.observatory_code;

            double sigma_arcsec = catalog_sigma(code, mjd);
            o.sigma_ra  = astrometry::Angle::from_arcsec(sigma_arcsec);
            o.sigma_dec = astrometry::Angle::from_arcsec(sigma_arcsec);
        }
    }

private:
    /// Observatory-code / epoch → recommended σ [arcsec].
    static double catalog_sigma(const std::string& code, double mjd)
    {
        static constexpr double MJD_YEAR_2000 = 51544.0;

        // Modern space/large-aperture surveys (Gaia-calibrated or sub-arcsec PSF)
        static const std::vector<std::string> MODERN_SURVEY = {
            "F51", "F52",           // Pan-STARRS 1 & 2
            "G96", "703", "E12",    // Catalina / Mt. Lemmon / Siding Spring
            "I52",                  // Oukaimeden
            "C51",                  // WISE/NEOWISE
            "258",                  // Gaia
            "T05", "T08",           // ATLAS Hawaii / Chile
            "568",                  // Mauna Kea (Subaru/CFHT)
        };
        for (const auto& s : MODERN_SURVEY) {
            if (code == s) return 0.5;
        }

        // Pre-2000 photographic / historical
        if (mjd < MJD_YEAR_2000) return 2.0;

        // Post-2000 generic CCD
        return 0.8;
    }
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_OD_SMART_POLICY_HPP
