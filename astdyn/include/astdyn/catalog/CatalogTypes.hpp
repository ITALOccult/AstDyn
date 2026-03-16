/**
 * @file CatalogTypes.hpp
 * @brief Core types for the astdyn::catalog stellar catalog module
 *
 * Defines the public data types used for stellar catalog queries.
 * The underlying implementation uses the Gaia DR3 catalog via IOC_GaiaLib.
 */

#ifndef ASTDYN_CATALOG_TYPES_HPP
#define ASTDYN_CATALOG_TYPES_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <optional>
#include <functional>
#include <tuple>

#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "src/core/frame_tags.hpp"

namespace astdyn::catalog {

using namespace astdyn::astrometry;

// ============================================================================
// Module-local constants
// ============================================================================
namespace constants {
    inline constexpr double AU_PER_PARSEC = 206264.806; ///< 1 pc = 206264.806 AU
}

// ============================================================================
// Star — Gaia DR3 stellar record
// ============================================================================

/**
 * @brief A single star record from the Gaia DR3 catalog.
 *
 * Positional data is in the ICRS frame (== GCRF for stellar objects).
 * Proper motions are given at the Gaia reference epoch J2016.0.
 */
struct Star {
    // --- Identification ---
    int64_t     source_id = 0;          ///< Gaia DR3 source identifier
    std::string common_name;            ///< IAU common name (if assigned)
    std::string hd_designation;         ///< Henry Draper catalog number
    std::string hip_designation;        ///< Hipparcos catalog number
    std::string tycho2_designation;     ///< Tycho-2 designation
    std::string sao_designation;        ///< SAO catalog number

    // --- Astrometry (ICRS, epoch J2016.0) ---
    RightAscension ra;                  ///< Right Ascension
    Declination    dec;                 ///< Declination
    Parallax       parallax;            ///< Parallax (mas)
    double         parallax_error_mas = 0.0;
    ProperMotion   pm_ra_cosdec;        ///< Proper motion in RA * cos(dec)
    ProperMotion   pm_dec;              ///< Proper motion in Dec
    double         pmra_error_mas_yr  = 0.0;
    double         pmdec_error_mas_yr = 0.0;

    // --- Photometry ---
    double g_mag  = 99.0;               ///< Gaia G-band magnitude
    double bp_mag = 99.0;               ///< Gaia BP-band magnitude
    double rp_mag = 99.0;               ///< Gaia RP-band magnitude
    double bp_rp  = 0.0;                ///< BP-RP color index

    // --- Quality ---
    double ruwe = 0.0;                  ///< Renormalised unit weight error
    double astrometric_excess_noise = 0.0;

    /// Predict star position at a specific TDB epoch (proper motion + annual parallax)
    [[nodiscard]] SkyCoord<core::GCRF> predict_at(
        time::EpochTDB target_time,
        const std::optional<Eigen::Vector3d>& observer_pos_ssb_m = std::nullopt) const;

    [[nodiscard]] bool has_parallax() const noexcept { return parallax.to_mas() > 0.0; }
    [[nodiscard]] bool has_proper_motion() const noexcept {
        return pm_ra_cosdec.to_mas_yr() != 0.0 || pm_dec.to_mas_yr() != 0.0;
    }
};

// ============================================================================
// Query parameter structs
// ============================================================================

/**
 * @brief Parameters for a circular cone search.
 */
struct ConeQuery {
    RightAscension ra;           ///< Center RA
    Declination    dec;          ///< Center Dec
    Angle          radius;       ///< Search radius
    double         max_magnitude = 18.0; ///< Faint limit (Gaia G)
    Parallax       min_parallax; ///< Optional parallax cut
};

/**
 * @brief A point on the sky, used to define a corridor path.
 */
struct SkyPoint {
    double ra_deg  = 0.0; ///< Right Ascension [deg]
    double dec_deg = 0.0; ///< Declination [deg]
};

/**
 * @brief Parameters for a corridor search along a polyline path.
 *
 * Useful for searching stars along the projected path of a solar system body.
 */
struct CorridorQuery {
    std::vector<SkyCoord<core::GCRF>> path; ///< Ordered sky positions
    Angle          width;               ///< Half-width of corridor
    double         max_magnitude = 18.0; ///< Faint magnitude limit (Gaia G)
    Parallax       min_parallax;        ///< Optional parallax cut
    size_t         max_results = 10000;  ///< Maximum number of returned stars
};

/**
 * @brief Chebyshev polynomial representing RA or Dec as a function of time.
 *
 * Covers the interval [t_start, t_end] (Julian Date, TDB).
 * The polynomial is evaluated as:
 *   f(t) = sum_k coeffs[k] * T_k(tau),  tau = 2*(t - t_start)/(t_end - t_start) - 1
 */
struct ChebyshevSegment {
    double t_start = 0.0;               ///< Segment start [JD TDB]
    double t_end   = 0.0;               ///< Segment end   [JD TDB]
    std::vector<double> ra_coeffs;      ///< Chebyshev coefficients for RA [deg]
    std::vector<double> dec_coeffs;     ///< Chebyshev coefficients for Dec [deg]
    std::vector<double> dist_coeffs;    ///< Chebyshev coefficients for Distance [AU]

    /**
     * @brief Evaluate polynomial at specific JD.
     * @param jd Julian Date [TDB]
     * @return tuple of RA (deg), Dec (deg), and Distance (AU).
     */
    [[nodiscard]] std::tuple<double, double, double> evaluate_all(double jd) const {
        auto [pos, vel] = evaluate_full(jd);
        return pos;
    }

    /**
     * @brief Evaluate polynomial and its derivative (velocity) at specific JD.
     * 
     * @param jd Julian Date [TDB]
     * @return pair of tuples: (RA deg, Dec deg, Dist AU), (vRA deg/day, vDec deg/day, vDist AU/day)
     */
    [[nodiscard]] std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>> 
    evaluate_full(double jd) const {
        if (ra_coeffs.empty() || dec_coeffs.empty()) {
            return {{0,0,0}, {0,0,0}};
        }

        double dt = t_end - t_start;
        double tau = (2.0 * jd - (t_end + t_start)) / dt;
        double dtau_dt = 2.0 / dt;

        auto eval = [&](const std::vector<double>& cs) -> std::pair<double, double> {
            if (cs.empty()) return {0.0, 0.0};
            if (cs.size() == 1) return {cs[0], 0.0};

            // Use Clenshaw-like recurrence for position and derivative
            // T_0 = 1, T_1 = x, T_n = 2xT_{n-1} - T_{n-2}
            // T'_n = n U_{n-1}
            // Faster: T_buf approach
            size_t n = cs.size();
            std::vector<double> T(n), U(n);
            T[0] = 1.0;
            U[0] = 0.0;
            if (n > 1) {
                T[1] = tau;
                U[1] = 1.0;
            }
            for (size_t k = 2; k < n; ++k) {
                T[k] = 2.0 * tau * T[k-1] - T[k-2];
                U[k] = 2.0 * tau * U[k-1] - U[k-2] + 2.0 * T[k-1];
            }

            double p = 0.0, v = 0.0;
            for (size_t k = 0; k < n; ++k) {
                p += cs[k] * T[k];
                v += cs[k] * U[k];
            }
            return {p, v * dtau_dt};
        };

        auto [ra, vra] = eval(ra_coeffs);
        auto [dec, vdec] = eval(dec_coeffs);
        auto [dist, vdist] = dist_coeffs.empty() ? std::make_pair(2.5, 0.0) : eval(dist_coeffs);

        return {{ra, dec, dist}, {vra, vdec, vdist}};
    }

    /// Backwards compatibility for RA/Dec only
    [[nodiscard]] std::pair<double, double> evaluate(double jd) const {
        auto [ra, dec, r] = evaluate_all(jd);
        return {ra, dec};
    }
};

/**
 * @brief Parameters for searching stars along a propagated orbit.
 *
 * The orbit is represented as a set of Chebyshev segments covering the full
 * time interval. Build this with make_orbit_query() in CatalogIntegration.hpp.
 */
struct OrbitQuery {
    time::EpochTDB t_start;             ///< Search window start
    time::EpochTDB t_end;               ///< Search window end
    std::vector<ChebyshevSegment> segments; ///< Chebyshev approximation of the orbit
    Angle          width;               ///< Corridor half-width
    double         max_magnitude = 18.0; ///< Faint magnitude limit (Gaia G)
    double         step_days     = 1.0;  ///< Discretisation step for orbit sampling [days]
};

// ============================================================================
// Progress callback
// ============================================================================

using ProgressCallback = std::function<void(double fraction, const std::string& message)>;

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_TYPES_HPP
