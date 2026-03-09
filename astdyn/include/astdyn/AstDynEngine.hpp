/**
 * @file AstDynEngine.hpp
 * @brief Main OrbFit engine for orbit determination
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 * 
 * This class replicates the main workflow of the original Fortran orbit fitting software program:
 * - Load observations
 * - Initial orbit determination (IOD)
 * - Differential correction (least squares fit)
 * - Orbit propagation
 * - Ephemeris generation
 * - Close approach analysis
 */

#ifndef ASTDYN_ORBFITENGINE_HPP
#define ASTDYN_ORBFITENGINE_HPP

#include "astdyn/core/physics_state.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/observations/ObservationManager.hpp"
#include "astdyn/orbit_determination/OrbitFitter.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/core/Enums.hpp"
#include "astdyn/close_approach/CloseApproach.hpp"
#include <memory>
#include <vector>
#include <string>
#include <iostream>

namespace astdyn {

/**
 * @brief Configuration for OrbFit engine
 */
struct AstDynConfig {
    // Propagation settings
    propagation::PropagatorSettings propagator_settings;
    
    // Integrator settings
    IntegratorType integrator_type = IntegratorType::RK4; ///< RK4, RKF78, AAS, etc.
    double initial_step_size = 0.1;                      ///< Initial step size [days]
    double tolerance = 1e-12;                            ///< Integration tolerance
    double aas_precision = 1e-4;                         ///< Precision metric for AAS integrator (step size control)

    // Ephemeris Configuration
    EphemerisType ephemeris_type = EphemerisType::Analytical; ///< Analytical, DE441
    std::string ephemeris_file;                               ///< Path to .bsp file (if DE441)
    std::string asteroid_ephemeris_file;                    ///< Path to asteroid .bsp file (optional)
    
    // Differential correction settings
    int max_iterations = 10;                 ///< Maximum DC iterations
    double convergence_threshold = 1e-6;     ///< Convergence threshold
    double outlier_sigma = 3.0;              ///< Outlier rejection threshold (legacy/default)
    double outlier_max_sigma = 10.0;         ///< Carpentry start sigma
    double outlier_min_sigma = 3.0;          ///< Carpentry end sigma
    
    // Close approach settings
    close_approach::CloseApproachSettings ca_settings;
    
    // Residual calculation settings
    bool aberration_correction = true;       ///< Apply annual aberration
    bool light_time_correction = true;       ///< Apply light-time correction
    
    // Output settings
    bool verbose = true;                     ///< Verbose output
    bool save_residuals = true;              ///< Save residual plots
    bool compute_ephemeris = true;           ///< Compute ephemeris
    
    // Time settings
    std::string eop_file = "";               ///< Path to IERS EOP file (finals.all)
    
    // Catalog Biases
    std::string catalog_bias_file = "";      ///< Path to catalog bias CSV file
};

/**
 * @brief Results from orbit determination
 */
struct OrbitDeterminationResult {
    physics::KeplerianStateTyped<core::ECLIPJ2000> orbit;    ///< Fitted orbital elements
    Eigen::MatrixXd covariance;              ///< Covariance matrix (6×6)
    std::vector<double> residuals_ra;        ///< RA residuals [arcsec]
    std::vector<double> residuals_dec;       ///< Dec residuals [arcsec]
    double rms_ra;                           ///< RMS of RA residuals [arcsec]
    double rms_dec;                          ///< RMS of Dec residuals [arcsec]
    double rms_total_arcsec;                 ///< Total RMS [arcsec]
    double chi_squared;                      ///< Chi-squared statistic
    int num_observations;                    ///< Number of observations used
    int num_rejected;                        ///< Number of rejected outliers
    int num_iterations;                      ///< Number of DC iterations
    bool converged;                          ///< Convergence flag
};

/**
 * @brief Apparent position results for occultations
 */
struct ApparentPlace {
    time::EpochTDB epoch;
    double ra;      ///< Topocentric RA [rad]
    double dec;     ///< Topocentric Dec [rad]
    double dist;    ///< Distance to observer [AU]
};

/**
 * @brief Main engine class for high-level orbit operations
 */
class AstDynEngine {
public:
    /**
     * @brief Construct engine with default configuration
     */
    AstDynEngine();
    
    /**
     * @brief Construct engine with specific configuration
     */
    explicit AstDynEngine(const AstDynConfig& config);
    
    // ========================================================================
    // Observation Management
    // ========================================================================
    
    /**
     * @brief Load observations from file (MPC format)
     */
    int load_observations(const std::string& filename);
    
    /**
     * @brief Add single observation
     */
    void add_observation(const observations::OpticalObservation& obs);
    
    /**
     * @brief Get all loaded observations
     */
    const std::vector<observations::OpticalObservation>& observations() const {
        return obs_context_.observations();
    }
    
    /**
     * @brief Clear all observations
     */
    void clear_observations() {
        obs_context_.clear();
    }
    
    // ========================================================================
    // Orbit Determination
    // ========================================================================
    
    /**
     * @brief Set the initial orbit for propagation
     */
    template <typename Frame>
    void set_initial_orbit(const physics::KeplerianStateTyped<Frame>& elements) {
        if constexpr (std::is_same_v<Frame, core::ECLIPJ2000>) {
            set_initial_orbit_ecl(elements);
        } else {
            auto cart = propagation::keplerian_to_cartesian(elements);
            auto pos_ecl = coordinates::ReferenceFrame::transform_pos<Frame, core::ECLIPJ2000>(cart.position);
            auto vel_ecl = coordinates::ReferenceFrame::transform_vel<Frame, core::ECLIPJ2000>(cart.position, cart.velocity);
            auto cart_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>(cart.epoch, pos_ecl, vel_ecl, cart.gm);
            set_initial_orbit_ecl(propagation::cartesian_to_keplerian<core::ECLIPJ2000>(cart_ecl));
        }
    }
    
    /**
     * @brief Perform initial orbit determination (IOD) from observations
     */
    physics::KeplerianStateTyped<core::ECLIPJ2000> initial_orbit_determination();
    
    /**
     * @brief Fit orbit to observations using differential correction
     */
    OrbitDeterminationResult fit_orbit();
    
    /**
     * @brief Get current best-fit orbit
     */
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& orbit() const {
        return current_orbit_;
    }
    
    /**
     * @brief Check if orbit is available
     */
    bool has_orbit() const {
        return has_orbit_;
    }
    
    // ========================================================================
    // Ephemeris Generation
    // ========================================================================
    
    std::vector<physics::CartesianStateTyped<core::GCRF>> compute_ephemeris(
        time::EpochTDB start_time,
        time::EpochTDB end_time,
        double step_days);
    
    physics::KeplerianStateTyped<core::ECLIPJ2000> propagate_to(time::EpochTDB target_time);
    
    // ========================================================================
    // Close Approach Analysis
    // ========================================================================
    
    std::vector<close_approach::CloseApproach> find_close_approaches(
        time::EpochTDB start_time,
        time::EpochTDB end_time);
    
    ApparentPlace compute_asteroid_apparent_place(time::EpochTDB t_occult, const std::string& observatory_code);

    ApparentPlace compute_star_apparent_place(
        double ra_j2000, double dec_j2000,
        double pm_ra, double pm_dec,
        time::EpochTDB t_occult, const std::string& observatory_code);
    
    double compute_moid(ephemeris::CelestialBody planet);
    
    // ========================================================================
    // Configuration and Settings
    // ========================================================================
    
    void set_config(const AstDynConfig& config) {
        config_ = config;
        update_propagator();
    }
    
    void load_config(const std::string& oop_file);
    
    const AstDynConfig& config() const {
        return config_;
    }
    
    void set_verbose(bool verbose) {
        config_.verbose = verbose;
    }
    
    // ========================================================================
    // Statistics and Diagnostics
    // ========================================================================
    
    /**
     * @brief Get shadow Hamiltonian drift from last AAS integration
     */
    double shadow_hamiltonian_drift() const;

    /**
     * @brief Get total number of force evaluations (NFE) from last integration
     */
    long total_force_evaluations() const;

    const OrbitDeterminationResult& last_result() const {
        return last_result_;
    }
    
    void print_orbit_summary(std::ostream& os = std::cout) const;
    void print_residuals_summary(std::ostream& os = std::cout) const;
    void export_orbit(const std::string& filename, const std::string& format = "oef");

private:
    void update_propagator();
    
    // Unified fit dispatcher
    template <typename Frame>
    OrbitDeterminationResult run_fit_in_frame();

    AstDynConfig config_;
    observations::ObservationManager obs_context_;
    
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<propagation::Propagator> propagator_;
    std::unique_ptr<close_approach::CloseApproachDetector> ca_detector_;

    physics::KeplerianStateTyped<core::ECLIPJ2000> current_orbit_;
    bool has_orbit_ = false;
    OrbitDeterminationResult last_result_;

    std::string loaded_ephemeris_file_ = "";
    bool ephemeris_loaded_ = false;

    void set_initial_orbit_ecl(const physics::KeplerianStateTyped<core::ECLIPJ2000>& elements);
};

// ============================================================================
// Template Implementations
// ============================================================================

template <typename Frame>
OrbitDeterminationResult AstDynEngine::run_fit_in_frame() {
    using namespace orbit_determination;
    
    OrbitFitter<Frame> fitter(ephemeris_, propagator_);
    fitter.set_corrections(config_.aberration_correction, config_.light_time_correction);
    
    DifferentialCorrectorSettings dc_settings;
    dc_settings.max_iterations = config_.max_iterations;
    dc_settings.convergence_tolerance = config_.convergence_threshold;
    dc_settings.outlier_sigma = config_.outlier_sigma;
    dc_settings.outlier_max_sigma = std::max(100.0, config_.outlier_sigma * 10.0); 
    dc_settings.outlier_min_sigma = config_.outlier_sigma;
    dc_settings.reject_outliers = true;
    dc_settings.verbose = config_.verbose;

    // Build initial state in Target Frame
    auto cart_legacy = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(current_orbit_);
    physics::CartesianStateTyped<Frame> initial_state;
    
    if constexpr (std::is_same_v<Frame, core::ECLIPJ2000>) {
        initial_state = cart_legacy;
    } else {
        auto pos_f = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, Frame>(cart_legacy.position, cart_legacy.epoch);
        auto vel_f = coordinates::ReferenceFrame::transform_vel<core::ECLIPJ2000, Frame>(cart_legacy.position, cart_legacy.velocity, cart_legacy.epoch);
        initial_state = physics::CartesianStateTyped<Frame>(cart_legacy.epoch, pos_f, vel_f, cart_legacy.gm);
    }

    auto result_dc = fitter.fit(obs_context_.observations(), initial_state, dc_settings);
    
    // Convert back to Ecliptic for storage
    physics::CartesianStateTyped<core::ECLIPJ2000> final_ecl;
    if constexpr (std::is_same_v<Frame, core::ECLIPJ2000>) {
        final_ecl = result_dc.final_state;
    } else {
        auto pos_e = coordinates::ReferenceFrame::transform_pos<Frame, core::ECLIPJ2000>(result_dc.final_state.position, result_dc.final_state.epoch);
        auto vel_e = coordinates::ReferenceFrame::transform_vel<Frame, core::ECLIPJ2000>(result_dc.final_state.position, result_dc.final_state.velocity, result_dc.final_state.epoch);
        final_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>(result_dc.final_state.epoch, pos_e, vel_e, result_dc.final_state.gm);
    }

    current_orbit_ = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(final_ecl);
    
    OrbitDeterminationResult result;
    result.orbit = current_orbit_;
    result.converged = result_dc.converged;
    result.rms_ra = result_dc.statistics.rms_ra.to_arcsec();
    result.rms_dec = result_dc.statistics.rms_dec.to_arcsec();
    result.rms_total_arcsec = result_dc.statistics.rms_total.to_arcsec();
    result.num_observations = result_dc.statistics.num_observations;
    result.num_iterations = result_dc.iterations;
    result.covariance = result_dc.covariance;
    result.chi_squared = result_dc.statistics.chi_squared;
    result.num_rejected = result_dc.statistics.num_outliers;
    
    for (const auto& res : result_dc.residuals) {
        if (!res.outlier) {
            result.residuals_ra.push_back(res.residual_ra.to_arcsec());
            result.residuals_dec.push_back(res.residual_dec.to_arcsec());
        }
    }
    
    return result;
}

} // namespace astdyn

#endif // ASTDYN_ORBFITENGINE_HPP
