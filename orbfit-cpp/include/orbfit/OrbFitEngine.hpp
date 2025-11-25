/**
 * @file OrbFitEngine.hpp
 * @brief Main OrbFit engine for orbit determination
 * @author OrbFit C++ Conversion Team
 * @date 2025-11-24
 * 
 * This class replicates the main workflow of the original Fortran OrbFit program:
 * - Load observations
 * - Initial orbit determination (IOD)
 * - Differential correction (least squares fit)
 * - Orbit propagation
 * - Ephemeris generation
 * - Close approach analysis
 */

#ifndef ORBFIT_ORBFITENGINE_HPP
#define ORBFIT_ORBFITENGINE_HPP

#include "orbfit/core/Types.hpp"
#include "orbfit/propagation/OrbitalElements.hpp"
#include "orbfit/propagation/Propagator.hpp"
#include "orbfit/observations/Observation.hpp"
#include "orbfit/observations/MPCReader.hpp"
#include "orbfit/orbit_determination/DifferentialCorrector.hpp"
#include "orbfit/orbit_determination/Residuals.hpp"
#include "orbfit/close_approach/CloseApproach.hpp"
#include "orbfit/ephemeris/PlanetaryEphemeris.hpp"
#include <memory>
#include <vector>
#include <string>
#include <iostream>

namespace orbfit {

/**
 * @brief Configuration for OrbFit engine
 */
struct OrbFitConfig {
    // Propagation settings
    propagation::PropagatorSettings propagator_settings;
    
    // Integrator settings
    double initial_step_size = 0.1;         ///< Initial step size [days]
    double tolerance = 1e-12;                ///< Integration tolerance
    
    // Differential correction settings
    int max_iterations = 10;                 ///< Maximum DC iterations
    double convergence_threshold = 1e-6;     ///< Convergence threshold
    double outlier_sigma = 3.0;              ///< Outlier rejection threshold (sigma)
    
    // Close approach settings
    close_approach::CloseApproachSettings ca_settings;
    
    // Output settings
    bool verbose = true;                     ///< Verbose output
    bool save_residuals = true;              ///< Save residual plots
    bool compute_ephemeris = true;           ///< Compute ephemeris
};

/**
 * @brief Results from orbit determination
 */
struct OrbitDeterminationResult {
    propagation::KeplerianElements orbit;    ///< Fitted orbital elements
    Eigen::MatrixXd covariance;              ///< Covariance matrix (6×6)
    std::vector<double> residuals_ra;        ///< RA residuals [arcsec]
    std::vector<double> residuals_dec;       ///< Dec residuals [arcsec]
    double rms_ra;                           ///< RMS of RA residuals [arcsec]
    double rms_dec;                          ///< RMS of Dec residuals [arcsec]
    double chi_squared;                      ///< Chi-squared statistic
    int num_observations;                    ///< Number of observations used
    int num_rejected;                        ///< Number of rejected outliers
    int num_iterations;                      ///< Number of DC iterations
    bool converged;                          ///< Convergence flag
};

/**
 * @brief Main OrbFit engine class
 * 
 * This class provides the complete OrbFit workflow:
 * 1. Load observations from MPC format
 * 2. Initial orbit determination (IOD)
 * 3. Differential correction (DC) to improve orbit
 * 4. Propagate orbit to generate ephemeris
 * 5. Analyze close approaches to planets
 * 
 * Example usage:
 * @code
 *   OrbFitEngine engine;
 *   
 *   // Load observations
 *   engine.load_observations("observations.txt");
 *   
 *   // Set initial orbit guess
 *   KeplerianElements initial_orbit = ...;
 *   engine.set_initial_orbit(initial_orbit);
 *   
 *   // Run orbit determination
 *   auto result = engine.fit_orbit();
 *   
 *   // Generate ephemeris
 *   auto ephemeris = engine.compute_ephemeris(start_mjd, end_mjd, step_days);
 *   
 *   // Check close approaches
 *   auto close_approaches = engine.find_close_approaches();
 * @endcode
 */
class OrbFitEngine {
public:
    /**
     * @brief Construct OrbFit engine with default configuration
     */
    OrbFitEngine();
    
    /**
     * @brief Construct OrbFit engine with custom configuration
     */
    explicit OrbFitEngine(const OrbFitConfig& config);
    
    // ========================================================================
    // Observation Management
    // ========================================================================
    
    /**
     * @brief Load observations from file (MPC format)
     * 
     * @param filename Path to observations file
     * @return Number of observations loaded
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
        return observations_;
    }
    
    /**
     * @brief Clear all observations
     */
    void clear_observations() {
        observations_.clear();
    }
    
    // ========================================================================
    // Orbit Determination
    // ========================================================================
    
    /**
     * @brief Set initial orbital elements (for refinement)
     * 
     * @param elements Initial orbital elements
     */
    void set_initial_orbit(const propagation::KeplerianElements& elements);
    
    /**
     * @brief Perform initial orbit determination (IOD) from observations
     * 
     * Determines initial orbit using Gauss method (minimum 3 observations).
     * 
     * @return Initial orbital elements
     * @throws std::runtime_error if insufficient observations
     */
    propagation::KeplerianElements initial_orbit_determination();
    
    /**
     * @brief Fit orbit to observations using differential correction
     * 
     * Performs iterative least-squares fit to refine orbit.
     * Requires initial orbit (from IOD or set_initial_orbit).
     * 
     * @return Orbit determination result with fitted orbit and statistics
     * @throws std::runtime_error if no initial orbit available
     */
    OrbitDeterminationResult fit_orbit();
    
    /**
     * @brief Get current best-fit orbit
     */
    const propagation::KeplerianElements& orbit() const {
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
    
    /**
     * @brief Compute ephemeris at multiple epochs
     * 
     * @param start_mjd Start epoch [MJD TDB]
     * @param end_mjd End epoch [MJD TDB]
     * @param step_days Step size [days]
     * @return Vector of Cartesian states at each epoch
     * @throws std::runtime_error if no orbit available
     */
    std::vector<propagation::CartesianElements> compute_ephemeris(
        double start_mjd,
        double end_mjd,
        double step_days);
    
    /**
     * @brief Propagate orbit to single epoch
     * 
     * @param target_mjd Target epoch [MJD TDB]
     * @return Orbital elements at target epoch
     */
    propagation::KeplerianElements propagate_to(double target_mjd);
    
    // ========================================================================
    // Close Approach Analysis
    // ========================================================================
    
    /**
     * @brief Find close approaches to planets
     * 
     * Searches for close approaches between start and end epochs.
     * 
     * @param start_mjd Start epoch [MJD TDB]
     * @param end_mjd End epoch [MJD TDB]
     * @return Vector of detected close approaches
     */
    std::vector<close_approach::CloseApproach> find_close_approaches(
        double start_mjd,
        double end_mjd);
    
    /**
     * @brief Compute MOID (Minimum Orbital Intersection Distance) with planet
     * 
     * @param planet Target planet
     * @return MOID [AU]
     */
    double compute_moid(ephemeris::CelestialBody planet);
    
    // ========================================================================
    // Configuration and Settings
    // ========================================================================
    
    /**
     * @brief Update configuration
     */
    void set_config(const OrbFitConfig& config) {
        config_ = config;
        update_propagator();
    }
    
    /**
     * @brief Get current configuration
     */
    const OrbFitConfig& config() const {
        return config_;
    }
    
    /**
     * @brief Enable/disable verbose output
     */
    void set_verbose(bool verbose) {
        config_.verbose = verbose;
    }
    
    // ========================================================================
    // Statistics and Diagnostics
    // ========================================================================
    
    /**
     * @brief Get last orbit determination result
     */
    const OrbitDeterminationResult& last_result() const {
        return last_result_;
    }
    
    /**
     * @brief Print orbit summary
     */
    void print_orbit_summary(std::ostream& os = std::cout) const;
    
    /**
     * @brief Print residuals summary
     */
    void print_residuals_summary(std::ostream& os = std::cout) const;
    
    /**
     * @brief Export orbit to file (various formats)
     * 
     * @param filename Output filename
     * @param format Format: "mpc", "oef", "json"
     */
    void export_orbit(const std::string& filename, 
                     const std::string& format = "oef");

private:
    // Update propagator with current configuration
    void update_propagator();
    
    // Validate observations (check timing, observatory codes, etc.)
    void validate_observations();
    
    // Configuration
    OrbFitConfig config_;
    
    // Observations
    std::vector<observations::OpticalObservation> observations_;
    
    // Current orbit
    propagation::KeplerianElements current_orbit_;
    bool has_orbit_ = false;
    
    // Last orbit determination result
    OrbitDeterminationResult last_result_;
    
    // Components
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<propagation::Propagator> propagator_;  // Shared because used by multiple components
    std::unique_ptr<orbit_determination::DifferentialCorrector> corrector_;
    std::unique_ptr<close_approach::CloseApproachDetector> ca_detector_;
};

} // namespace orbfit

#endif // ORBFIT_ORBFITENGINE_HPP
