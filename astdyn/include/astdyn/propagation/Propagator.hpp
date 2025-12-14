/**
 * @file Propagator.hpp
 * @brief Orbital propagation with n-body dynamics
 * 
 * This module provides high-level orbital propagation using
 * numerical integrators. It handles:
 * - Two-body dynamics (Keplerian)
 * - N-body perturbations (planets)
 * - Relativistic corrections (optional)
 * 
 * The propagator converts orbital elements to Cartesian state,
 * integrates the equations of motion, and converts back.
 */

#ifndef ASTDYN_PROPAGATOR_HPP
#define ASTDYN_PROPAGATOR_HPP

#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/AsteroidPerturbations.hpp"
#include <memory>

namespace astdyn::propagation {

/**
 * @brief Propagation settings
 */
struct PropagatorSettings {
    bool include_planets = true;        ///< Include planetary perturbations
    bool include_relativity = false;    ///< Include GR corrections
    bool include_moon = false;          ///< Include Moon separately
    bool include_asteroids = false;     ///< Include asteroid perturbations (AST17)
    
    // Planetary perturbations to include (if include_planets=true)
    bool perturb_mercury = false;
    bool perturb_venus = true;
    bool perturb_earth = true;
    bool perturb_mars = true;
    bool perturb_jupiter = true;
    bool perturb_saturn = true;
    bool perturb_uranus = false;
    bool perturb_neptune = false;
    
    double central_body_gm = constants::GMS; ///< Central body GM [AU³/day²] (heliocentric)
    
    // Relativity PPN parameters (Default: GR)
    double ppn_beta = 1.0;
    double ppn_gamma = 1.0;

    // Optional: Path to Asteroid SPK kernel (e.g. codes_300ast.bsp)
    // If empty, uses analytical approximation (AST17 constants)
    std::string asteroid_ephemeris_file = "";
};

/**
 * @brief Orbital propagator class
 * 
 * Propagates orbits using numerical integration of equations of motion.
 * Handles n-body perturbations from planets.
 * 
 * Example usage:
 * @code
 *   // Create ephemeris and integrator
 *   auto ephem = std::make_shared<PlanetaryEphemeris>();
 *   auto integrator = std::make_unique<RKF78Integrator>(0.1, 1e-12);
 *   
 *   // Create propagator
 *   Propagator prop(std::move(integrator), ephem);
 *   
 *   // Propagate orbit
 *   KeplerianElements initial = ...;
 *   KeplerianElements final = prop.propagate_keplerian(initial, target_mjd);
 * @endcode
 */
class Propagator {
public:
    /**
     * @brief Construct propagator
     * 
     * @param integrator Numerical integrator (RK4, RKF78, etc.)
     * @param ephemeris Planetary ephemeris for perturbations
     * @param settings Propagation settings
     */
    Propagator(std::unique_ptr<Integrator> integrator,
              std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris,
              const PropagatorSettings& settings = PropagatorSettings());
    
    /**
     * @brief Propagate Keplerian elements to target epoch
     * 
     * @param initial Initial Keplerian elements
     * @param target_mjd_tdb Target epoch (MJD TDB)
     * @return Keplerian elements at target epoch
     */
    KeplerianElements propagate_keplerian(const KeplerianElements& initial,
                                         double target_mjd_tdb);
    
    /**
     * @brief Propagate Cartesian state to target epoch
     * 
     * @param initial Initial Cartesian state
     * @param target_mjd_tdb Target epoch (MJD TDB)
     * @return Cartesian state at target epoch
     */
    CartesianElements propagate_cartesian(const CartesianElements& initial,
                                         double target_mjd_tdb);
    
    /**
     * @brief Propagate and generate ephemeris at multiple epochs
     * 
     * @param initial Initial state
     * @param epochs_mjd_tdb Target epochs
     * @return States at each epoch
     */
    std::vector<CartesianElements> propagate_ephemeris(
        const CartesianElements& initial,
        const std::vector<double>& epochs_mjd_tdb);
    
    /**
     * @brief Get integrator statistics from last propagation
     */
    const IntegrationStatistics& statistics() const {
        return integrator_->statistics();
    }
    
    /**
     * @brief Update propagation settings
     */
    void set_settings(const PropagatorSettings& settings) {
        settings_ = settings;
    }
    
    const PropagatorSettings& settings() const { return settings_; }
    PropagatorSettings& settings() { return settings_; }
    
    /**
     * @brief Compute accelerations for equations of motion
     * 
     * Computes d²r/dt² from gravitational forces.
     * Exposed for use in State Transition Matrix calculations.
     * 
     * @param t Time (MJD TDB)
     * @param state State vector [x, y, z, vx, vy, vz] in AU, AU/day
     * @return Derivative [vx, vy, vz, ax, ay, az]
     */
    Eigen::VectorXd compute_derivatives(double t, const Eigen::VectorXd& state);
    
private:
    
    /**
     * @brief Compute two-body acceleration (central body only)
     * 
     * @param position Position vector [AU]
     * @return Acceleration vector [AU/day²]
     */
    Eigen::Vector3d two_body_acceleration(const Eigen::Vector3d& position) const;
    
    /**
     * @brief Compute planetary perturbations (Heliocentric frame corrected)
     */
    Eigen::Vector3d planetary_perturbations(const Eigen::Vector3d& position, 
                                          double mjd_tdb,
                                          const Eigen::Vector3d& sun_pos_bary);
                                          
    /**
     * @brief Compute asteroid perturbations (Heliocentric frame corrected)
     */
    Eigen::Vector3d asteroid_perturbations(const Eigen::Vector3d& position, 
                                         double mjd_tdb,
                                         const Eigen::Vector3d& sun_pos_bary);
                                         
    // Relativistic correction (PPN)
    Eigen::Vector3d relativistic_correction(const Eigen::Vector3d& position, 
                                          const Eigen::Vector3d& velocity) const;
    
    std::unique_ptr<Integrator> integrator_;
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<ephemeris::AsteroidPerturbations> asteroids_;
    PropagatorSettings settings_;
};

/**
 * @brief Simple two-body propagator (analytical)
 * 
 * Uses Keplerian orbital mechanics for fast propagation
 * without perturbations. Useful for short arcs or preliminary orbits.
 */
class TwoBodyPropagator {
public:
    /**
     * @brief Propagate using Keplerian motion
     * 
     * @param initial Initial Keplerian elements
     * @param target_mjd_tdb Target epoch
     * @return Keplerian elements at target (only M changes)
     */
    static KeplerianElements propagate(const KeplerianElements& initial,
                                       double target_mjd_tdb);
    
    /**
     * @brief Compute mean anomaly at epoch from initial state
     * 
     * @param initial Initial elements
     * @param target_mjd_tdb Target epoch
     * @return Mean anomaly at target epoch [rad]
     */
    static double mean_anomaly_at_epoch(const KeplerianElements& initial,
                                        double target_mjd_tdb);
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATOR_HPP
