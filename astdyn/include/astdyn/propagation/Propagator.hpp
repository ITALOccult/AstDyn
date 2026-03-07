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
#include "astdyn/core/physics_state.hpp"
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

    // Frame Settings
    bool integrate_in_ecliptic = true; ///< If true, rotate all perturbations to Ecliptic J2000.
    
    // Non-Gravitational Forces (Yarkovsky)
    bool include_yarkovsky = false;
    double yarkovsky_a2 = 0.0; // AU/d^2 at 1 AU (Tangential Acceleration Parameter)
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
     * @brief Propagate Keplerian state
     * 
     * @param initial Initial Keplerian state (strongly-typed)
     * @param target_time Target epoch (time::EpochTDB)
     * @return Keplerian state at target epoch in the same reference frame
     */
    template <typename Frame>
    physics::KeplerianStateTyped<Frame> propagate_keplerian(
        const physics::KeplerianStateTyped<Frame>& initial,
        time::EpochTDB target_time);
    
    /**
     * @brief Propagate Cartesian state to target epoch
     * 
     * @param initial Initial Cartesian state (strongly-typed, SI units internally)
     * @param target_time Target epoch (time::EpochTDB)
     * @return Cartesian state at target epoch
     */
    template <typename Frame>
    physics::CartesianStateTyped<Frame> propagate_cartesian(
        const physics::CartesianStateTyped<Frame>& initial,
        time::EpochTDB target_time);
    
    template <typename Frame>
    std::vector<physics::CartesianStateTyped<Frame>> propagate_ephemeris(
        const physics::CartesianStateTyped<Frame>& initial,
        const std::vector<time::EpochTDB>& target_times);
    
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
    Eigen::VectorXd compute_derivatives(time::EpochTDB t, const Eigen::VectorXd& state);
    
    /**
     * @brief Raw core integration using untyped Eigen arrays.
     * 
     * State is strictly maintained in AU and AU/day to ensure integrator
     * precision and preventing floating point tolerance blowups seen with meters.
     */
    Eigen::VectorXd integrate_raw_au(const Eigen::VectorXd& y0_au, double t0_mjd, double tf_mjd);
    
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
                                          time::EpochTDB t,
                                          const Eigen::Vector3d& sun_pos_bary);
                                          
    /**
     * @brief Compute asteroid perturbations (Heliocentric frame corrected)
     */
    Eigen::Vector3d asteroid_perturbations(const Eigen::Vector3d& position, 
                                         time::EpochTDB t,
                                         const Eigen::Vector3d& sun_pos_bary);
                                         
    // Relativistic correction (PPN)
    Eigen::Vector3d relativistic_correction(const Eigen::Vector3d& position, 
                                          const Eigen::Vector3d& velocity) const;
    
    std::unique_ptr<Integrator> integrator_;
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<ephemeris::AsteroidPerturbations> asteroids_;
    PropagatorSettings settings_;
};

class TwoBodyPropagator {
public:
    /**
     * @brief Propagate using Keplerian motion
     * 
     * @param initial Initial Keplerian elements (typed)
     * @param target_time Target epoch
     * @return Keplerian elements at target (only M changes)
     */
    template <typename Frame>
    static physics::KeplerianStateTyped<Frame> propagate(
        const physics::KeplerianStateTyped<Frame>& initial,
        time::EpochTDB target_time)
    {
        double dt_days = target_time.mjd() - initial.epoch.mjd();
        double a_au = initial.a.to_au();
        double mu_au_d2 = initial.gm.to_au3_d2();
        
        // n = sqrt(mu / a^3) in rad/day
        double n = std::sqrt(mu_au_d2 / (a_au * a_au * a_au));
        
        double m_new = initial.M.to_rad() + n * dt_days;
        
        // Normalize M to [0, 2pi]
        m_new = std::fmod(m_new, constants::TWO_PI);
        if (m_new < 0) m_new += constants::TWO_PI;
        
        auto result = initial;
        result.epoch = target_time;
        result.M = astrometry::Angle::from_rad(m_new);
        return result;
    }
    
    /**
     * @brief Compute mean anomaly at epoch from initial state
     * 
     * @param initial Initial elements (typed)
     * @param target_time Target epoch
     * @return Mean anomaly at target epoch [rad]
     */
    template <typename Frame>
    static double mean_anomaly_at_epoch(
        const physics::KeplerianStateTyped<Frame>& initial,
        time::EpochTDB target_time)
    {
        return propagate(initial, target_time).M.to_rad();
    }
};

// ============================================================================
// Template Implementations (AstDyn 3.0 wrappers)
// ============================================================================

template <typename Frame>
physics::CartesianStateTyped<Frame> Propagator::propagate_cartesian(
    const physics::CartesianStateTyped<Frame>& initial,
    time::EpochTDB target_time) 
{
    // 1. Transform SI inputs strictly to AU, AU/day for numeric stability
    Eigen::VectorXd y0_au(6);
    y0_au.head<3>() = initial.position.to_eigen_si() / (constants::AU * 1000.0);
    y0_au.tail<3>() = initial.velocity.to_eigen_si() * 86400.0 / (constants::AU * 1000.0);
    
    // 2. Perform raw integration (no types, strict mathematical validity)
    double t0 = initial.epoch.mjd();
    double tf = target_time.mjd();
    Eigen::VectorXd yf_au = integrate_raw_au(y0_au, t0, tf);
    
    // 3. Immediately bring back into the safety of the typed 3.0 domain (SI outputs)
    Eigen::Vector3d pos_si = yf_au.head<3>() * (constants::AU * 1000.0);
    Eigen::Vector3d vel_si = yf_au.tail<3>() * (constants::AU * 1000.0) / 86400.0;
    
    return physics::CartesianStateTyped<Frame>::from_si(
        time::EpochTDB::from_mjd(tf),
        pos_si.x(), pos_si.y(), pos_si.z(),
        vel_si.x(), vel_si.y(), vel_si.z(),
        initial.gm.to_m3_s2()
    );
}

template <typename Frame>
std::vector<physics::CartesianStateTyped<Frame>> Propagator::propagate_ephemeris(
    const physics::CartesianStateTyped<Frame>& initial,
    const std::vector<time::EpochTDB>& target_times) 
{
    std::vector<physics::CartesianStateTyped<Frame>> results;
    results.reserve(target_times.size());
    for (const auto& t : target_times) {
        results.push_back(propagate_cartesian(initial, t));
    }
    return results;
}

template <typename Frame>
physics::KeplerianStateTyped<Frame> Propagator::propagate_keplerian(
    const physics::KeplerianStateTyped<Frame>& initial,
    time::EpochTDB target_time) 
{
    // 1. Convert to Cartesian (Type Safe)
    auto cart_initial = keplerian_to_cartesian(initial);
    
    // 2. Propagate Cartesian
    auto cart_final = propagate_cartesian(cart_initial, target_time);
    
    // 3. Convert back to Keplerian (Type Safe)
    return cartesian_to_keplerian<Frame>(cart_final);
}

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATOR_HPP
