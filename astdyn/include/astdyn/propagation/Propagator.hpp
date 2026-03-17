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
    bool include_moon = true;           ///< Include Moon separately
    bool include_asteroids = true;     ///< Include asteroid perturbations
    // Asteroid selection
    std::vector<int> include_asteroids_list = {};    ///< Specific asteroid numbers to include
    std::vector<int> exclude_asteroids_list = {};    ///< Specific asteroid numbers to exclude
    bool use_default_asteroid_set = true;          ///< If true, loads a predefined set of 17 massive asteroids (Pluto + 16)
    bool use_default_30_set = false;                ///< If true, loads the top 30 most massive asteroids (BC405 or similar)

    // Planetary perturbations to include (if include_planets=true)
    bool perturb_mercury = true;
    bool perturb_venus = true;
    bool perturb_earth = true;
    bool perturb_mars = true;
    bool perturb_jupiter = true;
    bool perturb_saturn = true;
    bool perturb_uranus = true;
    bool perturb_neptune = true;

    double central_body_gm = constants::GMS; ///< Central body GM [AU³/day²] (heliocentric)
    
    // Relativity PPN parameters (Default: GR)
    bool include_relativity = true;    ///< Include GR corrections (default true for precision)
    double ppn_beta = 1.0;
    double ppn_gamma = 1.0;

    // Harmonic corrections
    bool include_earth_j2 = true;     ///< Include Earth J2 perturbation
    bool include_sun_j2 = true;       ///< Include Sun J2 perturbation

    // Optional: Path to Asteroid SPK kernel (e.g. codes_300ast.bsp)
    // If empty, uses analytical approximation (AST17 constants)
    std::string asteroid_ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/sb441-n16.bsp";

    // Frame Settings
    bool integrate_in_ecliptic = true; ///< Forced true to match propagate_cartesian requirement.
    
    // Non-Gravitational Forces (Yarkovsky)
    bool include_yarkovsky = false;
    double yarkovsky_a2 = 0.0; // AU/d^2 at 1 AU (Tangential Acceleration Parameter)

    // Force Model Frame
    bool baricentric_integration = false; ///< Use SSB instead of Heliocentric origin
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
    Propagator(std::shared_ptr<Integrator> integrator,
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
     * @brief Get the underlying integrator
     */
    std::shared_ptr<Integrator> get_integrator() const { return integrator_; }
    
    /**
     * @brief Get the ephemeris provider
     */
    std::shared_ptr<ephemeris::PlanetaryEphemeris> get_ephemeris() const { return ephemeris_; }
    
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

    std::vector<Eigen::VectorXd> integrate_raw_au_batch(const Eigen::VectorXd& y0_au, double t0_mjd, const std::vector<double>& tf_mjds);
    
private:
    
    void update_force_cache(time::EpochTDB t);
    
    void setup_asteroid_perturbations();
    void setup_aas_parameters();

    Eigen::Vector3d compute_n_body_acceleration(const Eigen::Vector3d& position);
    
    Eigen::Vector3d compute_harmonic_acceleration(const Eigen::Vector3d& position, time::EpochTDB t);
    Eigen::Vector3d compute_earth_j2(const Eigen::Vector3d& pos, time::EpochTDB t);
    Eigen::Vector3d compute_sun_j2(const Eigen::Vector3d& pos);

    Eigen::Vector3d compute_non_gravitational_acceleration(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, time::EpochTDB t);

    Eigen::Vector3d relativistic_correction(const Eigen::Vector3d& position, 
                                           const Eigen::Vector3d& velocity) const;
    
    std::shared_ptr<Integrator> integrator_;
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<ephemeris::AsteroidPerturbations> asteroids_;
    PropagatorSettings settings_;
    Eigen::Matrix3d mat_ecl_; ///< Cached rotation matrix J2000 -> Ecliptic

    // Cache for force evaluations
    time::EpochTDB last_t_cache_;
    Eigen::Vector3d last_sun_pos_bary_cache_;
    struct PlanetStateCache {
        ephemeris::CelestialBody body;
        Eigen::Vector3d pos_bary_au;
    };
    PlanetStateCache planet_cache_[10]; // Fixed array to avoid allocations
    int num_planets_cached_ = 0;
    bool cache_valid_ = false;
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
// Template Implementations (AstDyn beta 0.9 wrappers)
// ============================================================================

template <typename Frame>
physics::CartesianStateTyped<Frame> Propagator::propagate_cartesian(
    const physics::CartesianStateTyped<Frame>& initial,
    time::EpochTDB target_time) 
{
    // === ENTRY POINT SI/AnyFrame -> AU/ECLIPJ2000 ===
    // 1. Transform to ECLIPJ2000 and AU/day for numeric stability
    auto initial_ecl = initial.template cast_frame<core::ECLIPJ2000>();
    Eigen::VectorXd y0_au = initial_ecl.to_eigen_au_aud();
    
    // 2. Perform raw integration in AU/day/Ecliptic/MJD
    Eigen::VectorXd yf_au = integrate_raw_au(y0_au, initial.epoch.mjd(), target_time.mjd());
    
    // === EXIT POINT AU/ECLIPJ2000 -> SI/AnyFrame ===
    // 3. Transform back to SI and the requested Frame
    auto final_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>::from_au_aud(
        target_time, yf_au, initial.gm);
    
    return final_ecl.template cast_frame<Frame>();
}

template <typename Frame>
std::vector<physics::CartesianStateTyped<Frame>> Propagator::propagate_ephemeris(
    const physics::CartesianStateTyped<Frame>& initial,
    const std::vector<time::EpochTDB>& target_times) 
{
    if (target_times.empty()) return {};

    // 1. Convert initial state to AU/Ecliptic
    auto initial_ecl = initial.template cast_frame<core::ECLIPJ2000>();
    Eigen::VectorXd y0_au = initial_ecl.to_eigen_au_aud();

    // 2. Prepare target MJD list
    std::vector<double> mjds;
    mjds.reserve(target_times.size());
    for (const auto& t : target_times) mjds.push_back(t.mjd());

    // 3. Perform batch integration
    auto results_au = integrate_raw_au_batch(y0_au, initial.epoch.mjd(), mjds);

    // 4. Convert back to requested frame and SI
    std::vector<physics::CartesianStateTyped<Frame>> results;
    results.reserve(target_times.size());
    for (size_t i = 0; i < target_times.size(); ++i) {
        auto state_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>::from_au_aud(
            target_times[i], results_au[i], initial.gm);
        results.push_back(state_ecl.template cast_frame<Frame>());
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
