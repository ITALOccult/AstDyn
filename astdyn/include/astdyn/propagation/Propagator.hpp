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
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/core/physics_state_au.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/AdamsIntegrator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/AsteroidPerturbations.hpp"
#include <memory>
#include <vector>

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
    bool integrate_in_ecliptic = false; ///< If true, rotate all perturbations to Ecliptic J2000.
    
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
using StateAU = physics::CartesianStateAU<core::ECLIPJ2000>;
using DerivativeAU = physics::DerivativeAU;

/**
 * @brief High-precision orbital propagator class using Type-Safe AU units.
 * 
 * All internal dynamics are computed in AU and AU/day to ensure 
 * numerical stability and prevent unit scaling errors.
 */
class Propagator {
public:
    Propagator(std::shared_ptr<Integrator<StateAU, DerivativeAU>> integrator,
               std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris,
               const PropagatorSettings& settings = PropagatorSettings());

    /**
     * @brief Propagate Cartesian state to target epoch.
     * PUBLIC SI INTERFACE: Converters to/from AU are handled internally.
     */
    template <typename Frame>
    physics::CartesianStateTyped<Frame> propagate_cartesian(
        const physics::CartesianStateTyped<Frame>& state,
        time::EpochTDB target_epoch) 
    {
        // 1. SI -> AU (Always forced to ECLIPJ2000 for integration)
        auto state_au = StateAU::from_si(state.template cast_frame<core::ECLIPJ2000>()); 
        
        // 2. Integration
        integrator_->reset_statistics();
        auto final_state_au = integrator_->integrate(
            [this](time::EpochTDB t, const StateAU& y) { return this->compute_derivatives_au(t, y); },
            state_au,
            state.epoch.mjd(),
            target_epoch.mjd()
        );
        
        // 3. AU -> SI (Back to user's Frame)
        return final_state_au.to_si().template cast_frame<Frame>();
    }

    template <typename Frame>
    physics::KeplerianStateTyped<Frame> propagate_keplerian(
        const physics::KeplerianStateTyped<Frame>& initial,
        time::EpochTDB target_time) 
    {
        auto cart_initial = keplerian_to_cartesian(initial);
        auto cart_final = propagate_cartesian(cart_initial, target_time);
        return cartesian_to_keplerian<Frame>(cart_final);
    }

    template <typename Frame>
    std::vector<physics::CartesianStateTyped<Frame>> propagate_ephemeris(
        const physics::CartesianStateTyped<Frame>& initial,
        const std::vector<time::EpochTDB>& target_times) 
    {
        std::vector<physics::CartesianStateTyped<Frame>> results;
        results.reserve(target_times.size());
        physics::CartesianStateTyped<Frame> current = initial;
        for (const auto& t : target_times) {
            current = (t.mjd() == current.epoch.mjd()) ? current : propagate_cartesian(current, t);
            results.push_back(current);
        }
        return results;
    }

    /**
     * @brief THE PHYSICS KERNEL: Compute state derivative in AU/day and AU/day^2.
     */
    DerivativeAU compute_derivatives_au(time::EpochTDB t, const StateAU& state);

    const IntegrationStatistics& statistics() const { return integrator_->statistics(); }
    void set_settings(const PropagatorSettings& s) { settings_ = s; }
    const PropagatorSettings& settings() const { return settings_; }
    PropagatorSettings& settings() { return settings_; }
    
    std::shared_ptr<ephemeris::PlanetaryEphemeris> get_ephemeris() const { return ephemeris_; }

private:
    std::shared_ptr<Integrator<StateAU, DerivativeAU>> integrator_;
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<ephemeris::AsteroidPerturbations> asteroids_;
    PropagatorSettings settings_;

    Eigen::Vector3d two_body_acceleration_au(const StateAU& state) const;
    Eigen::Vector3d planetary_perturbations_au(const StateAU& state, time::EpochTDB t);
};

/**
 * @brief Analytical 2-body propagation (Typed)
 */
class TwoBodyPropagator {
public:
    template <typename Frame>
    static physics::KeplerianStateTyped<Frame> propagate(
        const physics::KeplerianStateTyped<Frame>& initial, time::EpochTDB target_time) {
        double dt_days = target_time.mjd() - initial.epoch.mjd();
        double a_au = initial.a.to_au();
        double n = std::sqrt(initial.gm.to_au3_d2() / (a_au * a_au * a_au));
        double m_new = std::fmod(initial.M.to_rad() + n * dt_days, constants::TWO_PI);
        if (m_new < 0) m_new += constants::TWO_PI;
        auto res = initial;
        res.epoch = target_time;
        res.M = astrometry::Angle::from_rad(m_new);
        return res;
    }
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATOR_HPP
