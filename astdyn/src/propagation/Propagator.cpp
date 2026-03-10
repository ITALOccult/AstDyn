/**
 * @file Propagator.cpp
 * @brief Implementation of orbital propagation
 */

#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "src/core/frame_tags.hpp"
#include "src/propagation/kepler_propagator.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::propagation {
using namespace astdyn::constants;

using ephemeris::CelestialBody;

// ============================================================================
// Propagator Implementation
// ============================================================================

Propagator::Propagator(std::shared_ptr<Integrator<StateAU, DerivativeAU>> integrator,
                       std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris,
                       const PropagatorSettings& settings)
    : integrator_(std::move(integrator)),
      ephemeris_(std::move(ephemeris)),
      settings_(settings) {
    if (settings_.include_asteroids) {
        asteroids_ = std::make_shared<ephemeris::AsteroidPerturbations>();
        if (!settings_.asteroid_ephemeris_file.empty()) {
            asteroids_->loadSPK(settings_.asteroid_ephemeris_file);
        }
    }
}

DerivativeAU Propagator::compute_derivatives_au(time::EpochTDB t, const StateAU& state) {
    Eigen::Vector3d pos_au(state.x, state.y, state.z);
    Eigen::Vector3d vel_au(state.vx, state.vy, state.vz);

    Eigen::Vector3d acc_au = two_body_acceleration_au(state);

    if (settings_.include_planets) {
        acc_au += planetary_perturbations_au(state, t);
    }

    // Return type-safe derivative
    return DerivativeAU{
        physics::VelocityAUD::from_au_d(state.vx),
        physics::VelocityAUD::from_au_d(state.vy),
        physics::VelocityAUD::from_au_d(state.vz),
        physics::AccelerationAUD2::from_au_d2(acc_au.x()),
        physics::AccelerationAUD2::from_au_d2(acc_au.y()),
        physics::AccelerationAUD2::from_au_d2(acc_au.z())
    };
}

Eigen::Vector3d Propagator::two_body_acceleration_au(const StateAU& state) const {
    double r2 = state.x*state.x + state.y*state.y + state.z*state.z;
    double r3 = r2 * std::sqrt(r2);
    return -state.gm.to_au3_d2() * Eigen::Vector3d(state.x, state.y, state.z) / r3;
}

Eigen::Vector3d Propagator::planetary_perturbations_au(const StateAU& state, time::EpochTDB t) {
    Eigen::Vector3d perturbation = Eigen::Vector3d::Zero();
    Eigen::Vector3d ast_pos_au(state.x, state.y, state.z);
    
    auto add_planet = [&](ephemeris::CelestialBody planet, double gm_au) {
        auto p_pos_si = ephemeris_->getPosition(planet, t);
        auto p_pos_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(p_pos_si, t);
        Eigen::Vector3d p_pos_au = p_pos_ecl.to_eigen_si() / (constants::AU * 1000.0);
        
        Eigen::Vector3d delta = p_pos_au - ast_pos_au;
        double d3 = std::pow(delta.squaredNorm(), 1.5);
        double p3 = std::pow(p_pos_au.squaredNorm(), 1.5);
        
        perturbation += gm_au * (delta / d3 - p_pos_au / p3);
    };

    if (settings_.perturb_mercury) add_planet(CelestialBody::MERCURY, constants::GM_MERCURY_AU);
    if (settings_.perturb_venus)   add_planet(CelestialBody::VENUS,   constants::GM_VENUS_AU);
    if (settings_.perturb_earth)   add_planet(CelestialBody::EARTH,   constants::GM_EARTH_AU);
    if (settings_.perturb_mars)    add_planet(CelestialBody::MARS,    constants::GM_MARS_AU);
    if (settings_.perturb_jupiter) add_planet(CelestialBody::JUPITER, constants::GM_JUPITER_AU);
    if (settings_.perturb_saturn)  add_planet(CelestialBody::SATURN,  constants::GM_SATURN_AU);
    if (settings_.perturb_uranus)  add_planet(CelestialBody::URANUS,  constants::GM_URANUS_AU);
    if (settings_.perturb_neptune) add_planet(CelestialBody::NEPTUNE, constants::GM_NEPTUNE_AU);

    if (settings_.include_moon) {
        double moon_gm_au = physics::GravitationalParameter::from_km3_s2(constants::GM_MOON).to_au3_d2();
        add_planet(CelestialBody::MOON, moon_gm_au);
    }

    return perturbation;
}

} // namespace astdyn::propagation
