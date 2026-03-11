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

Propagator::Propagator(std::shared_ptr<Integrator> integrator,
                      std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris,
                      const PropagatorSettings& settings)
    : integrator_(std::move(integrator)),
      ephemeris_(std::move(ephemeris)),
      settings_(settings) {
    
    mat_ecl_ = coordinates::ReferenceFrame::j2000_to_ecliptic();

    if (settings_.include_asteroids) {
        asteroids_ = std::make_shared<ephemeris::AsteroidPerturbations>();
        if (!settings_.asteroid_ephemeris_file.empty()) {
            asteroids_->loadSPK(settings_.asteroid_ephemeris_file);
        }
    }
}

Eigen::VectorXd Propagator::compute_derivatives(time::EpochTDB t, const Eigen::VectorXd& state) {
    // State is strictly in AU and AU/day
    Eigen::Vector3d position = state.head<3>();
    Eigen::Vector3d velocity = state.tail<3>();
    
    // Get Sun Barycentric Position
    auto sun_pos_bary = ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(t);
    Eigen::Vector3d sun_pos_bary_vec = sun_pos_bary.to_eigen_si() / (constants::AU * 1000.0);
    
    // Rotate to Ecliptic if requested
    if (settings_.integrate_in_ecliptic) {
        sun_pos_bary_vec = mat_ecl_ * sun_pos_bary_vec;
    }
    
    // Compute total acceleration (in AU/day^2)
    Eigen::Vector3d acc = two_body_acceleration(position);
    
    if (settings_.include_planets) {
        acc += planetary_perturbations(position, t, sun_pos_bary_vec);
    }
    
    if (settings_.include_asteroids && asteroids_) {
        acc += asteroid_perturbations(position, t, sun_pos_bary_vec);
    }
    
    if (settings_.include_relativity) {
        acc += relativistic_correction(position, velocity);
    }
    
    // Yarkovsky Effect
    if (settings_.include_yarkovsky && std::abs(settings_.yarkovsky_a2) > 0.0) {
        double r_au = position.norm();
        double v_au_d = velocity.norm();
        if (r_au > 1e-4 && v_au_d > 1e-6) {
            // a2 is in AU/d^2.
            double acc_val = settings_.yarkovsky_a2 / (r_au * r_au);
            acc += acc_val * (velocity / v_au_d);
        }
    }
    
    // d[r, v]/dt = [v, a]
    Eigen::VectorXd xdot(6);
    xdot.head<3>() = velocity;
    xdot.tail<3>() = acc;
    
    return xdot;
}

Eigen::Vector3d Propagator::two_body_acceleration(const Eigen::Vector3d& position) const {
    double r = position.norm();
    double r3 = r * r * r;
    // central_body_gm is configured in AU^3/d^2
    double mu = settings_.central_body_gm;
    return -mu * position / r3;
}

Eigen::Vector3d Propagator::planetary_perturbations(const Eigen::Vector3d& position,
                                                  time::EpochTDB t,
                                                  const Eigen::Vector3d& sun_pos_bary) {
    Eigen::Vector3d perturbation = Eigen::Vector3d::Zero();
    auto provider = ephemeris::PlanetaryEphemeris::getProvider();
    const double au_m = constants::AU * 1000.0;

    auto add_planet_perturbation = [&]( CelestialBody planet, double planet_gm_au) {
        // Direct AU barycentric position
        auto p_pos_state = provider->getPosition(planet, t);
        Eigen::Vector3d p_pos_bary_au = p_pos_state.to_eigen_si() / (constants::AU * 1000.0);

        if (settings_.integrate_in_ecliptic) {
            p_pos_bary_au = mat_ecl_ * p_pos_bary_au;
        }
        
        Eigen::Vector3d planet_pos_helio_au = p_pos_bary_au - sun_pos_bary;
        Eigen::Vector3d delta = planet_pos_helio_au - position;
        
        const double d_norm = delta.norm();
        const double d3 = d_norm * d_norm * d_norm;
        const double p_dist = planet_pos_helio_au.norm();
        const double p_dist3 = p_dist * p_dist * p_dist;
        
        perturbation += planet_gm_au * (delta / d3 - planet_pos_helio_au / p_dist3);
    };
    
    if (settings_.perturb_mercury) add_planet_perturbation(CelestialBody::MERCURY, constants::GM_MERCURY_AU);
    if (settings_.perturb_venus)   add_planet_perturbation(CelestialBody::VENUS,   constants::GM_VENUS_AU);
    if (settings_.perturb_earth)   add_planet_perturbation(CelestialBody::EARTH,   constants::GM_EARTH_AU);
    if (settings_.perturb_mars)    add_planet_perturbation(CelestialBody::MARS,    constants::GM_MARS_AU);
    if (settings_.perturb_jupiter) add_planet_perturbation(CelestialBody::JUPITER, constants::GM_JUPITER_AU);
    if (settings_.perturb_saturn)  add_planet_perturbation(CelestialBody::SATURN,  constants::GM_SATURN_AU);
    if (settings_.perturb_uranus)  add_planet_perturbation(CelestialBody::URANUS,  constants::GM_URANUS_AU);
    if (settings_.perturb_neptune) add_planet_perturbation(CelestialBody::NEPTUNE, constants::GM_NEPTUNE_AU);
    
    if (settings_.include_moon) {
        double moon_gm_au = physics::GravitationalParameter::from_km3_s2(constants::GM_MOON).to_au3_d2();
        add_planet_perturbation(CelestialBody::MOON, moon_gm_au);
    }
    
    return perturbation;
}

Eigen::Vector3d Propagator::asteroid_perturbations(const Eigen::Vector3d& position,
                                                 time::EpochTDB t,
                                                 const Eigen::Vector3d& sun_pos_bary) {
    // Asteroids library assumes AU and returns AU/day^2.
    // Use the raw interface to avoid state creation overhead during integration.
    return asteroids_->computePerturbationRaw(position, t.mjd(), sun_pos_bary, settings_.integrate_in_ecliptic);
}

Eigen::Vector3d Propagator::relativistic_correction(const Eigen::Vector3d& position,
                                                   const Eigen::Vector3d& velocity) const {
    double r = position.norm();
    double v2 = velocity.squaredNorm();
    double c_au_d = physics::SpeedOfLight::to_au_d();
    double c2 = c_au_d * c_au_d;
    double mu = settings_.central_body_gm;
    
    double beta = settings_.ppn_beta;
    double gamma = settings_.ppn_gamma;
    
    double term1_coeff = 2.0 * (beta + gamma) * mu / r - gamma * v2;
    double rv_dot = position.dot(velocity);
    
    Eigen::Vector3d term1 = term1_coeff * position;
    Eigen::Vector3d term2 = 2.0 * (1.0 + gamma) * rv_dot * velocity;
    
    return (mu / (r * r * r * c2)) * (term1 + term2);
}

Eigen::VectorXd Propagator::integrate_raw_au(const Eigen::VectorXd& y0_au, double t0_mjd, double tf_mjd) {
    DerivativeFunction f = [this](double t_val, const Eigen::VectorXd& y) {
        return compute_derivatives(time::EpochTDB::from_mjd(t_val), y);
    };
    
    return integrator_->integrate(f, y0_au, t0_mjd, tf_mjd);
}

std::vector<Eigen::VectorXd> Propagator::integrate_raw_au_batch(const Eigen::VectorXd& y0_au, double t0_mjd, const std::vector<double>& tf_mjds) {
    DerivativeFunction f = [this](double t_val, const Eigen::VectorXd& y) {
        return compute_derivatives(time::EpochTDB::from_mjd(t_val), y);
    };
    
    return integrator_->integrate_at(f, y0_au, t0_mjd, tf_mjds);
}

} // namespace astdyn::propagation
