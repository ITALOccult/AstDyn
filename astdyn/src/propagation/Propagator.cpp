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

Propagator::Propagator(std::unique_ptr<Integrator> integrator,
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

Eigen::VectorXd Propagator::compute_derivatives(time::EpochTDB t, const Eigen::VectorXd& state) {
    // State is strictly in AU and AU/day
    Eigen::Vector3d position = state.head<3>();
    Eigen::Vector3d velocity = state.tail<3>();
    
    // Get Sun Barycentric Position (KM -> AU and rotate if needed)
    auto sun_pos_bary = ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(t);
    Eigen::Vector3d sun_pos_bary_vec = sun_pos_bary.to_eigen() / (constants::AU * 1000.0);
    
    // Rotate to Ecliptic if needed
    if (settings_.integrate_in_ecliptic) {
        sun_pos_bary_vec = coordinates::ReferenceFrame::j2000_to_ecliptic() * sun_pos_bary_vec;
    }
    
    // Compute total acceleration (in AU/day^2)
    Eigen::Vector3d acc = two_body_acceleration(position);
    
    if (settings_.include_planets) {
        acc += planetary_perturbations(position, t, sun_pos_bary_vec);
    }
    
    if (settings_.include_asteroids && asteroids_) {
        // Asteroids library assumes AU and returns AU/day^2 (or close to it based on previous scaling)
        // The asteroid perturbation function needs to handle the rotation if integrate_in_ecliptic is true.
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
    auto mat_ecl = coordinates::ReferenceFrame::j2000_to_ecliptic();

    auto add_planet_perturbation = [&]( CelestialBody planet, double planet_gm_au) {
        auto planet_state_bary = ephemeris_->getPosition(planet, t);
        Eigen::Vector3d p_pos_bary_au = planet_state_bary.to_eigen() / (constants::AU * 1000.0);
        
        if (settings_.integrate_in_ecliptic) {
            p_pos_bary_au = mat_ecl * p_pos_bary_au;
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
        double moon_gm_au = constants::GM_MOON * 1e9 * std::pow(86400.0, 2) / std::pow(constants::AU * 1000.0, 3);
        add_planet_perturbation(CelestialBody::MOON, moon_gm_au);
    }
    
    return perturbation;
}

Eigen::Vector3d Propagator::asteroid_perturbations(const Eigen::Vector3d& position,
                                                 time::EpochTDB t,
                                                 const Eigen::Vector3d& sun_pos_bary) {
    if (!asteroids_ || !settings_.include_asteroids) return Eigen::Vector3d::Zero();
    
    // The AsteroidPerturbations::computePerturbation expects SSB sun_pos and GCRF position.
    // However, our internal 'position' might be Ecliptic.
    // And asteroids_ internal ephemeris is GCRF.
    
    if (!settings_.integrate_in_ecliptic) {
         return asteroids_->computePerturbation(position, t.mjd(), sun_pos_bary);
    }
    
    // If in Ecliptic, rotate position back to Equatorial for the perturbation engine
    auto mat_eq = coordinates::ReferenceFrame::ecliptic_to_j2000();
    auto mat_ecl = coordinates::ReferenceFrame::j2000_to_ecliptic();
    
    Eigen::Vector3d pos_eq = mat_eq * position;
    Eigen::Vector3d sun_eq = mat_eq * sun_pos_bary;
    
    Eigen::Vector3d acc_eq = asteroids_->computePerturbation(pos_eq, t.mjd(), sun_eq);
    
    // Rotate result back to Ecliptic
    return mat_ecl * acc_eq;
}

Eigen::Vector3d Propagator::relativistic_correction(const Eigen::Vector3d& position,
                                                   const Eigen::Vector3d& velocity) const {
    double r = position.norm();
    double v2 = velocity.squaredNorm();
    double c_au_d = (constants::C_LIGHT * 1000.0) * 86400.0 / (constants::AU * 1000.0);
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

// Analytical TwoBodyPropagator
KeplerianElements TwoBodyPropagator::propagate(const KeplerianElements& initial,
                                               time::EpochTDB target_time) {
    using namespace astdyn::propagation;
    using ::astdyn::core::GCRF;
    using ::astdyn::types::TimedState;

    std::array<double, 6> arr = {
        initial.semi_major_axis, initial.eccentricity, initial.inclination,
        initial.longitude_ascending_node, initial.argument_perihelion, initial.mean_anomaly
    };
    const auto state_kep = astdyn::types::OrbitalState<GCRF, astdyn::types::KeplerianTag>(arr);
    
    KeplerPropagator<GCRF> engine(initial.gravitational_parameter);
    const auto result = engine.propagate(TimedState<GCRF, astdyn::types::KeplerianTag>(state_kep, initial.epoch), target_time);

    if (!result.has_value()) return initial;

    KeplerianElements final = initial;
    final.epoch = target_time;
    final.mean_anomaly = result->state.m_anomaly();
    
    return final;
}

double TwoBodyPropagator::mean_anomaly_at_epoch(const KeplerianElements& initial,
                                                time::EpochTDB target_time) {
    return propagate(initial, target_time).mean_anomaly;
}

} // namespace astdyn::propagation
