/**
 * @file Propagator.cpp
 * @brief Implementation of orbital propagation
 */

#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/core/Constants.hpp"
#include "src/utils/time_types.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "src/core/frame_tags.hpp"
#include "src/types/timed_state.hpp"
#include "src/propagation/kepler_propagator.hpp"
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

Eigen::VectorXd Propagator::compute_derivatives(utils::Instant t, const Eigen::VectorXd& state) {
    // State is in meters and meters/second
    Eigen::Vector3d position = state.head<3>();
    Eigen::Vector3d velocity = state.tail<3>();
    
    // Get Sun Barycentric Position
    auto sun_pos_bary = ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(t);
    Eigen::Vector3d sun_pos_bary_eigen = sun_pos_bary.to_eigen(); // in m
    
    // Compute total acceleration (in m/s^2)
    Eigen::Vector3d acceleration = two_body_acceleration(position);
    
    if (settings_.include_planets) {
        acceleration += planetary_perturbations(position, t, sun_pos_bary_eigen);
    }
    
    if (settings_.include_relativity) {
        acceleration += relativistic_correction(position, velocity);
    }
    
    if (settings_.include_asteroids && asteroids_) {
        // AsteroidPerturbations::computePerturbation might expect AU and return AU/day^2?
        // Let's assume for now it needs conversion if it hasn't been updated.
        // Actually, let's just make it work with what we have.
        // TODO: Update AsteroidPerturbations to use strong types too.
        double m_to_au = 1.0 / (constants::AU * 1000.0);
        double d_to_s = 86400.0;
        Eigen::Vector3d pos_au = position * m_to_au;
        Eigen::Vector3d sun_au = sun_pos_bary_eigen * m_to_au;
        Eigen::Vector3d acc_au_d2 = asteroids_->computePerturbation(pos_au, t.mjd.value, sun_au);
        acceleration += acc_au_d2 * (constants::AU * 1000.0 / (d_to_s * d_to_s));
    }
    
    // Yarkovsky Effect
    if (settings_.include_yarkovsky && std::abs(settings_.yarkovsky_a2) > 0.0) {
        double r_m = position.norm();
        double v_m = velocity.norm();
        if (r_m > 100.0 && v_m > 1e-6) {
            double r_au = r_m / (constants::AU * 1000.0);
            // a2 is in AU/d^2. Convert to m/s^2.
            double acc_val_m_s2 = (settings_.yarkovsky_a2 / (r_au * r_au)) * (constants::AU * 1000.0 / (86400.0 * 86400.0));
            acceleration += acc_val_m_s2 * (velocity / v_m);
        }
    }
    
    Eigen::VectorXd derivative(6);
    derivative.head<3>() = velocity;
    derivative.tail<3>() = acceleration;
    
    return derivative;
}

Eigen::Vector3d Propagator::two_body_acceleration(const Eigen::Vector3d& position) const {
    double r = position.norm();
    double r3 = r * r * r;
    // central_body_gm is in m^3/s^2 if we are working in Meters.
    // settings_.central_body_gm defaults to GMS (AU^3/d^2).
    double gm_m3_s2 = settings_.central_body_gm * (std::pow(constants::AU * 1000.0, 3) / std::pow(86400.0, 2));
    return -gm_m3_s2 * position / r3;
}

Eigen::Vector3d Propagator::planetary_perturbations(const Eigen::Vector3d& position,
                                                   utils::Instant t,
                                                   const Eigen::Vector3d& sun_pos_bary) {
    Eigen::Vector3d perturbation = Eigen::Vector3d::Zero();
    
    auto add_planet_perturbation = [&]( CelestialBody planet, double planet_gm_m3s2) {
        auto planet_pos_bary = ephemeris::PlanetaryEphemeris::getPosition(planet, t);
        Eigen::Vector3d planet_pos_helio = planet_pos_bary.to_eigen() - sun_pos_bary;
        Eigen::Vector3d delta = planet_pos_helio - position;
        
        const double d_norm = delta.norm();
        const double d3 = d_norm * d_norm * d_norm;
        const double p_dist = planet_pos_helio.norm();
        const double p_dist3 = p_dist * p_dist * p_dist;
        
        perturbation += planet_gm_m3s2 * (delta / d3 - planet_pos_helio / p_dist3);
    };
    
    // Convert GM to m^3/s^2
    auto to_si = [](double gm_au) { return gm_au * (std::pow(constants::AU * 1000.0, 3) / std::pow(86400.0, 2)); };

    if (settings_.perturb_mercury) add_planet_perturbation(CelestialBody::MERCURY, to_si(constants::GM_MERCURY_AU));
    if (settings_.perturb_venus)   add_planet_perturbation(CelestialBody::VENUS,   to_si(constants::GM_VENUS_AU));
    if (settings_.perturb_earth)   add_planet_perturbation(CelestialBody::EARTH,   to_si(constants::GM_EARTH_AU));
    if (settings_.perturb_mars)    add_planet_perturbation(CelestialBody::MARS,    to_si(constants::GM_MARS_AU));
    if (settings_.perturb_jupiter) add_planet_perturbation(CelestialBody::JUPITER, to_si(constants::GM_JUPITER_AU));
    if (settings_.perturb_saturn)  add_planet_perturbation(CelestialBody::SATURN,  to_si(constants::GM_SATURN_AU));
    if (settings_.perturb_uranus)  add_planet_perturbation(CelestialBody::URANUS,  to_si(constants::GM_URANUS_AU));
    if (settings_.perturb_neptune) add_planet_perturbation(CelestialBody::NEPTUNE, to_si(constants::GM_NEPTUNE_AU));
    
    if (settings_.include_moon) {
        // constants::GM_MOON is in km^3/s^2
        add_planet_perturbation(CelestialBody::MOON, constants::GM_MOON * 1e9);
    }
    
    return perturbation;
}

Eigen::Vector3d Propagator::relativistic_correction(const Eigen::Vector3d& position,
                                                   const Eigen::Vector3d& velocity) const {
    double r = position.norm();
    double v2 = velocity.squaredNorm();
    double c_ms = constants::C_LIGHT * 1000.0;
    double c2 = c_ms * c_ms;
    double mu = settings_.central_body_gm * (std::pow(constants::AU * 1000.0, 3) / std::pow(86400.0, 2));
    
    double beta = settings_.ppn_beta;
    double gamma = settings_.ppn_gamma;
    
    double term1_coeff = 2.0 * (beta + gamma) * mu / r - gamma * v2;
    double rv_dot = position.dot(velocity);
    
    Eigen::Vector3d term1 = term1_coeff * position;
    Eigen::Vector3d term2 = 2.0 * (1.0 + gamma) * rv_dot * velocity;
    
    return (mu / (r * r * r * c2)) * (term1 + term2);
}

CartesianElements Propagator::propagate_cartesian(const CartesianElements& initial,
                                                  utils::Instant target_time) {
    // initial.position and initial.velocity are types::Vector3 (m, m/s)
    Eigen::VectorXd y0(6);
    y0.head<3>() = initial.position.to_eigen();
    y0.tail<3>() = initial.velocity.to_eigen();
    
    // Integrate in absolute MJD
    double t0 = initial.epoch.mjd.value;
    double tf = target_time.mjd.value;
    
    DerivativeFunction f = [this](double t_val, const Eigen::VectorXd& y) {
        // Convert t_val (MJD) back to Instant
        auto t_inst = utils::Instant::from_tt(utils::ModifiedJulianDate(t_val));
        return compute_derivatives(t_inst, y);
    };
    
    Eigen::VectorXd yf = integrator_->integrate(f, y0, t0, tf);
    
    CartesianElements final;
    final.epoch = target_time;
    final.gravitational_parameter = initial.gravitational_parameter;
    final.position = types::Vector3<core::GCRF, core::Meter>(yf.head<3>());
    final.velocity = types::Vector3<core::GCRF, core::Meter>(yf.tail<3>());
    final.covariance = initial.covariance;
    
    return final;
}

KeplerianElements Propagator::propagate_keplerian(const KeplerianElements& initial,
                                                  utils::Instant target_time) {
    // Bridge to Cartesian
    CartesianElements cart = keplerian_to_cartesian(initial);
    cart.covariance = initial.covariance;
    
    CartesianElements cart_final = propagate_cartesian(cart, target_time);
    
    KeplerianElements final = cartesian_to_keplerian(cart_final);
    final.covariance = initial.covariance;
    
    return final;
}

std::vector<CartesianElements> Propagator::propagate_ephemeris(
    const CartesianElements& initial,
    const std::vector<utils::Instant>& target_times) {
    
    std::vector<CartesianElements> results;
    results.reserve(target_times.size());
    for (const auto& t : target_times) {
        results.push_back(propagate_cartesian(initial, t));
    }
    return results;
}

// Analytical TwoBodyPropagator
KeplerianElements TwoBodyPropagator::propagate(const KeplerianElements& initial,
                                               utils::Instant target_time) {
    using namespace astdyn::propagation;
    using ::astdyn::core::GCRF;
    using ::astdyn::types::TimedState;
    using ::astdyn::utils::Instant;

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
                                                utils::Instant target_time) {
    return propagate(initial, target_time).mean_anomaly;
}

} // namespace astdyn::propagation
