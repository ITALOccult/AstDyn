/**
 * @file Propagator.cpp
 * @brief Implementation of orbital propagation
 */

#include "orbfit/propagation/Propagator.hpp"
#include <cmath>

namespace orbfit::propagation {

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
    
    // Initialize asteroid perturbations if enabled
    if (settings_.include_asteroids) {
        asteroids_ = std::make_shared<ephemeris::AsteroidPerturbations>();
        asteroids_->loadDefaultAsteroids();
    }
}

Eigen::VectorXd Propagator::compute_derivatives(double t, const Eigen::VectorXd& state) {
    // State vector: [x, y, z, vx, vy, vz]
    Eigen::Vector3d position = state.head<3>();
    Eigen::Vector3d velocity = state.tail<3>();
    
    // Compute total acceleration
    Eigen::Vector3d acceleration = two_body_acceleration(position);
    
    if (settings_.include_planets) {
        acceleration += planetary_perturbations(position, t);
    }
    
    if (settings_.include_relativity) {
        acceleration += relativistic_correction(position, velocity);
    }
    
    if (settings_.include_asteroids && asteroids_) {
        acceleration += asteroid_perturbations(position, t);
    }
    
    // Derivative: [velocity, acceleration]
    Eigen::VectorXd derivative(6);
    derivative.head<3>() = velocity;
    derivative.tail<3>() = acceleration;
    
    return derivative;
}

Eigen::Vector3d Propagator::two_body_acceleration(const Eigen::Vector3d& position) const {
    double r = position.norm();
    double r3 = r * r * r;
    return -settings_.central_body_gm * position / r3;
}

Eigen::Vector3d Propagator::planetary_perturbations(const Eigen::Vector3d& position,
                                                   double mjd_tdb) {
    Eigen::Vector3d perturbation = Eigen::Vector3d::Zero();
    
    // Helper lambda for computing perturbation from one planet
    auto add_planet_perturbation = [&](CelestialBody planet, double planet_gm) {
        Eigen::Vector3d planet_pos = ephemeris::PlanetaryEphemeris::getPosition(planet, mjd_tdb);
        Eigen::Vector3d delta = planet_pos - position;
        double delta_norm = delta.norm();
        double delta3 = delta_norm * delta_norm * delta_norm;
        
        double planet_dist = planet_pos.norm();
        double planet_dist3 = planet_dist * planet_dist * planet_dist;
        
        // Indirect term: -GM_planet * r_planet / |r_planet|³
        // Direct term: GM_planet * (r_planet - r) / |r_planet - r|³
        perturbation += planet_gm * (delta / delta3 - planet_pos / planet_dist3);
    };
    
    // Add perturbations from each enabled planet
    // Convert GM from km³/s² to AU³/day²
    double conv = constants::GM_KM3S2_TO_AU3DAY2;
    
    if (settings_.perturb_mercury) {
        add_planet_perturbation(CelestialBody::MERCURY, constants::GM_MERCURY * conv);
    }
    if (settings_.perturb_venus) {
        add_planet_perturbation(CelestialBody::VENUS, constants::GM_VENUS * conv);
    }
    if (settings_.perturb_earth) {
        add_planet_perturbation(CelestialBody::EARTH, constants::GM_EARTH * conv);
    }
    if (settings_.perturb_mars) {
        add_planet_perturbation(CelestialBody::MARS, constants::GM_MARS * conv);
    }
    if (settings_.perturb_jupiter) {
        add_planet_perturbation(CelestialBody::JUPITER, constants::GM_JUPITER * conv);
    }
    if (settings_.perturb_saturn) {
        add_planet_perturbation(CelestialBody::SATURN, constants::GM_SATURN * conv);
    }
    if (settings_.perturb_uranus) {
        add_planet_perturbation(CelestialBody::URANUS, constants::GM_URANUS * conv);
    }
    if (settings_.perturb_neptune) {
        add_planet_perturbation(CelestialBody::NEPTUNE, constants::GM_NEPTUNE * conv);
    }
    
    if (settings_.include_moon) {
        add_planet_perturbation(CelestialBody::MOON, constants::GM_MOON * conv);
    }
    
    return perturbation;
}

Eigen::Vector3d Propagator::relativistic_correction(const Eigen::Vector3d& position,
                                                   const Eigen::Vector3d& velocity) const {
    // Schwarzschild metric correction (first order in 1/c²)
    // See Moyer (2003) "Formulation for Observed and Computed Values"
    
    double r = position.norm();
    double v2 = velocity.squaredNorm();
    double rdot = position.dot(velocity) / r;
    
    double c = constants::SPEED_OF_LIGHT_AU_PER_DAY;
    double c2 = c * c;
    double mu = settings_.central_body_gm;
    
    // PPN parameter β = γ = 1 (general relativity)
    // a_GR = (mu/r³) * [(4*mu/r - v²) * r + 4*rdot * v] / c²
    
    Eigen::Vector3d term1 = (4.0 * mu / r - v2) * position;
    Eigen::Vector3d term2 = 4.0 * rdot * velocity;
    
    return (mu / (r * r * r * c2)) * (term1 + term2);
}

Eigen::Vector3d Propagator::asteroid_perturbations(const Eigen::Vector3d& position,
                                                   double mjd_tdb) {
    if (!asteroids_) {
        return Eigen::Vector3d::Zero();
    }
    
    // Compute perturbation from all enabled asteroids
    // AsteroidPerturbations already handles direct + indirect terms
    return asteroids_->computePerturbation(position, mjd_tdb);
}

CartesianElements Propagator::propagate_cartesian(const CartesianElements& initial,
                                                  double target_mjd_tdb) {
    // Setup initial state vector [x, y, z, vx, vy, vz]
    Eigen::VectorXd y0(6);
    y0.head<3>() = initial.position;
    y0.tail<3>() = initial.velocity;
    
    // Temporarily update central body GM to match the orbit's GM
    // This ensures consistency between orbit definition and integration
    double original_gm = settings_.central_body_gm;
    settings_.central_body_gm = initial.gravitational_parameter;
    
    // Create derivative function
    DerivativeFunction f = [this](double t, const Eigen::VectorXd& y) {
        return compute_derivatives(t, y);
    };
    
    // Integrate
    Eigen::VectorXd yf = integrator_->integrate(f, y0, 
                                               initial.epoch_mjd_tdb, 
                                               target_mjd_tdb);
    
    // Restore original GM
    settings_.central_body_gm = original_gm;
    
    // Extract final state
    CartesianElements final;
    final.epoch_mjd_tdb = target_mjd_tdb;
    final.gravitational_parameter = initial.gravitational_parameter;
    final.position = yf.head<3>();
    final.velocity = yf.tail<3>();
    
    return final;
}

KeplerianElements Propagator::propagate_keplerian(const KeplerianElements& initial,
                                                  double target_mjd_tdb) {
    // Convert to Cartesian
    CartesianElements cart = keplerian_to_cartesian(initial);
    
    // Propagate
    CartesianElements cart_final = propagate_cartesian(cart, target_mjd_tdb);
    
    // Convert back to Keplerian
    return cartesian_to_keplerian(cart_final);
}

std::vector<CartesianElements> Propagator::propagate_ephemeris(
    const CartesianElements& initial,
    const std::vector<double>& epochs_mjd_tdb) {
    
    std::vector<CartesianElements> results;
    results.reserve(epochs_mjd_tdb.size());
    
    // Propagate to each epoch
    for (double epoch : epochs_mjd_tdb) {
        results.push_back(propagate_cartesian(initial, epoch));
    }
    
    return results;
}

// ============================================================================
// TwoBodyPropagator Implementation (Analytical)
// ============================================================================

KeplerianElements TwoBodyPropagator::propagate(const KeplerianElements& initial,
                                               double target_mjd_tdb) {
    KeplerianElements final = initial;
    final.epoch_mjd_tdb = target_mjd_tdb;
    
    // Only mean anomaly changes in two-body motion
    final.mean_anomaly = mean_anomaly_at_epoch(initial, target_mjd_tdb);
    
    return final;
}

double TwoBodyPropagator::mean_anomaly_at_epoch(const KeplerianElements& initial,
                                                double target_mjd_tdb) {
    double dt = target_mjd_tdb - initial.epoch_mjd_tdb;
    double n = initial.mean_motion();
    double M = initial.mean_anomaly + n * dt;
    
    // Normalize to [0, 2π)
    M = std::fmod(M, constants::TWO_PI);
    if (M < 0.0) M += constants::TWO_PI;
    
    return M;
}

} // namespace orbfit::propagation
