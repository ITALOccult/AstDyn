/**
 * @file Propagator.cpp
 * @brief Implementation of orbital propagation
 */

#include "astdyn/propagation/Propagator.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::propagation {

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
        
        // Load SPK if provided
        if (!settings_.asteroid_ephemeris_file.empty()) {
            asteroids_->loadSPK(settings_.asteroid_ephemeris_file);
        }
    }
}

Eigen::VectorXd Propagator::compute_derivatives(double t, const Eigen::VectorXd& state) {
    // State vector: [x, y, z, vx, vy, vz]
    Eigen::Vector3d position = state.head<3>();
    Eigen::Vector3d velocity = state.tail<3>();
    
    // Fetch Sun Barycentric Position once
    double jd_tdb = t + 2400000.5;
    Eigen::Vector3d sun_pos_bary = ephemeris::PlanetaryEphemeris::getPosition(
        ephemeris::CelestialBody::SUN, jd_tdb);
    
    // Compute total acceleration
    Eigen::Vector3d acceleration = two_body_acceleration(position);
    
    if (settings_.include_planets) {
        acceleration += planetary_perturbations(position, t, sun_pos_bary);
    }
    
    if (settings_.include_relativity) {
        acceleration += relativistic_correction(position, velocity);
    }
    
    if (settings_.include_asteroids && asteroids_) {
        acceleration += asteroid_perturbations(position, t, sun_pos_bary);
    }
    
    // Yarkovsky Effect (Non-Gravitational)
    if (settings_.include_yarkovsky && std::abs(settings_.yarkovsky_a2) > 0.0) {
        double r = position.norm();
        double v_norm = velocity.norm();
        if (r > 1e-6 && v_norm > 1e-9) { // Avoid division by zero
            // Tangential acceleration: a = (A2 / r^2) * (v / |v|)
            // A2 is AU/d^2 at 1 AU
            double r2 = r * r;
            acceleration += (settings_.yarkovsky_a2 / r2) * (velocity / v_norm);
        }
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
                                                   double mjd_tdb,
                                                   const Eigen::Vector3d& sun_pos_bary) {
    Eigen::Vector3d perturbation = Eigen::Vector3d::Zero();
    
    // Convert MJD to JD for ephemeris calls
    double jd_tdb = mjd_tdb + 2400000.5;
    
    // Helper lambda for computing perturbation from one planet
    auto add_planet_perturbation = [&](CelestialBody planet, double planet_gm) {
        Eigen::Vector3d planet_pos_bary = ephemeris::PlanetaryEphemeris::getPosition(planet, jd_tdb);
        
        // Convert to Heliocentric (Vector Sun->Planet)
        Eigen::Vector3d planet_pos_helio = planet_pos_bary - sun_pos_bary;
        
        Eigen::Vector3d delta = planet_pos_helio - position;
        double delta_norm = delta.norm();
        double delta3 = delta_norm * delta_norm * delta_norm;
        
        double planet_dist = planet_pos_helio.norm();
        double planet_dist3 = planet_dist * planet_dist * planet_dist;
        
        // Indirect term: -GM_planet * r_planet / |r_planet|³
        // Direct term: GM_planet * (r_planet - r) / |r_planet - r|³
        perturbation += planet_gm * (delta / delta3 - planet_pos_helio / planet_dist3);
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
    // double rdot = position.dot(velocity) / r; // Unused in PPN formulation (we use rv_dot)
    
    double c = constants::SPEED_OF_LIGHT_AU_PER_DAY;
    double c2 = c * c;
    double mu = settings_.central_body_gm;
    
    // PPN parameter formulation (General Relativity when β=1, γ=1)
    // a_GR = (mu/c²r³) * [ (2(β+γ)mu/r - γv²)r + 2(1+γ)(r·v)v ]
    
    double beta = settings_.ppn_beta;
    double gamma = settings_.ppn_gamma;
    
    // Term 1 coeff: 2(β+γ)μ/r - γv²
    double term1_coeff = 2.0 * (beta + gamma) * mu / r - gamma * v2;
    
    // Term 2 coeff: 2(1+γ)(r·v)
    // rdot = (r·v)/r, so (r·v) = rdot * r
    // But formula usually written as term along v uses (r·v).
    // Let's stick to using rdot * r for (r·v) if needed, or just rdot logic.
    // My previous code used rdot * velocity. 
    // Is previous code `4 * rdot * velocity` equivalent to `2(1+gamma)(r·v)v`?
    // If gamma=1, 2(2)(r·v)v = 4(r·v)v.
    // But acceleration has 1/r^3 factor.
    // (r·v) = r * rdot.
    // So 4 * r * rdot * v.
    // Divided by r^3 -> 4 * rdot * v / r^2.
    // Wait.
    // Previous code: return (mu / (r^3 * c2)) * ( ... + 4 * rdot * velocity )
    // Wait, rdot is dim L/T. Velocity L/T. 
    // 4 * rdot * velocity is dim L^2/T^2.
    // (4 mu/r - v^2) is L^2/T^2.
    // Matches.
    // So 4 * rdot * velocity is OK for inside bracket.
    // Let's check Term 2 in general form: 2(1+γ)(r·v)v
    // (r·v) = r * rdot.
    // So 2(1+γ) * r * rdot * v.
    // This has an extra 'r'.
    // The equation in Moyer/IERS is:
    // a = ... + 2(1+gamma) * (r . v) * v
    // But notice the common factor outside bracket: mu / (c^2 * r^3).
    // So inside term2 must have dimension L^2/T^2 * Length = L^3/T^2 ?? No.
    // Acceleration: L/T^2.
    // Outside: L^3/T^2 / (L^2/T^2 * L^3) = 1/L^2.
    // So inside must be L^3/T^2.
    // But (4mu/r - v^2) * r is (L^2/T^2) * L = L^3/T^2. Correct.
    // (r.v) * v is (L*L/T) * L/T = L^3/T^2. Correct.
    //
    // Previous code: `4.0 * rdot * velocity`.
    // rdot = r.v / r.
    // So rdot * v = (r.v/r) * v.
    // This is MISSING an 'r' factor compared to (r.v)*v ?
    // Or previous code assumed 1/r^2 outside?
    // Previous code: `return (mu / (r * r * r * c2)) * (term1 + term2);`
    // It has 1/r^3.
    // So inside must be L^3/T^2.
    // `4 * rdot * velocity` = 4 * (L/T) * (L/T) = L^2/T^2.
    // Dimension Mismatch in previous code??
    // Wait. `term1 = (4mu/r - v2) * position`. (L^2/T^2) * L = L^3/T^2. Correct.
    // `term2 = 4 * rdot * velocity`. L^2/T^2. INVALID.
    // It should have been `4 * rdot * r * velocity` or `4 * (r.v) * velocity`.
    // rdot = (r.v)/r.
    // So `4 * (r.v) * velocity` is correct.
    // But `4 * rdot * velocity` is NOT correct unless I multiply by r.
    // So previous code was BUGGY for term 2?
    // Let's check Moyer 2003 Eq 8-1.
    // Acceleration includes `4 * (r_dot_v) * v` inside.
    // Yes. `r_dot_v` is scalar product.
    // My previous code calculated `rdot = position.dot(velocity) / r;`.
    // So `rdot` is radial velocity.
    // `(r.v)` is `rdot * r`.
    // So I missed a factor of `r` in the second term.
    //
    // FIX: Using PPN formula correctly naturally fixes this.
    // Term 2 = 2(1+gamma) * (position.dot(velocity)) * velocity.
    
    double rv_dot = position.dot(velocity);
    
    Eigen::Vector3d term1 = term1_coeff * position;
    Eigen::Vector3d term2 = 2.0 * (1.0 + gamma) * rv_dot * velocity;
    
    return (mu / (r * r * r * c2)) * (term1 + term2);
}

Eigen::Vector3d Propagator::asteroid_perturbations(const Eigen::Vector3d& position,
                                                   double mjd_tdb,
                                                   const Eigen::Vector3d& sun_pos_bary) {
    if (!asteroids_) {
        return Eigen::Vector3d::Zero();
    }
    
    // Compute perturbation from all enabled asteroids
    // AsteroidPerturbations already handles direct + indirect terms
    return asteroids_->computePerturbation(position, mjd_tdb, sun_pos_bary);
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

} // namespace astdyn::propagation
