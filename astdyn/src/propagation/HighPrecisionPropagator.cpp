/**
 * @file HighPrecisionPropagator.cpp
 * @brief Implementation of HighPrecisionPropagator
 */

#include "astdyn/propagation/HighPrecisionPropagator.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>
#include <cmath>

namespace astdyn::propagation {

using namespace astdyn::ephemeris;
using namespace astdyn::constants;

HighPrecisionPropagator::HighPrecisionPropagator(const Config& config) 
    : config_(config) 
{
    // Initialize PlanetaryEphemeris wrapper
    planetary_ephemeris_ = std::make_shared<PlanetaryEphemeris>();

    // Load custom provider (DE441) if requested
    if (!config_.de441_path.empty()) {
        try {
            custom_provider_ = std::make_shared<DE441Provider>(config_.de441_path);
            PlanetaryEphemeris::setProvider(custom_provider_);
            // std::cout << "[HighPrecisionPropagator] Native DE441 loaded from: " << config_.de441_path << "\n";
        } catch (const std::exception& e) {
            std::cerr << "[HighPrecisionPropagator] Warning: Failed to load DE441 (" << e.what() 
                      << "). Using analytical fallback.\n";
        }
    }
}

HighPrecisionPropagator::~HighPrecisionPropagator() {
    // Reset global provider to avoid side effects if singleton employed (?)
    // PlanetaryEphemeris::setProvider(nullptr); 
    // Careful: PlanetaryEphemeris::setProvider is static global. 
    // If this class is used in multi-threaded env or multiple instances, this is tricky.
    // For now assuming single context usage or global config.
}

std::unique_ptr<Propagator> HighPrecisionPropagator::createPropagator() const {
    auto integrator = std::make_unique<RKF78Integrator>(config_.step_size, config_.tolerance);
    
    PropagatorSettings settings;
    settings.include_planets = config_.perturbations_planets;
    settings.include_relativity = config_.relativity;
    settings.include_asteroids = config_.perturbations_asteroids;
    settings.asteroid_ephemeris_file = config_.asteroid_ephemeris_file;
    // Enable all major perturbers if planets enabled
    if (settings.include_planets) {
        settings.perturb_jupiter = true;
        settings.perturb_saturn = true;
        settings.perturb_mars = true;
        settings.perturb_earth = true;
        settings.perturb_venus = true; 
        settings.perturb_mercury = false; // Usually negligible for main belt?
        settings.perturb_uranus = true;
        settings.perturb_neptune = true;
    }

    return std::make_unique<Propagator>(std::move(integrator), planetary_ephemeris_, settings);
}

// Helper: Ecliptic to Equatorial (J2000)
static Eigen::Vector3d ecl_to_eq(const Eigen::Vector3d& ecl) {
    double eps = 23.4392911 * M_PI / 180.0;
    double c = cos(eps);
    double s = sin(eps);
    return Eigen::Vector3d(ecl.x(), c*ecl.y() - s*ecl.z(), s*ecl.y() + c*ecl.z());
}

HighPrecisionPropagator::ObservationResult 
HighPrecisionPropagator::calculateGeocentricObservation(
    const KeplerianElements& initial_elements, 
    double target_jd_tdb
) {
    auto propagator = createPropagator();

    // 1. Initial State Setup
    // Convert input Keplerian (assumed Ecliptic J2000) to Cartesian ICRF (Equatorial)
    CartesianElements cart_ecl = keplerian_to_cartesian(initial_elements);
    
    CartesianElements cart_icrf;
    cart_icrf.epoch_mjd_tdb = initial_elements.epoch_mjd_tdb;
    cart_icrf.gravitational_parameter = initial_elements.gravitational_parameter;
    cart_icrf.position = ecl_to_eq(cart_ecl.position);
    cart_icrf.velocity = ecl_to_eq(cart_ecl.velocity);

    double target_mjd = target_jd_tdb - 2400000.5;

    // 2. Propagate to Target Approximation (Heliocentric)
    // Propagate asteroid center to target time
    CartesianElements state_ast = propagator->propagate_cartesian(cart_icrf, target_mjd);

    // 3. Earth Position at Target Time
    // PlanetaryEphemeris returns Ecliptic J2000 (even if DE441 wrapped)
    Eigen::Vector3d r_earth_ecl = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, target_jd_tdb);
    Eigen::Vector3d r_earth_eq = ecl_to_eq(r_earth_ecl);

    // 4. Light Time Correction Loop
    double t_emit_mjd = target_mjd;
    double lt_days = 0.0;
    Eigen::Vector3d r_ast_emit = state_ast.position; // Initial guess

    for (int i = 0; i < 3; ++i) {
        // Propagate asteroid to emission time (small step from target or full propagation)
        // Since Propagator might be stateful at the end, let's just re-call propagate
        // Efficiency note: Ideally we step back from state_ast, but Propagator interface is absolute.
        // If the propagator step allows, this will be fast. 
        // For MAX precision, we propagate from initial epoch (or closest checkpoint) to t_emit.
        // Assuming Propagator is robust:
        CartesianElements state_emit = propagator->propagate_cartesian(cart_icrf, t_emit_mjd);
        r_ast_emit = state_emit.position;

        // Geometric distance Earth(t_obs) to Asteroid(t_emit)
        double dist_au = (r_ast_emit - r_earth_eq).norm();
        lt_days = dist_au / SPEED_OF_LIGHT_AU_PER_DAY;
        
        // Update emission time
        t_emit_mjd = target_mjd - lt_days;
    }

    // 5. Final Geocentric Vector
    Eigen::Vector3d r_geo = r_ast_emit - r_earth_eq;
    double dist = r_geo.norm();
    
    // Convert to RA/DEC
    double ra = atan2(r_geo.y(), r_geo.x());
    double dec = asin(r_geo.z() / dist);

    // Normalize RA 0-360
    if (ra < 0) ra += 2.0 * M_PI;

    ObservationResult result;
    result.ra_deg = ra * 180.0 / M_PI;
    result.dec_deg = dec * 180.0 / M_PI;
    result.distance_au = dist;
    result.light_time_sec = lt_days * 86400.0;
    result.geocentric_position = r_geo;

    return result;
}

std::shared_ptr<EphemerisProvider> HighPrecisionPropagator::getEphemerisProvider() const {
    return custom_provider_;
}

} // namespace astdyn::propagation
