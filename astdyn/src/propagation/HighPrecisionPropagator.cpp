/**
 * @file HighPrecisionPropagator.cpp
 * @brief Implementation of HighPrecisionPropagator
 */

#include "astdyn/propagation/HighPrecisionPropagator.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "src/core/units.hpp"
#include "src/utils/time_types.hpp"
#include "src/types/orbital_state.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include <iostream>
#include <cmath>

namespace astdyn::propagation {

using namespace astdyn::ephemeris;
using namespace astdyn::constants;

HighPrecisionPropagator::HighPrecisionPropagator(const Config& config) 
    : config_(config) 
{
    planetary_ephemeris_ = std::make_shared<PlanetaryEphemeris>();
    if (!config_.de441_path.empty()) {
        try {
            custom_provider_ = std::make_shared<DE441Provider>(config_.de441_path);
            PlanetaryEphemeris::setProvider(custom_provider_);
        } catch (const std::exception& e) {
            std::cerr << "[HighPrecisionPropagator] Warning: Failed to load DE441 (" << e.what() << ").\n";
        }
    }
}

HighPrecisionPropagator::~HighPrecisionPropagator() {}

std::unique_ptr<Propagator> HighPrecisionPropagator::createPropagator() const {
    auto integrator = std::make_unique<RKF78Integrator>(config_.step_size, config_.tolerance);
    PropagatorSettings settings;
    settings.include_planets = config_.perturbations_planets;
    settings.include_relativity = config_.relativity;
    settings.include_asteroids = config_.perturbations_asteroids;
    settings.asteroid_ephemeris_file = config_.asteroid_ephemeris_file;
    
    if (settings.include_planets) {
        settings.perturb_jupiter = true;
        settings.perturb_saturn = true;
        settings.perturb_mars = true;
        settings.perturb_earth = true;
        settings.perturb_venus = true; 
        settings.perturb_mercury = true;
        settings.perturb_uranus = true;
        settings.perturb_neptune = true;
    }

    return std::make_unique<Propagator>(std::move(integrator), planetary_ephemeris_, settings);
}

HighPrecisionPropagator::ObservationResult 
HighPrecisionPropagator::calculateGeocentricObservation(
    const KeplerianElements& initial_elements, 
    utils::Instant target_time,
    InputFrame frame
) {
    if (!cached_propagator_) cached_propagator_ = createPropagator();
    
    CartesianElements cart_icrf = propagate_cartesian(initial_elements, target_time, frame);
    
    // Earth position at target time
    auto r_sun_ssb = PlanetaryEphemeris::getSunBarycentricPosition(target_time);
    auto r_earth_ssb = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, target_time);
    auto r_earth_helio = r_earth_ssb - r_sun_ssb;

    double lt_sec = 0.0;
    CartesianElements state_at_emit = cart_icrf;
    
    // Light-time correction loop
    for (int i = 0; i < 3; ++i) {
        auto t_emit = utils::Instant::from_tt(utils::ModifiedJulianDate(target_time.mjd.value - lt_sec / 86400.0));
        state_at_emit = cached_propagator_->propagate_cartesian(cart_icrf, t_emit);
        
        double dist = (state_at_emit.position - r_earth_helio).norm();
        lt_sec = dist / (constants::C_LIGHT * 1000.0);
    }

    auto r_geo = state_at_emit.position - r_earth_helio;
    double dist = r_geo.norm();
    
    double ra = atan2(r_geo.y, r_geo.x);
    double dec = asin(r_geo.z / dist);

    if (ra < 0) ra += astdyn::constants::TWO_PI;

    ObservationResult result;
    result.ra_deg = ra * astdyn::constants::RAD_TO_DEG;
    result.dec_deg = dec * astdyn::constants::RAD_TO_DEG;
    result.distance_au = dist / (constants::AU * 1000.0);
    result.light_time_sec = lt_sec;
    result.geocentric_position = r_geo;

    return result;
}

std::shared_ptr<EphemerisProvider> HighPrecisionPropagator::getEphemerisProvider() const {
    return custom_provider_;
}

void HighPrecisionPropagator::setPlanetaryEphemeris(std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephemeris) {
    planetary_ephemeris_ = ephemeris;
    cached_propagator_.reset();
}

CartesianElements HighPrecisionPropagator::propagate_cartesian(
    const KeplerianElements& initial_elements,
    utils::Instant target_time,
    InputFrame frame
) {
    if (!cached_propagator_) cached_propagator_ = createPropagator();
    
    CartesianElements cart_start = keplerian_to_cartesian(initial_elements);
    
    if (frame == InputFrame::ECLIPTIC) {
        // VSOP87 elements are usually in ecliptic frame. 
        // Need to ensure they are rotated to GCRF (Equatorial J2000) if requested.
        // For now assume keplerian_to_cartesian handles the requested frame or we convert here.
        // Actually, reference_frame.hpp transform_position handles this.
        auto pos_ecl = cart_start.position.to_eigen();
        auto vel_ecl = cart_start.velocity.to_eigen();
        auto pos_eq = astdyn::coordinates::ReferenceFrame::ecliptic_to_j2000() * pos_ecl;
        auto vel_eq = astdyn::coordinates::ReferenceFrame::ecliptic_to_j2000() * vel_ecl;
        cart_start.position = types::Vector3<core::GCRF, core::Meter>(pos_eq);
        cart_start.velocity = types::Vector3<core::GCRF, core::Meter>(vel_eq);
    }

    return cached_propagator_->propagate_cartesian(cart_start, target_time);
}

} // namespace astdyn::propagation
