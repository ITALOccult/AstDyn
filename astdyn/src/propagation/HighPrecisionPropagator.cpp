/**
 * @file HighPrecisionPropagator.cpp
 * @brief Implementation of HighPrecisionPropagator
 */

#include "astdyn/propagation/HighPrecisionPropagator.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
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
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements, 
    time::EpochTDB target_time
) {
    if (!cached_propagator_) cached_propagator_ = createPropagator();
    
    auto cart_icrf = propagate_cartesian(initial_elements, target_time);
    return calculateGeocentricObservation(cart_icrf, target_time);
}

HighPrecisionPropagator::ObservationResult 
HighPrecisionPropagator::calculateGeocentricObservation(
    const physics::CartesianStateTyped<core::GCRF>& cart_icrf, 
    time::EpochTDB target_time
) {
    if (!cached_propagator_) cached_propagator_ = createPropagator();
    
    // Earth position at target time
    auto target_time_inst = utils::Instant::from_tt(utils::ModifiedJulianDate(target_time.mjd()));
    auto r_sun_ssb = PlanetaryEphemeris::getSunBarycentricPosition(target_time_inst);
    auto r_earth_ssb = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, target_time_inst);
    auto r_earth_helio = r_earth_ssb - r_sun_ssb;

    double lt_sec = 0.0;
    physics::CartesianStateTyped<core::GCRF> state_at_emit = cart_icrf;
    
    // Light-time correction loop
    for (int i = 0; i < 3; ++i) {
        auto t_emit = time::EpochTDB::from_mjd(target_time.mjd() - lt_sec / 86400.0);
        state_at_emit = cached_propagator_->propagate_cartesian(cart_icrf, t_emit);
        
        double dist = (state_at_emit.position.to_eigen_si() - r_earth_helio.to_eigen()).norm();
        lt_sec = dist / (constants::C_LIGHT * 1000.0);
    }

    auto r_geo = state_at_emit.position.to_eigen_si() - r_earth_helio.to_eigen();
    double dist = r_geo.norm();
    
    double ra = atan2(r_geo.y(), r_geo.x());
    double dec = asin(r_geo.z() / dist);

    if (ra < 0) ra += astdyn::constants::TWO_PI;

    ObservationResult result;
    result.ra_deg = ra * astdyn::constants::RAD_TO_DEG;
    result.dec_deg = dec * astdyn::constants::RAD_TO_DEG;
    result.distance_au = dist / (constants::AU * 1000.0);
    result.light_time_sec = lt_sec;
    result.geocentric_position = types::Vector3<core::GCRF, core::Meter>(r_geo);

    return result;
}

std::shared_ptr<EphemerisProvider> HighPrecisionPropagator::getEphemerisProvider() const {
    return custom_provider_;
}

void HighPrecisionPropagator::setPlanetaryEphemeris(std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephemeris) {
    planetary_ephemeris_ = ephemeris;
    cached_propagator_.reset();
}

physics::CartesianStateTyped<core::GCRF> HighPrecisionPropagator::propagate_cartesian(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
    time::EpochTDB target_time
) {
    if (!cached_propagator_) cached_propagator_ = createPropagator();
    
    // Bridge to un-typed format for conversion math
    KeplerianElements kep_old;
    kep_old.epoch = utils::Instant::from_tt(utils::ModifiedJulianDate(initial_elements.epoch.mjd()));
    kep_old.semi_major_axis = initial_elements.a.to_au();
    kep_old.eccentricity = initial_elements.e;
    kep_old.inclination = initial_elements.i.to_rad();
    kep_old.longitude_ascending_node = initial_elements.node.to_rad();
    kep_old.argument_perihelion = initial_elements.omega.to_rad();
    kep_old.mean_anomaly = initial_elements.M.to_rad();
    kep_old.gravitational_parameter = initial_elements.gm.to_au3_d2();
    
    // Convert to un-typed Cartesian (which returns SI internally!)
    CartesianElements cart_old = keplerian_to_cartesian(kep_old);
    
    // Rotate to ICRF
    auto pos_ecl = cart_old.position.to_eigen();
    auto vel_ecl = cart_old.velocity.to_eigen();
    auto pos_eq = astdyn::coordinates::ReferenceFrame::ecliptic_to_j2000() * pos_ecl;
    auto vel_eq = astdyn::coordinates::ReferenceFrame::ecliptic_to_j2000() * vel_ecl;
    
    // Bring into AstDyn 3.0 Type Safety
    auto cart_typed = physics::CartesianStateTyped<core::GCRF>::from_si(
        initial_elements.epoch,
        pos_eq.x(), pos_eq.y(), pos_eq.z(),
        vel_eq.x(), vel_eq.y(), vel_eq.z(),
        constants::GM_SUN * 1e9
    );

    return cached_propagator_->propagate_cartesian(cart_typed, target_time);
}

physics::CartesianStateTyped<core::GCRF> HighPrecisionPropagator::propagate_cartesian(
    const physics::CartesianStateTyped<core::GCRF>& initial_elements,
    time::EpochTDB target_time
) {
    if (!cached_propagator_) cached_propagator_ = createPropagator();
    return cached_propagator_->propagate_cartesian(initial_elements, target_time);
}

} // namespace astdyn::propagation
