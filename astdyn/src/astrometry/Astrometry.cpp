#include "astdyn/astrometry/Astrometry.hpp"
#include <expected>
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/astrometry/AstrometricTypes.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/io/SPKReader.hpp"
#include "astdyn/ephemeris/AsteroidPerturbations.hpp"
#include <cmath>
#include <mutex>
#include <unordered_map>
#include <iostream>
#include <iomanip>

namespace astdyn::astrometry {

using namespace astdyn::constants;

static std::shared_ptr<ephemeris::DE441Provider> get_cached_provider(const std::string& path) {
    static std::mutex mtx;
    static std::unordered_map<std::string, std::shared_ptr<ephemeris::DE441Provider>> cache;
    if (path.empty()) return nullptr;
    std::lock_guard<std::mutex> lock(mtx);
    if (cache.find(path) == cache.end()) {
        cache[path] = std::make_shared<ephemeris::DE441Provider>(path);
    }
    return cache[path];
}

std::expected<AstrometricObservation, AstrometryError> AstrometryReducer::compute_observation(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
    const time::EpochTDB& t_elements,
    const time::EpochTDB& t_obs,
    const AstDynConfig& e_cfg,
    const AstrometricSettings& a_cfg) {
    
    // 1. Get Earth state at observation time (Cached Native DE441)
    if (e_cfg.ephemeris_file.empty()) return std::unexpected(AstrometryError::EphemerisUnavailable);
    auto de441 = get_cached_provider(e_cfg.ephemeris_file);
    if (!de441) return std::unexpected(AstrometryError::EphemerisUnavailable);
    
    // Sync PlanetaryEphemeris provider too
    ephemeris::PlanetaryEphemeris::setProvider(de441);

    // Get Earth state in GCRF and transform to Ecliptic J2000
    auto earth_p_eq = de441->getPosition(ephemeris::CelestialBody::EARTH, t_obs);
    auto earth_v_eq = de441->getVelocity(ephemeris::CelestialBody::EARTH, t_obs);
    
    auto earth_p_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(earth_p_eq);
    
    // 2. SUN position
    auto sun_p_eq = ephemeris::PlanetaryEphemeris::getPosition(ephemeris::CelestialBody::SUN, t_obs);
    auto sun_p_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(sun_p_eq);
    
    Eigen::Vector3d earth_pos_helio_ecl = earth_p_ecl.to_eigen_si() - sun_p_ecl.to_eigen_si();

    // 3. Light-time Iteration
    auto ast_pos_ecl = compute_light_time_corrected_pos(initial, t_elements, t_obs, earth_pos_helio_ecl, e_cfg);
    
    // DEBUG
    if (e_cfg.verbose) {
        std::cout << "    [DEBUG-ATR] Earth Helio (Ecl): " << (earth_pos_helio_ecl.norm() / (constants::AU*1000.0)) << " AU" << std::endl;
        std::cout << "    [DEBUG-ATR] Ast   Helio (Ecl): " << (ast_pos_ecl.norm() / (constants::AU*1000.0)) << " AU" << std::endl;
    }

    auto raw_rho_ecl = ast_pos_ecl - earth_pos_helio_ecl;

    // 4. Aberration
    // Aberration is usually computed in Equatorial frame using Velocity in GCRF
    // Let's transform rho to Equatorial using typesafe API
    auto rho_ecl_math = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(raw_rho_ecl.x(), raw_rho_ecl.y(), raw_rho_ecl.z());
    auto raw_rho_eq_math = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(rho_ecl_math);
    auto raw_rho_eq = raw_rho_eq_math.to_eigen_si();
    
    auto earth_v_gcrf = earth_v_eq.to_eigen_si();
    
    auto final_rho_eq = a_cfg.stellar_aberration ? apply_stellar_aberration(raw_rho_eq, earth_v_gcrf) : raw_rho_eq;

    // 5. Build Result
    return finalize_observation(final_rho_eq);
}

std::expected<AstrometricObservation, AstrometryError> AstrometryReducer::compute_observation_from_cartesian(
    const physics::CartesianStateTyped<core::GCRF>& initial,
    const time::EpochTDB& t_elements,
    const time::EpochTDB& t_obs,
    const AstDynConfig& e_cfg,
    const AstrometricSettings& a_cfg) {
    
    if (e_cfg.ephemeris_file.empty()) return std::unexpected(AstrometryError::EphemerisUnavailable);
    auto de441 = get_cached_provider(e_cfg.ephemeris_file);
    if (!de441) return std::unexpected(AstrometryError::EphemerisUnavailable);
    ephemeris::PlanetaryEphemeris::setProvider(de441);

    // Transform initial state to Ecliptic properly
    auto pos_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(initial.position);
    auto vel_ecl = coordinates::ReferenceFrame::transform_vel<core::GCRF, core::ECLIPJ2000>(initial.position, initial.velocity);
    
    // Bridge to Engine's expected Keplerian format (AstDynEngine currently fits/propagates in Ecliptic)
    // We'll use the proper typesafe bridge if we can or use the legacy conversion but with correct tags
    propagation::CartesianElements cart_ecl_old;
    cart_ecl_old.epoch = t_elements;
    cart_ecl_old.position = math::Vector3<core::GCRF, physics::Distance>::from_si(pos_ecl.x_si(), pos_ecl.y_si(), pos_ecl.z_si()); // Struct-enforced GCRF tag, but numbers are Ecliptic
    cart_ecl_old.velocity = math::Vector3<core::GCRF, physics::Velocity>::from_si(vel_ecl.x_si(), vel_ecl.y_si(), vel_ecl.z_si());
    cart_ecl_old.gravitational_parameter = initial.gm.to_m3_s2();
    
    auto kep_ecl_legacy = propagation::cartesian_to_keplerian(cart_ecl_old);
    
    auto initial_state_kep_ecl = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        t_elements,
        kep_ecl_legacy.semi_major_axis, kep_ecl_legacy.eccentricity,
        kep_ecl_legacy.inclination * constants::RAD_TO_DEG,
        kep_ecl_legacy.longitude_ascending_node * constants::RAD_TO_DEG,
        kep_ecl_legacy.argument_perihelion * constants::RAD_TO_DEG,
        kep_ecl_legacy.mean_anomaly * constants::RAD_TO_DEG,
        initial.gm
    );

    return compute_observation(initial_state_kep_ecl, t_elements, t_obs, e_cfg, a_cfg);
}

Eigen::Vector3d AstrometryReducer::compute_light_time_corrected_pos(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
    const time::EpochTDB& t_elements,
    const time::EpochTDB& t_obs,
    const Eigen::Vector3d& earth_pos_helio_ecl,
    const AstDynConfig& cfg) {
    
    AstDynEngine engine(cfg);
    engine.set_initial_orbit(initial);

    double tau_days = 0.0; 
    Eigen::Vector3d ast_p_ecl;
    
    for (int i = 0; i < 3; ++i) {
        auto el_kep_ecl = engine.propagate_to(time::EpochTDB::from_mjd(t_obs.mjd() - tau_days));
        
        // Convert Keplerian Ecliptic to Cartesian Ecliptic for distance check using type-safe API
        auto cart_ecl = propagation::keplerian_to_cartesian(el_kep_ecl);
        ast_p_ecl = cart_ecl.position.to_eigen_si();
        
        tau_days = (ast_p_ecl - earth_pos_helio_ecl).norm() / (constants::C_LIGHT * 1000.0 * 86400.0);
    }
    return ast_p_ecl;
}

Eigen::Vector3d AstrometryReducer::apply_stellar_aberration(
    const Eigen::Vector3d& rho_eq, const Eigen::Vector3d& earth_vel_eq) {
    double r = rho_eq.norm();
    Eigen::Vector3d u = rho_eq / r;
    Eigen::Vector3d v_c = earth_vel_eq / (C_LIGHT * 1000.0); // m/s to c
    
    double d_uv = u.dot(v_c);
    double beta_inv = std::sqrt(1.0 - v_c.squaredNorm());
    
    // Lorentz transform of direction vector (approx. SR aberration)
    return ((beta_inv * u + (1.0 + d_uv / (1.0 + beta_inv)) * v_c) / (1.0 + d_uv)) * r;
}

AstrometricObservation AstrometryReducer::finalize_observation(
    const Eigen::Vector3d& rho_eq) {
    double r = rho_eq.norm();
    double ra = std::atan2(rho_eq(1), rho_eq(0)); 
    if (ra < 0) ra += constants::TWO_PI;
    double dec = std::asin(rho_eq(2) / r);
    
    return { core::Radian(ra), core::Radian(dec), core::Meter(r) };
}

} // namespace astdyn::astrometry
