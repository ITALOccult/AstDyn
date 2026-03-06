#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
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

static propagation::KeplerianElements state_to_kep(
    const types::OrbitalState<core::ECLIPJ2000, types::KeplerianTag>& s, const utils::Instant& t) {
    propagation::KeplerianElements k;
    k.semi_major_axis = s.a(); k.eccentricity = s.e(); k.inclination = s.i();
    k.longitude_ascending_node = s.raan(); k.argument_perihelion = s.arg_peri();
    k.mean_anomaly = s.m_anomaly(); k.epoch = t; 
    k.gravitational_parameter = GMS;
    return k;
}

std::expected<AstrometricObservation, AstrometryError> AstrometryReducer::compute_observation(
    const types::OrbitalState<core::ECLIPJ2000, types::KeplerianTag>& initial,
    const utils::Instant& t_elements,
    const utils::Instant& t_obs,
    const AstDynConfig& e_cfg,
    const AstrometricSettings& a_cfg) {
    
    // 1. Get Earth state at observation time (Cached Native DE441)
    if (e_cfg.ephemeris_file.empty()) return std::unexpected(AstrometryError::EphemerisUnavailable);
    auto de441 = get_cached_provider(e_cfg.ephemeris_file);
    if (!de441) return std::unexpected(AstrometryError::EphemerisUnavailable);
    
    // Sync PlanetaryEphemeris provider too
    ephemeris::PlanetaryEphemeris::setProvider(de441);

    auto mat_ecl = coordinates::ReferenceFrame::j2000_to_ecliptic();
    auto earth_p_ssb = (mat_ecl * de441->getPosition(ephemeris::CelestialBody::EARTH, t_obs).to_eigen()).eval(); 
    auto earth_v_ssb = (mat_ecl * de441->getVelocity(ephemeris::CelestialBody::EARTH, t_obs).to_eigen()).eval();
    
    // 2. SUN position (to convert Earth to HELIOCENTRIC)
    auto sun_p_ssb = (mat_ecl * ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(t_obs).to_eigen()).eval();
    Eigen::Vector3d earth_pos = earth_p_ssb - sun_p_ssb;

    // 3. Light-time Iteration (Step 1)
    auto ast_pos = compute_light_time_corrected_pos(initial, t_elements, t_obs, earth_pos, e_cfg);
    
    // DEBUG
    std::cout << "    [DEBUG-ATR] Earth Helio (Ecl): " << (earth_pos.norm() / (constants::AU*1000.0)) << " AU" << std::endl;
    std::cout << "    [DEBUG-ATR] Ast   Helio (Ecl): " << (ast_pos.norm() / (constants::AU*1000.0)) << " AU" << std::endl;

    auto raw_rho = ast_pos - earth_pos;
    std::cout << "    [DEBUG-ATR] Geoc. Dist (Ecl): " << (raw_rho.norm() / (constants::AU*1000.0)) << " AU" << std::endl;

    // 4. Aberration (Step 2)
    auto aberr_rho = a_cfg.stellar_aberration ? apply_stellar_aberration(raw_rho, earth_v_ssb) : raw_rho;

    // 4. Frame Reduction (Step 3)
    auto final_rho = convert_frame_if_needed(aberr_rho, a_cfg);

    // 5. Build Result (Step 4)
    return finalize_observation(final_rho);
}

std::expected<AstrometricObservation, AstrometryError> AstrometryReducer::compute_observation_from_cartesian(
    const physics::CartesianStateTyped<core::GCRF>& initial,
    const utils::Instant& t_elements,
    const utils::Instant& t_obs,
    const AstDynConfig& e_cfg,
    const AstrometricSettings& a_cfg) {
    
    if (e_cfg.ephemeris_file.empty()) return std::unexpected(AstrometryError::EphemerisUnavailable);
    auto de441 = get_cached_provider(e_cfg.ephemeris_file);
    if (!de441) return std::unexpected(AstrometryError::EphemerisUnavailable);
    ephemeris::PlanetaryEphemeris::setProvider(de441);

    auto mat_ecl = coordinates::ReferenceFrame::j2000_to_ecliptic();
    auto earth_p_ssb = (mat_ecl * de441->getPosition(ephemeris::CelestialBody::EARTH, t_obs).to_eigen()).eval(); 
    auto earth_v_ssb = (mat_ecl * de441->getVelocity(ephemeris::CelestialBody::EARTH, t_obs).to_eigen()).eval();
    auto sun_p_ssb = (mat_ecl * ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(t_obs).to_eigen()).eval();
    Eigen::Vector3d earth_pos = earth_p_ssb - sun_p_ssb;

    // Convert to Keplerian (Ecliptic J2000)
    // First ensure we rotate from GCRF to Ecliptic
    auto pos_eq = initial.position.to_eigen_si();
    auto vel_eq = initial.velocity.to_eigen_si();
    auto mat_ecl_rot = coordinates::ReferenceFrame::j2000_to_ecliptic();
    
    propagation::CartesianElements cart_ecl;
    cart_ecl.epoch = t_elements;
    cart_ecl.position = types::Vector3<core::GCRF, core::Meter>(mat_ecl_rot * pos_eq);
    cart_ecl.velocity = types::Vector3<core::GCRF, core::Meter>(mat_ecl_rot * vel_eq);
    cart_ecl.gravitational_parameter = initial.gm.to_m3_s2();
    
    auto kep_ecl = propagation::cartesian_to_keplerian(cart_ecl);
    
    AstDynEngine engine(e_cfg);
    auto engine_init = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(kep_ecl.epoch.mjd.value),
        kep_ecl.semi_major_axis, kep_ecl.eccentricity,
        kep_ecl.inclination * constants::RAD_TO_DEG,
        kep_ecl.longitude_ascending_node * constants::RAD_TO_DEG,
        kep_ecl.argument_perihelion * constants::RAD_TO_DEG,
        kep_ecl.mean_anomaly * constants::RAD_TO_DEG,
        physics::GravitationalParameter::sun()
    );
    engine.set_initial_orbit(engine_init);

    auto ast_pos = compute_light_time_corrected_pos_internal(engine, t_obs, earth_pos);
    
    // DEBUG
    if (e_cfg.verbose) {
        std::cout << "    [DEBUG-ATR-V] Earth Helio (Ecl): " << (earth_pos.norm() / (constants::AU*1000.0)) << " AU" << std::endl;
        std::cout << "    [DEBUG-ATR-V] Ast   Helio (Ecl): " << (ast_pos.norm() / (constants::AU*1000.0)) << " AU" << std::endl;
    }

    auto raw_rho = ast_pos - earth_pos;
    if (e_cfg.verbose) {
        std::cout << "    [DEBUG-ATR-V] Geoc. Dist (Ecl): " << (raw_rho.norm() / (constants::AU*1000.0)) << " AU" << std::endl;
    }

    auto aberr_rho = a_cfg.stellar_aberration ? apply_stellar_aberration(raw_rho, earth_v_ssb) : raw_rho;
    auto final_rho = convert_frame_if_needed(aberr_rho, a_cfg);

    return finalize_observation(final_rho);
}

Eigen::Vector3d AstrometryReducer::compute_light_time_corrected_pos_internal(
    AstDynEngine& engine, const utils::Instant& t_obs, const Eigen::Vector3d& earth_pos) {
    
    // We work in Ecliptic J2000 for light-time
    auto mat_eq = coordinates::ReferenceFrame::ecliptic_to_j2000();
    auto mat_ecl = coordinates::ReferenceFrame::j2000_to_ecliptic();
    Eigen::Vector3d earth_pos_eq = mat_eq * earth_pos;
    
    double tau_days = 0.0; Eigen::Vector3d ast_p;
    for (int i = 0; i < 3; ++i) {
        auto el = engine.propagate_to(time::EpochTDB::from_mjd(t_obs.mjd.value - tau_days));
        
        // Convert typed ECLIP Keplerian to Equatorial Cartesian for LOS
        propagation::KeplerianElements kep_old;
        kep_old.epoch = utils::Instant::from_tt(utils::ModifiedJulianDate(el.epoch.mjd()));
        kep_old.semi_major_axis = el.a.to_au();
        kep_old.eccentricity = el.e;
        kep_old.inclination = el.i.to_rad();
        kep_old.longitude_ascending_node = el.node.to_rad();
        kep_old.argument_perihelion = el.omega.to_rad();
        kep_old.mean_anomaly = el.M.to_rad();
        kep_old.gravitational_parameter = el.gm.to_au3_d2();
        
        auto cart_ecl = propagation::keplerian_to_cartesian(kep_old);
        ast_p = mat_eq * cart_ecl.position.to_eigen(); // Equatorial
        tau_days = (ast_p - earth_pos_eq).norm() / (constants::C_LIGHT * 1000.0 * 86400.0);
    }
    return mat_ecl * ast_p;
}

Eigen::Vector3d AstrometryReducer::compute_light_time_corrected_pos(
    const types::OrbitalState<core::ECLIPJ2000, types::KeplerianTag>& initial,
    const utils::Instant& t_elements,
    const utils::Instant& t_obs,
    const Eigen::Vector3d& earth_pos,
    const AstDynConfig& cfg) {
    
    AstDynEngine engine(cfg);
    auto engine_init = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(t_elements.mjd.value),
        initial.a(), initial.e(), initial.i() * constants::RAD_TO_DEG,
        initial.raan() * constants::RAD_TO_DEG, initial.arg_peri() * constants::RAD_TO_DEG,
        initial.m_anomaly() * constants::RAD_TO_DEG,
        physics::GravitationalParameter::sun()
    );
    engine.set_initial_orbit(engine_init);

    double tau_days = 0.0; Eigen::Vector3d ast_p;
    auto mat_eq = coordinates::ReferenceFrame::ecliptic_to_j2000();
    auto mat_ecl = coordinates::ReferenceFrame::j2000_to_ecliptic();
    Eigen::Vector3d earth_pos_eq = mat_eq * earth_pos;

    for (int i = 0; i < 3; ++i) {
        auto el = engine.propagate_to(time::EpochTDB::from_mjd(t_obs.mjd.value - tau_days));
        
        // Convert typed ECLIP Keplerian to Equatorial Cartesian for LOS
        propagation::KeplerianElements kep_old;
        kep_old.epoch = utils::Instant::from_tt(utils::ModifiedJulianDate(el.epoch.mjd()));
        kep_old.semi_major_axis = el.a.to_au();
        kep_old.eccentricity = el.e;
        kep_old.inclination = el.i.to_rad();
        kep_old.longitude_ascending_node = el.node.to_rad();
        kep_old.argument_perihelion = el.omega.to_rad();
        kep_old.mean_anomaly = el.M.to_rad();
        kep_old.gravitational_parameter = el.gm.to_au3_d2();
        
        auto cart_ecl = propagation::keplerian_to_cartesian(kep_old);
        ast_p = mat_eq * cart_ecl.position.to_eigen(); // Equatorial
        tau_days = (ast_p - earth_pos_eq).norm() / (constants::C_LIGHT * 1000.0 * 86400.0);
    }
    return mat_ecl * ast_p;
}

Eigen::Vector3d AstrometryReducer::apply_stellar_aberration(
    const Eigen::Vector3d& rho, const Eigen::Vector3d& earth_vel) {
    double r = rho.norm();
    Eigen::Vector3d u = rho / r, v_c = earth_vel / (C_LIGHT * 1000.0);
    double d_uv = u.dot(v_c), b = std::sqrt(1.0 - v_c.squaredNorm());
    return ((b * u + (1.0 + d_uv / (1.0 + b)) * v_c) / (1.0 + d_uv)) * r;
}

Eigen::Vector3d AstrometryReducer::convert_frame_if_needed(
    const Eigen::Vector3d& vec,
    const AstrometricSettings& settings) {
    
    if (!settings.frame_conversion_to_equatorial) {
        std::cout << "[DEBUG-ATR] NO FRAME CONVERSION APPLIED!\n";
        return vec;
    }
    
    auto R = coordinates::ReferenceFrame::ecliptic_to_j2000();
    auto result = R * vec;
    
    std::cout << "[DEBUG-ATR] R:\n" << R << "\n";
    std::cout << "  Input  (Ecl): " << vec.transpose() << "\n";
    std::cout << "  Output (Eq) : " << result.transpose() << "\n";
    
    return result;
}

AstrometricObservation AstrometryReducer::finalize_observation(
    const Eigen::Vector3d& rho) {
    double r = rho.norm(), ra = std::atan2(rho.y(), rho.x()); if (ra < 0) ra += TWO_PI;
    return { core::Radian(ra), core::Radian(std::asin(rho.z() / r)), core::Meter(r) };
}

} // namespace astdyn::astrometry
