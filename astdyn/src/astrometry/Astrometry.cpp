#include "astdyn/astrometry/Astrometry.hpp"
#include <expected>
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/observations/ObservatoryDatabase.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/astrometry/AstrometricTypes.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/io/SPKReader.hpp"
#include "astdyn/ephemeris/AsteroidPerturbations.hpp"
#include "astdyn/astrometry/AstrometricCorrections.hpp"
#include "astdyn/astrometry/AstrometricCorrector.hpp"
#include <cmath>
#include <mutex>
#include <unordered_map>
#include <iostream>
#include <iomanip>

namespace astdyn::astrometry {

using namespace astdyn::constants;

std::shared_ptr<::astdyn::ephemeris::DE441Provider> AstrometryReducer::sync_ephemeris(const std::string& path) {
    if (path.empty()) return nullptr;
    static std::mutex mtx;
    static std::unordered_map<std::string, std::shared_ptr<::astdyn::ephemeris::DE441Provider>> cache;
    std::lock_guard lock(mtx);
    if (!cache.contains(path)) cache[path] = std::make_shared<::astdyn::ephemeris::DE441Provider>(path);
    ::astdyn::ephemeris::PlanetaryEphemeris::setGlobalProvider(cache[path]);
    return cache[path];
}

Eigen::Vector3d AstrometryReducer::compute_earth_helio(std::shared_ptr<::astdyn::ephemeris::DE441Provider> de441, const time::EpochTDB& t_obs) {
    auto p_earth = de441->getPosition(::astdyn::ephemeris::CelestialBody::EARTH, t_obs);
    auto p_sun = de441->getPosition(::astdyn::ephemeris::CelestialBody::SUN, t_obs);
    
    // Displacement in GCRF
    auto r_gcrf = math::Vector3<core::GCRF, physics::Distance>::from_si(
        p_earth.x_si() - p_sun.x_si(), 
        p_earth.y_si() - p_sun.y_si(), 
        p_earth.z_si() - p_sun.z_si()
    );
    
    return coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(r_gcrf).to_eigen_si();
}

std::expected<AstrometricObservation, AstrometryError> AstrometryReducer::compute_observation(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
    const time::EpochTDB& t_elements, const time::EpochTDB& t_obs,
    const AstDynConfig& e_cfg, const AstrometricSettings& a_cfg) 
{
    auto de441 = sync_ephemeris(e_cfg.ephemeris_file);
    if (!de441) return std::unexpected(AstrometryError::EphemerisUnavailable);
    
    Eigen::Vector3d earth_helio = compute_earth_helio(de441, t_obs);
    Eigen::Vector3d ast_pos = compute_light_time_corrected_pos(initial, t_elements, t_obs, earth_helio, e_cfg);
    
    Eigen::Vector3d diff = ast_pos - earth_helio;
    auto rho_ecl = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(diff.x(), diff.y(), diff.z());
    auto rho_eq = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(rho_ecl).to_eigen_si();
    
    Eigen::Vector3d q_sun = de441->getPosition(ephemeris::CelestialBody::SUN, t_obs).to_eigen_si() - 
                            de441->getPosition(ephemeris::CelestialBody::EARTH, t_obs).to_eigen_si();
    Eigen::Vector3d v_earth = de441->getVelocity(ephemeris::CelestialBody::EARTH, t_obs).to_eigen_si();
    rho_eq = AstrometricCorrector(a_cfg).apply(rho_eq, v_earth, q_sun);

    return finalize_observation(rho_eq);
}

std::expected<AstrometricObservation, AstrometryError> AstrometryReducer::compute_observation_from_cartesian(
    const physics::CartesianStateTyped<core::GCRF>& initial,
    const time::EpochTDB& t_elements, const time::EpochTDB& t_obs,
    const AstDynConfig& e_cfg, const AstrometricSettings& a_cfg) 
{
    if (e_cfg.ephemeris_file.empty()) return std::unexpected(AstrometryError::EphemerisUnavailable);
    
    auto pos_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(initial.position);
    auto vel_ecl = coordinates::ReferenceFrame::transform_vel<core::GCRF, core::ECLIPJ2000>(initial.position, initial.velocity);
    auto cart_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>(t_elements, pos_ecl, vel_ecl, initial.gm);
    auto initial_kep = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(cart_ecl);

    return compute_observation(initial_kep, t_elements, t_obs, e_cfg, a_cfg);
}

std::expected<AstrometricObservation, AstrometryError> AstrometryReducer::compute_topocentric_observation(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
    const time::EpochTDB& t_elements, const time::EpochTDB& t_obs,
    const std::string& obs_code, const AstDynConfig& engine_cfg) 
{
    auto de441 = AstrometryReducer::sync_ephemeris(engine_cfg.ephemeris_file);
    if (!de441) return std::unexpected(AstrometryError::EphemerisUnavailable);

    // Fetch observatory position
    Eigen::Vector3d obs_pos_gcrf = Eigen::Vector3d::Zero();
    if (!obs_code.empty() && obs_code != "500" && obs_code != "@ssb") {
        auto& db = observations::ObservatoryDatabase::getInstance();
        auto obs = db.getObservatory(obs_code);
        if (obs) obs_pos_gcrf = obs->getPositionGCRF(time::to_utc(t_obs)).to_eigen_si();
    }

    Eigen::Vector3d earth_bary = de441->getState(ephemeris::CelestialBody::EARTH, t_obs).position.to_eigen_si();
    Eigen::Vector3d sun_bary = de441->getState(ephemeris::CelestialBody::SUN, t_obs).position.to_eigen_si();
    Eigen::Vector3d observer_bary = earth_bary + obs_pos_gcrf;
    Eigen::Vector3d observer_helio_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(
        math::Vector3<core::GCRF, physics::Distance>::from_si(observer_bary.x() - sun_bary.x(), observer_bary.y() - sun_bary.y(), observer_bary.z() - sun_bary.z())
    ).to_eigen_si();

    Eigen::Vector3d ast_pos_helio = compute_light_time_corrected_pos(initial, t_elements, t_obs, observer_helio_ecl, engine_cfg);
    
    auto rho_helio_ecl = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(ast_pos_helio.x() - observer_helio_ecl.x(), ast_pos_helio.y() - observer_helio_ecl.y(), ast_pos_helio.z() - observer_helio_ecl.z());
    auto rho_eq = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(rho_helio_ecl).to_eigen_si();
    
    Eigen::Vector3d q_sun = sun_bary - (earth_bary + obs_pos_gcrf);
    Eigen::Vector3d v_earth = de441->getVelocity(ephemeris::CelestialBody::EARTH, t_obs).to_eigen_si();
    
    AstrometricSettings a_settings;
    a_settings.aberrazione_differenziale = engine_cfg.aberrazione_differenziale;
    a_settings.deflessione_relativistica = engine_cfg.deflessione_relativistica;
    
    rho_eq = AstrometricCorrector(a_settings).apply(rho_eq, v_earth, q_sun);

    return finalize_observation(rho_eq);
}

Eigen::Vector3d AstrometryReducer::compute_light_time_corrected_pos(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
    const time::EpochTDB& t_elements, const time::EpochTDB& t_obs,
    const Eigen::Vector3d& obs_helio, const AstDynConfig& cfg) 
{
    AstDynEngine engine(cfg);
    engine.set_initial_orbit(initial);
    double tau = 0.0; Eigen::Vector3d p_ast;
    
    // 5 iterations for high precision (usually 3 is enough for asteroids, 5 for comets/close approaches)
    for (int i = 0; i < 5; ++i) {
        time::EpochTDB t_emit = time::EpochTDB::from_mjd(t_obs.mjd() - tau);
        auto els = engine.propagate_to(t_emit);
        p_ast = propagation::keplerian_to_cartesian(els).position.to_eigen_si();
        tau = (p_ast - obs_helio).norm() / (C_LIGHT * SECONDS_PER_DAY);
    }
    return p_ast;
}

Eigen::Vector3d AstrometryReducer::aberrazione_differenziale(
    const Eigen::Vector3d& rho_eq, const Eigen::Vector3d& earth_vel_eq) {
    return ::astdyn::astrometry::aberrazione_differenziale(rho_eq, earth_vel_eq);
}

double AstrometryReducer::compute_cross_track_uncertainty(
    const Eigen::Matrix<double, 6, 6>& covariance_eq,
    const Eigen::Vector3d& rho_eq,
    const Eigen::Vector3d& velocity_eq,
    double star_ra_rad,
    double star_dec_rad)
{
    // 1. Star unit vector (shadow axis)
    Eigen::Vector3d s_hat(
        std::cos(star_dec_rad) * std::cos(star_ra_rad),
        std::cos(star_dec_rad) * std::sin(star_ra_rad),
        std::sin(star_dec_rad)
    );

    // 2. Build Fundamental Plane (FP) basis
    // u = East direction in FP orthogonal to s_hat
    Eigen::Vector3d u_hat = Eigen::Vector3d(-std::sin(star_ra_rad), std::cos(star_ra_rad), 0.0);
    // v = North/South direction in FP
    Eigen::Vector3d v_hat = s_hat.cross(u_hat);

    // 3. Project velocity into FP
    // Velocity of shadow in FP (assuming Earth center at rest for 1-sigma report)
    // velocity_eq is m/s, convert to km/s for covariance matching
    Eigen::Vector3d v_km_s = velocity_eq / 1000.0;
    double vx = v_km_s.dot(u_hat);
    double vy = v_km_s.dot(v_hat);
    
    // Shadow path direction in FP
    double v_mag_fp = std::sqrt(vx*vx + vy*vy);
    Eigen::Vector3d v_path_fp = (vx * u_hat + vy * v_hat) / v_mag_fp;
    
    // Perpendicular (cross-track) direction in FP
    Eigen::Vector3d cross_track_dir = s_hat.cross(v_path_fp).normalized();

    // 4. Project covariance matrix onto the cross-track direction
    // P_pos = top-left 3x3 of covariance (km^2)
    Eigen::Matrix3d P_pos = covariance_eq.block<3, 3>(0, 0);
    
    // sigma^2 = n^T * P * n
    double sigma2 = cross_track_dir.transpose() * P_pos * cross_track_dir;
    
    return std::sqrt(sigma2);
}

Eigen::Vector3d AstrometryReducer::deflessione_relativistica(
    const Eigen::Vector3d& rho_eq, const Eigen::Vector3d& earth_to_sun_eq) 
{
    return ::astdyn::astrometry::deflessione_relativistica(rho_eq, earth_to_sun_eq);
}

AstrometricObservation AstrometryReducer::finalize_observation(
    const Eigen::Vector3d& rho_eq) {
    double r = rho_eq.norm();
    double ra = std::atan2(rho_eq(1), rho_eq(0)); 
    if (ra < 0) ra += constants::TWO_PI;
    double dec = std::asin(rho_eq(2) / r);

    AstrometricObservation obs;
    obs.ra = RightAscension::from_rad(ra);
    obs.dec = Declination::from_rad(dec);
    obs.distance = physics::Distance::from_m(r);
    return obs;
}

} // namespace astdyn::astrometry
