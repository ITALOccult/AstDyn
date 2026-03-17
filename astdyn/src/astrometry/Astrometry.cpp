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
    auto p_earth_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(p_earth);
    auto p_sun_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(p_sun);
    return p_earth_ecl.to_eigen_si() - p_sun_ecl.to_eigen_si();
}

std::expected<AstrometricObservation, AstrometryError> AstrometryReducer::compute_observation(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial,
    const time::EpochTDB& t_elements, const time::EpochTDB& t_obs,
    const AstDynConfig& e_cfg, const AstrometricSettings& a_cfg) 
{
    auto de441 = AstrometryReducer::sync_ephemeris(e_cfg.ephemeris_file);
    if (!de441) return std::unexpected(AstrometryError::EphemerisUnavailable);
    
    Eigen::Vector3d earth_helio = AstrometryReducer::compute_earth_helio(de441, t_obs);
    Eigen::Vector3d ast_pos = compute_light_time_corrected_pos(initial, t_elements, t_obs, earth_helio, e_cfg);
    
    auto rho_ecl = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(ast_pos.x() - earth_helio.x(), ast_pos.y() - earth_helio.y(), ast_pos.z() - earth_helio.z());
    auto rho_eq = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(rho_ecl).to_eigen_si();
    
    Eigen::Vector3d q_sun = de441->getPosition(::astdyn::ephemeris::CelestialBody::SUN, t_obs).to_eigen_si() - de441->getPosition(::astdyn::ephemeris::CelestialBody::EARTH, t_obs).to_eigen_si();
    if (a_cfg.light_deflection) rho_eq = apply_light_deflection(rho_eq, q_sun);
    if (a_cfg.stellar_aberration) rho_eq = apply_stellar_aberration(rho_eq, de441->getVelocity(::astdyn::ephemeris::CelestialBody::EARTH, t_obs).to_eigen_si());

    return finalize_observation(rho_eq);
}

std::expected<AstrometricObservation, AstrometryError> AstrometryReducer::compute_observation_from_cartesian(
    const physics::CartesianStateTyped<core::GCRF>& initial,
    const time::EpochTDB& t_elements,
    const time::EpochTDB& t_obs,
    const AstDynConfig& e_cfg,
    const AstrometricSettings& a_cfg) {
    
    if (e_cfg.ephemeris_file.empty()) return std::unexpected(AstrometryError::EphemerisUnavailable);
    auto de441 = AstrometryReducer::sync_ephemeris(e_cfg.ephemeris_file);
    if (!de441) return std::unexpected(AstrometryError::EphemerisUnavailable);

    // Transform initial state to Ecliptic properly
    auto pos_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(initial.position);
    auto vel_ecl = coordinates::ReferenceFrame::transform_vel<core::GCRF, core::ECLIPJ2000>(initial.position, initial.velocity);
    
    // Bridge to Engine's expected Keplerian format (AstDynEngine currently fits/propagates in Ecliptic)
    propagation::CartesianElements cart_ecl_old;
    cart_ecl_old.epoch = t_elements;
    // Legacy CartesianElements expects GCRF tag in its type, but we pass the Ecliptic values
    // to match what the engine expects for its internal Ecliptic-based solver.
    cart_ecl_old.position = math::Vector3<core::GCRF, physics::Distance>::from_si(pos_ecl.x_si(), pos_ecl.y_si(), pos_ecl.z_si());
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
    if (engine_cfg.light_deflection) rho_eq = apply_light_deflection(rho_eq, q_sun);
    if (engine_cfg.aberration_correction) rho_eq = apply_stellar_aberration(rho_eq, de441->getVelocity(ephemeris::CelestialBody::EARTH, t_obs).to_eigen_si());

    return finalize_observation(rho_eq);
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
        
        auto cart_ecl = propagation::keplerian_to_cartesian(el_kep_ecl);
        ast_p_ecl = cart_ecl.position.to_eigen_si();
        
        tau_days = (ast_p_ecl - earth_pos_helio_ecl).norm() / (physics::SpeedOfLight::to_ms() * 86400.0);
    }
    return ast_p_ecl;
}

Eigen::Vector3d AstrometryReducer::apply_stellar_aberration(
    const Eigen::Vector3d& rho_eq, const Eigen::Vector3d& earth_vel_eq) {
    
    // We use the rigorous IAU 2000 definition:
    // p' = (p + V/c + (1 + p.V/c)/(1 + 1/gamma) * V/c) / (1 + p.V/c)
    // where p is the UNIT vector toward the object.
    
    double r = rho_eq.norm();
    Eigen::Vector3d p = rho_eq / r;
    Eigen::Vector3d v = earth_vel_eq / physics::SpeedOfLight::to_ms();
    
    double p_dot_v = p.dot(v);
    double v2 = v.squaredNorm();
    double inv_gamma = std::sqrt(1.0 - v2);
    
    // Exact Relativistic formula
    double denom = 1.0 + p_dot_v;
    Eigen::Vector3d p_prime = (inv_gamma * p + (1.0 + p_dot_v / (1.0 + inv_gamma)) * v) / denom;
    
    return p_prime * r;
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

Eigen::Vector3d AstrometryReducer::apply_light_deflection(
    const Eigen::Vector3d& rho_eq, const Eigen::Vector3d& earth_to_sun_eq) {
    
    // Standard Gravitational Deflection (Sun)
    // dU = 2*GM/(c^2*b) * ( (1 + cos psi) * (U x Un) )
    
    double r = rho_eq.norm();
    Eigen::Vector3d u = rho_eq / r;
    Eigen::Vector3d q = earth_to_sun_eq; // Vector Earth -> Sun
    double q_dist = q.norm();
    if (q_dist < 1000.0) return rho_eq; // Avoid NaNs if Earth is Sun (should not happen)
    
    Eigen::Vector3d e = q / q_dist;
    
    double u_dot_e = u.dot(e);
    
    // 2*mu/c^2 ~ 3 km
    double g = 2.0 * constants::GM_SUN / (std::pow(physics::SpeedOfLight::to_ms(), 2));
    
    // Formula from Kaplan / IAU
    double fac = g / q_dist * (1.0 + u_dot_e);
    // This is small, so we return the shifted vector
    return (u + fac * (e - u_dot_e * u)).normalized() * r;
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
