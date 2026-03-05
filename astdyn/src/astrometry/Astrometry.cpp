#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>

namespace astdyn::astrometry {

using namespace astdyn::constants;

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
    const utils::Instant& t_obs,
    const AstDynConfig& e_cfg,
    const AstrometricSettings& a_cfg) {
    
    // 1. Get Earth state at observation time (Native DE441)
    if (e_cfg.ephemeris_file.empty()) return std::unexpected(AstrometryError::EphemerisUnavailable);
    ephemeris::DE441Provider de441(e_cfg.ephemeris_file);
    auto earth_pos = de441.getPosition(ephemeris::CelestialBody::EARTH, t_obs).to_eigen(); 
    auto earth_vel = de441.getVelocity(ephemeris::CelestialBody::EARTH, t_obs).to_eigen();

    // 2. Light-time Iteration (Step 1)
    auto ast_pos = compute_light_time_corrected_pos(initial, t_obs, earth_pos, e_cfg);
    auto raw_rho = ast_pos - earth_pos;

    // 3. Aberration (Step 2)
    auto aberr_rho = a_cfg.stellar_aberration ? apply_stellar_aberration(raw_rho, earth_vel) : raw_rho;

    // 4. Frame Reduction (Step 3)
    auto final_rho = convert_frame_if_needed(aberr_rho, a_cfg);

    // 5. Build Result (Step 4)
    return finalize_observation(final_rho);
}

Eigen::Vector3d AstrometryReducer::compute_light_time_corrected_pos(
    const types::OrbitalState<core::ECLIPJ2000, types::KeplerianTag>& initial,
    const utils::Instant& t_obs,
    const Eigen::Vector3d& earth_pos,
    const AstDynConfig& cfg) {
    AstDynEngine engine(cfg);
    engine.set_initial_orbit(state_to_kep(initial, t_obs));
    double tau_days = 0.0; Eigen::Vector3d ast_p;
    for (int i = 0; i < 3; ++i) {
        auto el = engine.propagate_to(t_obs.mjd.value - tau_days);
        ast_p = propagation::keplerian_to_cartesian(el).position.to_eigen(); 
        tau_days = (ast_p - earth_pos).norm() / (C_LIGHT * 1000.0 * 86400.0);
    }
    return ast_p;
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
    
    if (!settings.frame_conversion_to_equatorial) return vec;
    // Rotation from Ecliptic to J2000 Equatorial
    return coordinates::ReferenceFrame::ecliptic_to_j2000() * vec;
}

AstrometricObservation AstrometryReducer::finalize_observation(
    const Eigen::Vector3d& rho) {
    
    double r = rho.norm();
    double ra = std::atan2(rho.y(), rho.x());
    double dec = std::asin(rho.z() / r);
    if (ra < 0) ra += 2.0 * PI;
    
    return { core::Radian(ra), core::Radian(dec), core::Meter(r) };
}

} // namespace astdyn::astrometry
