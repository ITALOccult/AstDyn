#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/io/MPCClient.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/io/SPKReader.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/astrometry/ChebyshevEphemerisManager.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/AstDynEngine.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace astdyn::astrometry {

namespace {

    double compute_squared_angular_distance(double t, const astdyn::catalog::ChebyshevSegment& segment, double target_ra, double target_dec) {
        auto [pos, vel] = segment.evaluate_full(t);
        double dra = std::get<0>(pos) - target_ra;
        double ddec = std::get<1>(pos) - target_dec;
        if (dra > 180.0) dra -= 360.0;
        if (dra < -180.0) dra += 360.0;
        dra *= std::cos(std::get<1>(pos) * astdyn::constants::DEG_TO_RAD);
        return dra*dra + ddec*ddec;
    }

    double find_tca(const astdyn::catalog::ChebyshevSegment& segment, const astdyn::catalog::Star& star, double t_start, double t_end) {
        auto target_ra = star.ra.to_deg();
        auto target_dec = star.dec.to_deg();
        auto dist_sq = [&](double t) { return compute_squared_angular_distance(t, segment, target_ra, target_dec); };

        double best_t = t_start + (t_end - t_start) / 2.0;
        double min_d2 = 1e18;
        int samples = 48; 
        for (int i = 0; i <= samples; ++i) {
            double t = t_start + (t_end - t_start) * i / samples;
            double d2 = dist_sq(t);
            if (d2 < min_d2) {
                min_d2 = d2;
                best_t = t;
            }
        }
        
        double step = (t_end - t_start) / samples / 2.0;
        for (int iter = 0; iter < 10; ++iter) {
            double t1 = best_t - step;
            double t2 = best_t + step;
            double d1 = (t1 >= t_start) ? dist_sq(t1) : 1e18;
            double d2 = (t2 <= t_end) ? dist_sq(t2) : 1e18;
            if (d1 < min_d2 && d1 < d2) {
                min_d2 = d1;
                best_t = t1;
            } else if (d2 < min_d2) {
                min_d2 = d2;
                best_t = t2;
            }
            step /= 2.0;
        }
        return best_t;
    }

} // namespace

OccultationParameters OccultationLogic::compute_parameters(
    const RightAscension& star_ra, const Declination& star_dec,
    const RightAscension& ast_ra, const Declination& ast_dec,
    const physics::Distance& ast_dist,
    const Angle& ast_dra_dt, const Angle& ast_ddec_dt,
    const physics::Velocity& ast_ddist_dt,
    const time::EpochTDB& t_ca,
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem)
{
    OccultationParameters params = compute_fundamental_plane_geometry(star_ra, star_dec, ast_ra, ast_dec, ast_dist);
    params.t_ca = t_ca;

    compute_shadow_velocity(params, star_dec, ast_dist, ast_dra_dt, ast_ddec_dt);
    
    params.total_apparent_rate = std::sqrt(ast_dra_dt.to_arcsec() * ast_dra_dt.to_arcsec() + 
                                           ast_ddec_dt.to_arcsec() * ast_ddec_dt.to_arcsec());
    
    if (ephem) {
        compute_sky_conditions(params, t_ca, ast_ra, ast_dec, ast_dist, star_ra, star_dec, ephem);
        compute_sub_asteroid_point(params, t_ca, ast_ra, ast_dec, ast_dist);
    }
    return params;
}

void OccultationLogic::compute_sub_asteroid_point(
    OccultationParameters& params,
    const time::EpochTDB& t_ca,
    const RightAscension& ast_ra, const Declination& ast_dec,
    const physics::Distance& ast_dist)
{
    double a = ast_ra.to_rad();
    double d = ast_dec.to_rad();
    Eigen::Vector3d r_ast_vec = ast_dist.to_m() * Eigen::Vector3d(std::cos(d) * std::cos(a), std::cos(d) * std::sin(a), std::sin(d));
    
    auto pos_ast_gcrf = math::Vector3<core::GCRF, physics::Distance>::from_si(r_ast_vec.x(), r_ast_vec.y(), r_ast_vec.z());
    auto pos_ast_itrf = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ITRF>(pos_ast_gcrf, t_ca);
    
    double r_itrf = pos_ast_itrf.norm().to_m();
    if (r_itrf > 1.0) {
        params.center_lat = Angle::from_rad(std::asin(pos_ast_itrf.z_si() / r_itrf));
        params.center_lon = Angle::from_rad(std::atan2(pos_ast_itrf.y_si(), pos_ast_itrf.x_si()));
    }
}

OccultationParameters OccultationLogic::compute_fundamental_plane_geometry(
    const RightAscension& star_ra, const Declination& star_dec,
    const RightAscension& ast_ra, const Declination& ast_dec,
    const physics::Distance& ast_dist)
{
    double as = star_ra.to_rad();
    double ds = star_dec.to_rad();
    double a = ast_ra.to_rad();
    double d = ast_dec.to_rad();

    Eigen::Vector3d k(std::cos(ds) * std::cos(as), std::cos(ds) * std::sin(as), std::sin(ds));
    Eigen::Vector3d r_ast = ast_dist.to_m() * Eigen::Vector3d(std::cos(d) * std::cos(a), std::cos(d) * std::sin(a), std::sin(d));
    
    // East/North Basis in FP
    Eigen::Vector3d east(-std::sin(as), std::cos(as), 0.0);
    Eigen::Vector3d north = k.cross(east);

    // Coordinate projection
    Eigen::Vector3d r_sub_gcrf = r_ast - (r_ast.dot(k)) * k;
    
    OccultationParameters params;
    params.xi_ca = physics::Distance::from_m(r_sub_gcrf.dot(east));
    params.eta_ca = physics::Distance::from_m(r_sub_gcrf.dot(north));
    params.impact_parameter = physics::Distance::from_m(r_sub_gcrf.norm());
    return params;
}

void OccultationLogic::compute_shadow_velocity(
    OccultationParameters& params,
    const Declination& star_dec,
    const physics::Distance& ast_dist,
    const Angle& ast_dra_dt,
    const Angle& ast_ddec_dt)
{
    double ds = star_dec.to_rad();
    double dra_s = ast_dra_dt.to_rad() / 3600.0;
    double ddec_s = ast_ddec_dt.to_rad() / 3600.0;
    double r_val = ast_dist.to_m();
    
    params.dxi_dt = physics::Velocity::from_ms(r_val * dra_s * std::cos(ds));
    params.deta_dt = physics::Velocity::from_ms(r_val * ddec_s);
    
    double vx = params.dxi_dt.to_ms();
    double vy = params.deta_dt.to_ms();
    params.shadow_velocity = physics::Velocity::from_ms(std::sqrt(vx*vx + vy*vy));
    params.position_angle = Angle::from_rad(std::atan2(vx, vy)).wrap_0_2pi();
    params.relative_velocity_mag = params.shadow_velocity;
}

void OccultationLogic::compute_sky_conditions(
    OccultationParameters& params,
    const time::EpochTDB& t_ca,
    const RightAscension& ast_ra, const Declination& ast_dec,
    const physics::Distance& ast_distance,
    const RightAscension& star_ra, const Declination& star_dec,
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem)
{
    // Sun Altitude
    auto sun_state = ephem->getState(ephemeris::CelestialBody::SUN, t_ca);
    Eigen::Vector3d n_sun = sun_state.position.to_eigen_si().normalized();
    Eigen::Vector3d n_ast = Eigen::Vector3d(
        std::cos(ast_dec.to_rad()) * std::cos(ast_ra.to_rad()),
        std::cos(ast_dec.to_rad()) * std::sin(ast_ra.to_rad()),
        std::sin(ast_dec.to_rad())
    );
    
    params.sun_altitude = Angle::from_rad(std::asin(std::clamp(n_ast.dot(n_sun), -1.0, 1.0)));
    params.is_daylight = params.sun_altitude.to_deg() > constants::SUN_ALTITUDE_LIMIT_DEG;

    // Moon Geometry
    auto moon_ssb = ephem->getState(ephemeris::CelestialBody::MOON, t_ca);
    auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t_ca);
    Eigen::Vector3d n_moon = (moon_ssb.position.to_eigen_si() - earth_ssb.position.to_eigen_si()).normalized();
    
    params.moon_altitude = Angle::from_rad(std::asin(std::clamp(n_ast.dot(n_moon), -1.0, 1.0)));
    
    double as = star_ra.to_rad();
    double ds = star_dec.to_rad();
    Eigen::Vector3d k_star(std::cos(ds) * std::cos(as), std::cos(ds) * std::sin(as), std::sin(ds));
    params.moon_dist = Angle::from_rad(std::acos(std::clamp(k_star.dot(n_moon), -1.0, 1.0)));
}

void OccultationLogic::process_day_window(
    std::vector<OccultationCandidate>& results,
    double day_jd,
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
    const OccultationConfig& config,
    AstDynEngine& engine,
    double t_start_jd, double t_end_jd)
{
    time::EpochTDB midnight = time::EpochTDB::from_jd(day_jd);
    try {
        auto segment = catalog::fit_chebyshev(initial_elements, midnight, 1.0, engine.config());
        auto [pos, vel] = segment.evaluate_full(midnight.jd());
        double dist_au = std::get<2>(pos);
        if (dist_au < 1e-6) return;

        Angle radius = Angle::from_rad(config.max_shadow_distance.to_m() / (dist_au * constants::AU * 1000.0));
        auto stars = catalog::find_stars_near_segment(catalog::GaiaDR3Catalog::instance(), segment, radius, config.max_mag_star);
        
        for (const auto& star : stars) {
            evaluate_candidate(results, segment, star, config, engine, t_start_jd, t_end_jd);
        }
    } catch (...) { }
}

void OccultationLogic::evaluate_candidate(
    std::vector<OccultationCandidate>& results,
    const catalog::ChebyshevSegment& segment,
    const catalog::Star& star,
    const OccultationConfig& config,
    AstDynEngine& engine,
    double t_start_jd, double t_end_jd)
{
    double t_ca_jd = find_tca(segment, star, segment.t_start, segment.t_end);
    if (t_ca_jd < t_start_jd || t_ca_jd > t_end_jd) return;

    time::EpochTDB t_ca = time::EpochTDB::from_jd(t_ca_jd);
    auto [pos, vel] = segment.evaluate_full(t_ca_jd);
    auto star_at_tca = star.predict_at(t_ca);
    
    // Rates from Chebyshev (deg/day) to arcsec/hr
    double dra_hr = (std::get<0>(vel) / 24.0) * 3600.0;
    double ddec_hr = (std::get<1>(vel) / 24.0) * 3600.0;

    OccultationParameters params = compute_parameters(
        star_at_tca.ra(), star_at_tca.dec(),
        RightAscension::from_deg(std::get<0>(pos)), Declination::from_deg(std::get<1>(pos)),
        physics::Distance::from_au(std::get<2>(pos)),
        Angle::from_arcsec(dra_hr), Angle::from_arcsec(ddec_hr),
        physics::Velocity::from_au_d(std::get<2>(vel)),
        t_ca, engine.getEphemeris());

    if (params.impact_parameter > config.max_shadow_distance) return;
    if (config.filter_daylight && params.is_daylight && params.sun_altitude.to_deg() > config.min_sun_altitude) return;

    results.push_back({"", star, params}); 
}

std::vector<OccultationCandidate> OccultationLogic::find_occultations(
    const std::string& asteroid_id,
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
    time::EpochTDB start,
    time::EpochTDB end,
    const OccultationConfig& config,
    AstDynEngine& engine,
    OccultationRefinementMode mode)
{
    std::vector<OccultationCandidate> results;
    std::vector<double> days;
    for (double d = std::floor(start.jd() - 0.5) + 0.5; d <= end.jd() + 0.5; d += 1.0) {
        if (d >= start.jd() - 0.5 && d <= end.jd() + 0.5) days.push_back(d);
    }

    #pragma omp parallel
    {
        std::vector<OccultationCandidate> thread_results;
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < days.size(); ++i) {
            process_day_window(thread_results, days[i], initial_elements, config, engine, start.jd(), end.jd());
        }
        #pragma omp critical
        results.insert(results.end(), thread_results.begin(), thread_results.end());
    }

    for (auto& cand : results) cand.asteroid_id = asteroid_id;
    std::sort(results.begin(), results.end(), [](const OccultationCandidate& a, const OccultationCandidate& b) {
        return a.params.t_ca.mjd() < b.params.t_ca.mjd();
    });

    return results;
}

std::vector<OccultationSystemCandidate> OccultationLogic::find_system_occultations(
    const std::vector<std::string>& body_ids,
    const std::string& bsp_path,
    time::EpochTDB start,
    time::EpochTDB end,
    const OccultationConfig& config,
    AstDynEngine& engine)
{
    std::vector<OccultationSystemCandidate> results;
    try {
        io::SPKReader spk(bsp_path);
        for (const auto& id : body_ids) {
            int naif_id = std::stoi(id);
            for (time::EpochTDB t = start; t <= end; t += time::TimeDuration::from_seconds(3600.0)) {
                double et = (t.mjd() - 51544.5) * 86400.0;
                auto state = spk.getState(naif_id, et);
                (void)state; // Placeholder search logic
            }
        }
    } catch (...) {}
    return results;
}

std::vector<OccultationCandidate> OccultationLogic::find_multi_asteroid_occultations(
    const std::vector<std::string>& asteroid_ids,
    ChebyshevEphemerisManager& manager,
    time::EpochTDB start,
    time::EpochTDB end,
    const OccultationConfig& config,
    AstDynEngine& engine)
{
    std::vector<OccultationCandidate> results;
    for (const auto& id : asteroid_ids) {
        if (!manager.has_body(id)) continue;
        
        // Use a 1-hour step for initial search
        for (time::EpochTDB t = start; t <= end; t += time::TimeDuration::from_seconds(3600.0)) {
            auto [pos, vel] = manager.evaluate_full(id, t);
            (void)pos; (void)vel; // Placeholder usage
            // logic to check stars in the area would go here
        }
    }
    return results;
}

void OccultationLogic::apply_uncertainty(
    OccultationParameters& params,
    const catalog::Star& star,
    const Eigen::Matrix<double, 6, 6>& covariance_t0,
    const physics::CartesianStateTyped<core::GCRF>& initial_state,
    AstDynEngine& engine)
{
    using namespace orbit_determination;
    StateTransitionMatrix<core::GCRF> stm_engine(engine.propagator());
    auto stm_res = stm_engine.compute(initial_state, params.t_ca);
    Eigen::Matrix<double, 6, 6> cov_tca = stm_res.map_covariance(covariance_t0);
    
    auto state_tca = stm_res.final_state;
    auto earth_tca = engine.getEphemeris()->getState(ephemeris::CelestialBody::EARTH, params.t_ca);
    Eigen::Vector3d rho_geo = state_tca.position.to_eigen_si() - earth_tca.position.to_eigen_si();

    params.cross_track_uncertainty = physics::Distance::from_km(
        AstrometryReducer::compute_cross_track_uncertainty(
            cov_tca, 
            rho_geo,
            state_tca.velocity.to_eigen_si(),
            star.ra.to_rad(), 
            star.dec.to_rad()
        )
    );
}

} // namespace astdyn::astrometry
