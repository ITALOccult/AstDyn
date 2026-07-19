#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/ephemeris/CelestialBody.hpp"
#include <utility>
#include <array>
#include <stdexcept>
#include "astdyn/coordinates/CelestialToTerrestrial.hpp"
#include "astdyn/propagation/StateTransitionTensor.hpp"
#include "astdyn/propagation/ForceField.hpp"
#include "astdyn/astrometry/StellarCovariance.hpp"
#include "astdyn/astrometry/OccultationCovariance.hpp"
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
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
#include "astdyn/astrometry/SPKChebyshevEphemeris.hpp"
#include "astdyn/astrometry/ClosestApproachFinder.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace astdyn::astrometry {

namespace {
/// Diagnostics for the occultation search: set ASTDYN_OCC_DEBUG=1 to enable.
inline bool occ_debug() {
    static const bool on = (std::getenv("ASTDYN_OCC_DEBUG") != nullptr);
    return on;
}
} // namespace


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
    // Kept for the occelmnt export, which reports them separately.
    params.geocentric_distance = ast_dist;
    params.d_ra_arcsec_hr  = ast_dra_dt.to_arcsec();
    params.d_dec_arcsec_hr = ast_ddec_dt.to_arcsec();

    // Heliocentric distance and phase angle, for the HG apparent magnitude.
    // Both the ephemeris and (ast_ra, ast_dec) are equatorial, so the geocentric
    // vector adds to the Earth's position without any rotation.
    if (ephem) {
        const auto earth = ephem->getState(ephemeris::CelestialBody::EARTH, t_ca);
        const double ca = ast_ra.to_rad(), cd = ast_dec.to_rad();
        const Eigen::Vector3d rho = ast_dist.to_m() *
            Eigen::Vector3d(std::cos(cd) * std::cos(ca),
                            std::cos(cd) * std::sin(ca), std::sin(cd));
        const Eigen::Vector3d r_helio = earth.position.to_eigen_si() + rho;
        params.heliocentric_distance = physics::Distance::from_m(r_helio.norm());
        // Phase angle: at the object, between the directions to Sun and to Earth.
        const double c = std::clamp((-r_helio).normalized().dot((-rho).normalized()),
                                    -1.0, 1.0);
        params.phase_angle = Angle::from_rad(std::acos(c));
    }
    
    if (ephem) {
        // Order matters: the shadow centre defines the point at which the Sun's
        // altitude is evaluated, so it must be computed first.
        compute_shadow_centre(params, t_ca, star_ra, star_dec);
        compute_sky_conditions(params, t_ca, ast_ra, ast_dec, ast_dist, star_ra, star_dec, ephem);
    }
    return params;
}

void OccultationLogic::compute_shadow_centre(
    OccultationParameters& params,
    const time::EpochTDB& t_ca,
    const RightAscension& star_ra, const Declination& star_dec)
{
    // The shadow centre is where the star-asteroid axis meets the Earth. It is
    // displaced from the sub-asteroid point by the impact parameter (xi, eta),
    // i.e. by up to several thousand km for a non-central event.
    const time::EpochUTC t   = time::to_utc(t_ca);
    const time::EpochTT  tt  = time::to_tt(t);
    const time::EpochUT1 ut1 = time::to_ut1(t);

    const double a = star_ra.to_rad(), d = star_dec.to_rad();
    const auto star_dir = math::Direction<core::GCRF>::from_xyz(
        std::cos(d) * std::cos(a), std::cos(d) * std::sin(a), std::sin(d));

    if (const auto sp = coordinates::shadow_point(params.xi_ca, params.eta_ca,
                                                  star_dir, tt, ut1)) {
        params.center_lat = sp->lat.angle();
        params.center_lon = sp->lon;
    }

    // Sub-star point: where the star is at the zenith, i.e. the fundamental-plane
    // reference Occult4 reports in <Earth>. It is simply the star direction
    // rotated into ITRF, and its latitude is GEOCENTRIC -- which is why it comes
    // out exactly equal to the star's apparent declination, as Occult4 publishes.
    {
        const auto k_itrf = coordinates::gcrf_to_itrf(tt, ut1) * star_dir;
        params.substar_lat = Angle::from_rad(std::asin(std::clamp(k_itrf.z(), -1.0, 1.0)));
        params.substar_lon = Angle::from_rad(std::atan2(k_itrf.y(), k_itrf.x()));
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
    // Geocentric direction to the Sun.
    //
    // NOTE: the ephemeris is HELIOCENTRIC, so getState(SUN) returns (0,0,0) by
    // construction. Normalising that yields (0,0,0) -- Eigen returns the vector
    // unchanged rather than NaN when the norm is zero -- which silently made
    // sun_altitude exactly 0 deg, is_daylight permanently true, and the daylight
    // filter reject every single event. The Sun as seen from the Earth is minus
    // the Earth's heliocentric position.
    const auto earth_state = ephem->getState(ephemeris::CelestialBody::EARTH, t_ca);
    const Eigen::Vector3d n_sun = (-earth_state.position.to_eigen_si()).normalized();

    Eigen::Vector3d n_ast = Eigen::Vector3d(
        std::cos(ast_dec.to_rad()) * std::cos(ast_ra.to_rad()),
        std::cos(ast_dec.to_rad()) * std::sin(ast_ra.to_rad()),
        std::sin(ast_dec.to_rad())
    );

    // Altitude of the Sun above the horizon AT THE SHADOW CENTRE. The previous
    // asin(n_ast . n_sun) was not an altitude at all: that dot product is the
    // cosine of the solar elongation, and an altitude depends on where on the
    // Earth the shadow falls, not on the direction of the asteroid.
    // params.center_lat is geodetic, so the local zenith is the ellipsoid normal.
    {
        const time::EpochUTC t_utc = time::to_utc(t_ca);
        const auto C = coordinates::gcrf_to_itrf(time::to_tt(t_utc), time::to_ut1(t_utc));
        const Eigen::Vector3d sun_itrf = C.matrix() * n_sun;
        const double la = params.center_lat.to_rad(), lo = params.center_lon.to_rad();
        const Eigen::Vector3d zenith(std::cos(la) * std::cos(lo),
                                     std::cos(la) * std::sin(lo),
                                     std::sin(la));
        params.sun_altitude = Angle::from_rad(
            std::asin(std::clamp(zenith.dot(sun_itrf), -1.0, 1.0)));
    }
    params.is_daylight = params.sun_altitude.to_deg() > constants::SUN_ALTITUDE_LIMIT_DEG;

    // Sub-solar point (geocentric), for the day/night terminator.
    {
        const time::EpochUTC t_utc = time::to_utc(t_ca);
        const auto C = coordinates::gcrf_to_itrf(time::to_tt(t_utc), time::to_ut1(t_utc));
        const Eigen::Vector3d s_itrf = C.matrix() * n_sun;
        params.subsolar_lat = Angle::from_rad(std::asin(std::clamp(s_itrf.z(), -1.0, 1.0)));
        params.subsolar_lon = Angle::from_rad(std::atan2(s_itrf.y(), s_itrf.x()));
    }

    // Moon Geometry
    auto moon_ssb = ephem->getState(ephemeris::CelestialBody::MOON, t_ca);
    auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t_ca);
    Eigen::Vector3d r_moon_earth = (moon_ssb.position.to_eigen_si() - earth_ssb.position.to_eigen_si());
    Eigen::Vector3d n_moon = r_moon_earth.normalized();
    
    params.moon_altitude = Angle::from_rad(std::asin(std::clamp(n_ast.dot(n_moon), -1.0, 1.0)));
    
    double as = star_ra.to_rad();
    double ds = star_dec.to_rad();
    Eigen::Vector3d k_star(std::cos(ds) * std::cos(as), std::cos(ds) * std::sin(as), std::sin(ds));
    params.moon_dist = Angle::from_rad(std::acos(std::clamp(k_star.dot(n_moon), -1.0, 1.0)));

    // Moon Phase (approximate from elongation)
    double cos_psi = std::clamp(n_sun.dot(n_moon), -1.0, 1.0);
    params.moon_phase = 0.5 * (1.0 - cos_psi);
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
    double t_start_jd, double t_end_jd,
    double diameter_km)
{
    // Use the high-precision ClosestApproachFinder
    auto [pos_center, vel_center] = segment.evaluate_full((segment.t_start + segment.t_end) / 2.0);
    double dist_au = std::get<2>(pos_center);
    Angle max_sep_angle = Angle::from_rad((config.max_shadow_distance.to_m() * 1.5) / (dist_au * constants::AU * 1000.0));

    auto ca_results = ClosestApproachFinder::find_in_segment(
        segment, star, 
        time::EpochTDB::from_jd(std::max(segment.t_start, t_start_jd)),
        time::EpochTDB::from_jd(std::min(segment.t_end, t_end_jd)),
        max_sep_angle,
        720); // 2-minute sampling for a 1-day segment is enough for detection

    if (occ_debug()) std::cerr << "[DBG]     star " << star.source_id << " G=" << star.g_mag
        << " -> " << ca_results.size() << " CA entro " << max_sep_angle.to_arcsec() << "\"\n";

    for (const auto& ca : ca_results) {
        // Map ClosestApproachResult to OccultationParameters
        auto [pos, vel] = segment.evaluate_full(ca.t_ca.jd());
        double dist_au_at_ca = std::get<2>(pos);
        
        physics::Distance impact_dist = physics::Distance::from_m(
            ca.separation.to_rad() * dist_au_at_ca * constants::AU * 1000.0);

        if (occ_debug()) std::cerr << "[DBG]       CA jd=" << ca.t_ca.jd() << " sep=" << ca.separation.to_arcsec()
            << "\" -> impatto=" << impact_dist.to_km() << " km (max=" << config.max_shadow_distance.to_km() << ")\n";
        if (impact_dist > config.max_shadow_distance) {
            if (occ_debug()) std::cerr << "[DBG]       SCARTATO: impatto oltre max_shadow_distance\n";
            continue;
        }

        auto star_at_tca = star.predict_at(ca.t_ca);
        
        double dra_hr = (std::get<0>(vel) / 24.0) * 3600.0;
        double ddec_hr = (std::get<1>(vel) / 24.0) * 3600.0;

        OccultationParameters params = compute_parameters(
            star_at_tca.ra(), star_at_tca.dec(),
            RightAscension::from_deg(std::get<0>(pos)), Declination::from_deg(std::get<1>(pos)),
            physics::Distance::from_au(std::get<2>(pos)),
            Angle::from_arcsec(dra_hr), Angle::from_arcsec(ddec_hr),
            physics::Velocity::from_au_d(std::get<2>(vel)),
            ca.t_ca, engine.getEphemeris());

        // Maximum duration: the full shadow width traversed at the shadow speed.
        // Without this max_duration stays at zero, and the moment min_duration_s
        // is configured above 0.1 the filter below rejects every event -- the
        // same trap the daylight filter had.
        if (diameter_km > 0.0 && params.shadow_velocity.to_ms() > 0.0) {
            params.max_duration = time::TimeDuration::from_seconds(
                diameter_km * 1000.0 / params.shadow_velocity.to_ms());
        }

        if (config.filter_daylight && params.is_daylight && params.sun_altitude.to_deg() > config.min_sun_altitude) {
            if (occ_debug()) std::cerr << "[DBG]       SCARTATO: giorno (sole alt="
                << params.sun_altitude.to_deg() << " > " << config.min_sun_altitude << ")\n";
            continue;
        }

        // --- NEW: Scientific & Proximity Filters ---

        // 1. Duration filter
        if (config.min_duration_s > 0.1 && params.max_duration.to_seconds() < config.min_duration_s) {
            if (occ_debug()) std::cerr << "[DBG]       SCARTATO: durata " << params.max_duration.to_seconds()
                << "s < " << config.min_duration_s << "s\n";
            continue;
        }

        // 2. Gaia Quality Filter (RUWE)
        if (star.ruwe > config.max_gaia_ruwe) {
            if (occ_debug()) std::cerr << "[DBG]       SCARTATO: RUWE " << star.ruwe
                << " > " << config.max_gaia_ruwe << "\n";
            continue;
        }

        // 3. Moon Filters
        if (params.moon_dist.to_deg() < config.min_moon_dist) {
            if (occ_debug()) std::cerr << "[DBG]       SCARTATO: Luna a " << params.moon_dist.to_deg()
                << " deg < " << config.min_moon_dist << "\n";
            continue;
        }
        if (params.moon_phase > config.max_moon_phase + 0.001) {
            if (occ_debug()) std::cerr << "[DBG]       SCARTATO: fase Luna " << params.moon_phase
                << " > " << config.max_moon_phase << "\n";
            continue;
        }

        // 4. Proximity filter (Great Circle Distance from center_lat/lon to observer)
        if (config.max_obs_dist_km > 1.0) {
            double phi1 = config.obs_lat * constants::DEG_TO_RAD;
            double phi2 = params.center_lat.to_rad();
            double dlam = (params.center_lon.to_deg() - config.obs_lon) * constants::DEG_TO_RAD;
            
            double d_rad = std::acos(std::sin(phi1)*std::sin(phi2) + std::cos(phi1)*std::cos(phi2)*std::cos(dlam));
            double d_km = d_rad * 6371.0; 
            
            if (d_km > config.max_obs_dist_km) continue;
        }

        results.push_back({"", star, params}); 
    }
}

void OccultationLogic::evaluate_candidate(
    const std::string& id,
    const catalog::Star& star,
    const catalog::ChebyshevSegment& segment,
    const std::vector<ClosestApproachResult>& candidates,
    std::vector<OccultationCandidate>& results,
    const OccultationConfig& config,
    AstDynEngine& engine)
{
    for (const auto& ca : candidates) {
        auto star_at_tca = star.predict_at(ca.t_ca);
        
        auto [pos, vel] = segment.evaluate_full(ca.t_ca.jd());
        double dist_au_at_ca = std::get<2>(pos);
        
        physics::Distance impact_dist = physics::Distance::from_m(
            ca.separation.to_rad() * dist_au_at_ca * constants::AU * 1000.0);

        if (impact_dist > config.max_shadow_distance) continue;

        double dra_hr = (std::get<0>(vel) / 24.0) * 3600.0;
        double ddec_hr = (std::get<1>(vel) / 24.0) * 3600.0;

        OccultationParameters params = compute_parameters(
            star_at_tca.ra(), star_at_tca.dec(),
            RightAscension::from_deg(std::get<0>(pos)), Declination::from_deg(std::get<1>(pos)),
            physics::Distance::from_au(std::get<2>(pos)),
            Angle::from_arcsec(dra_hr), Angle::from_arcsec(ddec_hr),
            physics::Velocity::from_au_d(std::get<2>(vel)),
            ca.t_ca, engine.getEphemeris());

        if (config.filter_daylight && params.is_daylight && params.sun_altitude.to_deg() > config.min_sun_altitude) {
            continue;
        }

        // --- Scientific Quality Filters ---

        // 1. Duration filter
        if (config.min_duration_s > 0.1 && params.max_duration.to_seconds() < config.min_duration_s) {
            continue;
        }

        // 2. Gaia Quality Filter (RUWE)
        if (star.ruwe > config.max_gaia_ruwe) {
            continue;
        }

        // 3. Moon Filters
        if (params.moon_dist.to_deg() < config.min_moon_dist) {
            continue;
        }
        if (params.moon_phase > config.max_moon_phase + 0.001) {
            continue;
        }

        // 4. Proximity filter
        if (config.max_obs_dist_km > 1.0) {
            double phi1 = config.obs_lat * constants::DEG_TO_RAD;
            double phi2 = params.center_lat.to_rad();
            double dlam = (params.center_lon.to_deg() - config.obs_lon) * constants::DEG_TO_RAD;
            
            double d_rad = std::acos(std::sin(phi1)*std::sin(phi2) + std::cos(phi1)*std::cos(phi2)*std::cos(dlam));
            double d_km = d_rad * 6371.0; 
            
            if (d_km > config.max_obs_dist_km) continue;
        }

        results.push_back({id, star, params}); 
    }
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
    std::vector<double> segment_starts;
    for (double jd = start.jd(); jd < end.jd(); jd += 1.0) {
        segment_starts.push_back(jd);
    }

    #pragma omp parallel
    {
        std::vector<OccultationCandidate> thread_results;
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < segment_starts.size(); ++i) {
            try {
                double jd = segment_starts[i];
                double segment_duration = std::min(jd + 1.0, end.jd()) - jd;
                time::EpochTDB t_center = time::EpochTDB::from_jd(jd + segment_duration / 2.0);
                
                auto segment = catalog::fit_chebyshev(initial_elements, t_center, segment_duration, engine.config());
                
                auto [pos, vel] = segment.evaluate_full(jd + segment_duration / 2.0);
                double dist_au = std::get<2>(pos);
                if (dist_au < 1e-6) continue;

                Angle radius = Angle::from_rad((config.max_shadow_distance.to_m() * 1.5) / (dist_au * constants::AU * 1000.0));
                auto stars = catalog::find_stars_near_segment(catalog::GaiaDR3Catalog::instance(), segment, radius, config.max_mag_star);
                
                for (const auto& star : stars) {
                    evaluate_candidate(thread_results, segment, star, config, engine, start.jd(), end.jd());
                    for (auto& cand : thread_results) {
                        if (cand.asteroid_id.empty()) cand.asteroid_id = asteroid_id;
                    }
                }
            } catch (...) {
                // Skip segment
            }
        }
        #pragma omp critical
        results.insert(results.end(), thread_results.begin(), thread_results.end());
    }

    // The asteroid_id is already assigned within the loop for each candidate
    // for (auto& cand : results) cand.asteroid_id = asteroid_id; 
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
        // Map to group occultations by the background star, as multiple bodies (e.g. Jovian moons)
        // might occult the same star during the same period.
        std::map<unsigned long long, OccultationSystemCandidate> star_map;

        for (const auto& id : body_ids) {
            int naif_id = std::stoi(id);
            // Convert SPK orbital data into daily Chebyshev segments for high-speed geometric search.
            SPKChebyshevEphemeris eph(naif_id, spk, start, end);
            
            // Use a more robust way to iterate over the time window: 
            // Iterate through all 1-day segments of the asteroid's ephemeris
            // and find stars for each.
            double jd = start.jd();
            while (jd < end.jd()) {
                try {
                    time::EpochTDB t = time::EpochTDB::from_jd(jd);
                    const auto& segment = eph.get_segment(t);
                    
                    // Evaluate path for the current period
                    auto [pos, vel] = segment.evaluate_full(std::clamp(jd, segment.t_start, segment.t_end));
                    double dist_au = std::get<2>(pos);
                    if (dist_au < 1e-6) {
                        jd = segment.t_end;
                        continue;
                    }

                    // Compute search radius in the sky based on maximum shadow distance (FP)
                    Angle radius = Angle::from_rad(config.max_shadow_distance.to_m() / (dist_au * constants::AU * 1000.0));
                    
                    // Bulk query against the Gaia database for stars near the body path
                    auto stars = catalog::find_stars_near_segment(catalog::GaiaDR3Catalog::instance(), segment, radius, config.max_mag_star);
                    
                    for (const auto& star : stars) {
                        std::vector<OccultationCandidate> temp_results;
                        // Refine the candidate to find the exact TCA and fundamental plane parameters
                        evaluate_candidate(temp_results, segment, star, config, engine, start.jd(), end.jd());
                        
                        for (const auto& cand : temp_results) {
                            // If this star wasn't seen before, initialize a new system candidate
                            if (star_map.find(star.source_id) == star_map.end()) {
                                star_map[star.source_id].star = star;
                            }
                            // Store the specific body interaction
                            BodyOccultation body;
                            body.name = id;
                            body.params = cand.params;
                            body.diameter = physics::Distance::from_km(100.0); // Default placeholder
                            star_map[star.source_id].bodies.push_back(body);
                        }
                    }
                    
                    // Advance to next segment
                    jd = segment.t_end;
                } catch (...) {
                    jd += 1.0; // Fallback
                }
            }
        }
        // Collect all unique star candidates found
        for (auto& pair : star_map) {
            results.push_back(pair.second);
        }
    } catch (...) {
        // Log or handle SPK load errors
    }
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
    // Iterate over the list of managed asteroids (pre-calculated with polynomials)
    for (const auto& id : asteroid_ids) {
        if (!manager.has_body(id)) {
            if (occ_debug()) std::cerr << "[DBG] '" << id
                << "' NON e' nel manager -> Horizons non ha restituito elementi\n";
            continue;
        }
        if (occ_debug()) std::cerr << "[DBG] '" << id << "' diametro=" << manager.get_diameter(id)
            << " km (filtro min=" << config.min_asteroid_diameter_km << ")\n";

        // --- NEW: Diameter Filter ---
        if (config.min_asteroid_diameter_km > 0.1 && manager.get_diameter(id) < config.min_asteroid_diameter_km) {
            if (occ_debug()) std::cerr << "[DBG]   SCARTATO: diametro sotto soglia\n";
            continue; 
        }
        
        const auto& eph = manager.get_ephemeris(id);
        
        // Use a more robust way to iterate over the time window: 
        // Iterate through all 1-day segments of the asteroid's ephemeris
        // and find stars for each.
        double jd = start.jd();
        while (jd < end.jd()) {
            try {
                time::EpochTDB t = time::EpochTDB::from_jd(jd);
                const auto& segment = eph.get_segment(t);
                
                // Identify stars near this segment
                auto [pos, vel] = segment.evaluate_full(std::clamp(jd, segment.t_start, segment.t_end));
                double dist_au = std::get<2>(pos);
                if (dist_au < 1e-6) {
                    jd = segment.t_end;
                    continue;
                }

                Angle radius = Angle::from_rad((config.max_shadow_distance.to_m() * 1.5) / (dist_au * constants::AU * 1000.0));

                if (occ_debug()) std::cerr << "[DBG]   segmento [" << segment.t_start << ", " << segment.t_end
                    << "] dist=" << dist_au << " AU  raggio_cono=" << radius.to_arcsec()
                    << "\"  mag_max=" << config.max_mag_star << "\n";

                auto stars = catalog::find_stars_near_segment(catalog::GaiaDR3Catalog::instance(), segment, radius, config.max_mag_star);

                if (occ_debug()) std::cerr << "[DBG]   stelle restituite dal catalogo: " << stars.size() << "\n";

                for (const auto& star : stars) {
                    evaluate_candidate(results, segment, star, config, engine,
                                       start.jd(), end.jd(), manager.get_diameter(id));
                    // Assign proper ID to results
                    if (!results.empty() && results.back().asteroid_id.empty()) {
                        results.back().asteroid_id = id;
                    }
                }
                
                // Advance to next segment
                jd = segment.t_end;
            } catch (const std::exception& e) {
                std::cerr << "[DBG] ECCEZIONE @ jd=" << jd << " : " << e.what() << "\n";
                jd += 1.0; // Fallback
            } catch (...) {
                std::cerr << "[DBG] ECCEZIONE SCONOSCIUTA @ jd=" << jd << "\n";
                jd += 1.0; // Fallback
            }
        }
    }
    return results;
}

double hg_magnitude(double h_mag, double g, double r_au, double d_au, Angle alpha) {
    if (h_mag <= 0.0 || r_au <= 0.0 || d_au <= 0.0) return -5.0;   // Occult4's sentinel
    const double t = std::tan(0.5 * alpha.to_rad());
    if (t < 0.0) return -5.0;
    const double phi1 = std::exp(-3.33 * std::pow(t, 0.63));
    const double phi2 = std::exp(-1.87 * std::pow(t, 1.22));
    const double f = (1.0 - g) * phi1 + g * phi2;
    if (f <= 0.0) return -5.0;
    return h_mag + 5.0 * std::log10(r_au * d_au) - 2.5 * std::log10(f);
}

void OccultationLogic::apply_uncertainty(
    OccultationParameters& params,
    const catalog::Star& star,
    const Eigen::Matrix<double, 6, 6>& covariance_t0,
    const physics::CartesianStateTyped<core::ECLIPJ2000>& initial_state,
    AstDynEngine& engine)
{
    // ---- 1. Transport: Phi AND Psi ------------------------------------------
    // The state transition MATRIX gives only Phi, which is the whole of the
    // linear theory. The tensor gives Psi as well, and Psi is what makes the
    // bias, the skewness and the index N computable at all.
    //
    // The perturbers are not optional here. Over the few weeks separating an
    // AstDyS epoch from a main-belt event, Kepler + J2 is adequate; over the
    // months separating it from a NEO event it is not, and the resulting error
    // lands squarely on the quantity being validated. The tensor already exposes
    // a per-step hook for exactly this.
    // Force model unificato: ForceField in ECLITTICA (integrate_in_ecliptic),
    // stessa fisica del vecchio PotentialModel di questa sezione: J2 solare OFF
    // (era model.j2=0), relativita' OFF, asteroidi OFF (solo perturbatori
    // planetari maggiori). Perturbatori e rotazione eq->ecl sono ora interni a
    // ForceField::n_body_perturbation; EphemerisRefresh e ASTDYN_PERTURBERS via.
    propagation::PropagatorSettings stt_settings = engine.config().propagator_settings;
    stt_settings.integrate_in_ecliptic   = true;
    stt_settings.include_planets         = true;
    stt_settings.include_moon            = true;
    stt_settings.include_sun_j2          = false;
    stt_settings.include_relativity      = false;
    stt_settings.include_asteroids       = false;
    stt_settings.baricentric_integration = false;
    auto stt_force = std::make_shared<propagation::ForceField>(
        stt_settings, engine.getEphemeris());

    // Rotazione eq->ecl per la geometria d'ombra a valle (osservatore, stella).
    const Matrix3d eq2ecl_frame = coordinates::ReferenceFrame::j2000_to_ecliptic();

    // Precisione dell'integratore come manopola misurabile (metrica di passo
    // AAS via ForceField::gradient_spectral_radius).
    double precision = 1e-4;
    if (const char* p = std::getenv("ASTDYN_STT_PRECISION")) precision = std::atof(p);

    propagation::StateTransitionTensor stt(stt_force, precision);
    const auto tr = stt.propagate(initial_state, params.t_ca);

    // Il numero grezzo, senza geometria di mezzo. Un test chiuso (STT contro
    // Keplero analitico, senza perturbatori) dava 0.1 km su 263 giorni, mentre
    // qui l'asse d'ombra sbaglia di decine di migliaia: uno dei due mente, e
    // l'unico modo di sapere quale e' confrontare la posizione propagata con
    // Horizons alla STESSA epoca, invece di dedurla dalla geometria a valle.
    if (occ_debug()) {
        const auto xi = initial_state.to_eigen_au_aud();
        const auto xf = tr.final_state.to_eigen_au_aud();
        std::cerr << "[scope] propagazione: da MJD(TDB) " << initial_state.epoch.mjd()
                  << " a " << params.t_ca.mjd() << "  (dt = "
                  << (params.t_ca.mjd() - initial_state.epoch.mjd()) << " d)\n";
        std::cerr << "[scope]   r0 = (" << xi(0) << ", " << xi(1) << ", " << xi(2)
                  << ") AU   |r0| = " << xi.head<3>().norm() << "\n";
        std::cerr << "[scope]   r1 = (" << xf(0) << ", " << xf(1) << ", " << xf(2)
                  << ") AU   |r1| = " << xf.head<3>().norm() << "\n";
        std::cerr << "[scope]   confronta r1 con Horizons a JD "
                  << (params.t_ca.mjd() + 2400000.5) << " (@sun, ecliptic)\n";
    }

    // ---- 2. Geometry at the event, in AU, ecliptic ---------------------------
    // PlanetaryEphemeris returns J2000 EQUATORIAL states -- its return type says
    // so, math::Vector3<core::GCRF, ...> -- while the tensor propagates in the
    // ecliptic. Calling to_eigen_au_aud() strips the tag, and once the numbers
    // are bare Eigen the compiler cannot help: subtracting an equatorial vector
    // from an ecliptic one is a perfectly valid subtraction of six doubles, and
    // it silently scaled the whole error ellipse by a factor of 1.6.
    const auto x_ast = tr.final_state.to_eigen_au_aud();
    const auto earth = engine.getEphemeris()->getState(ephemeris::CelestialBody::EARTH,
                                                       params.t_ca);
    const Vector3d r_ast = x_ast.head<3>();
    const Vector3d r_obs = eq2ecl_frame * earth.to_eigen_au_aud().head<3>();

    // ---- 2b. La STT concorda con l'effemeride? -------------------------------
    // La posizione dell'ombra in params viene dalla Chebyshev, cioe' da HORIZONS;
    // la covarianza qui sotto e' propagata dallo stato ASTDYS. Sono due soluzioni
    // orbitali diverse, e per 820987 differiscono di 37 km: un residuo misurato
    // contro l'una e normalizzato con l'altra porta dentro anche quel disaccordo.
    //
    // Il confronto va fatto nel PIANO FONDAMENTALE, non sulla distanza: l'errore
    // di una propagazione e' quasi tutto ALONG-TRACK, e l'along-track lascia la
    // distanza geocentrica quasi invariata. Confrontare le distanze e' cieco
    // proprio alla componente che conta.
    //
    // Questo blocco riporta il numero. Su main-belt a poche settimane e' piccolo;
    // su un NEO propagato mesi senza perturbatori planetari non lo sara'.

    // ---- 3. Plane-of-sky basis ----------------------------------------------
    // The star direction is equatorial, the dynamics is ecliptic, so the line of
    // sight must be rotated. Rotating the EQUATORIAL POLE too, and using it as
    // the roll reference, makes the resulting basis exactly (alpha*, delta) --
    // expressed in ecliptic components. That is what lets C_star, which Gaia
    // publishes in (alpha*, delta), be added without any further rotation.
    const auto star_ep = star.predict_at(params.t_ca);
    const double sa = star_ep.ra().to_rad(), sd = star_ep.dec().to_rad();
    const Vector3d s_eq(std::cos(sd) * std::cos(sa),
                        std::cos(sd) * std::sin(sa),
                        std::sin(sd));
    const Vector3d s_hat    = eq2ecl_frame * s_eq;
    const Vector3d pole_ecl = eq2ecl_frame * Vector3d::UnitZ();
    const auto basis = tangent_basis(s_hat, pole_ecl);

    // Asse d'ombra secondo la STT, nel piano fondamentale: la retta per
    // l'asteroide parallela alla direzione della stella, intersecata col piano
    // per il geocentro normale a s_hat.
    if (occ_debug() || true) {
        const Vector3d rho = r_ast - r_obs;
        const Vector3d axis = rho - rho.dot(s_hat) * s_hat;
        const double xi_stt  = axis.dot(basis[0]) * constants::AU;
        const double eta_stt = axis.dot(basis[1]) * constants::AU;
        const double dxi  = xi_stt  - params.xi_ca.to_km();
        const double deta = eta_stt - params.eta_ca.to_km();
        const double gap = std::hypot(dxi, deta);
        if (gap > 1.0) {
            std::cerr << "[scope] STT(AstDyS) vs effemeride(Horizons) nel piano fondamentale: "
                      << gap << " km   (dxi=" << dxi << " deta=" << deta << ")\n";
        }
    }

    // ---- Tempo-luce ---------------------------------------------------------
    // r_ast e' la posizione GEOMETRICA a t_ca: dove l'oggetto E'. Ma params
    // porta la posizione ASTROMETRICA, dedotta da un'effemeride osservata: dove
    // lo VEDI. Sono due cose diverse, e la differenza e' quanto l'asteroide
    // percorre mentre la sua luce arriva:
    //
    //   Phaethon a 36.6 km/s, tempo di volo 398 s (0.80 AU) -> 14.590 km
    //
    // Non e' un errore di propagazione, e si riconosce: non cresce con dt, e'
    // insensibile alla precisione dell'integratore, nessun perturbatore lo fa
    // crollare, e scala con la DISTANZA invece che col tempo.
    //
    // La posizione apparente e' quella all'istante di EMISSIONE:
    //   t_em = t_ca - |r(t_em) - r_obs| / c
    // Si risolve per iterazione. Qui r(t_em) e' approssimato linearmente con la
    // velocita' a t_ca: su ~400 s l'errore residuo e' sotto il chilometro,
    // ampiamente sufficiente per un confronto diagnostico.
    //
    // NOTA per la covarianza: l'ellisse NON ne risente. La sensibilita' della
    // posizione apparente allo stato iniziale contiene il termine
    // v * dt_em/dx0 ~ v/c ~ 1e-4, trascurabile su Phi. Era il metro a mentire,
    // non la covarianza.
    const Vector3d v_ast = x_ast.tail<3>();          // AU/giorno
    Vector3d r_app = r_ast;
    double lt_days = 0.0;
    for (int it = 0; it < 4; ++it) {
        lt_days = (r_app - r_obs).norm() / constants::SPEED_OF_LIGHT_AU_PER_DAY;
        r_app = r_ast - v_ast * lt_days;
    }

    // L'asse d'ombra secondo la STT, ora omogeneo con params.
    {
        const Vector3d rho  = r_app - r_obs;
        const Vector3d axis = rho - rho.dot(s_hat) * s_hat;
        const double xi_stt  = axis.dot(basis[0]) * constants::AU;   // km
        const double eta_stt = axis.dot(basis[1]) * constants::AU;
        const double gap = std::hypot(xi_stt - params.xi_ca.to_km(),
                                      eta_stt - params.eta_ca.to_km());
        if (occ_debug() || gap > 1.0) {
            const Vector3d ax_geo = (r_ast - r_obs)
                                  - (r_ast - r_obs).dot(s_hat) * s_hat;
            const double gap_geo = std::hypot(
                ax_geo.dot(basis[0]) * constants::AU - params.xi_ca.to_km(),
                ax_geo.dot(basis[1]) * constants::AU - params.eta_ca.to_km());
            std::cerr << "[scope] tempo-luce " << lt_days * 86400.0 << " s ("
                      << (r_app - r_obs).norm() << " AU)   |v| = "
                      << v_ast.norm() * constants::AU / 86400.0 << " km/s\n";
            std::cerr << "[scope] asse d'ombra: STT/AstDyS (" << xi_stt << ", "
                      << eta_stt << ") km  vs  effemeride ("
                      << params.xi_ca.to_km() << ", " << params.eta_ca.to_km()
                      << ") km  ->  scarto " << gap << " km"
                      << "   (senza correzione: " << gap_geo << " km)\n";
        }
    }

    // ---- 4. Composite map and moments ---------------------------------------
    const auto pm = projection_maps(r_ast, basis, r_obs);
    const auto cm = composite_map(pm, tr.phi, tr.psi);
    const auto m  = moments(cm, covariance_t0);

    // ---- 5. The star's own uncertainty (Eq. 17) ------------------------------
    // Independent of the orbit, so the plane-of-sky covariances simply add.
    Eigen::Matrix2d C_total = m.C_xi;
    if (star.astrometric_params_solved >= 5) {
        Eigen::Matrix<double, 5, 1> serr;
        serr << star.ra_error_mas, star.dec_error_mas, star.parallax_error_mas,
                star.pmra_error_mas_yr, star.pmdec_error_mas_yr;
        // All ten correlations, not a diagonal approximation: Gaia publishes
        // them because the ellipse is wrong without them.
        GaiaCorrelations gc;
        gc.ra_dec         = star.ra_dec_corr;
        gc.ra_parallax    = star.ra_parallax_corr;
        gc.ra_pmra        = star.ra_pmra_corr;
        gc.ra_pmdec       = star.ra_pmdec_corr;
        gc.dec_parallax   = star.dec_parallax_corr;
        gc.dec_pmra       = star.dec_pmra_corr;
        gc.dec_pmdec      = star.dec_pmdec_corr;
        gc.parallax_pmra  = star.parallax_pmra_corr;
        gc.parallax_pmdec = star.parallax_pmdec_corr;
        gc.pmra_pmdec     = star.pmra_pmdec_corr;

        // Julian year of the event, against Gaia's 2016.0 reference epoch. The
        // dt term is not a detail: ten years of proper-motion uncertainty
        // typically dwarfs the catalogue position error.
        const double dt_yr = dt_from_epoch(
            2000.0 + (params.t_ca.mjd() - constants::MJD2000) / 365.25);
        C_total += stellar_covariance(build_c5(serr, gc), dt_yr);
    }
    // A 2-parameter Gaia solution has no proper motion, hence no C_star: the
    // star's contribution is then unknown, not zero. Left out deliberately
    // rather than silently approximated.

    // ---- 6. Cross-track direction, in the plane of sky ------------------------
    // The shadow moves along (dxi/dt, deta/dt); cross-track is perpendicular to
    // it, and it is the only direction where the uncertainty displaces the path.
    const double vx = params.dxi_dt.to_ms(), vy = params.deta_dt.to_ms();
    const double vn = std::hypot(vx, vy);
    const Eigen::Vector2d n_hat = (vn > 0.0) ? Eigen::Vector2d(-vy / vn, vx / vn)
                                             : Eigen::Vector2d(1.0, 0.0);

    // ---- 7. Fill the parameters ---------------------------------------------
    const double rho_au = (r_ast - r_obs).norm();
    const double au_km  = constants::AU;

    // Cross-track sigma: angle -> distance at the object.
    const double sig_n_rad = std::sqrt(n_hat.transpose() * C_total * n_hat);
    params.cross_track_uncertainty = physics::Distance::from_km(sig_n_rad * rho_au * au_km);

    // Second-order bias along cross-track (Eq. 14). Identically zero in the
    // linear theory, so it is worth reporting on its own.
    params.bias_cross_track = physics::Distance::from_km(
        n_hat.dot(m.bias) * rho_au * au_km);

    // 1-sigma error ellipse.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(C_total);
    if (es.info() == Eigen::Success) {
        const double l0 = std::max(es.eigenvalues()(0), 0.0);
        const double l1 = std::max(es.eigenvalues()(1), 0.0);
        params.err_major = Angle::from_rad(std::sqrt(l1));
        params.err_minor = Angle::from_rad(std::sqrt(l0));
        // PA measured from north (delta) towards east (alpha*).
        const Eigen::Vector2d v = es.eigenvectors().col(1);
        params.err_pa = Angle::from_rad(std::atan2(v(0), v(1)));
    }

    // The index: available a priori, no Monte Carlo.
    params.nonlinearity_index = nonlinearity_index(cm, covariance_t0, n_hat);
}

} // namespace astdyn::astrometry
