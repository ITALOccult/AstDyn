/**
 * @file OccultationLogic.cpp
 * @brief Implementation of occultation physical parameters computation.
 */

#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/io/MPCClient.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/io/SPKReader.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/astrometry/ChebyshevEphemerisManager.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace astdyn::astrometry {

OccultationParameters OccultationLogic::compute_parameters(
    const RightAscension& star_ra, const Declination& star_dec,
    const RightAscension& ast_ra, const Declination& ast_dec,
    const physics::Distance& ast_dist,
    const Angle& ast_dra_dt, const Angle& ast_ddec_dt,
    const physics::Velocity& ast_ddist_dt,
    const time::EpochTDB& t_ca) 
{
    using namespace constants;

    double as = star_ra.to_rad();
    double ds = star_dec.to_rad();
    double a = ast_ra.to_rad();
    double d = ast_dec.to_rad();
    double rho = ast_dist.to_m();
    double ast_dra_dt_rad_s = ast_dra_dt.to_rad();
    double ast_ddec_dt_rad_s = ast_ddec_dt.to_rad();
    double ast_ddist_dt_m_s = ast_ddist_dt.to_ms();

    // 1. Besselian Basis (toward star)
    Eigen::Vector3d k(std::cos(as) * std::cos(ds), std::sin(as) * std::cos(ds), std::sin(ds));
    Eigen::Vector3d i(-std::sin(as), std::cos(as), 0.0);
    // Rigorous North vector formula: j = k x i
    Eigen::Vector3d j(-std::sin(ds) * std::cos(as), -std::sin(ds) * std::sin(as), std::cos(ds));

    // 2. Asteroid and its velocity in 3D (GCRF)
    Eigen::Vector3d r_ast_vec(
        rho * std::cos(a) * std::cos(d),
        rho * std::sin(a) * std::cos(d),
        rho * std::sin(d)
    );

    // Derivatives of r_ast vector
    double cos_a = std::cos(a);
    double sin_a = std::sin(a);
    double cos_d = std::cos(d);
    double sin_d = std::sin(d);

    Eigen::Vector3d v_ast_vec(
        ast_ddist_dt_m_s * cos_a * cos_d - rho * sin_a * cos_d * ast_dra_dt_rad_s - rho * cos_a * sin_d * ast_ddec_dt_rad_s,
        ast_ddist_dt_m_s * sin_a * cos_d + rho * cos_a * cos_d * ast_dra_dt_rad_s - rho * sin_a * sin_d * ast_ddec_dt_rad_s,
        ast_ddist_dt_m_s * sin_d + rho * cos_d * ast_ddec_dt_rad_s
    );

    // 3. Project to Fundamental Plane (xi, eta)
    double xi = r_ast_vec.dot(i);
    double eta = r_ast_vec.dot(j);
    double dxi = v_ast_vec.dot(i);
    double deta = v_ast_vec.dot(j);

    // 4. Closest Approach Analysis (Linear approximation on the plane)
    // Distance squared: f(t) = (xi + dxi*t)^2 + (eta + deta*t)^2
    // f'(t) = 2(xi + dxi*t)dxi + 2(eta + deta*t)deta = 0
    // xi*dxi + dxi^2*t + eta*deta + deta^2*t = 0
    // t_min = -(xi*dxi + eta*deta) / (dxi^2 + deta^2)
    
    double v2 = dxi * dxi + deta * deta;
    double t_ca_offset = 0.0;
    if (v2 > 1e-18) {
        t_ca_offset = -(xi * dxi + eta * deta) / v2;
    }

    double xi_ca_val = xi + dxi * t_ca_offset;
    double eta_ca_val = eta + deta * t_ca_offset;
    double b = std::sqrt(xi_ca_val * xi_ca_val + eta_ca_val * eta_ca_val);

    // Build Result - Preserve signs (Bug B Fix)
    OccultationParameters params;
    params.xi_ca = physics::Distance::from_m(xi_ca_val);
    params.eta_ca = physics::Distance::from_m(eta_ca_val);
    params.impact_parameter = physics::Distance::from_m(b);
    params.shadow_velocity = physics::Velocity::from_ms(std::sqrt(v2));
    params.dxi_dt = physics::Velocity::from_ms(dxi);
    params.deta_dt = physics::Velocity::from_ms(deta);
    params.t_ca = t_ca + time::TimeDuration::from_seconds(t_ca_offset);
    
    // Position angle of the track (direction of velocity on plane)
    // atan2(E, N)
    params.position_angle = Angle::from_rad(std::atan2(dxi, deta)).wrap_0_2pi();

    // Angular velocities (relative velocity components)
    // We use RA velocity scaled by cos(Dec) for the East-West component
    params.d_ra_cos_dec_per_sec = Angle::from_rad(ast_dra_dt_rad_s * std::cos(d));
    params.d_dec_per_sec = Angle::from_rad(ast_ddec_dt_rad_s);
    
    // Total apparent rate in arcsec/hr
    double total_rate_rad_s = std::sqrt(ast_dra_dt_rad_s * std::cos(d) * ast_dra_dt_rad_s * std::cos(d) + ast_ddec_dt_rad_s * ast_ddec_dt_rad_s);
    params.total_apparent_rate = total_rate_rad_s * 206265.0 * 3600.0;
    
    params.closest_approach_time_offset = time::TimeDuration::from_seconds(t_ca_offset);
    params.time_uncertainty = time::TimeDuration::zero();
    params.cross_track_uncertainty = physics::Distance::from_km(40.0); // Default for visualization

    // 5. Central Point & Duration Calculation
    // Sub-asteroid point on Earth (Central Point)
    // We assume spherical Earth for target point approximation
    double R_earth = 6371000.0;
    
    // Position of shadow axis at TCA (in FP basis i, j, k)
    // The "central point" is the point on Earth surface closest to the shadow axis.
    // Shadow axis at TCA is xi_ca*i + eta_ca*j + u*k
    // On sphere: r^2 = xi_ca^2 + eta_ca^2 + u^2 = R_earth^2
    // u = sqrt(R_earth^2 - b^2)
    
    Eigen::Vector3d r_sub_gcrf;
    if (b < R_earth) {
        double u = std::sqrt(R_earth * R_earth - b * b);
        r_sub_gcrf = xi_ca_val * i + eta_ca_val * j + u * k;
    } else {
        // Shadow misses Earth center, use closest point on surface
        r_sub_gcrf = (R_earth / b) * (xi_ca_val * i + eta_ca_val * j);
    }
    
    // Convert GCRF to Geodetic (Longitude/Latitude)
    // Simple GMST approximation: GMST = 280.4606 + 360.9856*(MJD - 51544.5)
    // We use a simplified Greenwich Sidereal Time calculation
    double mjd_utc = params.t_ca.mjd(); 
    double gmst_deg = 280.46061837 + 360.98564736629 * (mjd_utc - 51544.5);
    double gmst_rad = std::fmod(gmst_deg * M_PI / 180.0, 2.0 * M_PI);
    
    double x = r_sub_gcrf.x();
    double y = r_sub_gcrf.y();
    double z = r_sub_gcrf.z();
    
    double lon_rad = std::atan2(y, x) - gmst_rad;
    double lat_rad = std::asin(z / R_earth);
    
    params.center_lon = Angle::from_rad(lon_rad).wrap_pi();
    params.center_lat = Angle::from_rad(lat_rad);

    // Max Duration (Placeholder diameter: 100km if not known)
    double diam_km = 100.0; // Default placeholder
    double v_shadow = params.shadow_velocity.to_ms();
    if (v_shadow > 0.01) {
        params.max_duration = time::TimeDuration::from_seconds(diam_km * 1000.0 / v_shadow);
    } else {
        params.max_duration = time::TimeDuration::zero();
    }
    
    // 6. Observability Filters
    Eigen::Vector3d n_sub = r_sub_gcrf.normalized();

    // 6.1 SUN (Altitude at sub-asteroid point)
    auto sun_state = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::SUN, params.t_ca);
    Eigen::Vector3d n_sun = sun_state.position.to_eigen_si().normalized();
    double sun_alt_rad = std::asin(std::clamp(n_sub.dot(n_sun), -1.0, 1.0));
    params.sun_altitude = sun_alt_rad * 180.0 / M_PI;
    params.is_daylight = (params.sun_altitude > -0.83); 

    // 6.2 MOON
    auto moon_ssb = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::MOON, params.t_ca);
    auto earth_ssb = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, params.t_ca);
    Eigen::Vector3d n_moon = (moon_ssb.position.to_eigen_si() - earth_ssb.position.to_eigen_si()).normalized();
    double moon_alt_rad = std::asin(std::clamp(n_sub.dot(n_moon), -1.0, 1.0));
    params.moon_altitude = moon_alt_rad * 180.0 / M_PI;

    // 6.3 Moon Distance (from Star)
    params.moon_dist = std::acos(std::clamp(k.dot(n_moon), -1.0, 1.0)) * 180.0 / M_PI;

    // 6.4 Magnitude Drop (Simplified: m_comb - m_ast)
    double flux_ast_proxy = std::pow(10.0, -0.4 * 12.0); // Assume H=12
    double flux_star = std::pow(10.0, -0.4 * params.star_mag);
    params.mag_drop = -2.5 * std::log10(flux_star / (flux_star + flux_ast_proxy));

    return params;
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
    using namespace catalog;
    std::vector<OccultationCandidate> results;
    
    // 1. Initial State Initialization & Automatic Refinement
    physics::KeplerianStateTyped<core::ECLIPJ2000> working_elements = initial_elements;
    
    // Mode (b): Orbit Fitting with automatic download
    if (mode == OccultationRefinementMode::OrbitFitting) {
        io::MPCClient mpc;
        auto obs_res = mpc.fetch_observations(asteroid_id);
        if (obs_res) {
            engine.clear_observations();
            for (const auto& obs : *obs_res) engine.add_observation(obs);
            engine.set_initial_orbit(working_elements);
            auto fit_res = engine.fit_orbit();
            if (fit_res.converged) {
                working_elements = fit_res.orbit;
            }
        }
    }

    // 2. Daily Chebyshev Loop (Discovery Phase)
    double t_start_jd = start.jd();
    double t_end_jd = end.jd();
    auto& catalog_inst = GaiaDR3Catalog::instance(); 

    // Iterate through days, centered at midnight UTC
    for (double day_jd = std::floor(t_start_jd - 0.5) + 0.5; day_jd <= t_end_jd + 0.5; day_jd += 1.0) {
        time::EpochTDB midnight = time::EpochTDB::from_jd(day_jd);
        
        // (Step 1) Daily Chebyshev segment
        auto segment = catalog::fit_chebyshev(working_elements, midnight, 1.0, engine.config());
        
        // (Step 2) Query candidates in the daily corridor
        auto candidates = find_stars_near_segment(catalog_inst, segment, Angle::from_arcsec(300.0), config.max_mag_star);
        
        // (Step 3) Verification and Refinement for each found star
        for (const auto& star : candidates) {
            OccultationParameters params;
            
            // Depending on requested mode (a, b, c), choose the asteroid state source
            if (mode == OccultationRefinementMode::ChebyshevDaily) {
                // Option (a): Use Chebyshev polynomial evaluation (Default/Fast)
                double best_jd = day_jd;
                double min_dist_deg2 = 1e18;
                
                // 1. Coarse search (10 sec steps)
                for (double offset = -0.5; offset <= 0.5; offset += 10.0/86400.0) {
                    double t = day_jd + offset;
                    auto [ra_ast, dec_ast] = segment.evaluate(t);
                    double dra = (star.ra.to_deg() - ra_ast) * std::cos(dec_ast * M_PI / 180.0);
                    double ddec = star.dec.to_deg() - dec_ast;
                    double d2 = dra*dra + ddec*ddec;
                    if (d2 < min_dist_deg2) {
                        min_dist_deg2 = d2;
                        best_jd = t;
                    }
                }
                
                // 2. Compute path parameters at best_jd
                // 2. Compute path parameters at best_jd using analytical derivatives
                auto [pos, vel] = segment.evaluate_full(best_jd);
                auto [ra_final, dec_final, dist_final] = pos;
                auto [vra_day, vdec_day, vdist_day] = vel;
                
                std::cout << "[DEBUG] Refined CA | JD: " << best_jd << " | RA: " << ra_final 
                          << " | Dec: " << dec_final << " | vRA: " << vra_day << " deg/day\n";

                params = compute_parameters(
                    star.ra, star.dec,
                    RightAscension(Angle::from_deg(ra_final)), Declination(Angle::from_deg(dec_final)),
                    physics::Distance::from_au(dist_final), 
                    Angle::from_deg(vra_day / 86400.0), Angle::from_deg(vdec_day / 86400.0), 
                    physics::Velocity::from_au_d(vdist_day),
                    time::EpochTDB::from_jd(best_jd));
            } 
            else {
                // Propagation-based mode
                double best_jd = day_jd;
                double min_dist_deg2 = 1e18;
                for (double offset = -0.5; offset <= 0.5; offset += 60.0/86400.0) {
                    auto [ra_ast, dec_ast] = segment.evaluate(day_jd + offset);
                    double dra = (star.ra.to_deg() - ra_ast) * std::cos(dec_ast * M_PI / 180.0);
                    double ddec = star.dec.to_deg() - dec_ast;
                    if (dra*dra + ddec*ddec < min_dist_deg2) {
                        min_dist_deg2 = dra*dra + ddec*ddec;
                        best_jd = day_jd + offset;
                    }
                }

                auto obs = AstrometryReducer::compute_observation(
                    working_elements, working_elements.epoch, time::EpochTDB::from_jd(best_jd),
                    engine.config(), AstrometricSettings());
                
                if (obs) {
                    params = compute_parameters(
                        star.ra, star.dec,
                        RightAscension(Angle::from_rad((*obs).ra.value)), Declination(Angle::from_rad((*obs).dec.value)),
                        physics::Distance::from_m((*obs).distance.value),
                        Angle::zero(), Angle::zero(), physics::Velocity::zero(),
                        time::EpochTDB::from_jd(best_jd));
                }
            }

            // Diagnostic Log
            if (params.impact_parameter.to_km() < 200000.0) {
                 std::cout << "  [DEBUG] Star ID: " << star.source_id << " | Mag: " << star.g_mag 
                           << " | Impact: " << params.impact_parameter.to_km() << " km\n";
            }

            // Verification: check if impact parameter < 100,000 km (broad discovery)
            if (params.impact_parameter.to_km() < 100000.0) {
                params.star_id = std::to_string(star.source_id);
                params.star_mag = star.g_mag;
                OccultationCandidate cand;
                cand.asteroid_id = asteroid_id;
                cand.star = star;
                cand.params = params;
                results.push_back(cand);
            }
        }
    }
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
    using namespace catalog;
    std::vector<OccultationSystemCandidate> results;
    
    if (body_ids.empty()) return results;
    
    // 1. Initialize SPK Reader for the system
    io::SPKReader system_reader(bsp_path);
    auto& catalog_inst = GaiaDR3Catalog::instance();

    // 2. Daily loop triggered by the PRIMARY body (first in list)
    int primary_id = std::stoi(body_ids[0]);
    double t_start_jd = start.jd();
    double t_end_jd = end.jd();

    for (double day_jd = std::floor(t_start_jd - 0.5) + 0.5; day_jd <= t_end_jd + 0.5; day_jd += 1.0) {
        time::EpochTDB midnight = time::EpochTDB::from_jd(day_jd);
        double et_mid = (midnight.jd() - astdyn::constants::JD2000) * 86400.0;
        
        // Use primary body for discovery
        // Since we don't have Keplerian elements for a BSP body, we "simulate" a segment by sampling
        std::vector<physics::CartesianStateTyped<core::GCRF>> body_states;
        std::vector<physics::CartesianStateTyped<core::GCRF>> earth_states;
        
        // Sample +/- 12 hours
        for (double dt = -43200; dt <= 43200; dt += 3600) {
            double et = et_mid + dt;
            time::EpochTDB t = time::EpochTDB::from_jd(astdyn::constants::JD2000 + et / 86400.0);
            
            auto s_raw = system_reader.getState(primary_id, et);
            body_states.push_back(physics::CartesianStateTyped<core::GCRF>::from_si(
                t,
                s_raw[0]*1000.0, s_raw[1]*1000.0, s_raw[2]*1000.0, 
                s_raw[3]*1000.0, s_raw[4]*1000.0, s_raw[5]*1000.0));
            
            auto e_ssb = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t);
            earth_states.push_back(e_ssb); 
        }

        auto query = astdyn::catalog::make_orbit_query(body_states, earth_states, midnight - time::TimeDuration::from_hours(12), midnight + time::TimeDuration::from_hours(12), Angle::from_arcsec(300.0), config.max_mag_star);
        auto candidates = catalog_inst.query_orbit(query);

        for (const auto& star : candidates) {
            OccultationSystemCandidate system_res;
            system_res.star = star;
            
            // For each starch candidate, check ALL bodies in the system
            for (const auto& b_id_str : body_ids) {
                int b_id = std::stoi(b_id_str);
                
                // 1. Find TCA (Time of Closest Approach)
                double best_jd = day_jd;
                double min_dist_deg2 = 1e18;
                
                // Coarse search (+/- 12 hours, 1 minute steps)
                for (double offset = -0.5; offset <= 0.5; offset += 60.0/86400.0) {
                    double et = et_mid + offset * 86400.0;
                    auto s = system_reader.getState(b_id, et);
                    auto e = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, midnight + time::TimeDuration::from_seconds(offset * 86400.0));
                    
                    Eigen::Vector3d rho = (Eigen::Vector3d(s[0], s[1], s[2]) * 1000.0) - e.position.to_eigen_si();
                    double dist = rho.norm();
                    double ra_a = std::atan2(rho.y(), rho.x()) * 180.0 / M_PI;
                    if (ra_a < 0) ra_a += 360.0;
                    double dec_a = std::asin(rho.z() / dist) * 180.0 / M_PI;
                    
                    double dra = (star.ra.to_deg() - ra_a) * std::cos(dec_a * M_PI / 180.0);
                    double ddec = star.dec.to_deg() - dec_a;
                    double d2 = dra*dra + ddec*ddec;
                    if (d2 < min_dist_deg2) {
                        min_dist_deg2 = d2;
                        best_jd = day_jd + offset;
                    }
                }
                
                if (min_dist_deg2 < 0.01) { // within 0.1 deg
                    double et_ca = (best_jd - astdyn::constants::JD2000) * 86400.0;
                    auto earth_ca = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, time::EpochTDB::from_jd(best_jd));
                    
                    // Apparent State at TCA
                    auto s_ca = system_reader.getState(b_id, et_ca);
                    Eigen::Vector3d rho_ca = (Eigen::Vector3d(s_ca[0], s_ca[1], s_ca[2]) * 1000.0) - earth_ca.position.to_eigen_si();
                    rho_ca = AstrometryReducer::apply_stellar_aberration(rho_ca, earth_ca.velocity.to_eigen_si());
                    
                    auto ast_coord = SkyCoord<core::GCRF>::from_vector(math::Vector3<core::GCRF, physics::Distance>::from_si(rho_ca.x(), rho_ca.y(), rho_ca.z()));
                    
                    // Simple rate at TCA
                    auto s_next = system_reader.getState(b_id, et_ca + 1.0);
                    auto earth_next = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, time::EpochTDB::from_jd(best_jd + 1.0/86400.0));
                    Eigen::Vector3d rho_next = (Eigen::Vector3d(s_next[0], s_next[1], s_next[2]) * 1000.0) - earth_next.position.to_eigen_si();
                    rho_next = AstrometryReducer::apply_stellar_aberration(rho_next, earth_next.velocity.to_eigen_si());
                    auto ast_coord_next = SkyCoord<core::GCRF>::from_vector(math::Vector3<core::GCRF, physics::Distance>::from_si(rho_next.x(), rho_next.y(), rho_next.z()));
                    
                    auto params = compute_parameters(
                        star.ra, star.dec,
                        ast_coord.ra(), ast_coord.dec(),
                        physics::Distance::from_m(rho_ca.norm()),
                        ast_coord_next.ra() - ast_coord.ra(),
                        ast_coord_next.dec() - ast_coord.dec(),
                        physics::Velocity::from_ms(rho_next.norm() - rho_ca.norm()),
                        time::EpochTDB::from_jd(best_jd));
                    
                    if (params.impact_parameter.to_km() < 20000.0) {
                        BodyOccultation body_occ;
                        body_occ.name = b_id_str;
                        body_occ.params = params;
                        body_occ.diameter = physics::Distance::from_km(100.0); // Placeholder
                        system_res.bodies.push_back(body_occ);
                    }
                }
            }
            
            if (!system_res.bodies.empty()) {
                results.push_back(system_res);
            }
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
    
    // 1. Setup STM engine
    StateTransitionMatrix<core::GCRF> stm_engine(engine.propagator());
    
    // 2. Compute STM from t0 to TCA
    auto stm_res = stm_engine.compute(initial_state, params.t_ca);
    
    // 3. Map covariance to TCA
    Eigen::Matrix<double, 6, 6> cov_tca = stm_res.map_covariance(covariance_t0);
    
    // 4. Compute cross-track uncertainty in the B-plane
    auto state_tca = stm_res.final_state;
    auto earth_tca = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, params.t_ca);
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

std::vector<OccultationCandidate> OccultationLogic::find_multi_asteroid_occultations(
    const std::vector<std::string>& asteroid_ids,
    ChebyshevEphemerisManager& manager,
    time::EpochTDB start,
    time::EpochTDB end,
    const OccultationConfig& config,
    AstDynEngine& engine)
{
    using namespace catalog;
    std::vector<OccultationCandidate> results;
    auto& catalog_inst = GaiaDR3Catalog::instance();

    double t_start_jd = start.jd();
    double t_end_jd = end.jd();

    // Iterate through days
    for (double day_jd = std::floor(t_start_jd - 0.5) + 0.5; day_jd <= t_end_jd + 0.5; day_jd += 1.0) {
        time::EpochTDB midnight = time::EpochTDB::from_jd(day_jd);
        if (day_jd < t_start_jd - 0.5 || day_jd > t_end_jd + 0.5) continue;

        for (const auto& id : asteroid_ids) {
            if (!manager.has_body(id)) continue;
            
            const auto& ephem = manager.get_ephemeris(id);
            // Get the segment for this day
            try {
                auto segment = ephem.get_segment(midnight);
                
                // Query stars in the daily corridor for this asteroid
                auto stars = find_stars_near_segment(catalog_inst, segment, Angle::from_arcsec(300.0), config.max_mag_star);
                
                for (const auto& star : stars) {
                    // Refine each candidate
                    double best_jd = day_jd;
                    double min_dist_deg2 = 1e18;
                    
                    // 1. Coarse search
                    for (double offset = -0.5; offset <= 0.5; offset += 10.0/86400.0) {
                        double t = day_jd + offset;
                        if (t < t_start_jd || t > t_end_jd) continue;
                        auto [pos, vel] = segment.evaluate_full(t);
                        double dra = (star.ra.to_deg() - std::get<0>(pos)) * std::cos(std::get<1>(pos) * M_PI / 180.0);
                        double ddec = star.dec.to_deg() - std::get<1>(pos);
                        if (dra*dra + ddec*ddec < min_dist_deg2) {
                            min_dist_deg2 = dra*dra + ddec*ddec;
                            best_jd = t;
                        }
                    }
                    
                    if (min_dist_deg2 < 0.01) { // within 0.1 deg
                        // 2. Refine Star Position (Proper Motion & Parallax) at TCA
                        auto earth_state = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, time::EpochTDB::from_jd(best_jd));
                        auto star_refined = star.predict_at(time::EpochTDB::from_jd(best_jd), earth_state.position.to_eigen_si());

                        // 3. Compute final parameters
                        auto [pos_f, vel_f] = segment.evaluate_full(best_jd);
                        auto params = compute_parameters(
                            star_refined.ra(), star_refined.dec(),
                            RightAscension(Angle::from_deg(std::get<0>(pos_f))),
                            Declination(Angle::from_deg(std::get<1>(pos_f))),
                            physics::Distance::from_au(std::get<2>(pos_f)),
                            Angle::from_deg(std::get<0>(vel_f) / 86400.0),
                            Angle::from_deg(std::get<1>(vel_f) / 86400.0),
                            physics::Velocity::from_au_d(std::get<2>(vel_f)),
                            time::EpochTDB::from_jd(best_jd));
                        
                        // 4. Observability Filtering
                        bool filter_pass = true;
                        if (config.filter_daylight && params.is_daylight) filter_pass = false;
                        if (params.sun_altitude > config.min_sun_altitude) filter_pass = false;
                        if (params.moon_dist < config.min_moon_dist) filter_pass = false;
                        
                        if (filter_pass && params.impact_parameter.to_km() < 50000.0) {
                            params.star_id = std::to_string(star.source_id);
                            params.star_mag = star.g_mag;
                            OccultationCandidate cand;
                            cand.asteroid_id = id;
                            cand.star = star;
                            cand.params = params;
                            results.push_back(cand);
                        }
                    }
                }
            } catch (...) {
                continue;
            }
        }
    }
    return results;
}

} // namespace astdyn::astrometry
