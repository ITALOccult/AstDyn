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
#include <cmath>
#include <iostream>
#include <algorithm>

namespace astdyn::astrometry {

OccultationParameters OccultationLogic::compute_parameters(
    const RightAscension& star_ra, const Declination& star_dec,
    const RightAscension& ast_ra, const Declination& ast_dec,
    const physics::Distance& ast_dist,
    const Angle& ast_dra_dt, const Angle& ast_ddec_dt,
    const physics::Velocity& ast_ddist_dt) 
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
    double t_ca = 0.0;
    if (v2 > 1e-18) {
        t_ca = -(xi * dxi + eta * deta) / v2;
    }

    double xi_ca_val = xi + dxi * t_ca;
    double eta_ca_val = eta + deta * t_ca;
    double b = std::sqrt(xi_ca_val * xi_ca_val + eta_ca_val * eta_ca_val);

    // Build Result - Preserve signs (Bug B Fix)
    OccultationParameters params;
    params.xi_ca = physics::Distance::from_m(xi_ca_val);
    params.eta_ca = physics::Distance::from_m(eta_ca_val);
    params.impact_parameter = physics::Distance::from_m(b);
    params.shadow_velocity = physics::Velocity::from_ms(std::sqrt(v2));
    params.dxi_dt = physics::Velocity::from_ms(dxi);
    params.deta_dt = physics::Velocity::from_ms(deta);
    
    // Position angle of the track (direction of velocity on plane)
    // atan2(E, N)
    params.position_angle = Angle::from_rad(std::atan2(dxi, deta)).wrap_0_2pi();

    // Angular velocities (relative velocity components)
    // We use RA velocity scaled by cos(Dec) for the East-West component
    params.d_ra_cos_dec_per_sec = Angle::from_rad(ast_dra_dt_rad_s * std::cos(d));
    params.d_dec_per_sec = Angle::from_rad(ast_ddec_dt_rad_s);
    
    params.closest_approach_time_offset = time::TimeDuration::from_seconds(t_ca);
    params.time_uncertainty = time::TimeDuration::zero();
    params.cross_track_uncertainty = physics::Distance::zero();

    return params;
}

std::vector<OccultationCandidate> OccultationLogic::find_occultations(
    const std::string& asteroid_id,
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
    time::EpochTDB start,
    time::EpochTDB end,
    double max_mag,
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
        auto candidates = find_stars_near_segment(catalog_inst, segment, Angle::from_arcsec(600.0), max_mag);
        
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
                auto [ra_final, dec_final] = segment.evaluate(best_jd);
                
                // Redefine velocities for compute_parameters (numerical derivative)
                double dt_sec = 1.0;
                auto [ra2, dec2] = segment.evaluate(best_jd + dt_sec/86400.0);
                Angle dra_vel = Angle::from_deg(ra2 - ra_final) * (1.0 / dt_sec);
                Angle ddec_vel = Angle::from_deg(dec2 - dec_final) * (1.0 / dt_sec);
                
                params = compute_parameters(
                    star.ra, star.dec,
                    RightAscension(Angle::from_deg(ra_final)), Declination(Angle::from_deg(dec_final)),
                    physics::Distance::from_au(2.5), dra_vel, ddec_vel, physics::Velocity::zero());
                params.t_ca = time::EpochTDB::from_jd(best_jd);
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
                        Angle::zero(), Angle::zero(), physics::Velocity::zero());
                    params.t_ca = time::EpochTDB::from_jd(best_jd);
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
                results.push_back({star, params});
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
    double max_mag,
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

        auto query = astdyn::catalog::make_orbit_query(body_states, earth_states, midnight - time::TimeDuration::from_hours(12), midnight + time::TimeDuration::from_hours(12), Angle::from_arcsec(120.0), max_mag);
        auto candidates = catalog_inst.query_orbit(query);

        for (const auto& star : candidates) {
            OccultationSystemCandidate system_res;
            system_res.star = star;
            
            // For each starch candidate, check ALL bodies in the system
            for (const auto& b_id_str : body_ids) {
                int b_id = std::stoi(b_id_str);
                
                // Get Apparent State and compute parameters
                // (Iterative Light Time + Aberration)
                double lt = 0.0;
                auto earth_mid = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, midnight);
                Eigen::Vector3d rho_vec;
                
                for(int i=0; i<3; ++i) {
                    double et_emit = et_mid - lt;
                    auto s_emit = system_reader.getState(b_id, et_emit);
                    rho_vec = (Eigen::Vector3d(s_emit[0], s_emit[1], s_emit[2]) * 1000.0) - earth_mid.position.to_eigen_si();
                    lt = rho_vec.norm() / (astdyn::constants::C_LIGHT * 1000.0);
                }
                rho_vec = AstrometryReducer::apply_stellar_aberration(rho_vec, earth_mid.velocity.to_eigen_si());
                
                auto ast_coord = SkyCoord<core::GCRF>::from_vector(math::Vector3<core::GCRF, physics::Distance>::from_si(rho_vec.x(), rho_vec.y(), rho_vec.z()));
                
                // Simple rate (1s diff)
                double et1 = et_mid + 1.0;
                auto s1 = system_reader.getState(b_id, et1 - lt); 
                auto earth1 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, midnight + time::TimeDuration::from_seconds(1.0));
                Eigen::Vector3d rho1 = (Eigen::Vector3d(s1[0], s1[1], s1[2]) * 1000.0) - earth1.position.to_eigen_si();
                rho1 = AstrometryReducer::apply_stellar_aberration(rho1, earth1.velocity.to_eigen_si());
                auto ast_coord1 = SkyCoord<core::GCRF>::from_vector(math::Vector3<core::GCRF, physics::Distance>::from_si(rho1.x(), rho1.y(), rho1.z()));
                
                auto params = compute_parameters(
                    star.ra, star.dec,
                    ast_coord.ra(), ast_coord.dec(),
                    physics::Distance::from_m(rho_vec.norm()),
                    ast_coord1.ra() - ast_coord.ra(),
                    ast_coord1.dec() - ast_coord.dec(),
                    physics::Velocity::from_ms(rho1.norm() - rho_vec.norm()));
                
                if (params.impact_parameter.to_km() < 10000.0) {
                    BodyOccultation body_occ;
                    body_occ.name = b_id_str;
                    body_occ.params = params;
                    body_occ.params.t_ca = midnight + params.closest_approach_time_offset;
                    body_occ.diameter = physics::Distance::from_km(100.0); // Placeholder or lookup
                    system_res.bodies.push_back(body_occ);
                }
            }
            
            if (!system_res.bodies.empty()) {
                results.push_back(system_res);
            }
        }
    }
    return results;
}

} // namespace astdyn::astrometry
