/**
 * @file occultation_test_50936.cpp
 * @brief Refined Occultation Pipeline Test for asteroid 50936 (Nireus).
 * Computes full physical parameters using Apparent Places and compares with Occult4 XML.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/astrometry/AstrometricTypes.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;
using namespace astdyn::catalog;
using namespace astdyn::physics;
using namespace astdyn::astrometry;
using namespace astdyn::ephemeris;
using namespace astdyn::propagation;
using namespace astdyn::time;

/**
 * @brief Apply stellar aberration to a direction vector.
 */
SkyCoord<core::GCRF> get_apparent_asteroid(
    const physics::CartesianStateTyped<core::GCRF>& initial,
    time::EpochTDB t,
    propagation::Propagator& propagator,
    const std::shared_ptr<ephemeris::DE441Provider>& ephem) {
    
    auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t);
    auto sun_ssb = ephem->getState(ephemeris::CelestialBody::SUN, t);
    
    double tau = 0.0;
    Eigen::Vector3d p_ast_ssb;
    for (int i = 0; i < 3; ++i) {
        time::EpochTDB t_emit = time::EpochTDB::from_jd(t.jd() - tau);
        auto ast_helio = propagator.propagate_cartesian(initial, t_emit);
        auto sun_emit = ephem->getState(ephemeris::CelestialBody::SUN, t_emit);
        p_ast_ssb = ast_helio.position.to_eigen_si() + sun_emit.position.to_eigen_si();
        tau = (p_ast_ssb - earth_ssb.position.to_eigen_si()).norm() / (astdyn::constants::C_LIGHT * 1000.0 * 86400.0);
    }
    
    auto rho_geometric = p_ast_ssb - earth_ssb.position.to_eigen_si();
    // For occultation comparison, we need the geometric direction relative to the star.
    // Aberration is not applied here because it would also need to be applied to the star.
    return SkyCoord<core::GCRF>::from_rad(std::atan2(rho_geometric.y(), rho_geometric.x()), std::asin(rho_geometric.z() / rho_geometric.norm()));
}

SkyCoord<core::GCRF> get_apparent_star(const Star& star, time::EpochTDB t, const std::shared_ptr<ephemeris::DE441Provider>& ephem) {
    auto star_p = star.predict_at(t);
    auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t);
    auto sun_ssb = ephem->getState(ephemeris::CelestialBody::SUN, t);

    double ra_p = star_p.ra().to_rad();
    double dec_p = star_p.dec().to_rad();
    Eigen::Vector3d u_star(std::cos(dec_p)*std::cos(ra_p), std::cos(dec_p)*std::sin(ra_p), std::sin(dec_p));

    Eigen::Vector3d rho_geometric = u_star * 1e16; 
    auto rho_deflected = AstrometryReducer::apply_light_deflection(rho_geometric, sun_ssb.position.to_eigen_si() - earth_ssb.position.to_eigen_si());
    auto rho_apparent = AstrometryReducer::apply_stellar_aberration(rho_deflected, earth_ssb.velocity.to_eigen_si());

    return SkyCoord<core::GCRF>::from_rad(std::atan2(rho_apparent.y(), rho_apparent.x()), std::asin(rho_apparent.z() / rho_apparent.norm()));
}

int main() {
    std::string gaia_config_json = "{\"catalog_type\":\"online_esa\",\"timeout_seconds\":30}";
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";

    std::cout << "=== Precision Occultation Refinement: (50936) Nireus ===" << std::endl;
    
    auto ephem = std::make_shared<DE441Provider>(bsp_path);
    PlanetaryEphemeris::setProvider(ephem);
    GaiaDR3Catalog::initialize(gaia_config_json);

    EpochTDB t_start = to_tdb(EpochUTC::from_mjd(61112.041667)); // 01:00 UT
    EpochTDB t_end   = to_tdb(EpochUTC::from_mjd(61112.166667)); // 04:00 UT
    
    io::HorizonsClient horizons;
    auto initial_state = *horizons.query_vectors("50936", t_start, "@sun");

    PropagatorSettings settings;
    settings.include_planets = true;
    settings.integrate_in_ecliptic = false; // Stay in GCRF to match Horizons
    auto integrator = std::make_shared<RKF78Integrator>(0.001);
    Propagator propagator(integrator, std::make_shared<PlanetaryEphemeris>(), settings);

    // Use Official XML Star coordinates for TYC 5546-1581-1 event
    double xml_ra_deg = 204.735733; 
    double xml_dec_deg = -8.8142;
    
    Star target_star;
    target_star.ra = RightAscension(Angle::from_deg(xml_ra_deg));
    target_star.dec = Declination(Angle::from_deg(xml_dec_deg));
    target_star.pm_ra_cosdec = ProperMotion::from_mas_yr(0.0);
    target_star.pm_dec = ProperMotion::from_mas_yr(0.0);
    target_star.parallax = Parallax::from_mas(0.0);

    // 5. CATALOG CHECK
    std::cout << "\n--- CATALOG CHECK ---" << std::endl;
    std::cout << "Star: " << std::fixed << std::setprecision(6) << xml_ra_deg << " " << xml_dec_deg << std::endl;

    // 6. CALCOLO DELLE PREVISIONI
    EpochUTC t0_utc = EpochUTC::from_mjd(61112.115671); 
    EpochTDB t0_tdb = to_tdb(t0_utc);
    
    // 1. Get SSB states at t_obs
    auto sun_ssb = ephem->getState(ephemeris::CelestialBody::SUN, t0_tdb);
    auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t0_tdb);
    
    // 2. Iterative Light Time
    double lt = 0.0;
    Eigen::Vector3d rho_vec;
    for (int i=0; i<3; ++i) {
        auto ast_s = propagator.propagate_cartesian(initial_state, t0_tdb - time::TimeDuration::from_seconds(lt));
        rho_vec = ast_s.position.to_eigen_si() + sun_ssb.position.to_eigen_si() - earth_ssb.position.to_eigen_si();
        lt = rho_vec.norm() / (astdyn::constants::C_LIGHT * 1000.0);
    }
    
    // 3. Apply Stellar Aberration (Fixes ~6000km shift)
    rho_vec = AstrometryReducer::apply_stellar_aberration(rho_vec, earth_ssb.velocity.to_eigen_si());
    
    auto ast_coord_h = SkyCoord<core::GCRF>::from_vector(math::Vector3<core::GCRF, physics::Distance>::from_si(rho_vec.x(), rho_vec.y(), rho_vec.z()));

    // 4. Velocity/Rates at t_obs
    auto t1_tdb = t0_tdb + time::TimeDuration::from_seconds(1.0);
    auto sun_ssb1 = ephem->getState(ephemeris::CelestialBody::SUN, t1_tdb);
    auto earth_ssb1 = ephem->getState(ephemeris::CelestialBody::EARTH, t1_tdb);
    auto ast_s1 = propagator.propagate_cartesian(initial_state, t1_tdb - time::TimeDuration::from_seconds(lt));
    Eigen::Vector3d rho_vec1 = ast_s1.position.to_eigen_si() + sun_ssb1.position.to_eigen_si() - earth_ssb1.position.to_eigen_si();
    // Also aberrant velocity direction (approximate)
    rho_vec1 = AstrometryReducer::apply_stellar_aberration(rho_vec1, earth_ssb1.velocity.to_eigen_si());
    auto ast_coord_h1 = SkyCoord<core::GCRF>::from_vector(math::Vector3<core::GCRF, physics::Distance>::from_si(rho_vec1.x(), rho_vec1.y(), rho_vec1.z()));

    Angle dra = ast_coord_h1.ra() - ast_coord_h.ra();
    Angle ddec = ast_coord_h1.dec() - ast_coord_h.dec();
    physics::Velocity v_rad = physics::Velocity::from_ms(rho_vec1.norm() - rho_vec.norm());

    auto params_nominal = OccultationLogic::compute_parameters(
        target_star.ra, target_star.dec, 
        ast_coord_h.ra(), ast_coord_h.dec(), 
        physics::Distance::from_m(rho_vec.norm()), 
        dra, ddec, v_rad);

    std::cout << "RISULTATO AstDyn (Nominale) - Position Angle: " << params_nominal.position_angle.to_deg() << " deg" << std::endl;
    std::cout << "RISULTATO AstDyn (Nominale) - Shadow Velocity: " << params_nominal.shadow_velocity.to_km_s() << " km/s" << std::endl;
    std::cout << "RISULTATO AstDyn (Nominale) - Time Offset: " << params_nominal.closest_approach_time_offset.to_seconds() << " s" << std::endl;
    
    auto path_nominal = OccultationMapper::compute_path(params_nominal, target_star.ra, target_star.dec, physics::Distance::from_km(100.0), t0_utc, time::TimeDuration::from_seconds(1200.0));
    OccultationMapper::export_kml(path_nominal, "nireus_nominal_africa.kml");
    OccultationMapper::export_svg(path_nominal, "nireus_nominal_africa.svg");

    // --- TEST 2: OFFICIAL OCCULT4 PARAMETERS (USER DATA) ---
    std::cout << "\n--- TEST 2: OFFICIAL OCCULT4 PARAMETERS ---" << std::endl;
    OccultationParameters params_official;
    params_official.xi_ca = physics::Distance::from_km(96.0);
    params_official.eta_ca = physics::Distance::from_km(4596.0);
    params_official.impact_parameter = physics::Distance::from_km(4597.0); 
    params_official.shadow_velocity = physics::Velocity::from_km_s(7.31);
    params_official.position_angle = Angle::from_deg(271.0);
    params_official.dxi_dt = physics::Velocity::from_ms(params_official.shadow_velocity.to_ms() * std::sin(params_official.position_angle.to_rad()));
    params_official.deta_dt = physics::Velocity::from_ms(params_official.shadow_velocity.to_ms() * std::cos(params_official.position_angle.to_rad()));
    params_official.closest_approach_time_offset = time::TimeDuration::zero();
    params_official.time_uncertainty = time::TimeDuration::zero();
    params_official.cross_track_uncertainty = physics::Distance::zero();
    
    std::cout << "OFFICIAL Occult4 - Shadow Velocity: " << params_official.shadow_velocity.to_km_s() << " km/s" << std::endl;

    auto path_official = OccultationMapper::compute_path(params_official, target_star.ra, target_star.dec, physics::Distance::from_km(100.0), t0_utc, time::TimeDuration::from_seconds(1200.0));
    OccultationMapper::export_kml(path_official, "nireus_official_occult4.kml");

    // Map Global Comparison
    std::vector<OccultationPath> paths = {path_nominal, path_official};
    std::vector<std::string> labels = {"AstDyn Nominal (JPL#105)", "Occult4 Official (JPL#39)"};
    std::vector<std::string> colors = {"#06b6d4", "#ef4444"};
    OccultationMapper::export_global_svg(paths, labels, colors, "nireus_global_map.svg");

    // SINTESI CONFRONTO
    double d_xi = (params_nominal.xi_ca - params_official.xi_ca).to_km();
    double d_eta = (params_nominal.eta_ca - params_official.eta_ca).to_km();
    double d_dist = std::sqrt(d_xi*d_xi + d_eta*d_eta);
    double d_pa = (params_nominal.position_angle - params_official.position_angle).to_deg();

    std::cout << "\n--- SINTESI CONFRONTO ---" << std::endl;
    std::cout << "Differenza Centro (km): " << std::fixed << std::setprecision(2) << d_dist << " km" << std::endl;
    std::cout << "Differenza Direzione Traccia (deg): " << d_pa << " deg" << std::endl;
    
    if (!path_nominal.center_line.empty() && !path_official.center_line.empty()) {
        double lat_range_nom = std::abs(path_nominal.center_line.front().lat.to_deg() - path_nominal.center_line.back().lat.to_deg());
        double lon_range_nom = std::abs(path_nominal.center_line.front().lon.to_deg() - path_nominal.center_line.back().lon.to_deg());
        double lat_range_off = std::abs(path_official.center_line.front().lat.to_deg() - path_official.center_line.back().lat.to_deg());
        double lon_range_off = std::abs(path_official.center_line.front().lon.to_deg() - path_official.center_line.back().lon.to_deg());
        
        std::cout << "Lunghezza Traccia Nominale: ~" << std::sqrt(lat_range_nom*lat_range_nom + lon_range_nom*lon_range_nom) << " deg" << std::endl;
        std::cout << "Lunghezza Traccia Official: ~" << std::sqrt(lat_range_off*lat_range_off + lon_range_off*lon_range_off) << " deg" << std::endl;
    }

    std::cout << "\n--- CONCLUSIONE ---" << std::endl;
    std::cout << "Creati file KML e SVG per tutte le tracce." << std::endl;
    std::cout << "File SVG Globale con terre: nireus_global_map.svg" << std::endl;

    GaiaDR3Catalog::shutdown();
    return 0;
}
