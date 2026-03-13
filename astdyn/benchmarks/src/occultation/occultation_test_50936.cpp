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
    auto rho_deflected = AstrometryReducer::apply_light_deflection(rho_geometric, sun_ssb.position.to_eigen_si() - earth_ssb.position.to_eigen_si());
    auto rho_apparent = AstrometryReducer::apply_stellar_aberration(rho_deflected, earth_ssb.velocity.to_eigen_si());

    return SkyCoord<core::GCRF>::from_rad(std::atan2(rho_apparent.y(), rho_apparent.x()), std::asin(rho_apparent.z() / rho_apparent.norm()));
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

    // Use Official XML Star coordinates for direct comparison (Italy Track)
    double xml_ra_deg = 204.735733; 
    double xml_dec_deg = -8.8142;
    
    Star target_star;
    target_star.ra = RightAscension(Angle::from_deg(xml_ra_deg));
    target_star.dec = Declination(Angle::from_deg(xml_dec_deg));
    target_star.pm_ra_cosdec = ProperMotion::from_mas_yr(0.0);
    target_star.pm_dec = ProperMotion::from_mas_yr(0.0);
    target_star.parallax = Parallax::from_mas(0.0);

    // 5. CATALOG CHECK
    std::cout << "\n--- CATALOG CHECK (Gaia DR3) ---" << std::endl;
    auto& cat = GaiaDR3Catalog::instance();
    auto s_opt = cat.by_tycho2("5546-1581-1");
    if (s_opt) {
        const auto& s = *s_opt;
        std::cout << "Gaia DR3 Match:  RA " << std::fixed << std::setprecision(8) << s.ra.to_deg() 
                  << ", Dec " << s.dec.to_deg() << std::endl;
        std::cout << "Gaia Proper Motion: RA " << s.pm_ra_cosdec.to_arcsec_yr() << " arcsec/yr, Dec " << s.pm_dec.to_arcsec_yr() << " arcsec/yr" << std::endl;
    } else {
        std::cout << "No Gaia match found for TYC 5546-1581-1 in local catalog." << std::endl;
    }

    // 6. CALCOLO DELLE PREVISIONI
    TimeDuration h = TimeDuration::from_seconds(1.0);
    EpochTDB t_nominal = to_tdb(EpochUTC::from_mjd(61112.115671));
    auto ast_helio_check = propagator.propagate_cartesian(initial_state, t_nominal);
    auto sun_ssb_check = ephem->getState(ephemeris::CelestialBody::SUN, t_nominal);
    auto earth_ssb_check = ephem->getState(ephemeris::CelestialBody::EARTH, t_nominal);
    double dist_m = (ast_helio_check.position.to_eigen_si() + sun_ssb_check.position.to_eigen_si() - earth_ssb_check.position.to_eigen_si()).norm();

    // --- TEST 1: PREVISIONE NOMINALE ASTDYN (JPL#105) ---
    std::cout << "\n--- TEST 1: PREVISIONE NOMINALE ASTDYN (JPL#105) ---" << std::endl;
    EpochUTC t0_utc = EpochUTC::from_mjd(61112.115671); // 02:46:34 UT
    EpochTDB t0_tdb = to_tdb(t0_utc);
    
    auto ast_app_real = get_apparent_asteroid(initial_state, t0_tdb, propagator, ephem);
    auto ast_app_h = get_apparent_asteroid(initial_state, t0_tdb + h, propagator, ephem);
    
    double dra_real = (ast_app_h.ra().to_rad() - ast_app_real.ra().to_rad()) / 1.0;
    double ddec_real = (ast_app_h.dec().to_rad() - ast_app_real.dec().to_rad()) / 1.0;
    double dist_m_real = dist_m; // Use the already computed dist_m

    auto params_nominal = OccultationLogic::compute_parameters(
        target_star.ra, target_star.dec,
        ast_app_real.ra(), ast_app_real.dec(),
        dist_m_real, dra_real, ddec_real, 0.0
    );

    std::cout << "RISULTATO AstDyn (Nominale): Impact = " << params_nominal.impact_parameter_km << " km (AFRICA)" << std::endl;
    std::cout << "RISULTATO AstDyn (Nominale): Velocity = " << params_nominal.shadow_velocity_kms << " km/s" << std::endl;
    
    auto path_nominal = OccultationMapper::compute_path(params_nominal, target_star.ra, target_star.dec, t0_utc, 1200.0);
    OccultationMapper::export_kml(path_nominal, "nireus_nominal_africa.kml");
    OccultationMapper::export_svg(path_nominal, "nireus_nominal_africa.svg");

    // --- TEST 2: EMULAZIONE OCCULT4 (JPL#39) ---
    std::cout << "\n--- TEST 2: EMULAZIONE INPUT OCCULT4 (JPL#39) ---" << std::endl;
    
    double occ4_sep_north_deg = 3.52 / 3600.0; 
    auto ast_occ4 = SkyCoord<core::GCRF>::from_deg(xml_ra_deg, xml_dec_deg + occ4_sep_north_deg);
    auto star_occ4 = SkyCoord<core::GCRF>::from_deg(xml_ra_deg, xml_dec_deg);
    
    double dist_km = dist_m / 1000.0;
    double v_occ4 = 7.316; 
    double total_rad_s = (v_occ4 / dist_km);
    double pa_rad = 265.4 * astdyn::constants::DEG_TO_RAD;
    // PA = 0 is North (+Dec), PA = 90 is East (+RA)
    // So dra_cos_dec = sin(PA), ddec = cos(PA)
    double dra_occ4 = total_rad_s * std::sin(pa_rad);
    double ddec_occ4 = total_rad_s * std::cos(pa_rad);

    auto params_occ4 = OccultationLogic::compute_from_sky(star_occ4, ast_occ4, dist_m, dra_occ4, ddec_occ4);
    params_occ4.time_uncertainty_sec = 0.8;

    std::cout << "RISULTATO AstDyn (Emulazione): Impact = " << params_occ4.impact_parameter_km << " km (ITALIA)" << std::endl;
    std::cout << "  - xi_ca  = " << params_occ4.xi_ca_km << " km" << std::endl;
    std::cout << "  - eta_ca = " << params_occ4.eta_ca_km << " km" << std::endl;
    
    auto path_occ4 = OccultationMapper::compute_path(params_occ4, star_occ4.ra(), star_occ4.dec(), t0_utc, 600.0);
    OccultationMapper::export_kml(path_occ4, "nireus_emulated_italy.kml");
    OccultationMapper::export_svg(path_occ4, "nireus_emulated_italy.svg");

    // NEW: Global Map with Earth background
    std::cout << "\n--- GENERAZIONE MAPPA GLOBALE (CON TERRE) ---" << std::endl;
    
    // Create a 3rd 'Shifted' path for comparison (e.g., user manual correction)
    auto params_occ4_shifted = params_occ4;
    params_occ4_shifted.eta_ca_km -= 1000.0; // Shift south by 1000km
    auto path_occ4_shifted = OccultationMapper::compute_path(params_occ4_shifted, star_occ4.ra(), star_occ4.dec(), t0_utc, 600.0);

    std::vector<OccultationPath> global_paths = {path_nominal, path_occ4, path_occ4_shifted};
    std::vector<std::string> global_labels = {
        "AstDyn Nominal (JPL#105 - Africa)", 
        "Occult4 Emulated (JPL#39 - Italy)",
        "Occult4 Shifted (-1000km Eta)"
    };
    std::vector<std::string> global_colors = {"#06b6d4", "#ef4444", "#fbbf24"};
    
    OccultationMapper::export_global_svg(global_paths, global_labels, global_colors, "nireus_global_map.svg");
    std::cout << "File SVG Globale: nireus_global_map.svg" << std::endl;

    std::cout << "\n--- CONCLUSIONE ---" << std::endl;
    std::cout << "Creati file KML e SVG per tutte le tracce." << std::endl;
    std::cout << "File SVG Globale con terre: nireus_global_map.svg" << std::endl;

    GaiaDR3Catalog::shutdown();
    return 0;
}
