#include <astdyn/AstDyn.hpp>
#include <astdyn/api/OrbitFitAPI.hpp>
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    astdyn::initialize();

    std::string eq1_file = "../data/34713.eq1";
    std::string rwo_file = "../data/34713.rwo";
    
    // Target: 2026-Jan-10 00:00:00 UTC
    double jpl_ra = 85.96320833;
    double jpl_dec = 26.51625;
    
    // JPL State Vector (ICRF EQUATORIAL) -> AU, AU/d
    // Source: JPL Horizons, Target Item: 34713, Coord: ICRF/J2000
    double km_to_au = 1.0 / 149597870.7;
    double day_to_sec = 86400.0;
    
    // Values extracted for ICRF (Example)
    // Position: -0.2114755106 AU, 2.5259468962 AU, 1.2058444983 AU
    Eigen::Vector3d jpl_pos(-0.2114755106, 2.5259468962, 1.2058444983);
    Eigen::Vector3d jpl_vel(-0.0101344445, -0.0006399134, -0.0002877543); // AU/d (Approx from state)

    // Time conversion
    double target_mjd_tdb = astdyn::time::utc_to_tdb(61050.5);
    double target_jd_tdb = astdyn::time::mjd_to_jd(target_mjd_tdb);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "=== Asteroid 34713 Verification ===\n";
    std::cout << "Target Date: 2026-01-10 00:00:00 UTC (MJD TDB: " << target_mjd_tdb << ")\n";
    std::cout << "JPL Truth RA:  " << jpl_ra << ", Dec: " << jpl_dec << "\n";
    std::cout << "JPL Truth Pos: " << jpl_pos.transpose() << " AU\n";

    // Standardize Obliquity
    // coordinates::ReferenceFrame is now using 23.439291 deg (IAU 2000)

    propagation::HighPrecisionPropagator::Config config;
    config.de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    propagation::HighPrecisionPropagator propagator(config);

    // --- CASE A: NO FIT (ASTDYS ELEMENTS) ---
    std::cout << "\n--- Case A: Nominal Propagation (AstDys Elements, No Fit) ---\n";
    auto initial_equ = api::OrbitFitAPI::parse_eq1(eq1_file);
    auto initial_kep_ecl = propagation::equinoctial_to_keplerian(initial_equ);

    auto res_a = propagator.calculateGeocentricObservation(initial_kep_ecl, target_jd_tdb, propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC);
    auto state_a = propagator.propagate_cartesian(initial_kep_ecl, target_mjd_tdb, propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC);
    
    std::cout << "Nominal RA:  " << res_a.ra_deg << " deg\n";
    std::cout << "Nominal Dec: " << res_a.dec_deg << " deg\n";
    std::cout << "Diff vs JPL RA:  " << (res_a.ra_deg - jpl_ra) * 3600.0 << " arcsec\n";
    std::cout << "Diff vs JPL Dec: " << (res_a.dec_deg - jpl_dec) * 3600.0 << " arcsec\n";
    std::cout << "Nominal Pos: " << state_a.position.transpose() << " AU\n";
    std::cout << "Diff Pos Dist: " << (state_a.position - jpl_pos).norm() * constants::AU << " km\n";

    // --- CASE B: WITH FIT ---
    std::cout << "\n--- Case B: Differential Correction (Fitting Observations) ---\n";
    auto fit_result = api::OrbitFitAPI::run_fit(eq1_file, rwo_file, "", true, config.de441_path);
    
    if (fit_result.success) {
        std::cout << "\nFit success. RMS RA: " << fit_result.rms_ra << ", RMS Dec: " << fit_result.rms_dec << " arcsec\n";
        auto res_b = propagator.calculateGeocentricObservation(fit_result.fitted_orbit, target_jd_tdb, propagation::HighPrecisionPropagator::InputFrame::EQUATORIAL);
        auto state_b = propagator.propagate_cartesian(fit_result.fitted_orbit, target_mjd_tdb, propagation::HighPrecisionPropagator::InputFrame::EQUATORIAL);
        
        std::cout << "Fitted RA:   " << res_b.ra_deg << " deg\n";
        std::cout << "Fitted Dec:  " << res_b.dec_deg << " deg\n";
        std::cout << "Diff vs JPL RA:  " << (res_b.ra_deg - jpl_ra) * 3600.0 << " arcsec\n";
        std::cout << "Diff vs JPL Dec: " << (res_b.dec_deg - jpl_dec) * 3600.0 << " arcsec\n";
        std::cout << "Fitted Pos:  " << state_b.position.transpose() << " AU\n";
        std::cout << "Diff Pos Dist: " << (state_b.position - jpl_pos).norm() * constants::AU << " km\n";
    }

    astdyn::shutdown();
    return 0;
}
