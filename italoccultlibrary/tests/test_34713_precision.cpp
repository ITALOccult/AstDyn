#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include "astdyn_wrapper.h"
#include "orbital_conversions.h"
#include <astdyn/AstDyn.hpp>
#include <astdyn/api/OrbitFitAPI.hpp>
#include <astdyn/propagation/HighPrecisionPropagator.hpp>

using namespace ioccultcalc;

int main() {
    astdyn::initialize();
    
    std::cout << "================================================================================" << std::endl;
    std::cout << "PRECISION VALIDATION FOR ASTEROID 34713 (2001 PZ28)" << std::endl;
    std::cout << "================================================================================" << std::endl;

    // --- REFERENCE DATA (JPL Horizons Truth) ---
    // Target: 2026-Jan-10 00:00:00 UTC (MJD 61050.0 TDB)
    // Source: JPL Horizons API (Observer Table, 500@399)
    const double target_jd_tdb = 2461050.5;
    const double target_mjd_tdb = 61050.0;
    
    // RA: 05 43 51.17 -> (5 + 43/60.0 + 51.17/3600.0) * 15.0 = 85.963208333
    // Dec: +26 30 58.5 -> 26 + 30/60.0 + 58.5/3600.0 = 26.51625
    const double jpl_ra_truth = 85.963208333;
    const double jpl_dec_truth = 26.51625;

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Target UTC: 2026-01-10 00:00:00 (MJD TDB: " << target_mjd_tdb << ")" << std::endl;
    std::cout << "JPL Truth RA:  " << jpl_ra_truth << " deg" << std::endl;
    std::cout << "JPL Truth Dec: " << jpl_dec_truth << " deg" << std::endl << std::endl;

    // --- PROPAGATOR SETUP ---
    astdyn::propagation::HighPrecisionPropagator::Config prop_cfg;
    prop_cfg.de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    prop_cfg.tolerance = 1e-13;
    astdyn::propagation::HighPrecisionPropagator hp_prop(prop_cfg);

    // Helper for reporting
    auto report = [&](const std::string& label, const astdyn::propagation::KeplerianElements& elements, 
                      astdyn::propagation::HighPrecisionPropagator::InputFrame frame,
                      const astdyn::OrbitContext& ctx) {
        auto res = hp_prop.calculateGeocentricObservation(elements, target_jd_tdb, frame);
        
        std::cout << "[" << label << "]" << std::endl;
        std::cout << "  Context: " << ctx.toString() << std::endl;
        std::cout << "  RA:  " << std::setw(15) << res.ra_deg << " (Diff: " << std::setw(10) << (res.ra_deg - jpl_ra_truth) * 3600.0 << " arcsec)" << std::endl;
        std::cout << "  Dec: " << std::setw(15) << res.dec_deg << " (Diff: " << std::setw(10) << (res.dec_deg - jpl_dec_truth) * 3600.0 << " arcsec)" << std::endl;
        std::cout << "  Dist: " << res.distance_au << " AU" << std::endl << std::endl;
    };

    // --- CASE A: NO FIT (ASTDYS ELEMENTS - CONVERTED TO OSCULATING) ---
    std::string eq1_path = "/Users/michelebigi/Documents/Develop/ASTDYN/ITALOccultLibrary/astdyn/data/34713.eq1";
    auto equ_a = astdyn::api::OrbitFitAPI::parse_eq1(eq1_path);
    auto kep_a_osc = astdyn::api::OrbitFitAPI::convert_mean_equinoctial_to_osculating(equ_a);

    astdyn::OrbitContext ctx_a_osc;
    ctx_a_osc.frame = astdyn::ReferenceFrame::EQUATORIAL_J2000;
    ctx_a_osc.model = astdyn::OrbitModel::OSCULATING;

    report("CASE A: AstDys Nominal (Mean -> Osculating)", kep_a_osc, astdyn::propagation::HighPrecisionPropagator::InputFrame::EQUATORIAL, ctx_a_osc);

    // --- CASE C: JPL ELEMENTS (Moved up for faster feedback) ---
    astdyn::propagation::KeplerianElements kep_c;
    // Values from JPL Horizons elements fetch for 61050.0
    kep_c.semi_major_axis = 4.195820769563343E+08 / astdyn::constants::AU;
    kep_c.eccentricity = 0.1732535539331025;
    kep_c.inclination = 7.978063854954303 * M_PI / 180.0;
    kep_c.longitude_ascending_node = 289.3398783780134 * M_PI / 180.0;
    kep_c.argument_perihelion = 64.67935047633746 * M_PI / 180.0;
    kep_c.mean_anomaly = 80.19889798428873 * M_PI / 180.0;
    kep_c.epoch_mjd_tdb = 61050.0;
    kep_c.gravitational_parameter = astdyn::constants::GMS;

    astdyn::OrbitContext ctx_c;
    ctx_c.frame = astdyn::ReferenceFrame::ECLIPTIC_J2000;
    ctx_c.model = astdyn::OrbitModel::OSCULATING;

    report("CASE C: JPL Elements", kep_c, astdyn::propagation::HighPrecisionPropagator::InputFrame::ECLIPTIC, ctx_c);

    // --- CASE B: FITTED ---
    std::string rwo_path = "/Users/michelebigi/Documents/Develop/ASTDYN/ITALOccultLibrary/astdyn/data/34713.rwo";
    auto fit_res = astdyn::api::OrbitFitAPI::run_fit(eq1_path, rwo_path, "", true, prop_cfg.de441_path);
    if (fit_res.success) {
        astdyn::OrbitContext ctx_b;
        ctx_b.frame = astdyn::ReferenceFrame::EQUATORIAL_J2000;
        ctx_b.model = astdyn::OrbitModel::OSCULATING;
        report("CASE B: Fitted Propagation", fit_res.fitted_orbit, astdyn::propagation::HighPrecisionPropagator::InputFrame::EQUATORIAL, ctx_b);
    }


    std::cout << "================================================================================" << std::endl;
    return 0;
}
