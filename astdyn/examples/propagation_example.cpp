/**
 * @file propagation_example.cpp
 * @brief Rigorous propagation using Keplerian elements for verification
 */

#include <astdyn/propagation/HighPrecisionPropagator.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/time/TimeScale.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::constants;

int main() {
    std::cout << std::fixed << std::setprecision(10);
    
    // 1. Initial State provided by user (Ecliptic J2000)
    // Epoch: 2458315.5 JD TDB
    double epoch_jd_tdb = 2458315.5;
    double epoch_mjd_tdb = epoch_jd_tdb - 2400000.5;

    // Elements:
    // EC = 0.1748137420083936
    // QR = 2.313247539763603 (Perihelion distance AU)
    // OM = 289.4707834155099 (deg)
    // W  = 63.81392291128233 (deg)
    // IN = 7.967893011038914 (deg)
    // TP = 2458949.8019304527
    
    double e = 0.1748137420083936;
    double q_au = 2.313247539763603;
    double a_au = q_au / (1.0 - e);
    double tp_jd_tdb = 2458949.8019304527;
    
    // Mean motion n = sqrt(GMS / a^3) in rad/day
    double n = std::sqrt(GMS / (a_au * a_au * a_au));
    double m_rad = n * (epoch_jd_tdb - tp_jd_tdb);
    
    KeplerianElements elements;
    elements.epoch_mjd_tdb = epoch_mjd_tdb;
    elements.semi_major_axis = a_au;
    elements.eccentricity = e;
    elements.inclination = 7.967893011038914 * DEG_TO_RAD;
    elements.longitude_ascending_node = 289.4707834155099 * DEG_TO_RAD;
    elements.argument_perihelion = 63.81392291128233 * DEG_TO_RAD;
    elements.mean_anomaly = m_rad;
    elements.gravitational_parameter = GMS;

    // 2. Propagator Configuration
    HighPrecisionPropagator::Config config;
    config.de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    config.perturbations_planets = true;
    config.perturbations_asteroids = true;
    config.relativity = true;
    config.tolerance = 1e-14;

    HighPrecisionPropagator propagator(config);

    // 3. Target Date: 2026-Jan-10 01:56:34 UTC
    double target_fraction = (1.0 + 56.0/60.0 + 34.0/3600.0) / 24.0;
    double target_mjd_utc = time::calendar_to_mjd(2026, 1, 10, target_fraction);
    double target_jd_tdb = time::mjd_to_jd(time::utc_to_tdb(target_mjd_utc));

    std::cout << "Target Date:  2026-01-10 01:56:34 UTC\n";
    std::cout << "Target JD TDB: " << target_jd_tdb << "\n\n";
    
    // 4. Calculate Geocentric Position
    auto result = propagator.calculateGeocentricObservation(elements, target_jd_tdb);

    // 5. Display Results
    std::cout << "--- AstDyn Rigorous Position (with DE441) ---\n";
    std::cout << "RA  (deg): " << result.ra_deg << " deg\n";
    std::cout << "Dec (deg): " << result.dec_deg << " deg\n";
    std::cout << "Distance:  " << result.distance_au << " AU\n";
    
    // Conversion to HMS/DMS
    int ra_h = static_cast<int>(result.ra_deg / 15.0);
    int ra_m = static_cast<int>((result.ra_deg / 15.0 - ra_h) * 60.0);
    double ra_s = ((result.ra_deg / 15.0 - ra_h) * 60.0 - ra_m) * 60.0;
    
    double abs_dec = std::abs(result.dec_deg);
    int dec_d = static_cast<int>(abs_dec);
    int dec_m = static_cast<int>((abs_dec - dec_d) * 60.0);
    double dec_s = ((abs_dec - dec_d) * 60.0 - dec_m) * 60.0;
    char dec_sign = (result.dec_deg >= 0) ? '+' : '-';

    std::cout << "\nComparison with JPL Horizons:\n";
    std::cout << "           AstDyn          JPL (approx)\n";
    std::cout << "RA  (HMS): " << std::setw(2) << std::setfill('0') << ra_h << " " 
              << std::setw(2) << std::setfill('0') << ra_m << " " 
              << std::fixed << std::setprecision(2) << std::setw(5) << ra_s << "   05 43 47.27\n";
    std::cout << "Dec (DMS): " << dec_sign << std::setw(2) << std::setfill('0') << dec_d << " " 
              << std::setw(2) << std::setfill('0') << dec_m << " " 
              << std::fixed << std::setprecision(1) << std::setw(4) << dec_s << "   +26 30 42.0\n";
    std::cout << "Dist (AU): " << std::setprecision(10) << result.distance_au << "   1.87556697\n";

    return 0;
}
