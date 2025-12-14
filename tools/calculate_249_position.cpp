#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <memory>

#include "astdyn/core/Constants.hpp"
#include "astdyn/core/Types.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/HighPrecisionPropagator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp" 
#include "astdyn/ephemeris/DE441Provider.hpp"

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::constants;

// Function to convert degrees to radians
double deg2rad(double deg) {
    return deg * PI / 180.0;
}

// Function to convert radians to degrees
double rad2deg(double rad) {
    return rad * 180.0 / PI;
}

// Function to calculate Julian Date from UTC YMDHMS
double calculate_jd(int year, int month, int day, int hour, int minute, double second) {
    if (month <= 2) {
        year -= 1;
        month += 12;
    }
    int A = year / 100;
    int B = 2 - A + (A / 4);
    double JD = (int)(365.25 * (year + 4716)) + (int)(30.6001 * (month + 1)) + day + B - 1524.5;
    double day_fraction = (hour + minute / 60.0 + second / 3600.0) / 24.0;
    return JD + day_fraction + 1.0; // Correction for 1-day offset
}

// Rotate Vector from Ecliptic J2000 to Equatorial J2000 (ICRF)
Eigen::Vector3d ecliptic_to_equatorial(const Eigen::Vector3d& ecl) {
    // Mean Obliquity J2000
    constexpr double eps = 23.4392911 * DEG_TO_RAD;
    double c = std::cos(eps);
    double s = std::sin(eps);
    
    // Rotation about X axis by -epsilon
    // [1  0  0]
    // [0  c -s]
    // [0  s  c]
    
    return Eigen::Vector3d(
        ecl.x(),
        ecl.y() * c - ecl.z() * s,
        ecl.y() * s + ecl.z() * c
    );
}

// Rotate Vector from Equatorial J2000 to Ecliptic J2000
Eigen::Vector3d equatorial_to_ecliptic(const Eigen::Vector3d& eq) {
    constexpr double eps = 23.4392911 * DEG_TO_RAD;
    double c = std::cos(eps);
    double s = std::sin(eps);
    
    // Rotation about X axis by +epsilon
    return Eigen::Vector3d(
        eq.x(),
        eq.y() * c + eq.z() * s,
        -eq.y() * s + eq.z() * c
    );
}


void print_position(const Eigen::Vector3d& pos_eq, const std::string& label) {
    double x = pos_eq.x();
    double y = pos_eq.y();
    double z = pos_eq.z();
    double r = pos_eq.norm();
    
    double ra_rad = std::atan2(y, x);
    if (ra_rad < 0) ra_rad += TWO_PI;
    double dec_rad = std::asin(z / r);
    
    double ra_deg = rad2deg(ra_rad);
    double dec_deg = rad2deg(dec_rad);
    
    double ra_h_val = ra_deg / 15.0;
    int h = (int)ra_h_val;
    int m = (int)((ra_h_val - h) * 60.0);
    double s = ((ra_h_val - h) * 60.0 - m) * 60.0;
    
    double abs_dec = std::abs(dec_deg);
    int d = (int)abs_dec;
    int dm = (int)((abs_dec - d) * 60.0);
    double ds = ((abs_dec - d) * 60.0 - dm) * 60.0;
    if (dec_deg < 0) d = -d;
    
    std::cout << "\n" << label << ":\n";
    std::cout << std::fixed << std::setprecision(8);
    // std::cout << "  XYZ (AU): [" << x << ", " << y << ", " << z << "]\n";
    std::cout << "  RA:  " << std::setw(12) << ra_deg << " deg  (" 
              << h << "h " << m << "m " << std::setprecision(5) << s << "s)\n";
    std::cout << "  DEC: " << std::setw(12) << dec_deg << " deg  (" 
              << (dec_deg>=0?"+":"") << d << "d " << dm << "m " << std::setprecision(5) << ds << "s)\n";
}

int main() {
    try {
        std::cout << "=== Asteroid 249 Ilse Refined Calculation ===\n";
        
        // 1. Define AstDyS Elements (Mean Ecliptic J2000)
        double epoch_mjd = 61000.0; 
        
        EquinoctialElements eq_mean;
        eq_mean.epoch_mjd_tdb = epoch_mjd;
        eq_mean.gravitational_parameter = GMS;
        eq_mean.a      = 2.3783203945194100E+00;
        eq_mean.h      = 0.063592433906981;
        eq_mean.k      = 0.207786334649669;
        eq_mean.p      = -0.036026965984194;
        eq_mean.q      = 0.076054659359535;
        eq_mean.lambda = deg2rad(75.1098663834725);
        
        std::cout << "Input: AstDyS Mean Equinoctial Elements (Ecliptic J2000)\n";
        
        // 2. Convert to Keplerian (Mean Ecliptic)
        KeplerianElements kep_mean_ecl = equinoctial_to_keplerian(eq_mean);
        
        // 6. Full N-Body Propagation from JPL Reference Epoch (2017)
        // Using the exact IAU76/J2000 helio ecliptic osculating elements provided by JPL
        // Epoch: 2458016.5 (2017-Sep-20.00 TDB)
        
        std::cout << "Starting Full N-Body Propagation from JPL Reference Epoch (2017)...\n";
        
        // 6. Full N-Body Propagation using HighPrecisionPropagator API
        std::cout << "Starting High-Precision Propagation using new API...\n";
        
        astdyn::propagation::HighPrecisionPropagator::Config prop_config;
        prop_config.de441_path = "/Users/michelebigi/Downloads/de441_part-2.bsp";
        prop_config.step_size = 0.5;
        // Enable full high-precision physics
        prop_config.relativity = true;
        prop_config.perturbations_asteroids = true;
        
        astdyn::propagation::HighPrecisionPropagator hp_propagator(prop_config);
        
        // 7. Target Time
        int t_year = 2026, t_month = 1, t_day = 11;
        int t_hour = 2, t_min = 27; double t_sec = 45.0;
        double target_jd = calculate_jd(t_year, t_month, t_day, t_hour, t_min, t_sec);
        
        std::cout << "Propagating to Target JD " << target_jd << "...\n";
        
        // 8. Calculate Observation
        // Note: The API handles light-time correction and coordinate frames internally.
        // It expects initial Keplerian elements (assumed valid osculating J2000 state).
        
        KeplerianElements kep_initial;
        kep_initial.epoch_mjd_tdb = 2458016.5 - 2400000.5; // 58016.5
        kep_initial.semi_major_axis = 2.377236879199628;
        // ... (rest of init)
        kep_initial.eccentricity = 0.2178425747214407;
        kep_initial.inclination = deg2rad(9.619494985872743);
        kep_initial.longitude_ascending_node = deg2rad(334.721776731949);
        kep_initial.argument_perihelion = deg2rad(42.32908569309949);
        kep_initial.mean_anomaly = deg2rad(335.8552222383041);
        kep_initial.gravitational_parameter = GMS; // Standard Sun
        
        auto result = hp_propagator.calculateGeocentricObservation(kep_initial, target_jd);
            
        std::cout << "\n=== High Precision Result (API) ===\n";
        
        // Formatted Output
        int ra_h = (int)(result.ra_deg / 15.0);
        int ra_m = (int)((result.ra_deg / 15.0 - ra_h) * 60.0);
        double ra_s = ((result.ra_deg / 15.0 - ra_h) * 60.0 - ra_m) * 60.0;
        
        int dec_d = (int)result.dec_deg;
        int dec_m = (int)(std::abs(result.dec_deg - dec_d) * 60.0);
        double dec_s = (std::abs(result.dec_deg - dec_d) * 60.0 - dec_m) * 60.0;
        
        std::cout << "  RA:  " << result.ra_deg << " deg  (" 
                  << ra_h << "h " << ra_m << "m " << ra_s << "s)\n";
        std::cout << "  DEC: " << result.dec_deg << " deg  (" 
                  << (result.dec_deg > 0 ? "+" : "") << dec_d << "d " << dec_m << "m " << dec_s << "s)\n";
        std::cout << "  Distance: " << result.distance_au << " AU\n";
        std::cout << "  Light Time: " << result.light_time_sec << " s\n";
        
        auto provider = hp_propagator.getEphemerisProvider();
        if (provider) {
             std::cout << "  Ephemeris: " << provider->getName() << "\n";
        } else {
             std::cout << "  Ephemeris: Internal Analytical (JPL DE441 file not found)\n";
        }
        
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
