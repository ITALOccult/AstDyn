
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <memory>

#include "astdyn/ephemeris/AsteroidFitter.hpp"
#include "astdyn/ephemeris/AsteroidFitConfig.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "astdyn/core/Constants.hpp"

using namespace astdyn;
using namespace astdyn;
using namespace astdyn::ephemeris;

struct SkyCoord {
    double ra_deg;
    double dec_deg;
};

// J2000 Truth
SkyCoord get_jpl_truth() {
    // 2026-Mar-27 00:00 UT
    // RA: 02 47 44.11 = 41.9337917 deg
    // Dec: +14 19 02.2 = 14.3172778 deg
    return {41.9337917, 14.3172778};
}

SkyCoord cartesian_to_radec_deg(const Eigen::Vector3d& pos) {
    double x = pos.x();
    double y = pos.y();
    double z = pos.z();
    double r = std::sqrt(x*x + y*y + z*z);
    
    double dec = std::asin(z / r);
    double ra = std::atan2(y, x);
    if (ra < 0) ra += 2.0 * constants::PI;
    
    return {ra * constants::RAD_TO_DEG, dec * constants::RAD_TO_DEG};
}

double angular_separation_arcsec(const SkyCoord& c1, const SkyCoord& c2) {
    double r1 = c1.ra_deg * constants::DEG_TO_RAD;
    double d1 = c1.dec_deg * constants::DEG_TO_RAD;
    double r2 = c2.ra_deg * constants::DEG_TO_RAD;
    double d2 = c2.dec_deg * constants::DEG_TO_RAD;

    double dra = (r1 - r2) * std::cos(d1);
    double ddec = d1 - d2;
    double dist_rad = std::sqrt(dra*dra + ddec*ddec);
    return dist_rad * 3600.0 * constants::RAD_TO_DEG;
}

void print_header() {
    std::cout << "\n------------------------------------------------------------" << std::endl;
    std::cout << "| Scenario      | RA (deg)    | Dec (deg)   | Error (\")  | Time (s) |" << std::endl;
    std::cout << "|---------------|-------------|-------------|------------|----------|" << std::endl;
    std::cout.flush();
}

void print_row(const std::string& name, const SkyCoord& c, double err, double time) {
    std::cout << "| " << std::left << std::setw(13) << name << " | "
              << std::fixed << std::setprecision(5) << std::setw(11) << c.ra_deg << " | "
              << std::setw(11) << c.dec_deg << " | "
              << std::setprecision(3) << std::setw(10) << (err >= 0 ? std::to_string(err) : "   -   ") << " | "
              << std::setprecision(3) << std::setw(8) << (time >= 0 ? std::to_string(time) : "   -   ") << " |" << std::endl;
    std::cout.flush(); // Ensure output is visible immediately
}

#include "astdyn/ephemeris/DE441Provider.hpp"

// Helper to get Earth State
static std::shared_ptr<ephemeris::DE441Provider> g_provider;

void init_ephemeris(const std::string& bsp_path) {
    if (std::filesystem::exists(bsp_path)) {
        g_provider = std::make_shared<ephemeris::DE441Provider>(bsp_path);
        ephemeris::PlanetaryEphemeris::setProvider(g_provider);
        std::cout << "Loaded DE441 for Earth position: " << bsp_path << std::endl;
    } else {
        std::cerr << "WARNING: DE441 not found at " << bsp_path << ". Using Analytical (Low Precision)." << std::endl;
        ephemeris::PlanetaryEphemeris::setProvider(nullptr);
    }
}

Eigen::Vector3d get_earth_position(double mjd_tdb) {
    double jd = time::mjd_to_jd(mjd_tdb);
    auto earth_state = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, jd);
    auto sun_state = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::SUN, jd);
    return earth_state.position() - sun_state.position();
}

double measure_execution(std::function<ephemeris::AsteroidFitResult()> func, const std::string& name, const SkyCoord& truth, double target_mjd) {
    auto start = std::chrono::high_resolution_clock::now();
    ephemeris::AsteroidFitResult res = func();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    if (res.fitted_positions.empty()) {
        std::cerr << "!! " << name << " FAILED: No positions generated." << std::endl;
        print_row(name, {0,0}, -1.0, elapsed.count());
        return elapsed.count();
    }
    
    Eigen::Vector3d earth_pos = get_earth_position(target_mjd);
    Eigen::Vector3d ast_helio = res.fitted_positions.back();
    Eigen::Vector3d rho = ast_helio - earth_pos; // Geocentric Vector
    
    SkyCoord c = cartesian_to_radec_deg(rho);
    double err = angular_separation_arcsec(truth, c);
    print_row(name, c, err, elapsed.count());
    
    if (!res.success) {
        std::cout << "   (Note: " << name << " fit reported failure)" << std::endl;
    }
    return elapsed.count();
}

int main() {
    double target_mjd = time::calendar_to_mjd(2026, 3, 27, 0.0);
    SkyCoord truth = get_jpl_truth();
    
    std::cout << "=== Asteroid 79450 Benchmark ===" << std::endl;
    std::cout << "Target Date: 2026-03-27 00:00 UTC (MJD " << target_mjd << ")" << std::endl;
    
    print_header();
    // Print Truth
    std::cout << "| JPL Horizons  | " 
              << std::fixed << std::setprecision(5) << std::setw(11) << truth.ra_deg << " | "
              << std::setw(11) << truth.dec_deg << " | "
              << "    -       |    -     |" << std::endl;

    std::cout << "Loading DE441 Ephemeris: /Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp..." << std::endl;
    init_ephemeris("/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp");

    // 1. No Fit
    measure_execution([&]() {
        io::parsers::OrbFitEQ1Parser parser;
        // Assuming running from build directory where 'data' link exists
        auto eq1 = parser.parse("data/79450.eq1");
        
        AsteroidFitConfig cfg;
        cfg.orbit.a = eq1.semi_major_axis;
        cfg.orbit.e = eq1.eccentricity;
        cfg.orbit.i = eq1.inclination;
        cfg.orbit.Omega = eq1.longitude_asc_node;
        cfg.orbit.omega = eq1.argument_perihelion;
        cfg.orbit.M = eq1.mean_anomaly;
        cfg.orbit.epoch_mjd_tdb = eq1.epoch_mjd_tdb;
        
        cfg.mjd_observations = { target_mjd };
        cfg.output_equatorial = true;
        
        // Full Physics (No Fit needs Engine too for accurate propagation of initial state)
        cfg.engine_config.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
        cfg.engine_config.include_planets = true;
        cfg.engine_config.asteroid_mode = "none"; // Disable for speed debug
        cfg.engine_config.integrator_type = "RKF78";
        
        return AsteroidFitter::fitFromConfig(cfg);
    }, "No Fit", truth, target_mjd);

    // 2. With Fit (High Precision Configuration via JSON-like struct)
    std::cout << "Running Fit (High Precision Config)..." << std::endl;
    measure_execution([&]() {
        AsteroidFitConfig cfg;
        cfg.eq1_file = "data/79450.eq1";
        cfg.rwo_file = "data/79450.rwo";
        cfg.mjd_observations = { target_mjd };
        cfg.output_equatorial = true;
        
        // Full Physics Configuration matching propagate_79148.cpp
        cfg.engine_config.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
        cfg.engine_config.include_planets = true;
        cfg.engine_config.asteroid_mode = "major_17"; // Enable Major 17 Asteroids
        cfg.engine_config.integrator_type = "RKF78";
        cfg.engine_config.initial_step_size = 0.5;
        cfg.engine_config.tolerance = 1e-13;
        
        return AsteroidFitter::fitFromConfig(cfg);
    }, "With Fit", truth, target_mjd);

    std::cout << "------------------------------------------------------------" << std::endl;
    return 0;
}
