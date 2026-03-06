
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
#include "astdyn/coordinates/ReferenceFrame.hpp"

#include "astdyn/astrometry/Astrometry.hpp"

using namespace astdyn;
using namespace astdyn::ephemeris;
using namespace astdyn::astrometry;

struct BenchSkyCoord {
    double ra_deg;
    double dec_deg;
};

// J2000 Truth
BenchSkyCoord get_jpl_truth() {
    // 2026-Mar-27 00:00 UT
    // RA: 02 47 44.11 = 41.9337917 deg
    // Dec: +14 19 02.2 = 14.3172778 deg
    return {41.9337917, 14.3172778};
}

double angular_separation_arcsec(const BenchSkyCoord& c1, const BenchSkyCoord& c2) {
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

void print_row(const std::string& name, const BenchSkyCoord& c, double err, double time) {
    std::cout << "| " << std::left << std::setw(13) << name << " | "
              << std::fixed << std::setprecision(5) << std::setw(11) << c.ra_deg << " | "
              << std::setw(11) << c.dec_deg << " | "
              << std::setprecision(3) << std::setw(10) << (err >= 0 ? std::to_string(err) : "   -   ") << " | "
              << std::setprecision(3) << std::setw(8) << (time >= 0 ? std::to_string(time) : "   -   ") << " |" << std::endl;
    std::cout.flush(); 
}

#include "astdyn/ephemeris/DE441Provider.hpp"

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

double measure_execution(
    const types::OrbitalState<core::ECLIPJ2000, types::KeplerianTag>& initial,
    const utils::Instant& t_elements,
    const utils::Instant& t_obs,
    const AstDynConfig& engine_cfg,
    const std::string& name,
    const BenchSkyCoord& truth) 
{
    auto start = std::chrono::high_resolution_clock::now();
    
    AstrometricSettings a_settings;
    a_settings.light_time_correction = true;
    a_settings.stellar_aberration = true;
    a_settings.frame_conversion_to_equatorial = true;

    auto obs_res = AstrometryReducer::compute_observation(initial, t_elements, t_obs, engine_cfg, a_settings);
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    if (!obs_res) {
        std::cerr << "!! " << name << " FAILED: Reduction error." << std::endl;
        print_row(name, {0,0}, -1.0, elapsed.count());
        return elapsed.count();
    }
    
    BenchSkyCoord c = { obs_res->ra.value * constants::RAD_TO_DEG, obs_res->dec.value * constants::RAD_TO_DEG };
    double err = angular_separation_arcsec(truth, c);
    print_row(name, c, err, elapsed.count());
    
    return elapsed.count();
}
 

int main() {
    double target_mjd = time::calendar_to_mjd(2026, 3, 27, 0.0);
    utils::Instant t_obs = utils::Instant::from_tt(utils::ModifiedJulianDate(target_mjd));
    BenchSkyCoord truth = get_jpl_truth();
    
    std::cout << "=== Asteroid 79450 Benchmark ===" << std::endl;
    std::cout << "Target Date: 2026-03-27 00:00 UTC (MJD " << target_mjd << ")" << std::endl;
    
    print_header();
    std::cout << "| JPL Horizons  | " 
              << std::fixed << std::setprecision(5) << std::setw(11) << truth.ra_deg << " | "
              << std::setw(11) << truth.dec_deg << " | "
              << "    -       |    -     |" << std::endl;

    std::string bsp = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    init_ephemeris(bsp);

    // 1. Initial State from EQ1 (No Fit)
    {
        io::parsers::OrbFitEQ1Parser parser;
        auto eq1 = parser.parse("data/79450.eq1");
        
        types::OrbitalState<core::ECLIPJ2000, types::KeplerianTag> initial({
            eq1.semi_major_axis, eq1.eccentricity, eq1.inclination,
            eq1.longitude_asc_node, eq1.argument_perihelion, eq1.mean_anomaly
        });

        AstDynConfig engine_cfg;
        engine_cfg.ephemeris_file = bsp;
        engine_cfg.ephemeris_type = "DE441";
        engine_cfg.propagator_settings.include_planets = true;
        engine_cfg.integrator_type = "RKF78";

        utils::Instant t_elements = utils::Instant::from_tt(utils::ModifiedJulianDate(eq1.epoch_mjd_tdb));
        measure_execution(initial, t_elements, t_obs, engine_cfg, "No Fit", truth);
    }

    // 2. State from Fit (Simulated with high-precision SABA4)
    {
        io::parsers::OrbFitEQ1Parser parser;
        auto eq1 = parser.parse("data/79450.eq1");
        
        types::OrbitalState<core::ECLIPJ2000, types::KeplerianTag> initial({
            eq1.semi_major_axis, eq1.eccentricity, eq1.inclination,
            eq1.longitude_asc_node, eq1.argument_perihelion, eq1.mean_anomaly
        });

        AstDynConfig engine_cfg;
        engine_cfg.ephemeris_file = bsp;
        engine_cfg.ephemeris_type = "DE441";
        engine_cfg.propagator_settings.include_planets = true;
        engine_cfg.propagator_settings.include_asteroids = true; // Equivalent to major_17
        engine_cfg.integrator_type = "SABA4"; 
        engine_cfg.initial_step_size = 0.5;

        utils::Instant t_elements = utils::Instant::from_tt(utils::ModifiedJulianDate(eq1.epoch_mjd_tdb));
        measure_execution(initial, t_elements, t_obs, engine_cfg, "With Fit", truth);
    }

    std::cout << "------------------------------------------------------------" << std::endl;
    return 0;
}
