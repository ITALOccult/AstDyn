// AsteroidFitConfig.hpp
// Configuration structure for AsteroidFitter, supporting local files and URLs.

#ifndef ASTDYN_EPHEMERIS_ASTEROIDFITCONFIG_HPP
#define ASTDYN_EPHEMERIS_ASTEROIDFITCONFIG_HPP

#include <cstring>
#include <string>
#include <vector>
#include <optional>
#include <fstream>
#include <stdexcept>
#include <nlohmann/json.hpp>

namespace astdyn {
namespace ephemeris {

/**
 * @brief Simple public Keplerian elements used for JSON configuration.
 *
 * Public members match the test expectations (a, e, i, Omega, omega, M).
 */
struct SimpleKeplerianElements {
    double a = 0.0;      // semi-major axis (AU)
    double e = 0.0;      // eccentricity
    double i = 0.0;      // inclination (rad)
    double Omega = 0.0;  // longitude of ascending node (rad)
    double omega = 0.0;  // argument of periapsis (rad)
    double M = 0.0;      // mean anomaly (rad)
    double epoch_mjd_tdb = 0.0; // Epoch of elements (MJD TDB)
};


/**
 * @brief Configuration for the propagation engine.
 */
struct EngineConfig {
    std::string ephemeris_file;
    std::string asteroid_ephemeris_file;
    std::string catalog_bias_file; // path to CSV file
    
    // Asteroid Model Mode: "none", "major_17" (Internal/Analytical), "full_300" (External SPK)
    std::string asteroid_mode = "major_17"; 
    
    bool include_planets = true;
    bool include_relativity = false; // GR
    bool include_moon = true;
    
    // Detailed Planetary Perturbations (if include_planets is true)
    bool perturb_mercury = true;
    bool perturb_venus = true;
    bool perturb_earth = true;
    bool perturb_mars = true;
    bool perturb_jupiter = true;
    bool perturb_saturn = true;
    bool perturb_uranus = true;
    bool perturb_neptune = true;
    bool perturb_pluto = false;

    std::string integrator_type = "RKF78";
    double initial_step_size = 0.5;
    double tolerance = 1e-12;

    // Non-Gravitational Forces
    bool include_yarkovsky = false;
    double yarkovsky_a2 = 0.0;
    
    // Differential Correction (Carpentry)
    double outlier_max_sigma = 10.0;
    double outlier_min_sigma = 3.0;
};

/**
 * @brief Configuration for the high‑precision asteroid fitting workflow.
 */
struct AsteroidFitConfig {
    // Local file paths (optional). Empty string means not provided.
    std::string eq1_file;
    std::string rwo_file;
    std::string oop_file; // optional configuration for AstDyn engine

    // URLs for remote resources (optional). If set, the file will be fetched.
    std::string eq1_url;
    std::string rwo_url;
    
    // Output file (optional) for saving results
    std::string output_file;
    
    // Engine configuration (overrides defaults and oop_file if specified)
    EngineConfig engine_config;

    // Orbital elements to use when files are not supplied.
    SimpleKeplerianElements orbit;

    // Observation epochs (MJD, UTC). Used when no .rwo file is supplied.
    std::vector<double> mjd_observations;

    // Output frame: true = Equatorial J2000, false = Ecliptic J2000.
    bool output_equatorial = true;
};

// JSON Parser Implementation
inline AsteroidFitConfig loadAsteroidFitConfig(const std::string& config_file) {
    AsteroidFitConfig config;
    std::ifstream f(config_file);
    if (!f.is_open()) {
        throw std::runtime_error("Could not open config file: " + config_file);
    }
    
    nlohmann::json j;
    f >> j;
    
    // Parse File Paths
    if (j.contains("eq1_file")) config.eq1_file = j["eq1_file"];
    if (j.contains("rwo_file")) config.rwo_file = j["rwo_file"];
    if (j.contains("oop_file")) config.oop_file = j["oop_file"];
    if (j.contains("output_file")) config.output_file = j["output_file"];
    
    // Parse Observations
    if (j.contains("mjd_observations")) {
        config.mjd_observations = j["mjd_observations"].get<std::vector<double>>();
    }
    
    // Parse Initial Orbit (if provided directly)
    if (j.contains("initial_orbit")) {
        auto& orb = j["initial_orbit"];
        config.orbit.a = orb.value("a", 0.0);
        config.orbit.e = orb.value("e", 0.0);
        config.orbit.i = orb.value("i", 0.0);
        config.orbit.Omega = orb.value("Omega", 0.0);
        config.orbit.omega = orb.value("omega", 0.0);
        config.orbit.M = orb.value("M", 0.0);
        config.orbit.epoch_mjd_tdb = orb.value("epoch_mjd_tdb", 0.0);
    }
    
    // Parse Output Options
    config.output_equatorial = j.value("output_equatorial", true);
    
    // Parse Engine Config
    if (j.contains("engine_config")) {
        auto& ec = j["engine_config"];
        config.engine_config.ephemeris_file = ec.value("ephemeris_file", "");
        config.engine_config.asteroid_ephemeris_file = ec.value("asteroid_ephemeris_file", "");
        config.engine_config.catalog_bias_file = ec.value("catalog_bias_file", "");
        config.engine_config.asteroid_mode = ec.value("asteroid_mode", "major_17");
        
        config.engine_config.include_planets = ec.value("include_planets", true);
        config.engine_config.include_relativity = ec.value("include_relativity", false);
        config.engine_config.include_moon = ec.value("include_moon", true);
        
        config.engine_config.perturb_mercury = ec.value("perturb_mercury", true);
        config.engine_config.perturb_venus = ec.value("perturb_venus", true);
        config.engine_config.perturb_earth = ec.value("perturb_earth", true);
        config.engine_config.perturb_mars = ec.value("perturb_mars", true);
        config.engine_config.perturb_jupiter = ec.value("perturb_jupiter", true);
        config.engine_config.perturb_saturn = ec.value("perturb_saturn", true);
        config.engine_config.perturb_uranus = ec.value("perturb_uranus", true);
        config.engine_config.perturb_neptune = ec.value("perturb_neptune", true);
        config.engine_config.perturb_pluto = ec.value("perturb_pluto", false);

        config.engine_config.include_yarkovsky = ec.value("include_yarkovsky", false);
        config.engine_config.yarkovsky_a2 = ec.value("yarkovsky_a2", 0.0);
        
        config.engine_config.outlier_max_sigma = ec.value("outlier_max_sigma", 10.0);
        config.engine_config.outlier_min_sigma = ec.value("outlier_min_sigma", 3.0);

        config.engine_config.integrator_type = ec.value("integrator_type", "RKF78");
        config.engine_config.initial_step_size = ec.value("initial_step_size", 0.5);
        config.engine_config.tolerance = ec.value("tolerance", 1e-12);
    }
    
    return config;
}

} // namespace ephemeris
} // namespace astdyn

#endif // ASTDYN_EPHEMERIS_ASTEROIDFITCONFIG_HPP
