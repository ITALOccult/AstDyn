// AsteroidFitter.hpp
// High‑precision fitting utility that combines orbital elements, RWO observations
// and the existing AstDyn engine to produce a best‑fit orbit.
// It uses PositionCalculator for fast state propagation and OrbitFitAPI for the
// differential correction. All standard corrections (light‑time, stellar
// aberration, relativistic terms, EOP/UT1 handling) are applied via the AstDyn
// engine – this class merely orchestrates the workflow.

#ifndef ASTDYN_EPHEMERIS_ASTEROIDFITTER_HPP
#define ASTDYN_EPHEMERIS_ASTEROIDFITTER_HPP

#include <cstring>
#include <vector>
#include <Eigen/Dense>
#include <filesystem>
#include <cstdlib>
#include <stdexcept>

#include "astdyn/api/OrbitFitAPI.hpp"
#include "astdyn/observations/RWOReader.hpp"
#include "astdyn/ephemeris/PositionCalculator.hpp"
#include "astdyn/ephemeris/AsteroidFitConfig.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>

namespace astdyn {
namespace ephemeris {

/**
 * @brief Result of a high‑precision asteroid fit.
 */
struct AsteroidFitResult {
    bool success;
    std::string message;
    // Final fitted Keplerian elements (Equatorial J2000)
    propagation::KeplerianElements fitted_orbit;
    // RMS of residuals (arcsec)
    double rms_ra;
    double rms_dec;
    // Number of observations used after outlier rejection
    int num_observations;
    // Full list of heliocentric positions (AU) for each observation epoch
    std::vector<Eigen::Vector3d> fitted_positions;
};

/**
 * @brief High‑precision fitting helper.
 */
class AsteroidFitter {
public:
    /**
     * @brief Fit using .eq1 and .rwo files.
     */
    static AsteroidFitResult fit(const std::string& eq1_file,
                                 const std::string& rwo_file,
                                 const std::string& oop_file = "",
                                 bool verbose = true) {
        auto api_result = api::OrbitFitAPI::run_fit(eq1_file, rwo_file, oop_file, verbose);
        AsteroidFitResult result;
        result.success = api_result.success;
        result.message = api_result.message;
        result.fitted_orbit = api_result.fitted_orbit;
        result.rms_ra = api_result.rms_ra;
        result.rms_dec = api_result.rms_dec;
        result.num_observations = api_result.num_observations;
        if (!api_result.success) {
            return result;
        }
        // Load observation epochs (optical observations)
        std::vector<astdyn::observations::OpticalObservation> records =
            astdyn::observations::RWOReader::readFile(rwo_file);
        std::vector<double> mjd_obs;
        mjd_obs.reserve(records.size());
        for (const auto& rec : records) {
            mjd_obs.push_back(rec.mjd_utc);
        }
        // Build PositionCalculator element struct from fitted orbit (Equatorial)
        KeplerianElements calc_elem;
        calc_elem.a = result.fitted_orbit.semi_major_axis;
        calc_elem.e = result.fitted_orbit.eccentricity;
        calc_elem.i = result.fitted_orbit.inclination;
        calc_elem.Omega = result.fitted_orbit.longitude_ascending_node;
        calc_elem.omega = result.fitted_orbit.argument_perihelion;
        calc_elem.M = result.fitted_orbit.mean_anomaly;
        double epoch_mjd = mjd_obs.empty() ? 0.0 : mjd_obs.front();
        result.fitted_positions.reserve(mjd_obs.size());
        for (double mjd : mjd_obs) {
            Eigen::Vector3d pos = PositionCalculator::computePosition(
                calc_elem, mjd, epoch_mjd, true /* outputEquatorial */);
            result.fitted_positions.push_back(pos);
        }
        return result;
    }

    /**
     * @brief Compute positions from an in‑memory orbit configuration.
     */
    static AsteroidFitResult computeFromMemory(const astdyn::ephemeris::SimpleKeplerianElements& orbit,
                                               const std::vector<double>& mjd_observations,
                                               bool output_equatorial = true) {
        AsteroidFitResult result;
        result.success = true;
        result.message = "Positions computed from in‑memory orbital elements.";
        // Populate fitted_orbit (propagation namespace) from SimpleKeplerianElements
        result.fitted_orbit.semi_major_axis = orbit.a;
        result.fitted_orbit.eccentricity = orbit.e;
        result.fitted_orbit.inclination = orbit.i;
        result.fitted_orbit.longitude_ascending_node = orbit.Omega;
        result.fitted_orbit.argument_perihelion = orbit.omega;
        result.fitted_orbit.mean_anomaly = orbit.M;
        result.rms_ra = 0.0;
        result.rms_dec = 0.0;
        result.num_observations = static_cast<int>(mjd_observations.size());
        result.fitted_positions.reserve(mjd_observations.size());
        // Build PositionCalculator element struct
        KeplerianElements calc_elem;
        calc_elem.a = orbit.a;
        calc_elem.e = orbit.e;
        calc_elem.i = orbit.i;
        calc_elem.Omega = orbit.Omega;
        calc_elem.omega = orbit.omega;
        calc_elem.M = orbit.M;
        // Use the explicit epoch from the orbital elements
        double epoch_mjd = orbit.epoch_mjd_tdb;
        // Fallback or warning could be added if epoch is 0, but 0 is valid MJD (though unlikely for asteroid)
        if (epoch_mjd == 0.0 && !mjd_observations.empty()) {
             // If not set, maybe fallback to first observation? 
             // Ideally we should trust the user provided epoch.
             // Let's keep it simple: assume user sets it.
        }
        for (double mjd : mjd_observations) {
            Eigen::Vector3d pos = PositionCalculator::computePosition(
                calc_elem, mjd, epoch_mjd, output_equatorial);
            result.fitted_positions.push_back(pos);
        }
        return result;
    }

    /**
     * @brief Generate a temporary .oop configuration file from EngineConfig.
     */
    static std::string generateTempConfig(const AsteroidFitConfig& cfg) {
        std::string temp_path = "temp_astdyn_fit.oop";
        std::ofstream f(temp_path);
        if (!f) return "";
        
        const auto& ec = cfg.engine_config;
        
        // Ephemeris
        if (!ec.ephemeris_file.empty()) {
            f << "ephemeris.type = DE441\n";
            f << "ephemeris.file = " << ec.ephemeris_file << "\n";
        } else {
            f << "ephemeris.type = Analytical\n";
        }
        
        if (!ec.asteroid_ephemeris_file.empty()) {
            f << "ephemeris.asteroid_file = " << ec.asteroid_ephemeris_file << "\n";
        }
        
        // Global Perturbations
        f << "perturb.planets = " << (ec.include_planets ? "true" : "false") << "\n";
        f << "perturb.relativity = " << (ec.include_relativity ? "true" : "false") << "\n";
        f << "perturb.moon = " << (ec.include_moon ? "true" : "false") << "\n";
        
        // Non-Gravitational
        f << "perturb.yarkovsky = " << (ec.include_yarkovsky ? "true" : "false") << "\n";
        f << "perturb.yarkovsky_a2 = " << ec.yarkovsky_a2 << "\n";
        
        if (!ec.catalog_bias_file.empty()) {
            f << "catalog.bias_file = " << ec.catalog_bias_file << "\n";
        }
        
        // Carpentry
        f << "diffcorr.outlier_max = " << ec.outlier_max_sigma << "\n";
        f << "diffcorr.outlier_min = " << ec.outlier_min_sigma << "\n";
        
        // Detailed Planets
        if (ec.include_planets) {
             f << "perturb.mercury = " << (ec.perturb_mercury ? "true" : "false") << "\n";
             f << "perturb.venus = " << (ec.perturb_venus ? "true" : "false") << "\n";
             f << "perturb.earth = " << (ec.perturb_earth ? "true" : "false") << "\n";
             f << "perturb.mars = " << (ec.perturb_mars ? "true" : "false") << "\n";
             f << "perturb.jupiter = " << (ec.perturb_jupiter ? "true" : "false") << "\n";
             f << "perturb.saturn = " << (ec.perturb_saturn ? "true" : "false") << "\n";
             f << "perturb.uranus = " << (ec.perturb_uranus ? "true" : "false") << "\n";
             f << "perturb.neptune = " << (ec.perturb_neptune ? "true" : "false") << "\n";
        }
        
        // Asteroid Mode
        if (ec.asteroid_mode == "none") {
            f << "perturb.asteroids = false\n";
        } else if (ec.asteroid_mode == "major_17") {
            f << "perturb.asteroids = true\n";
        } else if (ec.asteroid_mode == "full_300") {
            f << "perturb.asteroids = true\n";
        }
        
        // Integrator
        f << "integrator.type = " << ec.integrator_type << "\n";
        f << "integrator.step_size = " << ec.initial_step_size << "\n";
        f << "integrator.tolerance = " << ec.tolerance << "\n";
        
        f.close();
        return temp_path;
    }

    /**
     * @brief Fit from a JSON configuration (file‑based or in‑memory).
     */
    static AsteroidFitResult fitFromConfig(const AsteroidFitConfig& cfg) {
        if (!cfg.eq1_file.empty() && !cfg.rwo_file.empty()) {
            // Determine which config to use for the FIT
            std::string oop_to_use = cfg.oop_file;
            std::string temp_oop = "";
            bool using_temp = false;
            
            // Generate temp config if engine_config provided
            if (!cfg.engine_config.ephemeris_file.empty() || cfg.engine_config.include_planets) {
                 temp_oop = generateTempConfig(cfg);
                 if (!temp_oop.empty()) {
                     oop_to_use = temp_oop;
                     using_temp = true;
                     // std::cout << "[DEBUG] Generated temporary config: " << temp_oop << std::endl;
                 }
            }

            AsteroidFitResult result = fit(cfg.eq1_file, cfg.rwo_file, oop_to_use);
            
            // Cleanup temp file
            if (using_temp && std::filesystem::exists(temp_oop)) {
                std::filesystem::remove(temp_oop);
            }
            
            // Post-fit propagation if requested and fit success
            if (result.success && !cfg.mjd_observations.empty()) {
                result.fitted_positions.clear();
                
                // Initialize Engine for propagation
                AstDynEngine engine;
                if (!cfg.oop_file.empty()) {
                    engine.load_config(cfg.oop_file);
                }
                
                // Override/Set configuration from EngineConfig if provided
                if (!cfg.engine_config.ephemeris_file.empty() || !cfg.engine_config.asteroid_ephemeris_file.empty() || cfg.engine_config.include_planets) {
                    AstDynConfig conf = engine.config();
                    if (!cfg.engine_config.ephemeris_file.empty()) {
                        conf.ephemeris_type = "DE441";
                        conf.ephemeris_file = cfg.engine_config.ephemeris_file;
                    }
                    
                    // PLANETS & RELATIVITY
                    conf.propagator_settings.include_planets = cfg.engine_config.include_planets;
                    conf.propagator_settings.include_relativity = cfg.engine_config.include_relativity;
                    conf.propagator_settings.include_moon = cfg.engine_config.include_moon;

                    conf.propagator_settings.perturb_mercury = cfg.engine_config.perturb_mercury;
                    conf.propagator_settings.perturb_venus = cfg.engine_config.perturb_venus;
                    conf.propagator_settings.perturb_earth = cfg.engine_config.perturb_earth;
                    conf.propagator_settings.perturb_mars = cfg.engine_config.perturb_mars;
                    conf.propagator_settings.perturb_jupiter = cfg.engine_config.perturb_jupiter;
                    conf.propagator_settings.perturb_saturn = cfg.engine_config.perturb_saturn;
                    conf.propagator_settings.perturb_uranus = cfg.engine_config.perturb_uranus;
                    conf.propagator_settings.perturb_neptune = cfg.engine_config.perturb_neptune;
                    
                    // ASTEROID MODE
                    if (cfg.engine_config.asteroid_mode == "none") {
                        conf.propagator_settings.include_asteroids = false;
                    } else if (cfg.engine_config.asteroid_mode == "major_17") {
                        conf.propagator_settings.include_asteroids = true;
                        conf.propagator_settings.asteroid_ephemeris_file = ""; // Use Internal AST17
                    } else if (cfg.engine_config.asteroid_mode == "full_300") {
                        conf.propagator_settings.include_asteroids = true;
                        conf.propagator_settings.asteroid_ephemeris_file = cfg.engine_config.asteroid_ephemeris_file;
                        // Also sync to config struct top level if needed
                        conf.asteroid_ephemeris_file = cfg.engine_config.asteroid_ephemeris_file;
                    }

                    conf.integrator_type = cfg.engine_config.integrator_type;
                    conf.initial_step_size = cfg.engine_config.initial_step_size;
                    conf.tolerance = cfg.engine_config.tolerance;
                    
                    engine.set_config(conf);
                }

                // Set the fitted orbit
                engine.set_initial_orbit(result.fitted_orbit); 
                
                for (double target_mjd : cfg.mjd_observations) {
                     // Propagate to target MJD
                     auto state_at_target = engine.propagate_to(target_mjd);
                     
                     // Convert to PositionCalculator format
                     KeplerianElements calc_elem;
                     calc_elem.a = state_at_target.semi_major_axis;
                     calc_elem.e = state_at_target.eccentricity;
                     calc_elem.i = state_at_target.inclination;
                     calc_elem.Omega = state_at_target.longitude_ascending_node;
                     calc_elem.omega = state_at_target.argument_perihelion;
                     calc_elem.M = state_at_target.mean_anomaly;
                     
                     calc_elem.equatorial = true;
                     
                     Eigen::Vector3d pos = PositionCalculator::computePosition(
                        calc_elem, target_mjd, target_mjd, cfg.output_equatorial);
                     
                     result.fitted_positions.push_back(pos);
                }
            }
            return result;
        }

        // Check if High Precision Engine is requested
        if (!cfg.engine_config.ephemeris_file.empty() || !cfg.engine_config.asteroid_ephemeris_file.empty()) {
             AsteroidFitResult result;
             result.success = true; // Propagation is successful by default
             result.message = "High-precision propagation (No Fit).";
             result.num_observations = 0;
             result.rms_ra = 0.0;
             result.rms_dec = 0.0;
             
             // Setup Engine
             AstDynEngine engine;
             if (!cfg.oop_file.empty()) {
                 engine.load_config(cfg.oop_file);
             }
             
             AstDynConfig conf = engine.config();
             if (!cfg.engine_config.ephemeris_file.empty()) {
                 conf.ephemeris_type = "DE441";
                 conf.ephemeris_file = cfg.engine_config.ephemeris_file;
             }
             
             // PLANETS & RELATIVITY
             conf.propagator_settings.include_planets = cfg.engine_config.include_planets; 
             conf.propagator_settings.include_relativity = cfg.engine_config.include_relativity;
             conf.propagator_settings.include_moon = cfg.engine_config.include_moon;

             conf.propagator_settings.perturb_mercury = cfg.engine_config.perturb_mercury;
             conf.propagator_settings.perturb_venus = cfg.engine_config.perturb_venus;
             conf.propagator_settings.perturb_earth = cfg.engine_config.perturb_earth;
             conf.propagator_settings.perturb_mars = cfg.engine_config.perturb_mars;
             conf.propagator_settings.perturb_jupiter = cfg.engine_config.perturb_jupiter;
             conf.propagator_settings.perturb_saturn = cfg.engine_config.perturb_saturn;
             conf.propagator_settings.perturb_uranus = cfg.engine_config.perturb_uranus;
             conf.propagator_settings.perturb_neptune = cfg.engine_config.perturb_neptune;
             
             // ASTEROID MODE
             if (cfg.engine_config.asteroid_mode == "none") {
                 conf.propagator_settings.include_asteroids = false;
             } else if (cfg.engine_config.asteroid_mode == "major_17") {
                 conf.propagator_settings.include_asteroids = true;
                 conf.propagator_settings.asteroid_ephemeris_file = ""; // Use Internal AST17
             } else if (cfg.engine_config.asteroid_mode == "full_300") {
                 conf.propagator_settings.include_asteroids = true;
                 conf.propagator_settings.asteroid_ephemeris_file = cfg.engine_config.asteroid_ephemeris_file;
                 conf.asteroid_ephemeris_file = cfg.engine_config.asteroid_ephemeris_file;
             }

             conf.integrator_type = cfg.engine_config.integrator_type;
             conf.initial_step_size = cfg.engine_config.initial_step_size;
             conf.tolerance = cfg.engine_config.tolerance;
             engine.set_config(conf);
             
             // Set Initial Orbit (Convert Simple -> Keplerian)
             // Note: SimpleKeplerianElements are usually Ecliptic. 
             // AstDynEngine usually takes Equatorial. 
             // We must transform it!
             
             // 1. Fill Keplerian (Ecliptic)
             propagation::KeplerianElements kep_ecl;
             kep_ecl.semi_major_axis = cfg.orbit.a;
             kep_ecl.eccentricity = cfg.orbit.e;
             kep_ecl.inclination = cfg.orbit.i;
             kep_ecl.longitude_ascending_node = cfg.orbit.Omega;
             kep_ecl.argument_perihelion = cfg.orbit.omega;
             kep_ecl.mean_anomaly = cfg.orbit.M;
             kep_ecl.epoch_mjd_tdb = cfg.orbit.epoch_mjd_tdb;
             kep_ecl.gravitational_parameter = constants::GMS;
             
             // 2. Transform to Equatorial J2000
             auto cart_ecl = propagation::keplerian_to_cartesian(kep_ecl);
             coordinates::CartesianState state_ecl(cart_ecl.position, cart_ecl.velocity, cart_ecl.gravitational_parameter);
             
             auto state_eq = coordinates::ReferenceFrame::transform_state(state_ecl, coordinates::FrameType::ECLIPTIC, coordinates::FrameType::J2000);
             
             propagation::CartesianElements ce;
             ce.epoch_mjd_tdb = kep_ecl.epoch_mjd_tdb;
             ce.position = state_eq.position();
             ce.velocity = state_eq.velocity();
             ce.gravitational_parameter = state_eq.mu();
             propagation::KeplerianElements kep_eq = propagation::cartesian_to_keplerian(ce);
             
             engine.set_initial_orbit(kep_eq);
             
             // 3. Propagate
             for (double target_mjd : cfg.mjd_observations) {
                 auto state_at_target = engine.propagate_to(target_mjd);
                 
                 KeplerianElements calc_elem;
                 calc_elem.a = state_at_target.semi_major_axis;
                 calc_elem.e = state_at_target.eccentricity;
                 calc_elem.i = state_at_target.inclination;
                 calc_elem.Omega = state_at_target.longitude_ascending_node;
                 calc_elem.omega = state_at_target.argument_perihelion;
                 calc_elem.M = state_at_target.mean_anomaly;
                 calc_elem.equatorial = true; // Engine output is Equatorial
                 
                 Eigen::Vector3d pos = PositionCalculator::computePosition(
                    calc_elem, target_mjd, target_mjd, cfg.output_equatorial);
                 result.fitted_positions.push_back(pos);
             }
             return result;
        }

        return computeFromMemory(cfg.orbit, cfg.mjd_observations, cfg.output_equatorial);
    }
    
    /**
     * @brief Deprecated overload: fit directly from SimpleKeplerianElements.
     */
    static AsteroidFitResult fit(const astdyn::ephemeris::SimpleKeplerianElements& orbit,
                                 const std::vector<double>& mjd_observations,
                                 bool output_equatorial = true) {
        return computeFromMemory(orbit, mjd_observations, output_equatorial);
    }
};

} // namespace ephemeris
} // namespace astdyn

#endif // ASTDYN_EPHEMERIS_ASTEROIDFITTER_HPP
