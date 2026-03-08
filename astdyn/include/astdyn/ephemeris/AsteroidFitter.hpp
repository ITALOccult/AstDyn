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
#include "astdyn/time/epoch.hpp"
#include "src/types/vectors.hpp"
#include "src/core/frame_tags.hpp"
#include "src/core/units.hpp"
#include <iostream>
#include <fstream>

namespace astdyn {
namespace ephemeris {

/**
 * @brief Result of a high‑precision asteroid fit.
 */
struct AsteroidFitResult {
    bool success;
    std::string message;
    // Final fitted Keplerian elements (Equatorial J2000)
    physics::KeplerianStateTyped<core::GCRF> fitted_orbit;
    // RMS of residuals (arcsec)
    double rms_ra;
    double rms_dec;
    // Number of observations used after outlier rejection
    int num_observations;
    // Full list of heliocentric positions for each observation epoch
    std::vector<math::Vector3<core::GCRF, physics::Distance>> fitted_positions;
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
        
        if (!api_result.success || !api_result.fitted_state) {
            return result;
        }

        // Populate typed fitted_orbit
        const auto& state = *api_result.fitted_state;
        auto eq1_data = api::OrbitFitAPI::parse_eq1(eq1_file);
        
        result.fitted_orbit = physics::KeplerianStateTyped<core::GCRF>::from_traditional(
            eq1_data.epoch,
            state.a(),
            state.e(),
            state.i() * constants::RAD_TO_DEG,
            state.raan() * constants::RAD_TO_DEG,
            state.arg_peri() * constants::RAD_TO_DEG,
            state.m_anomaly() * constants::RAD_TO_DEG
        );
        
        result.rms_ra = api_result.rms_ra.value / 1000.0; // mas -> arcsec
        result.rms_dec = api_result.rms_dec.value / 1000.0;
        result.num_observations = api_result.num_observations;
        // Load observation epochs (optical observations)
        std::vector<astdyn::observations::OpticalObservation> records =
            astdyn::observations::RWOReader::readFile(rwo_file);
        std::vector<time::EpochUTC> t_obs;
        t_obs.reserve(records.size());
        for (const auto& rec : records) {
            t_obs.push_back(rec.time);
        }
        // Build PositionCalculator element struct from fitted orbit (Equatorial)
        KeplerianElements calc_elem;
        calc_elem.a = result.fitted_orbit.a.to_au();
        calc_elem.e = result.fitted_orbit.e;
        calc_elem.i = result.fitted_orbit.i.to_rad();
        calc_elem.Omega = result.fitted_orbit.node.to_rad();
        calc_elem.omega = result.fitted_orbit.omega.to_rad();
        calc_elem.M = result.fitted_orbit.M.to_rad();
        calc_elem.epoch = t_obs.empty() ? time::EpochTDB::from_mjd(0.0) : time::EpochTDB::from_mjd(t_obs.front().mjd());
        calc_elem.equatorial = true;
        result.fitted_positions.reserve(t_obs.size());
        for (const auto& t : t_obs) {
            auto pos = PositionCalculator::computePosition(
                calc_elem, time::EpochTDB::from_mjd(t.mjd()), true /* outputEquatorial */);
            result.fitted_positions.push_back(pos);
        }
        return result;
    }

    /**
     * @brief Compute positions from an in‑memory orbit configuration.
     */
    static AsteroidFitResult computeFromMemory(const astdyn::ephemeris::SimpleKeplerianElements& orbit,
                                               const std::vector<time::EpochUTC>& t_observations,
                                               bool output_equatorial = true) {
        AsteroidFitResult result;
        result.success = true;
        result.message = "Positions computed from in‑memory orbital elements.";
        // Populate typed fitted_orbit (GCRF assumed if from memory for now?)
        result.fitted_orbit = physics::KeplerianStateTyped<core::GCRF>::from_traditional(
            orbit.epoch, orbit.a, orbit.e, 
            orbit.i * constants::RAD_TO_DEG, orbit.Omega * constants::RAD_TO_DEG, 
            orbit.omega * constants::RAD_TO_DEG, orbit.M * constants::RAD_TO_DEG
        );
        result.rms_ra = 0.0;
        result.rms_dec = 0.0;
        result.num_observations = static_cast<int>(t_observations.size());
        result.fitted_positions.reserve(t_observations.size());
        // Build PositionCalculator element struct
        KeplerianElements calc_elem;
        calc_elem.a = orbit.a;
        calc_elem.e = orbit.e;
        calc_elem.i = orbit.i;
        calc_elem.Omega = orbit.Omega;
        calc_elem.omega = orbit.omega;
        calc_elem.M = orbit.M;
        calc_elem.epoch = orbit.epoch;
        
        for (const auto& t : t_observations) {
            auto pos = PositionCalculator::computePosition(
                calc_elem, time::EpochTDB::from_mjd(t.mjd()), output_equatorial);
            result.fitted_positions.push_back(pos);
        }
        return result;
    }

    /**
     * @brief Generate a temporary .json configuration file from EngineConfig.
     */
    static std::string generateTempConfig(const AsteroidFitConfig& cfg) {
        std::string temp_path = "temp_astdyn_fit.json";
        std::ofstream f(temp_path);
        if (!f) return "";
        
        const auto& ec = cfg.engine_config;
        nlohmann::json j;
        
        // Ephemeris
        if (!ec.ephemeris_file.empty()) {
            j["ephemeris"]["type"] = "DE441";
            j["ephemeris"]["file"] = ec.ephemeris_file;
        } else {
            j["ephemeris"]["type"] = "Analytical";
        }
        
        if (!ec.asteroid_ephemeris_file.empty()) {
            j["ephemeris"]["asteroid_file"] = ec.asteroid_ephemeris_file;
        }
        
        // Perturbations
        j["perturb"]["planets"] = ec.include_planets;
        j["perturb"]["relativity"] = ec.include_relativity;
        j["perturb"]["moon"] = ec.include_moon;
        
        // Non-Gravitational
        j["perturb"]["yarkovsky"] = ec.include_yarkovsky;
        j["perturb"]["yarkovsky_a2"] = ec.yarkovsky_a2;
        
        if (!ec.catalog_bias_file.empty()) {
            j["catalog"]["bias_file"] = ec.catalog_bias_file;
        }
        
        // Differential Correction
        j["diffcorr"]["outlier_max"] = ec.outlier_max_sigma;
        j["diffcorr"]["outlier_min"] = ec.outlier_min_sigma;
        
        // Detailed Planets
        if (ec.include_planets) {
             j["perturb"]["mercury"] = ec.perturb_mercury;
             j["perturb"]["venus"] = ec.perturb_venus;
             j["perturb"]["earth"] = ec.perturb_earth;
             j["perturb"]["mars"] = ec.perturb_mars;
             j["perturb"]["jupiter"] = ec.perturb_jupiter;
             j["perturb"]["saturn"] = ec.perturb_saturn;
             j["perturb"]["uranus"] = ec.perturb_uranus;
             j["perturb"]["neptune"] = ec.perturb_neptune;
        }
        
        // Asteroid Mode
        if (ec.asteroid_mode == "none") {
            j["perturb"]["asteroids"] = false;
        } else if (ec.asteroid_mode == "major_17") {
            j["perturb"]["asteroids"] = true;
        } else if (ec.asteroid_mode == "full_300") {
            j["perturb"]["asteroids"] = true;
        }
        
        // Integrator
        j["integrator"]["type"] = ec.integrator_type;
        j["integrator"]["step_size"] = ec.initial_step_size;
        j["integrator"]["tolerance"] = ec.tolerance;
        
        f << j.dump(4);
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
            if (result.success && !cfg.t_observations.empty()) {
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
                
                for (const auto& target_t : cfg.t_observations) {
                     auto state_at_target = engine.propagate_to(time::EpochTDB::from_mjd(target_t.mjd()));
                     Eigen::Vector3d target_pos;                    // Convert to PositionCalculator format
                     KeplerianElements calc_elem;
                     calc_elem.a = state_at_target.a.to_au();
                     calc_elem.e = state_at_target.e;
                     calc_elem.i = state_at_target.i.to_rad();
                     calc_elem.Omega = state_at_target.node.to_rad();
                     calc_elem.omega = state_at_target.omega.to_rad();
                     calc_elem.M = state_at_target.M.to_rad();
                     calc_elem.epoch = state_at_target.epoch;
                     calc_elem.equatorial = true;
                     
                     auto pos = PositionCalculator::computePosition(
                        calc_elem, time::EpochTDB::from_mjd(target_t.mjd()), cfg.output_equatorial);
                     
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
             
             auto initial_ecl = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
                 cfg.orbit.epoch, cfg.orbit.a, cfg.orbit.e, 
                 cfg.orbit.i * constants::RAD_TO_DEG, cfg.orbit.Omega * constants::RAD_TO_DEG, 
                 cfg.orbit.omega * constants::RAD_TO_DEG, cfg.orbit.M * constants::RAD_TO_DEG
             );
             auto cart_ecl = propagation::keplerian_to_cartesian(initial_ecl);
             auto pos_eq = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(cart_ecl.position);
             auto vel_eq = coordinates::ReferenceFrame::transform_vel<core::ECLIPJ2000, core::GCRF>(cart_ecl.position, cart_ecl.velocity);
             auto cart_eq = physics::CartesianStateTyped<core::GCRF>(cart_ecl.epoch, pos_eq, vel_eq, cart_ecl.gm);
             
             // AstDynEngine requires ECLIPJ2000 for its internal integrator consistency
             auto pos_ecl_back = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(cart_eq.position);
             auto vel_ecl_back = coordinates::ReferenceFrame::transform_vel<core::GCRF, core::ECLIPJ2000>(cart_eq.position, cart_eq.velocity);
             auto cart_ecl_back = physics::CartesianStateTyped<core::ECLIPJ2000>(cart_eq.epoch, pos_ecl_back, vel_ecl_back, cart_eq.gm);

             engine.set_initial_orbit(propagation::cartesian_to_keplerian<core::ECLIPJ2000>(cart_ecl_back));
             
             // 3. Propagate
             for (const auto& target_t : cfg.t_observations) {
                 auto state_at_target = engine.propagate_to(time::EpochTDB::from_mjd(target_t.mjd()));
                 
                  KeplerianElements calc_elem;
                  calc_elem.a = state_at_target.a.to_au();
                  calc_elem.e = state_at_target.e;
                  calc_elem.i = state_at_target.i.to_rad();
                  calc_elem.Omega = state_at_target.node.to_rad();
                  calc_elem.omega = state_at_target.omega.to_rad();
                  calc_elem.M = state_at_target.M.to_rad();
                  calc_elem.epoch = state_at_target.epoch;
                 calc_elem.equatorial = true; // Engine output is Equatorial
                 
                 auto pos = PositionCalculator::computePosition(
                    calc_elem, time::EpochTDB::from_mjd(target_t.mjd()), cfg.output_equatorial);
                 result.fitted_positions.push_back(pos);
             }
             return result;
        }

        return computeFromMemory(cfg.orbit, cfg.t_observations, cfg.output_equatorial);
    }
    
    /**
     * @brief Deprecated overload: fit directly from SimpleKeplerianElements.
     */
    static AsteroidFitResult fit(const astdyn::ephemeris::SimpleKeplerianElements& orbit,
                                 const std::vector<time::EpochUTC>& t_observations,
                                 bool output_equatorial = true) {
        return computeFromMemory(orbit, t_observations, output_equatorial);
    }
};

} // namespace ephemeris
} // namespace astdyn

#endif // ASTDYN_EPHEMERIS_ASTEROIDFITTER_HPP
