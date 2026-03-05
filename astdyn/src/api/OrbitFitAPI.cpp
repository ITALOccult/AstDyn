/**
 * @file OrbitFitAPI.cpp
 * @brief Implementation of OrbitFitAPI
 * @author ITALOccult AstDyn Team
 * @date 2025-12-13
 */

#include "astdyn/api/OrbitFitAPI.hpp"
#include "astdyn/observations/RWOReader.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include "src/core/frame_tags.hpp"
#include "src/types/orbital_state.hpp"
#include "src/coordinates/state_conversions.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace astdyn {
namespace api {

using namespace astdyn::propagation;

std::string OrbitFitAPI::ltrim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t");
    return (start == std::string::npos) ? "" : s.substr(start);
}

propagation::EquinoctialElements OrbitFitAPI::parse_eq1(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    propagation::EquinoctialElements equ;
    equ.epoch = utils::Instant::from_tt(utils::ModifiedJulianDate(0.0));
    std::string line;
    bool found_equ = false;
    bool found_mjd = false;

    while (std::getline(file, line)) {
        std::string trimmed = ltrim(line);
        if (trimmed.empty() || trimmed[0] == '!') continue;

        if (trimmed.substr(0, 3) == "EQU") {
            std::istringstream iss(trimmed.substr(3));
            iss >> equ.a >> equ.h >> equ.k >> equ.p >> equ.q >> equ.lambda;
            // Convert mean longitude from degrees to radians
            equ.lambda *= constants::PI / 180.0;
            found_equ = true;
        }
        else if (trimmed.substr(0, 3) == "MJD") {
            std::istringstream iss(trimmed.substr(3));
            double mjd_val;
            iss >> mjd_val;
            equ.epoch = utils::Instant::from_tt(utils::ModifiedJulianDate(mjd_val));
            found_mjd = true;
        }

        if (found_equ && found_mjd) break;
    }

    if (!found_equ || !found_mjd) {
        throw std::runtime_error("Could not parse equinoctial elements from file: " + filepath);
    }

    // Standardize to new GMS unit (AU^3/day^2) for internal consistency
    equ.gravitational_parameter = constants::GMS; 
    return equ;
}

types::OrbitalState<core::GCRF, types::KeplerianTag> 
OrbitFitAPI::prepare_initial_state(const propagation::EquinoctialElements& mean_equ) {
    // 1. Convert Equinoctial (Mean, Ecliptic) to Keplerian (Mean, Ecliptic)
    auto kep_ecl = propagation::equinoctial_to_keplerian(mean_equ);
    
    // 2. Convert to Cartesian (Ecliptic) for rotation
    auto cart_ecl = propagation::keplerian_to_cartesian(kep_ecl);
    
    // 3. Rotate Cartesian: Ecliptic -> Equatorial (GCRF)
    auto mat = coordinates::ReferenceFrame::ecliptic_to_j2000();
    auto final_pos = mat * cart_ecl.position.to_eigen();
    auto final_vel = mat * cart_ecl.velocity.to_eigen();
    
    // 4. Create a legacy Cartesian state in Equatorial frame to use existing conversion logic
    propagation::CartesianElements cart_eq;
    cart_eq.epoch = mean_equ.epoch;
    cart_eq.position = types::Vector3<core::GCRF, core::Meter>(final_pos);
    cart_eq.velocity = types::Vector3<core::GCRF, core::Meter>(final_vel);
    cart_eq.gravitational_parameter = kep_ecl.gravitational_parameter;
    
    // 5. Convert back to Keplerian (Equatorial / GCRF)
    auto kep_eq = propagation::cartesian_to_keplerian(cart_eq);
    
    // 6. Return typed GCRF OrbitalState
    return types::OrbitalState<core::GCRF, types::KeplerianTag>(std::array<double, 6>{
        kep_eq.semi_major_axis,
        kep_eq.eccentricity,
        kep_eq.inclination,
        kep_eq.longitude_ascending_node,
        kep_eq.argument_perihelion,
        kep_eq.mean_anomaly
    });
}

OrbitFitResult OrbitFitAPI::run_fit(
    const std::string& eq1_file,
    const std::string& rwo_file,
    const std::string& config_file,
    bool verbose,
    const std::string& de441_path
) {
    OrbitFitResult result;
    result.success = false;

    try {
        if (verbose) std::cout << "[OrbitFitAPI] Loading elements from " << eq1_file << "...\n";
        auto equ = parse_eq1(eq1_file);
        auto initial_state = OrbitFitAPI::prepare_initial_state(equ);

        auto obs_list = observations::RWOReader::readFile(rwo_file);
        if (verbose) std::cout << "   Loaded " << obs_list.size() << " observations\n";

        AstDynEngine engine;
        if (!de441_path.empty()) {
            AstDynConfig cfg = engine.config();
            cfg.ephemeris_type = "DE441";
            cfg.ephemeris_file = de441_path;
            engine.set_config(cfg);
        }
        
        // CTFYH: Ensure planetary perturbations are enabled by default for precision
        auto config = engine.config();
        config.propagator_settings.include_planets = true;
        config.integrator_type = "RKF78";
        config.verbose = verbose;
        engine.set_config(config);

        for (const auto& obs : obs_list) engine.add_observation(obs);

        // Convert OrbitalState to legacy KeplerianElements for current engine compatibility
        propagation::KeplerianElements legacy_init;
        legacy_init.semi_major_axis = initial_state.a();
        legacy_init.eccentricity = initial_state.e();
        legacy_init.inclination = initial_state.i();
        legacy_init.longitude_ascending_node = initial_state.raan();
        legacy_init.argument_perihelion = initial_state.arg_peri();
        legacy_init.mean_anomaly = initial_state.m_anomaly();
        legacy_init.epoch = equ.epoch;
        
        engine.set_initial_orbit(legacy_init);
        auto engine_result = engine.fit_orbit();
        
        result.success = true;
        result.converged = engine_result.converged;
        
        // Pack into OrbitalState Result
        result.fitted_state = types::OrbitalState<core::GCRF, types::KeplerianTag>(std::array<double, 6>{
            engine_result.orbit.semi_major_axis,
            engine_result.orbit.eccentricity,
            engine_result.orbit.inclination,
            engine_result.orbit.longitude_ascending_node,
            engine_result.orbit.argument_perihelion,
            engine_result.orbit.mean_anomaly
        });
        
        result.rms_ra = core::MilliArcSecond(engine_result.rms_ra * 1000.0);
        result.rms_dec = core::MilliArcSecond(engine_result.rms_dec * 1000.0);
        result.num_observations = engine_result.num_observations;
        result.iterations = engine_result.num_iterations;

        // Comparison against initial guess
        result.delta_a_km = std::abs(engine_result.orbit.semi_major_axis - initial_state.a()) * constants::AU;
        result.delta_e = std::abs(engine_result.orbit.eccentricity - initial_state.e());

    } catch (const std::exception& e) {
        result.success = false;
        result.message = e.what();
    }
    return result;
}

} // namespace api
} // namespace astdyn
