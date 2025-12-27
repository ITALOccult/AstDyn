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
            iss >> equ.epoch_mjd_tdb;
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

propagation::KeplerianElements OrbitFitAPI::convert_mean_equinoctial_to_osculating(
    const propagation::EquinoctialElements& mean_equ
) {
    // 1. Convert Equinoctial to Keplerian (Ecliptic)
    // Note: Assuming input is "mean" in the sense of averaged elements, but we treat them
    // geometrically for frame transformation.
    auto keplerian_ecliptic = propagation::equinoctial_to_keplerian(mean_equ);
    
    // 2. To Cartesian (Ecliptic J2000)
    auto cart_ecliptic = propagation::keplerian_to_cartesian(keplerian_ecliptic);
    
    // 3. Transform Frame: Ecliptic J2000 -> Equatorial J2000 (ICRF)
    coordinates::CartesianState state_ecliptic(
        cart_ecliptic.position,
        cart_ecliptic.velocity,
        cart_ecliptic.gravitational_parameter
    );
    
    auto state_equatorial = coordinates::ReferenceFrame::transform_state(
        state_ecliptic,
        coordinates::FrameType::ECLIPTIC,
        coordinates::FrameType::J2000
    );
    
    // 4. Convert back to Keplerian (Equatorial / Osculating)
    propagation::CartesianElements cart_equatorial;
    cart_equatorial.epoch_mjd_tdb = cart_ecliptic.epoch_mjd_tdb;
    cart_equatorial.position = state_equatorial.position();
    cart_equatorial.velocity = state_equatorial.velocity();
    cart_equatorial.gravitational_parameter = state_equatorial.mu();
    
    return propagation::cartesian_to_keplerian(cart_equatorial);
}

static Eigen::Vector3d ecl_to_eq_std(const Eigen::Vector3d& ecl) {
    return coordinates::ReferenceFrame::ecliptic_to_j2000() * ecl;
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
        if (verbose) std::cout << "[OrbitFitAPI] Loading OrbFit elements from " << eq1_file << "...\n";
        
        // 1. Parse Initial Orbit
        auto equ = parse_eq1(eq1_file);
        auto orbfit_orbit_ecliptic = propagation::equinoctial_to_keplerian(equ);

        // 2. Transform Frame: Ecliptic J2000 -> Equatorial J2000
        auto cart_ecliptic_elem = propagation::keplerian_to_cartesian(orbfit_orbit_ecliptic);

        OrbitContext ctx_in;
        ctx_in.frame = ReferenceFrame::ECLIPTIC_J2000;
        ctx_in.center = OriginCenter::SUN;
        ctx_in.model = OrbitModel::MEAN;
        ctx_in.format = ElementFormat::KEPLERIAN;

        propagation::CartesianElements cart_equatorial_elem;
        cart_equatorial_elem.epoch_mjd_tdb = cart_ecliptic_elem.epoch_mjd_tdb;
        cart_equatorial_elem.position = coordinates::ReferenceFrame::ecliptic_to_j2000() * cart_ecliptic_elem.position;
        cart_equatorial_elem.velocity = coordinates::ReferenceFrame::ecliptic_to_j2000() * cart_ecliptic_elem.velocity;
        cart_equatorial_elem.gravitational_parameter = cart_ecliptic_elem.gravitational_parameter;
        
        auto initial_orbit_equatorial = propagation::cartesian_to_keplerian(cart_equatorial_elem);
        
        OrbitContext ctx_fit;
        ctx_fit.frame = ReferenceFrame::EQUATORIAL_J2000;
        ctx_fit.center = OriginCenter::SUN;
        ctx_fit.model = OrbitModel::OSCULATING;
        ctx_fit.format = ElementFormat::KEPLERIAN;

        if (verbose) {
            std::cout << "   Initial Context: " << ctx_in.toString() << "\n";
            std::cout << "   Target Context:  " << ctx_fit.toString() << "\n";
            std::cout << "   Initial Orbit (Ecliptic): a=" << orbfit_orbit_ecliptic.semi_major_axis 
                      << " e=" << orbfit_orbit_ecliptic.eccentricity 
                      << " i=" << orbfit_orbit_ecliptic.inclination * constants::RAD_TO_DEG << "deg\n";
            std::cout << "   Initial Orbit (Transformed Equatorial): i=" 
                      << initial_orbit_equatorial.inclination * constants::RAD_TO_DEG << "deg\n";
        }

        // 3. Load Observations
        if (verbose) std::cout << "[OrbitFitAPI] Loading observations from " << rwo_file << "...\n";
        auto obs_list = observations::RWOReader::readFile(rwo_file);
        if (verbose) std::cout << "   Loaded " << obs_list.size() << " observations\n";

        // 4. Setup Engine
        AstDynEngine engine;
        
        if (!config_file.empty()) {
             std::ifstream cfg_check(config_file);
             if (cfg_check.good()) {
                 if (verbose) std::cout << "   Loading configuration from " << config_file << "...\n";
                 engine.load_config(config_file);
             } else {
                 if (verbose) std::cout << "   Warning: Config file " << config_file << " not found. Using defaults + planetary perturbations.\n";
             }
        }
        
        if (!de441_path.empty()) {
            if (verbose) std::cout << "   Overriding DE441 path: " << de441_path << "\n";
            AstDynConfig cfg = engine.config();
            cfg.ephemeris_type = "DE441";
            cfg.ephemeris_file = de441_path;
            engine.set_config(cfg);
        }
        
        // Ensure robust configuration if loaded one is minimal/missing
        AstDynConfig config = engine.config();
        
        // Force Planetary Perturbations (Sun+Planets)
        config.propagator_settings.include_planets = true;
        config.propagator_settings.perturb_mercury = true;
        config.propagator_settings.perturb_venus = true;
        config.propagator_settings.perturb_earth = true;
        config.propagator_settings.perturb_mars = true;
        config.propagator_settings.perturb_jupiter = true;
        config.propagator_settings.perturb_saturn = true;
        config.propagator_settings.perturb_uranus = true;
        config.propagator_settings.perturb_neptune = true;
        
        // Force Asteroid Perturbations (AST17 model) - DISABLED for stability defaults
        // This is critical for high precision over long arcs
        // config.propagator_settings.include_asteroids = true;
        
        // Force RKF78 for speed/accuracy balance
        config.integrator_type = "RKF78";
        config.initial_step_size = 0.5; // days
        config.tolerance = 1e-12;
        
        config.verbose = verbose;
        engine.set_config(config);
        
        engine.clear_observations();
        for (const auto& obs : obs_list) {
            engine.add_observation(obs);
        }
        
        // 5. Propagate to Meaningful Epoch (Mean observation time)
        // Helps convergence if the initial epoch is far from observations
        double sum_mjd = 0.0;
        for (const auto& obs : obs_list) sum_mjd += obs.mjd_utc;
        double mean_epoch = obs_list.empty() ? initial_orbit_equatorial.epoch_mjd_tdb : (sum_mjd / obs_list.size());
        
        engine.set_initial_orbit(initial_orbit_equatorial);
        auto start_orbit = engine.propagate_to(mean_epoch);
        engine.set_initial_orbit(start_orbit);
        
        // 6. Run Fit
        if (verbose) std::cout << "[OrbitFitAPI] Running differential correction...\n";
        auto engine_result = engine.fit_orbit();
        
        // 7. Populate Result
        result.success = true;
        result.message = "Fit completed";
        result.converged = engine_result.converged;
        result.fitted_orbit = engine_result.orbit;
        result.rms_ra = engine_result.rms_ra;
        result.rms_dec = engine_result.rms_dec;
        result.iterations = engine_result.num_iterations;
        result.num_observations = engine_result.num_observations;
        result.num_outliers = engine_result.num_rejected;
        
        // Comparison (Propagate initial orbit to result epoch)
        engine.set_initial_orbit(initial_orbit_equatorial);
        auto initial_at_result_epoch = engine.propagate_to(result.fitted_orbit.epoch_mjd_tdb);
        
        result.delta_a_km = (result.fitted_orbit.semi_major_axis - initial_at_result_epoch.semi_major_axis) * constants::AU;
        result.delta_e = result.fitted_orbit.eccentricity - initial_at_result_epoch.eccentricity;
        result.delta_i_arcsec = (result.fitted_orbit.inclination - initial_at_result_epoch.inclination) * constants::RAD_TO_ARCSEC;

    } catch (const std::exception& e) {
        result.success = false;
        result.message = e.what();
    }

    return result;
}

} // namespace api
} // namespace astdyn
