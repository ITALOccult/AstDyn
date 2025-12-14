/**
 * @file simple_pompeja_fit.cpp
 * @brief Simple example showing the complete orbit fitting workflow
 * 
 * This demonstrates the essential steps:
 * 1. Parse OrbFit .eq1 file (equinoctial elements)
 * 2. Load .rwo observations  
 * 3. Run differential correction
 * 4. Display results
 */

#include <astdyn/AstDynEngine.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/coordinates/CartesianState.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <numeric>
#include <algorithm>

using namespace astdyn;

// Helper to trim leading whitespace
std::string ltrim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t");
    return (start == std::string::npos) ? "" : s.substr(start);
}

// Parse OrbFit .eq1 file into propagation::EquinoctialElements
propagation::EquinoctialElements parse_eq1_file(const std::string& filepath) {
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
        throw std::runtime_error("Could not parse equinoctial elements from file");
    }

    equ.gravitational_parameter = constants::GMS;
    return equ;
}

int main(int argc, char** argv) {
    try {
        std::string eq1_file = "17030_astdys.eq1";
        std::string rwo_file = "17030_astdys.rwo";
        std::string oop_file = "17030.oop"; // Use the existing config file


        // Step 1: Parse OrbFit elements
        std::cout << "1. Loading OrbFit elements from " << eq1_file << "...\n";
        auto equ = parse_eq1_file(eq1_file);
        
        // Convert to Keplerian (Ecliptic J2000)
        auto orbfit_orbit_ecliptic = propagation::equinoctial_to_keplerian(equ);

        std::cout << "   Epoch: MJD " << std::fixed << std::setprecision(1) << orbfit_orbit_ecliptic.epoch_mjd_tdb << "\n";
        std::cout << "   a = " << std::setprecision(6) << orbfit_orbit_ecliptic.semi_major_axis << " AU\n";
        std::cout << "   e = " << std::setprecision(6) << orbfit_orbit_ecliptic.eccentricity << "\n";
        std::cout << "   i = " << std::setprecision(6) << orbfit_orbit_ecliptic.inclination * 180.0 / constants::PI << "° (Ecliptic)\n\n";

        // Convert Ecliptic J2000 -> Equatorial J2000
        std::cout << "   Converting Frame: Ecliptic J2000 -> Equatorial J2000\n";
        
        // 1. Keplerian (Ecliptic) -> Cartesian (Ecliptic)
        auto cart_ecliptic_elem = propagation::keplerian_to_cartesian(orbfit_orbit_ecliptic);
        
        // Convert to CartesianState for rotation
        coordinates::CartesianState state_ecliptic(
            cart_ecliptic_elem.position,
            cart_ecliptic_elem.velocity,
            cart_ecliptic_elem.gravitational_parameter
        );
        
        // 2. Rotate Cartesian State: Ecliptic -> Equatorial
        auto state_equatorial = coordinates::ReferenceFrame::transform_state(
            state_ecliptic,
            coordinates::FrameType::ECLIPTIC,
            coordinates::FrameType::J2000
        );
        
        // Convert back to CartesianElements
        propagation::CartesianElements cart_equatorial_elem;
        cart_equatorial_elem.epoch_mjd_tdb = cart_ecliptic_elem.epoch_mjd_tdb;
        cart_equatorial_elem.position = state_equatorial.position();
        cart_equatorial_elem.velocity = state_equatorial.velocity();
        cart_equatorial_elem.gravitational_parameter = state_equatorial.mu();
        
        // 3. Cartesian (Equatorial) -> Keplerian (Equatorial)
        auto orbfit_orbit_equatorial = propagation::cartesian_to_keplerian(cart_equatorial_elem);
        
        std::cout << "   Transformed Orbit (Equatorial):\n";
        std::cout << "   i = " << std::setprecision(6) << orbfit_orbit_equatorial.inclination * 180.0 / constants::PI << "°\n";
        std::cout << "   Ω = " << std::setprecision(6) << orbfit_orbit_equatorial.longitude_ascending_node * 180.0 / constants::PI << "°\n\n";

        // Step 2: Load observations
        std::cout << "2. Loading observations from " << rwo_file << "...\n";
        auto obs_list = observations::RWOReader::readFile(rwo_file); 


        // Step 3: Setup AstDyn engine
        std::cout << "3. Setting up AstDyn engine...\n";
        AstDynEngine engine;
        
        std::ifstream oop_check(oop_file);
        if (oop_check.good()) {
            std::cout << "   Loading configuration from " << oop_file << "...\n";
            engine.load_config(oop_file);
        } else {
            std::cout << "   Warning: " << oop_file << " not found, using default configuration.\n";
            std::cout << "   Forcing planetary perturbations ON (Sun+Planets).\n";
            AstDynConfig config = engine.config();
            config.propagator_settings.include_planets = true;
            config.propagator_settings.perturb_mercury = true;
            config.propagator_settings.perturb_venus = true;
            config.propagator_settings.perturb_earth = true;
            config.propagator_settings.perturb_mars = true;
            config.propagator_settings.perturb_jupiter = true;
            config.propagator_settings.perturb_saturn = true;
            config.propagator_settings.perturb_uranus = true;
            config.propagator_settings.perturb_neptune = true;
            engine.set_config(config);
        }

        // Force RKF78 Integrator for speed (RK4 is too slow for 100 years of data)
        AstDynConfig eco_config = engine.config();
        eco_config.integrator_type = "RKF78";
        eco_config.initial_step_size = 1.0; 
        eco_config.tolerance = 1e-12;
        engine.set_config(eco_config);
        
        // Add observations
        engine.clear_observations();
        for (const auto& obs : obs_list) {
            engine.add_observation(obs);
        }

        // Calculate mean epoch
        double sum_mjd = 0.0;
        for (const auto& obs : obs_list) {
            sum_mjd += obs.mjd_utc; 
        }
        double mean_epoch = obs_list.empty() ? orbfit_orbit_equatorial.epoch_mjd_tdb : (sum_mjd / obs_list.size());
        
        std::cout << "   Current Working Directory: " << "(current path)" << "\n";
        
        // Debug: Check Earth Ephemeris
        double check_epoch = 2451545.0; // J2000
        try {
            auto earth = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, check_epoch);
            std::cout << "   DEBUG: Earth at J2000: " << earth.position().transpose() << " AU\n";
            std::cout << "   DEBUG: Earth Radius: " << earth.position().norm() << " AU\n";
            if (earth.position().norm() < 0.1) std::cout << "   CRITICAL WARNING: Earth position near zero! Ephemeris likely not loaded.\n";
        } catch (const std::exception& e) {
             std::cout << "   CRITICAL ERROR: Ephemeris lookup failed: " << e.what() << "\n";
        }

        // Set initial orbit (Equatorial)
        engine.set_initial_orbit(orbfit_orbit_equatorial);

        std::cout << "4. Propagating orbit to mean epoch...\n";
        auto initial_state_at_mean = engine.propagate_to(mean_epoch);
        engine.set_initial_orbit(initial_state_at_mean);
        std::cout << "   Done.\n\n";

        // Step 4: Run differential correction
        std::cout << "5. Running differential correction...\n";
        
        AstDynConfig config = engine.config();
        config.max_iterations = 20;
        config.convergence_threshold = 1e-6;
        config.outlier_sigma = 3.0;
        config.verbose = true;
        engine.set_config(config);

        auto result = engine.fit_orbit();

        // Step 5: Display results
        std::cout << "\n╔═══════════════════════════════════════════════════════════╗\n";
        std::cout << "║                    RESULTS                                ║\n";
        std::cout << "╚═══════════════════════════════════════════════════════════╝\n\n";

        std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
        std::cout << "Iterations: " << result.num_iterations << "\n";
        std::cout << "RMS RA: " << std::fixed << std::setprecision(3) << result.rms_ra << " arcsec\n";
        std::cout << "RMS Dec: " << std::fixed << std::setprecision(3) << result.rms_dec << " arcsec\n";
        std::cout << "Observations used: " << result.num_observations << " / " << obs_list.size() << "\n";
        std::cout << "Outliers rejected: " << result.num_rejected << "\n\n";

        // Fitted orbit
        auto fitted_orbit = result.orbit;
        std::cout << "Fitted orbit at MJD " << std::setprecision(2) << fitted_orbit.epoch_mjd_tdb << ":\n";
        std::cout << "  a = " << std::setprecision(6) << fitted_orbit.semi_major_axis << " AU\n";
        std::cout << "  e = " << std::setprecision(6) << fitted_orbit.eccentricity << "\n";
        std::cout << "  i = " << std::setprecision(4) << fitted_orbit.inclination * 180.0 / constants::PI << "°\n";
        std::cout << "  Ω = " << std::setprecision(4) << fitted_orbit.longitude_ascending_node * 180.0 / constants::PI << "°\n";
        std::cout << "  ω = " << std::setprecision(4) << fitted_orbit.argument_perihelion * 180.0 / constants::PI << "°\n";
        std::cout << "  M = " << std::setprecision(4) << fitted_orbit.mean_anomaly * 180.0 / constants::PI << "°\n\n";

        // Comparison
        engine.set_initial_orbit(orbfit_orbit_equatorial);
        auto orbfit_at_result_epoch = engine.propagate_to(fitted_orbit.epoch_mjd_tdb);
        
        double delta_a_km = (fitted_orbit.semi_major_axis - orbfit_at_result_epoch.semi_major_axis) * constants::AU / 1000.0;
        double delta_e = fitted_orbit.eccentricity - orbfit_at_result_epoch.eccentricity;
        double delta_i_arcsec = (fitted_orbit.inclination - orbfit_at_result_epoch.inclination) * 180.0 * 3600.0 / constants::PI;

        std::cout << "Comparison with OrbFit (at result epoch):\n";
        std::cout << "  Δa = " << std::fixed << std::setprecision(2) << delta_a_km << " km\n";
        std::cout << "  Δe = " << std::scientific << std::setprecision(3) << delta_e << "\n";
        std::cout << "  Δi = " << std::fixed << std::setprecision(2) << delta_i_arcsec << " arcsec\n\n";

        std::cout << "╔═══════════════════════════════════════════════════════════╗\n";
        std::cout << "║                COMPLETED SUCCESSFULLY                     ║\n";
        std::cout << "╚═══════════════════════════════════════════════════════════╝\n";

        return result.converged ? 0 : 1;

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
