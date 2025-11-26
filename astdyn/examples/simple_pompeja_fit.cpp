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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace astdyn;

// Structure for equinoctial elements
struct EquinoctialElements {
    double epoch_mjd_tdb;
    double semi_major_axis;  // a
    double h;                // e·sin(ϖ)
    double k;                // e·cos(ϖ)
    double p;                // tan(i/2)·sin(Ω)
    double q;                // tan(i/2)·cos(Ω)
    double mean_longitude;   // λ (radians)
};

// Parse OrbFit .eq1 file
EquinoctialElements parse_eq1_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }

    EquinoctialElements equ;
    std::string line;
    bool found_equ = false;
    bool found_mjd = false;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '!') continue;

        if (line.substr(0, 3) == "EQU") {
            std::istringstream iss(line.substr(3));
            iss >> equ.semi_major_axis >> equ.h >> equ.k >> equ.p >> equ.q >> equ.mean_longitude;
            // Convert mean longitude from degrees to radians
            equ.mean_longitude *= constants::PI / 180.0;
            found_equ = true;
        }
        else if (line.substr(0, 3) == "MJD") {
            std::istringstream iss(line.substr(3));
            iss >> equ.epoch_mjd_tdb;
            found_mjd = true;
        }

        if (found_equ && found_mjd) break;
    }

    if (!found_equ || !found_mjd) {
        throw std::runtime_error("Could not parse equinoctial elements from file");
    }

    return equ;
}

// Convert equinoctial to Keplerian
orbit::KeplerianOrbit equinoctial_to_keplerian(const EquinoctialElements& equ) {
    orbit::KeplerianOrbit kep;

    kep.semi_major_axis = equ.semi_major_axis;
    kep.epoch_mjd_tdb = equ.epoch_mjd_tdb;

    // Eccentricity
    kep.eccentricity = std::sqrt(equ.h * equ.h + equ.k * equ.k);

    // Inclination
    double tan_half_i = std::sqrt(equ.p * equ.p + equ.q * equ.q);
    kep.inclination = 2.0 * std::atan(tan_half_i);

    // Longitude of ascending node
    kep.longitude_ascending_node = std::atan2(equ.p, equ.q);
    if (kep.longitude_ascending_node < 0) {
        kep.longitude_ascending_node += 2.0 * constants::PI;
    }

    // Argument of perihelion
    double omega_plus_Omega = std::atan2(equ.h, equ.k);
    if (omega_plus_Omega < 0) {
        omega_plus_Omega += 2.0 * constants::PI;
    }
    kep.argument_perihelion = omega_plus_Omega - kep.longitude_ascending_node;
    if (kep.argument_perihelion < 0) {
        kep.argument_perihelion += 2.0 * constants::PI;
    }

    // Mean anomaly
    kep.mean_anomaly = equ.mean_longitude - omega_plus_Omega;
    while (kep.mean_anomaly < 0) kep.mean_anomaly += 2.0 * constants::PI;
    while (kep.mean_anomaly >= 2.0 * constants::PI) kep.mean_anomaly -= 2.0 * constants::PI;

    return kep;
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 4) {
            std::cerr << "Usage: " << argv[0] << " <observations.rwo> <elements.eq1> <config.oop>\n\n";
            std::cerr << "Example:\n";
            std::cerr << "  " << argv[0] << " astdyn/tools/203_astdys_recent100.rwo "
                      << "astdyn/tools/203_astdys_latest.eq1 astdyn/tools/203.oop\n";
            return 1;
        }

        std::string rwo_file = argv[1];
        std::string eq1_file = argv[2];
        std::string oop_file = argv[3];

        std::cout << "╔═══════════════════════════════════════════════════════════╗\n";
        std::cout << "║          Simple Pompeja Orbit Fitting Example            ║\n";
        std::cout << "╚═══════════════════════════════════════════════════════════╝\n\n";

        // Step 1: Parse OrbFit elements
        std::cout << "1. Loading OrbFit elements from " << eq1_file << "...\n";
        auto equ = parse_eq1_file(eq1_file);
        auto orbfit_orbit = equinoctial_to_keplerian(equ);

        std::cout << "   Epoch: MJD " << std::fixed << std::setprecision(1) << orbfit_orbit.epoch_mjd_tdb << "\n";
        std::cout << "   a = " << std::setprecision(6) << orbfit_orbit.semi_major_axis << " AU\n";
        std::cout << "   e = " << std::setprecision(6) << orbfit_orbit.eccentricity << "\n";
        std::cout << "   i = " << std::setprecision(6) << orbfit_orbit.inclination * 180.0 / constants::PI << "°\n\n";

        // Step 2: Load observations
        std::cout << "2. Loading observations from " << rwo_file << "...\n";
        observations::RWOReader reader;
        auto obs_list = reader.read_file(rwo_file);
        std::cout << "   Loaded " << obs_list.size() << " observations\n\n";

        // Step 3: Setup AstDyn engine
        std::cout << "3. Setting up AstDyn engine...\n";
        AstDynEngine engine(oop_file);
        engine.set_observations(obs_list);

        // Get mean observation epoch
        double mean_epoch = engine.get_mean_observation_epoch();
        std::cout << "   Mean observation epoch: MJD " << std::fixed << std::setprecision(2) << mean_epoch << "\n\n";

        // Propagate OrbFit orbit to mean epoch
        std::cout << "4. Propagating orbit to mean epoch...\n";
        auto initial_state = engine.propagate_orbit(orbfit_orbit, mean_epoch);
        engine.set_initial_state(initial_state);
        std::cout << "   Done.\n\n";

        // Step 4: Run differential correction
        std::cout << "5. Running differential correction...\n";
        orbit_determination::DifferentialCorrectorSettings settings;
        settings.max_iterations = 20;
        settings.convergence_tolerance = 1e-6;
        settings.reject_outliers = true;
        settings.outlier_sigma = 3.0;

        auto result = engine.run_differential_correction(settings);

        // Step 5: Display results
        std::cout << "\n╔═══════════════════════════════════════════════════════════╗\n";
        std::cout << "║                    RESULTS                                ║\n";
        std::cout << "╚═══════════════════════════════════════════════════════════╝\n\n";

        std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
        std::cout << "Iterations: " << result.iterations << "\n";
        std::cout << "RMS: " << std::fixed << std::setprecision(3) << result.rms_arcsec << " arcsec\n";
        std::cout << "Observations used: " << result.num_observations_used << " / " << obs_list.size() << "\n";
        std::cout << "Outliers rejected: " << (obs_list.size() - result.num_observations_used) << "\n\n";

        // Fitted orbit
        auto fitted_orbit = result.final_orbit;
        std::cout << "Fitted orbit at MJD " << std::setprecision(2) << fitted_orbit.epoch_mjd_tdb << ":\n";
        std::cout << "  a = " << std::setprecision(6) << fitted_orbit.semi_major_axis << " AU\n";
        std::cout << "  e = " << std::setprecision(6) << fitted_orbit.eccentricity << "\n";
        std::cout << "  i = " << std::setprecision(4) << fitted_orbit.inclination * 180.0 / constants::PI << "°\n";
        std::cout << "  Ω = " << std::setprecision(4) << fitted_orbit.longitude_ascending_node * 180.0 / constants::PI << "°\n";
        std::cout << "  ω = " << std::setprecision(4) << fitted_orbit.argument_perihelion * 180.0 / constants::PI << "°\n";
        std::cout << "  M = " << std::setprecision(4) << fitted_orbit.mean_anomaly * 180.0 / constants::PI << "°\n\n";

        // Comparison (propagate OrbFit to same epoch)
        auto orbfit_at_mean = engine.propagate_orbit(orbfit_orbit, mean_epoch);
        double delta_a_km = (fitted_orbit.semi_major_axis - orbfit_at_mean.semi_major_axis) * constants::AU / 1000.0;
        double delta_e = fitted_orbit.eccentricity - orbfit_at_mean.eccentricity;
        double delta_i_arcsec = (fitted_orbit.inclination - orbfit_at_mean.inclination) * 180.0 * 3600.0 / constants::PI;

        std::cout << "Comparison with OrbFit (at same epoch):\n";
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
