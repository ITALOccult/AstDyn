/**
 * @file example_astdys_fitter.cpp
 * @brief Example showing how to use AstDysOrbitFitter for simple orbit fitting
 * 
 * This example demonstrates the complete workflow in just a few lines:
 * 1. Load observations from .rwo file
 * 2. Load initial orbit from .eq1 file
 * 3. Run differential correction
 * 4. Compare with reference orbit
 */

#include <astdyn/io/AstDysOrbitFitter.hpp>
#include <iostream>
#include <iomanip>

int main(int argc, char* argv[]) {
    try {
        // Check arguments
        if (argc < 4) {
            std::cerr << "Usage: " << argv[0] << " <observations.rwo> <elements.eq1> <config.json>\n";
            std::cerr << "\nExample:\n";
            std::cerr << "  " << argv[0] << " astdyn/tools/203_astdys_recent100.rwo "
                      << "astdyn/tools/203_astdys_latest.eq1 "
                      << "astdyn/tools/astdyn.json\n";
            return 1;
        }

        std::string rwo_file = argv[1];
        std::string eq1_file = argv[2];
        std::string config_file = argv[3];

        std::cout << "╔═══════════════════════════════════════════════════════════╗\n";
        std::cout << "║      AstDyS Orbit Fitter - Simple Example                ║\n";
        std::cout << "╚═══════════════════════════════════════════════════════════╝\n\n";

        // Create fitter
        astdyn::io::AstDysOrbitFitter fitter;
        fitter.set_verbose(true);

        // Step 1: Load observations
        std::cout << "1. Loading observations...\n";
        fitter.set_observations_file(rwo_file);
        std::cout << "\n";

        // Step 2: Load initial orbital elements
        std::cout << "2. Loading orbital elements...\n";
        fitter.set_elements_file(eq1_file);
        std::cout << "\n";

        // Step 3: Load configuration
        std::cout << "3. Loading configuration...\n";
        fitter.set_config_file(config_file);
        std::cout << "\n";

        // Step 4: Run differential correction
        std::cout << "4. Running differential correction...\n";
        auto result = fitter.fit();
        std::cout << "\n";

        // Step 5: Display results
        std::cout << "╔═══════════════════════════════════════════════════════════╗\n";
        std::cout << "║                    FIT RESULTS                            ║\n";
        std::cout << "╚═══════════════════════════════════════════════════════════╝\n\n";

        std::cout << "Object: " << result.object_name << "\n";
        std::cout << "Observations: " << result.num_observations_used 
                  << " / " << result.num_observations_loaded << " used\n";
        std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
        std::cout << "Iterations: " << result.num_iterations << "\n";
        std::cout << "Outliers: " << result.num_outliers << "\n\n";

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "RMS RA:  " << result.rms_ra << " arcsec\n";
        std::cout << "RMS Dec: " << result.rms_dec << " arcsec\n";
        std::cout << "Chi²:    " << result.chi_squared << "\n\n";

        // Fitted orbital elements
        std::cout << "Fitted Orbital Elements:\n";
        std::cout << "  a = " << std::setprecision(6) << result.fitted_orbit.semi_major_axis << " AU\n";
        std::cout << "  e = " << std::setprecision(6) << result.fitted_orbit.eccentricity << "\n";
        std::cout << "  i = " << std::setprecision(4) << result.fitted_orbit.inclination * 180.0 / M_PI << "°\n";
        std::cout << "  Ω = " << std::setprecision(4) << result.fitted_orbit.longitude_ascending_node * 180.0 / M_PI << "°\n";
        std::cout << "  ω = " << std::setprecision(4) << result.fitted_orbit.argument_perihelion * 180.0 / M_PI << "°\n";
        std::cout << "  M = " << std::setprecision(4) << result.fitted_orbit.mean_anomaly * 180.0 / M_PI << "°\n";
        std::cout << "  Epoch: MJD " << std::setprecision(6) << result.fitted_orbit.epoch_mjd_tdb << " TDB\n\n";

        // Comparison with reference orbit
        if (result.delta_a_km.has_value()) {
            std::cout << "Comparison with Reference Orbit:\n";
            std::cout << "  Δa = " << std::setprecision(2) << result.delta_a_km.value() << " km\n";
            std::cout << "  Δe = " << std::setprecision(6) << result.delta_e.value() << "\n";
            std::cout << "  Δi = " << std::setprecision(2) << result.delta_i_arcsec.value() << " arcsec\n";
        }

        std::cout << "\n╔═══════════════════════════════════════════════════════════╗\n";
        std::cout << "║                 FIT COMPLETED                             ║\n";
        std::cout << "╚═══════════════════════════════════════════════════════════╝\n";

        return result.converged ? 0 : 1;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
