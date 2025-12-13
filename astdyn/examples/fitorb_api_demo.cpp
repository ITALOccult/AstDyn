/**
 * @file fitorb_api_demo.cpp
 * @brief Demonstration of the high-level OrbitFitAPI
 * 
 * Uses the simplified run_fit() function to perform the orbit fitting
 * workflow that previously required ~200 lines of code.
 */

#include "astdyn/api/OrbitFitAPI.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>
#include <iomanip>

int main(int argc, char** argv) {
    std::string eq1_file = "203_astdys.eq1";
    std::string rwo_file = "203.rwo";
    std::string oop_file = "203.oop"; // Optional configuration

    std::cout << "=== Running OrbitFitAPI Demo ===\n";
    std::cout << "Orbit: " << eq1_file << "\n";
    std::cout << "Obs:   " << rwo_file << "\n\n";

    // One-shot execution
    auto result = astdyn::api::OrbitFitAPI::run_fit(eq1_file, rwo_file, oop_file, true);

    if (result.success) {
        std::cout << "\n=== FIT SUCCESS ===\n";
        std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
        std::cout << "Iterations: " << result.iterations << "\n";
        std::cout << "RMS RA: " << std::fixed << std::setprecision(3) << result.rms_ra << " arcsec\n";
        std::cout << "RMS Dec: " << result.rms_dec << " arcsec\n";
        std::cout << "Outliers: " << result.num_outliers << "/" << result.num_observations << "\n";
        
        std::cout << "\nChange from Initial Orbit:\n";
        std::cout << "  da: " << result.delta_a_km << " km\n";
        std::cout << "  di: " << result.delta_i_arcsec << " arcsec\n";

        return 0;
    } else {
        std::cerr << "\n=== FIT FAILED ===\n";
        std::cerr << "Error: " << result.message << "\n";
        return 1;
    }
}
