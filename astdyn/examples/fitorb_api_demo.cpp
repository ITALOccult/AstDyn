/**
 * @file fitorb_api_demo.cpp
 * @brief Demonstration of the high-level OrbFitAPI
 * 
 * Uses the simplified run_fit() function to perform the orbit fitting
 * workflow that previously required ~200 lines of code.
 */

#include "astdyn/orbit_determination/OrbFitAPI.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>
#include <iomanip>
#include <filesystem>

namespace {

std::string resolve_input_path(const std::string& input_path) {
    if (input_path.empty()) {
        return input_path;
    }

    const std::filesystem::path path(input_path);
    if (path.is_absolute() || std::filesystem::exists(path)) {
        return input_path;
    }

    const std::filesystem::path candidate_in_tools = std::filesystem::path("tools") / path;
    if (std::filesystem::exists(candidate_in_tools)) {
        return candidate_in_tools.string();
    }

    const std::filesystem::path candidate_in_repo_tools = std::filesystem::path("astdyn") / "tools" / path;
    if (std::filesystem::exists(candidate_in_repo_tools)) {
        return candidate_in_repo_tools.string();
    }

    return input_path;
}

} // namespace

int main(int argc, char** argv) {
    std::string eq1_file = resolve_input_path((argc > 1) ? argv[1] : "203_astdys.eq1");
    std::string rwo_file = resolve_input_path((argc > 2) ? argv[2] : "203.rwo");
    std::string oop_file = (argc > 3) ? argv[3] : ""; // Optional configuration

    std::cout << "=== Running OrbFitAPI Demo ===\n";
    std::cout << "Orbit: " << eq1_file << "\n";
    std::cout << "Obs:   " << rwo_file << "\n\n";

    // One-shot execution
    auto result = astdyn::orbit_determination::OrbFitAPI::run_fit(eq1_file, rwo_file, oop_file, true);

    if (result.success) {
        std::cout << "\n=== FIT SUCCESS ===\n";
        std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
        std::cout << "Iterations: " << result.iterations << "\n";
        std::cout << "RMS RA: " << std::fixed << std::setprecision(3) << result.rms_ra.value / 1000.0 << " arcsec\n";
        std::cout << "RMS Dec: " << result.rms_dec.value / 1000.0 << " arcsec\n";
        std::cout << "Outliers: " << result.num_outliers << "/" << result.num_observations << "\n";
        
        std::cout << "\nChange from Initial Orbit:\n";
        std::cout << "  da: " << result.delta_a_km << " km\n";
        // result.delta_i_arcsec removed as it's not in the relocation plan to keep structure clean

        return 0;
    } else {
        std::cerr << "\n=== FIT FAILED ===\n";
        std::cerr << "Error: " << result.message << "\n";
        return 1;
    }
}
