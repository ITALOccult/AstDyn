/**
 * @file benchmark_json_api.cpp
 * @brief Benchmark for the high-level AsteroidFitter JSON API
 */

#include "astdyn/ephemeris/AsteroidFitter.hpp"
#include "astdyn/ephemeris/AsteroidFitConfig.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <filesystem>

using namespace astdyn::ephemeris;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config.json>" << std::endl;
        return 1;
    }

    std::string config_path = argv[1];
    if (!std::filesystem::exists(config_path)) {
        std::cerr << "Error: Config file not found: " << config_path << std::endl;
        return 1;
    }

    std::cout << "╔════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║       ASTDYN HIGH-LEVEL JSON API BENCHMARK                 ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════╝" << std::endl << std::endl;

    try {
        std::cout << "[1/2] Loading Configuration: " << config_path << "..." << std::endl;
        auto config = loadAsteroidFitConfig(config_path);

        std::cout << "[2/2] Running Complete Workflow (Fit + Propagation)..." << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        AsteroidFitResult result = AsteroidFitter::fitFromConfig(config);
        auto end = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double> elapsed = end - start;

        if (result.success) {
            std::cout << std::endl << "✅ BENCHMARK SUCCESSFUL" << std::endl;
            std::cout << "────────────────────────────────────────────────────────────" << std::endl;
            std::cout << std::left << std::setw(25) << "Total Execution Time:" << std::fixed << std::setprecision(3) << elapsed.count() << " seconds" << std::endl;
            std::cout << std::left << std::setw(25) << "RMS RA:" << result.rms_ra << " arcsec" << std::endl;
            std::cout << std::left << std::setw(25) << "RMS Dec:" << result.rms_dec << " arcsec" << std::endl;
            std::cout << std::left << std::setw(25) << "Observations Used:" << result.num_observations << std::endl;
            std::cout << std::left << std::setw(25) << "Positions Generated:" << result.fitted_positions.size() << std::endl;
            std::cout << "────────────────────────────────────────────────────────────" << std::endl;
            
            if (!result.fitted_positions.empty()) {
                std::cout << "Last Geocentric Position (AU):" << std::endl;
                std::cout << "  X: " << result.fitted_positions.back().x() << std::endl;
                std::cout << "  Y: " << result.fitted_positions.back().y() << std::endl;
                std::cout << "  Z: " << result.fitted_positions.back().z() << std::endl;
            }
        } else {
            std::cerr << std::endl << "❌ BENCHMARK FAILED" << std::endl;
            std::cerr << "Message: " << result.message << std::endl;
            return 1;
        }

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
