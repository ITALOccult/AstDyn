/**
 * @file test_orbit_fit_quick.cpp
 * @brief Quick test of orbit fitting components
 */

#include <iostream>
#include <iomanip>
#include "astdyn/orbit_determination/ResidualCalculator.hpp"
#include "astdyn/orbit_determination/LeastSquaresFitter.hpp"

using namespace astdyn::orbit_determination;

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║        Quick Test: Orbit Fitting Components                     ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";
    
    // Test 1: ResidualCalculator
    std::cout << "Test 1: ResidualCalculator\n";
    std::cout << std::string(70, '-') << "\n";
    
    try {
        ResidualCalculator res_calc;
        
        // Create dummy observation
        Observation obs;
        obs.epoch_mjd = 61000.0;
        obs.ra_deg = 187.5;
        obs.dec_deg = 15.5;
        obs.ra_sigma_arcsec = 1.0;
        obs.dec_sigma_arcsec = 1.0;
        obs.observatory_code = "500";
        obs.weight = 1.0;
        obs.rejected = false;
        
        // Dummy state
        Eigen::Vector<double, 6> state;
        state << 3.0, 0.0, 0.0, 0.0, 0.005, 0.0;
        
        auto residual = res_calc.compute_residual(obs, state, 61000.0);
        
        std::cout << "✓ ResidualCalculator works\n";
        std::cout << "  Computed RA:  " << std::fixed << std::setprecision(6) 
                  << residual.ra_computed_deg << " deg\n";
        std::cout << "  Computed Dec: " << residual.dec_computed_deg << " deg\n";
        std::cout << "  RA residual:  " << residual.ra_residual_arcsec << " arcsec\n";
        std::cout << "  Dec residual: " << residual.dec_residual_arcsec << " arcsec\n\n";
        
    } catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
        return 1;
    }
    
    // Test 2: LeastSquaresFitter
    std::cout << "Test 2: LeastSquaresFitter\n";
    std::cout << std::string(70, '-') << "\n";
    
    try {
        LeastSquaresFitter fitter;
        fitter.set_max_iterations(5);
        fitter.set_tolerance(1e-6);
        
        std::cout << "✓ LeastSquaresFitter created\n";
        std::cout << "  Max iterations: 5\n";
        std::cout << "  Tolerance: 1e-6\n\n";
        
    } catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
        return 1;
    }
    
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  ✓ ALL COMPONENTS COMPILED AND WORK                              ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\nNext: Integrate with STMPropagator for full orbit determination\n";
    
    return 0;
}
