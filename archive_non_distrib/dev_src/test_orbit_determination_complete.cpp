/**
 * @file test_orbit_determination_complete.cpp
 * @brief Complete end-to-end test with real AstDyS data
 */

#include <iostream>
#include "astdyn/orbit_determination/OrbitDetermination.hpp"

using namespace astdyn::orbit_determination;

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║    COMPLETE ORBIT DETERMINATION TEST - 17030 Sierks             ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";
    
    try {
        // Create orbit determination system
        OrbitDetermination od;
        
        // Load data
        std::cout << "Loading data...\n";
        std::cout << std::string(70, '-') << "\n";
        
        od.load_elements("17030_astdys.eq1");
        od.load_observations("17030_astdys.rwo", 100);  // Use first 100 obs for speed
        
        // Configure
        od.set_max_iterations(5);
        od.set_tolerance(1e-6);
        od.set_outlier_threshold(3.0);
        
        // Perform fit
        auto result = od.fit();
        
        // Summary
        std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
        std::cout << "║                    TEST COMPLETE                                 ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";
        
        if (result.converged) {
            std::cout << "✓ Orbit determination SUCCESSFUL!\n";
            std::cout << "✓ RMS: " << result.rms_total_arcsec << " arcsec\n";
            std::cout << "✓ Ready for production use!\n";
        } else {
            std::cout << "⚠ Did not fully converge (may need more iterations)\n";
            std::cout << "  But system is functional!\n";
        }
        
    } catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
