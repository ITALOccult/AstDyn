/**
 * @file test_astdys_parser.cpp
 * @brief Test AstDyS RWO parser with real data
 */

#include <iostream>
#include <iomanip>
#include "astdyn/io/parsers/AstDysRWOParser.hpp"

using namespace astdyn::io::parsers;

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║        Test AstDyS RWO Parser (MPC Extended Format)             ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";
    
    try {
        AstDysRWOParser parser;
        auto observations = parser.parse("17030_astdys.rwo");
        
        std::cout << "✓ Successfully parsed 17030_astdys.rwo\n\n";
        std::cout << "Total observations: " << observations.size() << "\n\n";
        
        // Show first 10
        std::cout << "FIRST 10 OBSERVATIONS:\n";
        std::cout << std::string(80, '-') << "\n";
        std::cout << std::setw(12) << "MJD"
                  << std::setw(15) << "RA (deg)"
                  << std::setw(15) << "Dec (deg)"
                  << std::setw(8) << "Mag"
                  << std::setw(8) << "Code"
                  << "\n";
        std::cout << std::string(80, '-') << "\n";
        
        for (size_t i = 0; i < std::min(size_t(10), observations.size()); ++i) {
            const auto& obs = observations[i];
            std::cout << std::fixed << std::setprecision(5)
                      << std::setw(12) << obs.mjd_utc
                      << std::setprecision(6)
                      << std::setw(15) << obs.ra * 180/M_PI
                      << std::setw(15) << obs.dec * 180/M_PI
                      << std::setprecision(1)
                      << std::setw(8) << obs.mag
                      << std::setw(8) << obs.obs_code
                      << "\n";
        }
        
        // Statistics
        double min_mjd = observations[0].mjd_utc;
        double max_mjd = observations[0].mjd_utc;
        for (const auto& obs : observations) {
            min_mjd = std::min(min_mjd, obs.mjd_utc);
            max_mjd = std::max(max_mjd, obs.mjd_utc);
        }
        
        std::cout << "\nSTATISTICS:\n";
        std::cout << std::string(80, '-') << "\n";
        std::cout << "  First obs: MJD " << std::setprecision(5) << min_mjd << "\n";
        std::cout << "  Last obs:  MJD " << max_mjd << "\n";
        std::cout << "  Timespan:  " << std::setprecision(1) 
                  << (max_mjd - min_mjd) / 365.25 << " years\n";
        
        std::cout << "\n✓ Parser works correctly!\n";
        
    } catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
