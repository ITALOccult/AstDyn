/**
 * @file test_parsers.cpp
 * @brief Test .eq1 and .rwo parsers
 */

#include <iostream>
#include <iomanip>
#include "astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "astdyn/io/parsers/OrbFitRWOParser.hpp"

using namespace astdyn::io::parsers;

int main() {
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║          Parser Test (EQ1 + RWO)                           ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    // Test 1: Create sample .eq1 file
    std::cout << "Creating sample .eq1 file...\n";
    {
        std::ofstream eq1("test_asteroid.eq1");
        eq1 << "! Object 17030\n";
        eq1 << " EQU  3.17553  0.0  0.0  0.0  0.0  180.0\n";
        eq1 << " MJD  61000.0\n";
        eq1 << " MAG  14.5  0.15\n";
        eq1.close();
    }
    
    // Test 2: Parse .eq1
    std::cout << "\nTesting EQ1 Parser:\n";
    std::cout << std::string(60, '-') << "\n";
    
    try {
        OrbFitEQ1Parser eq1_parser;
        auto elements = eq1_parser.parse("test_asteroid.eq1");
        
        std::cout << "✓ Successfully parsed .eq1 file\n";
        std::cout << "  Object: " << elements.object_name << "\n";
        std::cout << "  Epoch (MJD): " << std::fixed << std::setprecision(1) 
                  << elements.epoch_mjd_tdb << "\n";
        std::cout << "  a = " << std::setprecision(6) << elements.semi_major_axis << " AU\n";
        std::cout << "  e = " << elements.eccentricity << "\n";
        std::cout << "  i = " << elements.inclination * 180/M_PI << " deg\n";
        std::cout << "  Ω = " << elements.longitude_asc_node * 180/M_PI << " deg\n";
        std::cout << "  ω = " << elements.argument_perihelion * 180/M_PI << " deg\n";
        std::cout << "  M = " << elements.mean_anomaly * 180/M_PI << " deg\n";
        std::cout << "  H = " << elements.magnitude << "\n";
        
    } catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
        return 1;
    }
    
    // Test 3: Create sample .rwo file
    std::cout << "\nCreating sample .rwo file...\n";
    {
        std::ofstream rwo("test_observations.rwo");
        rwo << "! Test observations for 17030\n";
        rwo << "17030 61000.0 12:30:45.6 +15:30:20.5 14.5 500\n";
        rwo << "17030 61001.0 12:31:10.2 +15:28:15.3 14.5 500\n";
        rwo << "17030 61002.0 12:31:34.8 +15:26:10.1 14.5 500\n";
        rwo.close();
    }
    
    // Test 4: Parse .rwo
    std::cout << "\nTesting RWO Parser:\n";
    std::cout << std::string(60, '-') << "\n";
    
    try {
        OrbFitRWOParser rwo_parser;
        auto observations = rwo_parser.parse("test_observations.rwo");
        
        std::cout << "✓ Successfully parsed .rwo file\n";
        std::cout << "  Number of observations: " << observations.size() << "\n\n";
        
        std::cout << std::setw(12) << "MJD"
                  << std::setw(15) << "RA (deg)"
                  << std::setw(15) << "Dec (deg)"
                  << std::setw(10) << "Mag"
                  << std::setw(10) << "Code"
                  << "\n";
        std::cout << std::string(60, '-') << "\n";
        
        for (const auto& obs : observations) {
            std::cout << std::fixed << std::setprecision(1)
                      << std::setw(12) << obs.mjd_utc
                      << std::setprecision(6)
                      << std::setw(15) << obs.ra * 180/M_PI
                      << std::setw(15) << obs.dec * 180/M_PI
                      << std::setprecision(1)
                      << std::setw(10) << obs.mag
                      << std::setw(10) << obs.obs_code
                      << "\n";
        }
        
    } catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
        return 1;
    }
    
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  ✓ ALL PARSER TESTS PASSED                                 ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";
    
    // Cleanup
    std::remove("test_asteroid.eq1");
    std::remove("test_observations.rwo");
    
    return 0;
}
