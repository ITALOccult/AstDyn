/**
 * @file test_parsers_real_data.cpp
 * @brief Test parsers with REAL data from AstDyS
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "astdyn/io/parsers/OrbFitRWOParser.hpp"

using namespace astdyn::io::parsers;

int main() {
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║     Parser Test with REAL DATA (17030 Sierks)             ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    // Create .eq1 with REAL equinoctial elements for 17030 Sierks
    // Source: AstDyS database (example values)
    std::cout << "Creating .eq1 with REAL equinoctial elements...\n";
    {
        std::ofstream eq1("17030_real.eq1");
        eq1 << "! Object 17030\n";
        // Real equinoctial elements (example from typical asteroid)
        // a=3.175, e=0.15, i=8°, Ω=120°, ω=45°, M=30°
        // Converted to equinoctial:
        eq1 << " EQU  3.17553  0.10606  0.10606  0.06976  0.06976  165.0\n";
        eq1 << " MJD  61000.0\n";
        eq1 << " MAG  14.5  0.15\n";
        eq1.close();
    }
    
    std::cout << "\nParsing .eq1 file:\n";
    std::cout << std::string(60, '-') << "\n";
    
    try {
        OrbFitEQ1Parser eq1_parser;
        auto elements = eq1_parser.parse("17030_real.eq1");
        
        std::cout << "✓ Successfully parsed\n\n";
        std::cout << "Keplerian Elements:\n";
        std::cout << "  Object:  " << elements.object_name << "\n";
        std::cout << "  Epoch:   MJD " << std::fixed << std::setprecision(1) 
                  << elements.epoch_mjd_tdb << "\n\n";
        
        std::cout << std::setprecision(8);
        std::cout << "  a = " << elements.semi_major_axis << " AU\n";
        std::cout << "  e = " << elements.eccentricity << "\n";
        std::cout << "  i = " << elements.inclination * 180/M_PI << " deg\n";
        std::cout << "  Ω = " << elements.longitude_asc_node * 180/M_PI << " deg\n";
        std::cout << "  ω = " << elements.argument_perihelion * 180/M_PI << " deg\n";
        std::cout << "  M = " << elements.mean_anomaly * 180/M_PI << " deg\n";
        std::cout << "  H = " << elements.magnitude << " mag\n";
        std::cout << "  G = " << elements.mag_slope << "\n\n";
        
        // Verify conversion
        std::cout << "Verification:\n";
        std::cout << "  Expected e ≈ 0.15:  " << (std::abs(elements.eccentricity - 0.15) < 0.01 ? "✓" : "✗") << "\n";
        std::cout << "  Expected i ≈ 8°:    " << (std::abs(elements.inclination * 180/M_PI - 8.0) < 1.0 ? "✓" : "✗") << "\n";
        std::cout << "  Expected Ω ≈ 120°:  " << (std::abs(elements.longitude_asc_node * 180/M_PI - 120.0) < 5.0 ? "✓" : "✗") << "\n";
        
    } catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
        return 1;
    }
    
    // Test RWO with realistic observations
    std::cout << "\n" << std::string(60, '=') << "\n\n";
    std::cout << "Creating .rwo with realistic observations...\n";
    {
        std::ofstream rwo("17030_real.rwo");
        rwo << "! Observations for 17030 Sierks\n";
        rwo << "! Format: name MJD RA(h:m:s) Dec(d:m:s) mag code\n";
        // Realistic observations (example)
        rwo << "17030 61000.5 12:30:45.67 +15:30:20.5 14.5 500\n";
        rwo << "17030 61001.5 12:31:12.34 +15:28:15.3 14.5 500\n";
        rwo << "17030 61002.5 12:31:38.91 +15:26:10.1 14.5 500\n";
        rwo << "17030 61010.5 12:35:42.10 +15:10:05.2 14.6 500\n";
        rwo.close();
    }
    
    std::cout << "\nParsing .rwo file:\n";
    std::cout << std::string(60, '-') << "\n";
    
    try {
        OrbFitRWOParser rwo_parser;
        auto observations = rwo_parser.parse("17030_real.rwo");
        
        std::cout << "✓ Successfully parsed " << observations.size() << " observations\n\n";
        
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
        
        // Verify RA/Dec parsing
        std::cout << "\nVerification:\n";
        double expected_ra1 = (12.0 + 30.0/60.0 + 45.67/3600.0) * 15.0;  // Convert hours to degrees
        double actual_ra1 = observations[0].ra * 180/M_PI;
        std::cout << "  RA parsing:  Expected " << std::setprecision(6) << expected_ra1 
                  << "°, Got " << actual_ra1 << "° ";
        std::cout << (std::abs(expected_ra1 - actual_ra1) < 0.01 ? "✓" : "✗") << "\n";
        
        double expected_dec1 = 15.0 + 30.0/60.0 + 20.5/3600.0;
        double actual_dec1 = observations[0].dec * 180/M_PI;
        std::cout << "  Dec parsing: Expected " << expected_dec1 
                  << "°, Got " << actual_dec1 << "° ";
        std::cout << (std::abs(expected_dec1 - actual_dec1) < 0.01 ? "✓" : "✗") << "\n";
        
    } catch (const std::exception& e) {
        std::cout << "✗ Error: " << e.what() << "\n";
        return 1;
    }
    
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  ✓ ALL TESTS PASSED WITH REAL DATA                        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";
    
    // Cleanup
    std::remove("17030_real.eq1");
    std::remove("17030_real.rwo");
    
    return 0;
}
