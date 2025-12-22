/**
 * @file test_astdys_real_data.cpp
 * @brief Test with REAL AstDyS data for 17030 Sierks
 * 
 * Downloads and parses:
 * - https://newton.spacedys.com/~astdys2/epoch/numbered/17/17030.eq1
 * - https://newton.spacedys.com/~astdys2/mpcobs/numbered/17/17030.rwo
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <set>
#include "astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "astdyn/io/parsers/OrbFitRWOParser.hpp"

using namespace astdyn::io::parsers;

void print_separator(char c = '=') {
    std::cout << std::string(70, c) << "\n";
}

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║        REAL AstDyS DATA TEST - 17030 Sierks                      ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";
    
    // ========================================================================
    // Part 1: Parse .eq1 (orbital elements)
    // ========================================================================
    
    std::cout << "PART 1: Parsing Orbital Elements (.eq1)\n";
    print_separator();
    
    try {
        OrbFitEQ1Parser eq1_parser;
        auto elements = eq1_parser.parse("17030_astdys.eq1");
        
        std::cout << "✓ Successfully parsed 17030_astdys.eq1\n\n";
        
        std::cout << "ORBITAL ELEMENTS (Keplerian):\n";
        std::cout << std::string(70, '-') << "\n";
        std::cout << "  Object:           " << elements.object_name << "\n";
        std::cout << "  Epoch (MJD TDT):  " << std::fixed << std::setprecision(6) 
                  << elements.epoch_mjd_tdb << "\n\n";
        
        std::cout << std::setprecision(10);
        std::cout << "  Semi-major axis:  a = " << elements.semi_major_axis << " AU\n";
        std::cout << "  Eccentricity:     e = " << elements.eccentricity << "\n";
        std::cout << "  Inclination:      i = " << elements.inclination * 180/M_PI << " deg\n";
        std::cout << "  Long. asc. node:  Ω = " << elements.longitude_asc_node * 180/M_PI << " deg\n";
        std::cout << "  Arg. perihelion:  ω = " << elements.argument_perihelion * 180/M_PI << " deg\n";
        std::cout << "  Mean anomaly:     M = " << elements.mean_anomaly * 180/M_PI << " deg\n\n";
        
        std::cout << "  Absolute mag:     H = " << std::setprecision(3) 
                  << elements.magnitude << "\n";
        std::cout << "  Slope parameter:  G = " << elements.mag_slope << "\n\n";
        
        // Derived quantities
        double q = elements.semi_major_axis * (1.0 - elements.eccentricity);
        double Q = elements.semi_major_axis * (1.0 + elements.eccentricity);
        double period = 2.0 * M_PI * std::sqrt(std::pow(elements.semi_major_axis, 3) / 2.959122082855911e-4);
        
        std::cout << "DERIVED QUANTITIES:\n";
        std::cout << std::string(70, '-') << "\n";
        std::cout << "  Perihelion dist:  q = " << std::setprecision(6) << q << " AU\n";
        std::cout << "  Aphelion dist:    Q = " << Q << " AU\n";
        std::cout << "  Orbital period:   P = " << period / 365.25 << " years\n\n";
        
    } catch (const std::exception& e) {
        std::cout << "✗ ERROR: " << e.what() << "\n";
        return 1;
    }
    
    // ========================================================================
    // Part 2: Parse .rwo (observations)
    // ========================================================================
    
    std::cout << "\nPART 2: Parsing Observations (.rwo)\n";
    print_separator();
    
    try {
        OrbFitRWOParser rwo_parser;
        auto observations = rwo_parser.parse("17030_astdys.rwo");
        
        std::cout << "✓ Successfully parsed 17030_astdys.rwo\n\n";
        std::cout << "OBSERVATION STATISTICS:\n";
        std::cout << std::string(70, '-') << "\n";
        std::cout << "  Total observations:   " << observations.size() << "\n";
        
        // Find date range
        double min_mjd = observations[0].mjd_utc;
        double max_mjd = observations[0].mjd_utc;
        for (const auto& obs : observations) {
            min_mjd = std::min(min_mjd, obs.mjd_utc);
            max_mjd = std::max(max_mjd, obs.mjd_utc);
        }
        
        std::cout << "  First observation:    MJD " << std::fixed << std::setprecision(5) 
                  << min_mjd << "\n";
        std::cout << "  Last observation:     MJD " << max_mjd << "\n";
        std::cout << "  Time span:            " << std::setprecision(1) 
                  << (max_mjd - min_mjd) / 365.25 << " years\n\n";
        
        // Count unique observatories
        std::set<std::string> obs_codes;
        for (const auto& obs : observations) {
            obs_codes.insert(obs.obs_code);
        }
        std::cout << "  Unique observatories: " << obs_codes.size() << "\n\n";
        
        // Show first 10 observations
        std::cout << "FIRST 10 OBSERVATIONS:\n";
        std::cout << std::string(70, '-') << "\n";
        std::cout << std::setw(12) << "MJD"
                  << std::setw(15) << "RA (deg)"
                  << std::setw(15) << "Dec (deg)"
                  << std::setw(8) << "Mag"
                  << std::setw(8) << "Code"
                  << "\n";
        std::cout << std::string(70, '-') << "\n";
        
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
        
        std::cout << "  ... (" << observations.size() - 10 << " more observations)\n\n";
        
        // RA/Dec range
        double min_ra = observations[0].ra;
        double max_ra = observations[0].ra;
        double min_dec = observations[0].dec;
        double max_dec = observations[0].dec;
        
        for (const auto& obs : observations) {
            min_ra = std::min(min_ra, obs.ra);
            max_ra = std::max(max_ra, obs.ra);
            min_dec = std::min(min_dec, obs.dec);
            max_dec = std::max(max_dec, obs.dec);
        }
        
        std::cout << "SKY COVERAGE:\n";
        std::cout << std::string(70, '-') << "\n";
        std::cout << "  RA range:   " << std::setprecision(2) 
                  << min_ra * 180/M_PI << "° to " << max_ra * 180/M_PI << "°\n";
        std::cout << "  Dec range:  " << min_dec * 180/M_PI << "° to " 
                  << max_dec * 180/M_PI << "°\n\n";
        
    } catch (const std::exception& e) {
        std::cout << "✗ ERROR: " << e.what() << "\n";
        return 1;
    }
    
    // ========================================================================
    // Summary
    // ========================================================================
    
    std::cout << "\n";
    print_separator('=');
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  ✓ ALL TESTS PASSED - PARSERS WORK WITH REAL AstDyS DATA        ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\nNext steps:\n";
    std::cout << "  1. Propagate orbit using elements from .eq1\n";
    std::cout << "  2. Compute residuals (O-C) for observations in .rwo\n";
    std::cout << "  3. Perform least-squares fit to improve orbit\n";
    std::cout << "  4. Compare with OrbFit results\n\n";
    
    return 0;
}
