/**
 * @file example_configurable_parsers.cpp
 * @brief Demonstrates configurable parser system for orbit fitting
 * 
 * This example shows how to:
 * 1. Use ParserFactory with auto-detection (default)
 * 2. Explicitly specify parser format
 * 3. Perform orbit fitting with parsed data
 */

#include <iostream>
#include <iomanip>
#include "astdyn/io/IParser.hpp"

using namespace astdyn::io;

void print_orbital_elements(const IOrbitParser::OrbitalElements& elem) {
    std::cout << "\n=== Orbital Elements ===" << std::endl;
    std::cout << "Object: " << elem.object_name << std::endl;
    std::cout << "Epoch: MJD " << std::fixed << std::setprecision(6) 
              << elem.epoch_mjd_tdb << " TDB" << std::endl;
    std::cout << "\nKeplerian Elements:" << std::endl;
    std::cout << "  a = " << std::setprecision(10) << elem.semi_major_axis << " AU" << std::endl;
    std::cout << "  e = " << std::setprecision(8) << elem.eccentricity << std::endl;
    std::cout << "  i = " << std::setprecision(6) << elem.inclination * 180.0 / M_PI << "°" << std::endl;
    std::cout << "  Ω = " << elem.longitude_asc_node * 180.0 / M_PI << "°" << std::endl;
    std::cout << "  ω = " << elem.argument_perihelion * 180.0 / M_PI << "°" << std::endl;
    std::cout << "  M = " << elem.mean_anomaly * 180.0 / M_PI << "°" << std::endl;
    
    if (elem.magnitude > 0) {
        std::cout << "\nMagnitude: H = " << elem.magnitude 
                  << ", G = " << elem.mag_slope << std::endl;
    }
}

void print_observations_summary(const std::vector<IObservationParser::OpticalObservation>& obs) {
    std::cout << "\n=== Observations Summary ===" << std::endl;
    std::cout << "Total observations: " << obs.size() << std::endl;
    
    if (!obs.empty()) {
        std::cout << "First observation: MJD " << std::fixed << std::setprecision(5) 
                  << obs.front().mjd_utc << " UTC" << std::endl;
        std::cout << "Last observation:  MJD " << obs.back().mjd_utc << " UTC" << std::endl;
        std::cout << "Time span: " << std::setprecision(2) 
                  << (obs.back().mjd_utc - obs.front().mjd_utc) << " days" << std::endl;
        
        // Show first observation details
        std::cout << "\nFirst observation details:" << std::endl;
        std::cout << "  Object: " << obs.front().object_name << std::endl;
        std::cout << "  RA:  " << std::setprecision(6) << obs.front().ra * 180.0 / M_PI << "°" << std::endl;
        std::cout << "  DEC: " << obs.front().dec * 180.0 / M_PI << "°" << std::endl;
        std::cout << "  Mag: " << std::setprecision(2) << obs.front().mag << std::endl;
        std::cout << "  Observatory: " << obs.front().obs_code << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <orbit_file> [observation_file] [format_hint]" << std::endl;
        std::cerr << "\nExamples:" << std::endl;
        std::cerr << "  " << argv[0] << " data.eq1                    # Auto-detect format" << std::endl;
        std::cerr << "  " << argv[0] << " data.eq1 obs.rwo           # Parse both files" << std::endl;
        std::cerr << "  " << argv[0] << " data.txt auto eq1          # Explicit format hint" << std::endl;
        return 1;
    }

    std::string orbit_file = argv[1];
    std::string obs_file = (argc > 2) ? argv[2] : "";
    std::string format_hint = (argc > 3) ? argv[3] : "auto";

    try {
        // ========== EXAMPLE 1: Auto-detection (default) ==========
        std::cout << "==========================================\n";
        std::cout << "EXAMPLE 1: Auto-detection parser\n";
        std::cout << "==========================================\n";
        
        auto orbit_parser = ParserFactory::create_orbit_parser(orbit_file);
        std::cout << "Using parser: " << orbit_parser->name() << std::endl;
        
        auto elements = orbit_parser->parse(orbit_file);
        print_orbital_elements(elements);

        // ========== EXAMPLE 2: Parse observations (if provided) ==========
        if (!obs_file.empty()) {
            std::cout << "\n==========================================\n";
            std::cout << "EXAMPLE 2: Parse observations\n";
            std::cout << "==========================================\n";
            
            auto obs_parser = ParserFactory::create_observation_parser(obs_file);
            std::cout << "Using parser: " << obs_parser->name() << std::endl;
            
            // Parse first 10 observations
            auto observations = obs_parser->parse(obs_file, 10);
            print_observations_summary(observations);
            
            // Parse all observations
            std::cout << "\nParsing all observations..." << std::endl;
            observations = obs_parser->parse(obs_file, 0);
            print_observations_summary(observations);
        }

        // ========== EXAMPLE 3: Explicit format specification ==========
        if (format_hint != "auto") {
            std::cout << "\n==========================================\n";
            std::cout << "EXAMPLE 3: Explicit format (" << format_hint << ")\n";
            std::cout << "==========================================\n";
            
            auto parser_explicit = ParserFactory::create_orbit_parser(orbit_file, format_hint);
            std::cout << "Using parser: " << parser_explicit->name() << std::endl;
            
            auto elem_explicit = parser_explicit->parse(orbit_file);
            print_orbital_elements(elem_explicit);
        }

        // ========== EXAMPLE 4: Check parser capabilities ==========
        std::cout << "\n==========================================\n";
        std::cout << "EXAMPLE 4: Parser capabilities\n";
        std::cout << "==========================================\n";
        
        auto test_parser = ParserFactory::create_orbit_parser(orbit_file);
        std::cout << "Parser name: " << test_parser->name() << std::endl;
        std::cout << "Can handle '" << orbit_file << "': " 
                  << (test_parser->can_handle(orbit_file) ? "YES" : "NO") << std::endl;
        std::cout << "Can handle 'test.mpc': " 
                  << (test_parser->can_handle("test.mpc") ? "YES" : "NO") << std::endl;

        std::cout << "\n✅ All examples completed successfully!" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
