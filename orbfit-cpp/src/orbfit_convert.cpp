/**
 * @file orbfit_convert.cpp
 * @brief Tool to convert OrbFit Fortran configuration to C++ format
 * 
 * Usage:
 *   orbfit_convert --input path/to/object.oef [options]
 * 
 * Converts:
 * - MEAN elements → OSCULATING elements
 * - .rwo files → validated observations
 * - .opt files → modern configuration
 */

#include <orbfit/io/OrbFitConfig.hpp>
#include <orbfit/OrbFitEngine.hpp>
#include <iostream>
#include <iomanip>

using namespace orbfit;
using namespace orbfit::config;

void print_usage() {
    std::cout << R"(
OrbFit Configuration Converter
===============================

Usage:
  orbfit_convert --config <path/to/config>
  orbfit_convert --oef <file.oef> --obs <file.obs>
  
Options:
  --config <dir>       Directory with object.oef, object.obs, object.opt
  --object <name>      Object name (default: from filename)
  --oef <file>         Input OEF file
  --obs <file>         Input observation file (.obs or .rwo)
  --output <dir>       Output directory (default: current)
  --convert-mean       Convert MEAN → OSCULATING elements
  --verbose           Print detailed information

Examples:
  # Load complete configuration
  orbfit_convert --config ./data --object ceres
  
  # Convert mean elements
  orbfit_convert --oef ceres_mean.oef --convert-mean --output ./output
  
  # Process observations
  orbfit_convert --oef ceres.oef --obs ceres.rwo --output ./results

)" << std::endl;
}

struct ProgramOptions {
    std::string config_dir;
    std::string object_name;
    std::string oef_file;
    std::string obs_file;
    std::string output_dir = ".";
    bool convert_mean = false;
    bool verbose = false;
};

ProgramOptions parse_args(int argc, char** argv) {
    ProgramOptions opts;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            print_usage();
            exit(0);
        } else if (arg == "--config" && i + 1 < argc) {
            opts.config_dir = argv[++i];
        } else if (arg == "--object" && i + 1 < argc) {
            opts.object_name = argv[++i];
        } else if (arg == "--oef" && i + 1 < argc) {
            opts.oef_file = argv[++i];
        } else if (arg == "--obs" && i + 1 < argc) {
            opts.obs_file = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            opts.output_dir = argv[++i];
        } else if (arg == "--convert-mean") {
            opts.convert_mean = true;
        } else if (arg == "--verbose" || arg == "-v") {
            opts.verbose = true;
        }
    }
    
    return opts;
}

void print_elements(const propagation::KeplerianElements& elem, const std::string& label) {
    std::cout << "\n=== " << label << " ===\n";
    std::cout << std::fixed << std::setprecision(9);
    std::cout << "Epoch:  " << elem.epoch_mjd_tdb << " MJD TDB\n";
    std::cout << "a:      " << elem.semi_major_axis << " AU\n";
    std::cout << "e:      " << elem.eccentricity << "\n";
    std::cout << std::setprecision(6);
    std::cout << "i:      " << (elem.inclination * constants::RAD_TO_DEG) << " deg\n";
    std::cout << "Ω:      " << (elem.longitude_ascending_node * constants::RAD_TO_DEG) << " deg\n";
    std::cout << "ω:      " << (elem.argument_perihelion * constants::RAD_TO_DEG) << " deg\n";
    std::cout << "M:      " << (elem.mean_anomaly * constants::RAD_TO_DEG) << " deg\n";
}

void compare_elements(const propagation::KeplerianElements& mean_elem,
                     const propagation::KeplerianElements& osc_elem) {
    std::cout << "\n=== Differences (OSCULATING - MEAN) ===\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Δa:     " << (osc_elem.semi_major_axis - mean_elem.semi_major_axis) << " AU\n";
    std::cout << "Δe:     " << (osc_elem.eccentricity - mean_elem.eccentricity) << "\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Δi:     " << ((osc_elem.inclination - mean_elem.inclination) * constants::RAD_TO_DEG) << " deg\n";
    std::cout << "ΔΩ:     " << ((osc_elem.longitude_ascending_node - mean_elem.longitude_ascending_node) * constants::RAD_TO_DEG) << " deg\n";
    std::cout << "Δω:     " << ((osc_elem.argument_perihelion - mean_elem.argument_perihelion) * constants::RAD_TO_DEG) << " deg\n";
    std::cout << "ΔM:     " << ((osc_elem.mean_anomaly - mean_elem.mean_anomaly) * constants::RAD_TO_DEG) << " deg\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    auto opts = parse_args(argc, argv);
    
    std::cout << "OrbFit Configuration Converter\n";
    std::cout << "==============================\n\n";
    
    try {
        OrbFitConfigManager config_mgr;
        
        // Load configuration
        if (!opts.config_dir.empty()) {
            std::cout << "Loading configuration from: " << opts.config_dir << "\n";
            std::cout << "Object name: " << opts.object_name << "\n\n";
            
            if (!config_mgr.loadConfiguration(opts.config_dir, opts.object_name)) {
                std::cerr << "Error: Failed to load configuration\n";
                return 1;
            }
        } else if (!opts.oef_file.empty()) {
            std::cout << "Loading OEF file: " << opts.oef_file << "\n\n";
            
            // Load individual files
            auto oef_data = OEFFileHandler::read(opts.oef_file);
            
            std::cout << "Object: " << oef_data.object_name << "\n";
            std::cout << "Format: " << oef_data.element_format << "\n";
            std::cout << "Type:   " << (oef_data.element_type == OrbitalElementSubType::MEAN ? "MEAN" : "OSCULATING") << "\n";
            
            // Print original elements
            print_elements(oef_data.keplerian, "Original Elements");
            
            // Convert if MEAN
            if (oef_data.element_type == OrbitalElementSubType::MEAN || opts.convert_mean) {
                auto osc_elem = OEFFileHandler::meanToOsculating(oef_data.keplerian);
                print_elements(osc_elem, "Osculating Elements");
                compare_elements(oef_data.keplerian, osc_elem);
                
                // Save converted
                OrbitalElementFile oef_osc = oef_data;
                oef_osc.keplerian = osc_elem;
                oef_osc.element_type = OrbitalElementSubType::OSCULATING;
                
                std::string output_file = opts.output_dir + "/" + oef_data.object_name + "_osc.oef";
                OEFFileHandler::write(output_file, oef_osc);
                std::cout << "\nConverted elements saved to: " << output_file << "\n";
            }
            
            return 0;
        } else {
            std::cerr << "Error: Must specify --config or --oef\n";
            print_usage();
            return 1;
        }
        
        // Get elements (auto-converts MEAN → OSCULATING)
        auto elements = config_mgr.getOsculatingElements();
        auto original = config_mgr.getOriginalElements();
        
        std::cout << "Element type: " << 
            (config_mgr.getElementType() == OrbitalElementSubType::MEAN ? "MEAN" : "OSCULATING") << "\n";
        
        print_elements(original, "Original Elements");
        
        if (config_mgr.getElementType() == OrbitalElementSubType::MEAN) {
            print_elements(elements, "Converted to Osculating");
            compare_elements(original, elements);
        }
        
        // Process observations
        auto observations = config_mgr.getValidObservations();
        std::cout << "\nObservations: " << observations.size() << " loaded\n";
        
        if (!observations.empty()) {
            std::cout << "Time span: " << observations.front().mjd_utc 
                     << " - " << observations.back().mjd_utc << " MJD\n";
        }
        
        // Save to output directory
        if (!opts.output_dir.empty()) {
            config_mgr.saveConfiguration(opts.output_dir, opts.object_name);
            std::cout << "\nConfiguration saved to: " << opts.output_dir << "\n";
        }
        
        std::cout << "\nConversion complete!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
