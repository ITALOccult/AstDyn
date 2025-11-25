/**
 * @file orbfit_main.cpp
 * @brief Main OrbFit program - orbit determination from observations
 * @author OrbFit C++ Conversion Team
 * @date 2025-11-24
 * 
 * This program replicates the functionality of the original Fortran OrbFit:
 * - Reads optical observations (MPC format)
 * - Performs orbit determination via differential correction
 * - Generates ephemeris
 * - Analyzes close approaches
 * - Exports results
 */

#include "orbfit/OrbFitEngine.hpp"
#include "orbfit/core/Constants.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

using namespace orbfit;

// ============================================================================
// Command Line Interface
// ============================================================================

struct ProgramOptions {
    std::string obs_file;                  ///< Observations file
    std::string output_prefix = "orbfit";  ///< Output file prefix
    
    // Initial orbit (if provided)
    bool have_initial_orbit = false;
    double epoch_mjd = 60000.0;
    double a = 2.5;
    double e = 0.1;
    double i_deg = 10.0;
    double omega_deg = 80.0;
    double Omega_deg = 45.0;
    double M_deg = 30.0;
    
    // Settings
    bool compute_ephemeris = true;
    double ephem_start = 0.0;
    double ephem_end = 0.0;
    double ephem_step = 1.0;
    
    bool find_close_approaches = true;
    double ca_start = 0.0;
    double ca_end = 0.0;
    
    bool compute_moid = false;
    int moid_planet = 3;  // Earth
    
    bool verbose = true;
    int max_iterations = 10;
    double tolerance = 1e-12;
};

void print_usage(const char* program_name) {
    std::cout << "\nUSAGE: " << program_name << " [options] <observations_file>\n\n";
    std::cout << "OPTIONS:\n";
    std::cout << "  -o PREFIX       Output file prefix (default: orbfit)\n";
    std::cout << "  -v              Verbose output (default: on)\n";
    std::cout << "  -q              Quiet mode\n";
    std::cout << "  -i N            Maximum iterations (default: 10)\n";
    std::cout << "  -t TOL          Integration tolerance (default: 1e-12)\n";
    std::cout << "\n";
    std::cout << "INITIAL ORBIT (optional):\n";
    std::cout << "  -epoch MJD      Epoch (MJD TDB)\n";
    std::cout << "  -a AU           Semi-major axis (AU)\n";
    std::cout << "  -e              Eccentricity\n";
    std::cout << "  -i DEG          Inclination (degrees)\n";
    std::cout << "  -omega DEG      Argument of perihelion (degrees)\n";
    std::cout << "  -Omega DEG      Longitude of ascending node (degrees)\n";
    std::cout << "  -M DEG          Mean anomaly (degrees)\n";
    std::cout << "\n";
    std::cout << "EPHEMERIS:\n";
    std::cout << "  --ephem START END STEP   Generate ephemeris (MJD, days)\n";
    std::cout << "  --no-ephem               Skip ephemeris generation\n";
    std::cout << "\n";
    std::cout << "CLOSE APPROACHES:\n";
    std::cout << "  --ca START END           Search close approaches (MJD)\n";
    std::cout << "  --no-ca                  Skip close approach search\n";
    std::cout << "  --moid PLANET            Compute MOID with planet (1=Merc, 3=Earth, etc)\n";
    std::cout << "\n";
    std::cout << "EXAMPLES:\n";
    std::cout << "  " << program_name << " observations.txt\n";
    std::cout << "  " << program_name << " -o asteroid99 -i 15 obs.txt\n";
    std::cout << "  " << program_name << " --ephem 60000 61000 10 obs.txt\n";
    std::cout << "  " << program_name << " --moid 3 obs.txt   # MOID with Earth\n";
    std::cout << "\n";
}

bool parse_args(int argc, char** argv, ProgramOptions& opts) {
    if (argc < 2) {
        return false;
    }
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            return false;
        }
        else if (arg == "-o" && i + 1 < argc) {
            opts.output_prefix = argv[++i];
        }
        else if (arg == "-v") {
            opts.verbose = true;
        }
        else if (arg == "-q") {
            opts.verbose = false;
        }
        else if (arg == "-i" && i + 1 < argc) {
            opts.max_iterations = std::stoi(argv[++i]);
        }
        else if (arg == "-t" && i + 1 < argc) {
            opts.tolerance = std::stod(argv[++i]);
        }
        else if (arg == "-epoch" && i + 1 < argc) {
            opts.epoch_mjd = std::stod(argv[++i]);
            opts.have_initial_orbit = true;
        }
        else if (arg == "-a" && i + 1 < argc) {
            opts.a = std::stod(argv[++i]);
            opts.have_initial_orbit = true;
        }
        else if (arg == "-e" && i + 1 < argc) {
            opts.e = std::stod(argv[++i]);
            opts.have_initial_orbit = true;
        }
        else if (arg == "-i" && i + 1 < argc) {
            opts.i_deg = std::stod(argv[++i]);
            opts.have_initial_orbit = true;
        }
        else if (arg == "-omega" && i + 1 < argc) {
            opts.omega_deg = std::stod(argv[++i]);
            opts.have_initial_orbit = true;
        }
        else if (arg == "-Omega" && i + 1 < argc) {
            opts.Omega_deg = std::stod(argv[++i]);
            opts.have_initial_orbit = true;
        }
        else if (arg == "-M" && i + 1 < argc) {
            opts.M_deg = std::stod(argv[++i]);
            opts.have_initial_orbit = true;
        }
        else if (arg == "--ephem" && i + 3 < argc) {
            opts.ephem_start = std::stod(argv[++i]);
            opts.ephem_end = std::stod(argv[++i]);
            opts.ephem_step = std::stod(argv[++i]);
            opts.compute_ephemeris = true;
        }
        else if (arg == "--no-ephem") {
            opts.compute_ephemeris = false;
        }
        else if (arg == "--ca" && i + 2 < argc) {
            opts.ca_start = std::stod(argv[++i]);
            opts.ca_end = std::stod(argv[++i]);
            opts.find_close_approaches = true;
        }
        else if (arg == "--no-ca") {
            opts.find_close_approaches = false;
        }
        else if (arg == "--moid" && i + 1 < argc) {
            opts.moid_planet = std::stoi(argv[++i]);
            opts.compute_moid = true;
        }
        else if (arg[0] != '-') {
            opts.obs_file = arg;
        }
    }
    
    return !opts.obs_file.empty();
}

// ============================================================================
// Main Program
// ============================================================================

int main(int argc, char** argv) {
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                         ORBFIT C++                            ║\n";
    std::cout << "║          Orbit Determination from Observations                ║\n";
    std::cout << "║                                                               ║\n";
    std::cout << "║  Conversion from Fortran OrbFit (Milani & Gronchi)           ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════╝\n\n";
    
    // Parse command line
    ProgramOptions opts;
    if (!parse_args(argc, argv, opts)) {
        print_usage(argv[0]);
        return 1;
    }
    
    try {
        // ====================================================================
        // 1. Setup Configuration
        // ====================================================================
        
        OrbFitConfig config;
        config.verbose = opts.verbose;
        config.max_iterations = opts.max_iterations;
        config.tolerance = opts.tolerance;
        config.propagator_settings.include_planets = true;
        config.propagator_settings.include_asteroids = false;  // Can enable if needed
        config.propagator_settings.perturb_jupiter = true;
        config.propagator_settings.perturb_saturn = true;
        config.propagator_settings.perturb_earth = true;
        config.propagator_settings.perturb_venus = true;
        
        OrbFitEngine engine(config);
        
        // ====================================================================
        // 2. Load Observations
        // ====================================================================
        
        std::cout << "═══════════════════════════════════════════════════════════\n";
        std::cout << "STEP 1: Loading Observations\n";
        std::cout << "═══════════════════════════════════════════════════════════\n\n";
        
        int num_obs = engine.load_observations(opts.obs_file);
        if (num_obs == 0) {
            std::cerr << "ERROR: No observations loaded\n";
            return 1;
        }
        
        // ====================================================================
        // 3. Set Initial Orbit
        // ====================================================================
        
        std::cout << "\n═══════════════════════════════════════════════════════════\n";
        std::cout << "STEP 2: Initial Orbit\n";
        std::cout << "═══════════════════════════════════════════════════════════\n\n";
        
        if (opts.have_initial_orbit) {
            propagation::KeplerianElements initial;
            initial.epoch_mjd_tdb = opts.epoch_mjd;
            initial.semi_major_axis = opts.a;
            initial.eccentricity = opts.e;
            initial.inclination = opts.i_deg * constants::DEG_TO_RAD;
            initial.longitude_ascending_node = opts.Omega_deg * constants::DEG_TO_RAD;
            initial.argument_perihelion = opts.omega_deg * constants::DEG_TO_RAD;
            initial.mean_anomaly = opts.M_deg * constants::DEG_TO_RAD;
            initial.gravitational_parameter = constants::GMS;
            
            engine.set_initial_orbit(initial);
        } else {
            std::cout << "No initial orbit provided - attempting IOD...\n";
            try {
                auto initial = engine.initial_orbit_determination();
                engine.set_initial_orbit(initial);
            } catch (const std::exception& e) {
                std::cerr << "ERROR: IOD failed: " << e.what() << "\n";
                std::cerr << "Please provide initial orbit via command line options\n";
                return 1;
            }
        }
        
        // ====================================================================
        // 4. Orbit Determination (Differential Correction)
        // ====================================================================
        
        std::cout << "\n═══════════════════════════════════════════════════════════\n";
        std::cout << "STEP 3: Orbit Determination\n";
        std::cout << "═══════════════════════════════════════════════════════════\n";
        
        auto result = engine.fit_orbit();
        
        // Export orbit
        std::string orbit_file = opts.output_prefix + ".oef";
        engine.export_orbit(orbit_file, "oef");
        
        // ====================================================================
        // 5. Generate Ephemeris (optional)
        // ====================================================================
        
        if (opts.compute_ephemeris) {
            std::cout << "\n═══════════════════════════════════════════════════════════\n";
            std::cout << "STEP 4: Ephemeris Generation\n";
            std::cout << "═══════════════════════════════════════════════════════════\n\n";
            
            // Use observation timespan if not specified
            if (opts.ephem_start == 0.0) {
                const auto& obs = engine.observations();
                opts.ephem_start = obs.front().mjd_utc;
                opts.ephem_end = obs.back().mjd_utc;
            }
            
            auto ephemeris = engine.compute_ephemeris(
                opts.ephem_start,
                opts.ephem_end,
                opts.ephem_step);
            
            // Save ephemeris to file
            std::string ephem_file = opts.output_prefix + "_ephemeris.txt";
            std::ofstream file(ephem_file);
            file << "# OrbFit Ephemeris\n";
            file << "# MJD_TDB   X(AU)   Y(AU)   Z(AU)   VX(AU/day)   VY(AU/day)   VZ(AU/day)\n";
            file << std::fixed << std::setprecision(9);
            
            for (const auto& state : ephemeris) {
                file << state.epoch_mjd_tdb << "  "
                     << state.position[0] << "  "
                     << state.position[1] << "  "
                     << state.position[2] << "  "
                     << state.velocity[0] << "  "
                     << state.velocity[1] << "  "
                     << state.velocity[2] << "\n";
            }
            
            std::cout << "Ephemeris saved to: " << ephem_file << "\n";
        }
        
        // ====================================================================
        // 6. Close Approach Analysis (optional)
        // ====================================================================
        
        if (opts.find_close_approaches) {
            std::cout << "\n═══════════════════════════════════════════════════════════\n";
            std::cout << "STEP 5: Close Approach Analysis\n";
            std::cout << "═══════════════════════════════════════════════════════════\n\n";
            
            // Use extended timespan if not specified
            if (opts.ca_start == 0.0) {
                const auto& obs = engine.observations();
                double span = obs.back().mjd_utc - obs.front().mjd_utc;
                opts.ca_start = obs.front().mjd_utc - span;
                opts.ca_end = obs.back().mjd_utc + span;
            }
            
            auto approaches = engine.find_close_approaches(
                opts.ca_start,
                opts.ca_end);
            
            if (!approaches.empty()) {
                std::string ca_file = opts.output_prefix + "_close_approaches.txt";
                std::ofstream file(ca_file);
                file << "# Close Approaches\n";
                file << "# Body  MJD_TDB  Distance(AU)  Distance(Radii)  V_rel(km/s)\n";
                file << std::fixed << std::setprecision(6);
                
                for (const auto& ca : approaches) {
                    file << static_cast<int>(ca.body) << "  "
                         << ca.time_mjd_tdb << "  "
                         << ca.distance_au << "  "
                         << ca.distance_in_radii << "  "
                         << ca.relative_velocity_kms << "\n";
                }
                
                std::cout << "Close approaches saved to: " << ca_file << "\n";
            }
        }
        
        // ====================================================================
        // 7. MOID Computation (optional)
        // ====================================================================
        
        if (opts.compute_moid) {
            std::cout << "\n═══════════════════════════════════════════════════════════\n";
            std::cout << "STEP 6: MOID Computation\n";
            std::cout << "═══════════════════════════════════════════════════════════\n\n";
            
            auto planet = static_cast<ephemeris::CelestialBody>(opts.moid_planet);
            double moid = engine.compute_moid(planet);
            
            std::cout << "MOID = " << std::fixed << std::setprecision(6) 
                     << moid << " AU\n";
        }
        
        // ====================================================================
        // Summary
        // ====================================================================
        
        std::cout << "\n═══════════════════════════════════════════════════════════\n";
        std::cout << "SUMMARY\n";
        std::cout << "═══════════════════════════════════════════════════════════\n\n";
        
        std::cout << "✓ Observations processed: " << result.num_observations << "\n";
        std::cout << "✓ Outliers rejected: " << result.num_rejected << "\n";
        std::cout << "✓ RMS residuals: " 
                 << std::fixed << std::setprecision(3)
                 << result.rms_ra << " × " << result.rms_dec << " arcsec\n";
        std::cout << "✓ Iterations: " << result.num_iterations << "\n";
        std::cout << "✓ Converged: " << (result.converged ? "Yes" : "No") << "\n";
        std::cout << "✓ Orbit saved: " << orbit_file << "\n";
        
        std::cout << "\n═══════════════════════════════════════════════════════════\n";
        std::cout << "OrbFit completed successfully!\n";
        std::cout << "═══════════════════════════════════════════════════════════\n\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << "\n\n";
        return 1;
    }
}
