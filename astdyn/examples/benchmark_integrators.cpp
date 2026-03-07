/**
 * @file benchmark_integrators.cpp
 * @brief Benchmark tool for comparing different integrators and force models against JPL Horizons truth.
 */

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <filesystem>

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/GaussIntegrator.hpp"
#include "astdyn/propagation/RadauIntegrator.hpp"
#include "astdyn/propagation/saba4_integrator.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/core/physics_state.hpp"

using namespace astdyn;
using namespace astdyn::io;
using namespace astdyn::propagation;
using namespace astdyn::astrometry;
using namespace astdyn::physics;

struct BenchmarkResult {
    std::string integrator;
    std::string model;
    double error_mas; // milliarcseconds
    double execution_time_ms;
};

// --- Time Offset Parser ---
double parse_offset_days(const std::string& offset) {
    if (offset.empty()) return 30.0;
    
    char unit = offset.back();
    double value = std::stod(offset.substr(0, offset.size() - 1));
    
    if (unit == 'M' || unit == 'm') return value * 30.4375; // Average month
    if (unit == 'G' || unit == 'g' || unit == 'D' || unit == 'd') return value;
    if (unit == 'A' || unit == 'a' || unit == 'Y' || unit == 'y') return value * 365.25;
    
    return std::stod(offset); // Assume days if no unit
}

void print_help(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "Options:\n"
              << "  -t, --target <name>      Asteroid name/number (default: 79518)\n"
              << "  -e, --epoch <jd>         Initial epoch in JD (default: current)\n"
              << "  -o, --offset <string>    Offset: +1M (Month), +10G (Day), +1A (Year) (default: +1M)\n"
              << "  -out, --output <file>    CSV output file (default: benchmark_results.csv)\n"
              << "  --bsp <path>             Path to planetary BSP file\n"
              << "  --absp <path>            Path to asteroid BSP file\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    std::string target = "79518";
    std::string offset_str = "+0G";
    std::string output_file = "benchmark_results.csv";
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    std::string absp_path = "";
    double initial_jd = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-t" || arg == "--target") && i + 1 < argc) target = argv[++i];
        else if ((arg == "-o" || arg == "--offset") && i + 1 < argc) offset_str = argv[++i];
        else if ((arg == "-out" || arg == "--output") && i + 1 < argc) output_file = argv[++i];
        else if ((arg == "-e" || arg == "--epoch") && i + 1 < argc) initial_jd = std::stod(argv[++i]);
        else if (arg == "--bsp" && i + 1 < argc) bsp_path = argv[++i];
        else if (arg == "--absp" && i + 1 < argc) absp_path = argv[++i];
        else if (arg == "--help") { print_help(argv[0]); return 0; }
    }

    if (offset_str[0] == '+') offset_str = offset_str.substr(1);
    double offset_days = parse_offset_days(offset_str);

    // 1. Setup Horizons Client
    HorizonsClient horizons;
    time::EpochTT start_time;
    if (initial_jd > 0) {
        start_time = time::EpochTT::from_mjd(initial_jd - 2400000.5);
    } else {
        start_time = time::EpochTT::from_mjd(time::now(TimeScale::TT));
    }
    
    time::EpochTT target_time_tt = time::EpochTT::from_mjd(start_time.mjd() + offset_days);
    time::EpochTDB start_time_tdb = astdyn::time::to_tdb(start_time);
    time::EpochTDB target_time_tdb = astdyn::time::to_tdb(target_time_tt);

    std::cout << "--- Benchmark Settings ---" << std::endl;
    std::cout << "Target:   " << target << std::endl;
    std::cout << "Start:    MJD " << std::fixed << std::setprecision(6) << start_time_tdb.mjd() << std::endl;
    std::cout << "Target:   MJD " << target_time_tdb.mjd() << " (+" << offset_days << " days)" << std::endl;
    std::cout << "BSP:      " << bsp_path << std::endl;

    // 2. Fetch Initial State (GCRF) from Horizons
    std::cout << "Fetching initial vectors from Horizons..." << std::endl;
    auto orbit_res = horizons.query_vectors(target, start_time_tdb, "@sun");
    if (!orbit_res) {
        std::cerr << "Error fetching vectors: " << (int)orbit_res.error() << std::endl;
        return 1;
    }
    const physics::CartesianStateTyped<core::GCRF> initial_vectors_eq = *orbit_res;

    // 3. Transform to ECLIPTIC (AstDyn native propagation frame) using typesafe API
    auto pos_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(initial_vectors_eq.position);
    auto vel_ecl = coordinates::ReferenceFrame::transform_vel<core::GCRF, core::ECLIPJ2000>(initial_vectors_eq.position, initial_vectors_eq.velocity);
    
    const physics::CartesianStateTyped<core::ECLIPJ2000> initial_vectors_ecl(
        start_time_tdb, pos_ecl, vel_ecl, initial_vectors_eq.gm
    );

    // 4. Fetch Truth Observation from Horizons
    std::cout << "Fetching truth observation from Horizons (target=" << target << ", mjd=" << target_time_tdb.mjd() << ")..." << std::endl;
    auto truth_res = horizons.query_observation(target, target_time_tdb, "500@399"); // Geocentric
    std::cout << "  [DEBUG-BENCH] query_observation returned." << std::endl;
    if (!truth_res) {
        std::cerr << "Error fetching truth: " << (int)truth_res.error() << std::endl;
        return 1;
    }
    auto truth_obs = *truth_res;

    // 5. Build Keplerian state (ECLIPTIC) - used for standard APIs
    // Conversion to Keplerian is still a bit manual in current header but let's use the provided bridge
    propagation::CartesianElements ce_old;
    ce_old.epoch = start_time_tdb;
    ce_old.position = types::Vector3<core::GCRF, core::Meter>(pos_ecl.to_eigen_si()); // bridge
    ce_old.velocity = types::Vector3<core::GCRF, core::Meter>(vel_ecl.to_eigen_si());
    ce_old.gravitational_parameter = initial_vectors_eq.gm.to_m3_s2();
    
    auto kep_old = propagation::cartesian_to_keplerian(ce_old);
    auto initial_state_kep = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        start_time_tdb, kep_old.semi_major_axis, kep_old.eccentricity,
        kep_old.inclination * constants::RAD_TO_DEG,
        kep_old.longitude_ascending_node * constants::RAD_TO_DEG,
        kep_old.argument_perihelion * constants::RAD_TO_DEG,
        kep_old.mean_anomaly * constants::RAD_TO_DEG,
        initial_vectors_eq.gm
    );

    std::cout << "\nInitial Keplerian Elements (Ecliptic J2000):" << std::endl;
    std::cout << "  a  = " << initial_state_kep.a.to_au() << " AU" << std::endl;
    std::cout << "  e  = " << initial_state_kep.e << std::endl;
    std::cout << "  i  = " << initial_state_kep.i.to_deg() << " deg" << std::endl;
    std::cout << "  OM = " << initial_state_kep.node.to_deg() << " deg" << std::endl;
    std::cout << "  W  = " << initial_state_kep.omega.to_deg() << " deg" << std::endl;
    std::cout << "  MA = " << initial_state_kep.M.to_deg() << " deg" << std::endl;

    std::cout << "\nTruth Observation (Equatorial J2000):" << std::endl;
    std::cout << "  RA  = " << truth_obs.ra.value * constants::RAD_TO_DEG << " deg" << std::endl;
    std::cout << "  DEC = " << truth_obs.dec.value * constants::RAD_TO_DEG << " deg" << std::endl;
    std::cout << "  Dist= " << truth_obs.distance.value / (constants::AU * 1000.0) << " AU" << std::endl;

    // 6. Define Test Matrix
    struct ModelDef {
        std::string name;
        PropagatorSettings settings;
    };

    std::vector<ModelDef> models = {
        {"no-force", [](){
            PropagatorSettings s;
            s.include_planets = false;
            s.include_relativity = false;
            s.include_asteroids = false;
            return s;
        }()},
        {"Full", [&](){
            PropagatorSettings s;
            s.include_planets = true;
            s.include_relativity = true;
            s.include_asteroids = true;
            s.asteroid_ephemeris_file = absp_path;
            return s;
        }()}
    };

    std::vector<std::string> integrators = {"RKF78", "GAUSS", "SABA4", "RADAU"};
    std::vector<BenchmarkResult> results;

    std::cout << "\nStarting Benchmark Loop...\n" << std::endl;
    std::cout << std::left << std::setw(12) << "Integrator" 
              << std::setw(12) << "Model" 
              << std::right << std::setw(18) << "Error (mas)" 
              << std::setw(3) << " " 
              << std::setw(15) << "Time (ms)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (const auto& int_name : integrators) {
        std::cout << ">> Integrator: " << int_name << " ..." << std::endl;
        for (const auto& model : models) {
            std::cout << "   Model: " << model.name << " ..." << std::endl;
            AstDynConfig cfg;
            cfg.propagator_settings = model.settings;
            cfg.integrator_type = int_name;
            cfg.ephemeris_file = bsp_path;
            cfg.initial_step_size = (int_name == "RK4" || int_name == "SABA4") ? 0.1 : 1.0;
            cfg.tolerance = 1e-13;
            cfg.ephemeris_type = "DE441";
            cfg.verbose = false;

            AstrometricSettings a_settings;
            a_settings.light_time_correction = true;
            a_settings.stellar_aberration = false;
            a_settings.frame_conversion_to_equatorial = true;

            // Run Propagate & Compare
            auto start = std::chrono::high_resolution_clock::now();
            auto obs_res = AstrometryReducer::compute_observation(initial_state_kep, start_time_tdb, target_time_tdb, cfg, a_settings);
            auto end = std::chrono::high_resolution_clock::now();
            double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

            auto calc_error = [&](auto res) -> double {
                if (!res) return -1.0;
                double ra_diff = (res->ra.value - truth_obs.ra.value);
                while (ra_diff > M_PI) ra_diff -= 2.0 * M_PI;
                while (ra_diff < -M_PI) ra_diff += 2.0 * M_PI;
                double dra = ra_diff * std::cos(truth_obs.dec.value);
                double ddec = res->dec.value - truth_obs.dec.value;
                return std::sqrt(dra*dra + ddec*ddec) * constants::RAD_TO_DEG * 3600.0 * 1000.0;
            };

            double error_mas = calc_error(obs_res);
            results.push_back({int_name, model.name, error_mas, elapsed_ms});

            std::cout << std::left << std::setw(12) << int_name 
                      << std::left << std::setw(12) << model.name 
                      << std::right << std::fixed << std::setprecision(3) << std::setw(18) << error_mas 
                      << std::setw(3) << " " 
                      << std::fixed << std::setprecision(2) << std::setw(15) << elapsed_ms << std::endl;
            
            // Diagnostic Run from Vectors (should be much better)
            if (model.name == "Full") {
                 auto obs_vec = AstrometryReducer::compute_observation_from_cartesian(initial_vectors_eq, start_time_tdb, target_time_tdb, cfg, a_settings);
                 if (obs_vec) {
                     double error_vec = calc_error(obs_vec);
                     std::cout << "    [DIAG] Calc:  RA=" << std::fixed << std::setprecision(8) << obs_vec->ra.value*constants::RAD_TO_DEG 
                               << ", DEC=" << obs_vec->dec.value*constants::RAD_TO_DEG << " deg" << std::endl;
                     std::cout << "    [DIAG] Truth: RA=" << std::fixed << std::setprecision(8) << truth_obs.ra.value*constants::RAD_TO_DEG 
                               << ", DEC=" << truth_obs.dec.value*constants::RAD_TO_DEG << " deg" << std::endl;
                     std::cout << "    [DIAG] From Vectors Error: " << std::fixed << std::setprecision(3) << error_vec << " mas" << std::endl;
                 }
            }
        }
    }

    // 7. Save to CSV
    std::ofstream csv(output_file);
    csv << "Integrator,Model,Error_mas,Time_ms\n";
    for (const auto& r : results) {
        csv << r.integrator << "," << r.model << "," << r.error_mas << "," << r.execution_time_ms << "\n";
    }
    csv.close();

    std::cout << "\n✅ Benchmark completed. Results saved to " << output_file << std::endl;

    return 0;
}
