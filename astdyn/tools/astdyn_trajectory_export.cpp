/**
 * @file astdyn_trajectory_export.cpp
 * @brief Command-line tool for exporting asteroid trajectories with high-precision physics.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/propagation/GRKNIntegrator.hpp"
#include "astdyn/propagation/GaussIntegrator.hpp"
#include "astdyn/propagation/saba4_integrator.hpp"
#include "astdyn/propagation/RadauIntegrator.hpp"
#include <boost/program_options.hpp>
#include <chrono>
#include <atomic>
#include <mutex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <thread>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::io;
using namespace astdyn::physics;
using namespace astdyn::constants;
namespace po = boost::program_options;

std::string error_to_string(HorizonsError err) {
    switch (err) {
        case HorizonsError::NetworkError: return "NetworkError";
        case HorizonsError::InvalidResponse: return "InvalidResponse";
        case HorizonsError::TargetNotFound: return "TargetNotFound";
        case HorizonsError::ParsingError: return "ParsingError";
        case HorizonsError::Timeout: return "Timeout";
        default: return "UnknownError";
    }
}

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    po::options_description desc("astdyn_trajectory_export: Export trajectory to CSV\nUsage: astdyn_trajectory_export --asteroid <id1> <id2> ... [options]");
    desc.add_options()
        ("help,h", "produce help message")
        ("asteroid", po::value<std::vector<std::string>>()->multitoken(), "asteroid ID(s) (number or designation)")
        ("t0", po::value<double>()->default_value(60310.0), "initial epoch (MJD TDB)")
        ("tf", po::value<double>(), "final epoch (MJD TDB)")
        ("step", po::value<double>()->default_value(30.0), "output step (days)")
        ("integrator", po::value<std::string>()->default_value("AAS"), "AAS | RKF78 | GL8 | GRKN | SABA4 | IAS15")
        ("tolerance", po::value<double>()->default_value(1e-4), "tolerance (for RKF78/GL8) or precision (for AAS)")
        ("forces", po::value<std::string>(), "full | twobody (shortcut to toggle all)")
        ("output", po::value<std::string>(), "output CSV file (or prefix for multiple asteroids)")
        ("sun-j2", po::value<bool>()->default_value(true), "include Sun J2")
        ("earth-j2", po::value<bool>()->default_value(true), "include Earth J2")
        ("relativity", po::value<bool>()->default_value(true), "include General Relativity")
        ("asteroids", po::value<bool>()->default_value(true), "include asteroid perturbations")
        ("asteroid-set", po::value<std::string>()->default_value("17"), "17 (AstDyn set) | 30 (BC405 set)")
        ("ephem", po::value<std::string>(), "path to planetary ephemeris BSP")
        ("stm", "include State Transition Matrix integration")
    ;

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const std::exception& e) {
        std::cerr << "Command line error: " << e.what() << "\n";
        return 1;
    }

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    if (!vm.count("asteroid") || !vm.count("tf")) {
        std::cerr << "Error: --asteroid and --tf are mandatory.\n";
        return 1;
    }

    auto asteroid_ids = vm["asteroid"].as<std::vector<std::string>>();
    double t0_mjd = vm["t0"].as<double>();
    double tf_mjd = vm["tf"].as<double>();
    double step_days = vm["step"].as<double>();
    std::string integrator_name = vm["integrator"].as<std::string>();
    double tol = vm["tolerance"].as<double>();

    // --- setup Ephemeris (Global) ---
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    if (vm.count("ephem")) {
        bsp_path = vm["ephem"].as<std::string>();
    } else if (!std::filesystem::exists(bsp_path)) {
        bsp_path = "de441.bsp"; 
    }
    
    std::shared_ptr<ephemeris::DE441Provider> de441;
    try {
        de441 = std::make_shared<ephemeris::DE441Provider>(bsp_path);
        ephemeris::PlanetaryEphemeris::setGlobalProvider(de441);
    } catch (...) {
        std::cerr << "Error: Ephemeris setup failed.\n"; return 1;
    }
    auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>(de441);

    PropagatorSettings settings;
    if (vm.count("forces") && vm["forces"].as<std::string>() == "twobody") {
        settings.include_planets = settings.include_moon = settings.include_relativity = settings.include_asteroids = settings.include_sun_j2 = settings.include_earth_j2 = false;
    } else {
        settings.include_planets = settings.include_moon = true;
        settings.include_relativity = vm["relativity"].as<bool>();
        settings.include_asteroids = vm["asteroids"].as<bool>();
        settings.include_sun_j2 = vm["sun-j2"].as<bool>();
        settings.include_earth_j2 = vm["earth-j2"].as<bool>();
        if (settings.include_asteroids) {
            std::string set = vm["asteroid-set"].as<std::string>();
            settings.use_default_asteroid_set = (set == "17");
            settings.use_default_30_set = (set == "30");
            settings.asteroid_ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/sb441-n16.bsp";
            if (!std::filesystem::exists(settings.asteroid_ephemeris_file)) settings.asteroid_ephemeris_file = "sb441-n16.bsp";
        }
        settings.baricentric_integration = true;
    }

    std::cout << "[trajectory_export] Fetching initial states for " << asteroid_ids.size() << " asteroids using Horizons (serial)...\n";
    struct AsteroidJob {
        std::string id;
        Eigen::VectorXd y0;
        double t0;
    };
    std::vector<AsteroidJob> jobs;

    for (const auto& id : asteroid_ids) {
        time::EpochTDB t0 = time::EpochTDB::from_mjd(t0_mjd);
        HorizonsClient horizons_local;
        
        bool success = false;
        for (int retry = 0; retry < 3; ++retry) {
            auto res = horizons_local.query_vectors(id, t0, "@0");
            if (res) {
                jobs.push_back({id, res->to_eigen_au_aud(), t0_mjd});
                success = true;
                break;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
        
        if (!success) {
            std::cerr << "Warning: Failed to fetch state for " << id << "\n";
        }
    }

    std::cout << "[trajectory_export] Processing " << jobs.size() << " asteroids in parallel...\n";

    long total_planned = 0;
    for (const auto& job : jobs) {
        total_planned += (long)((tf_mjd - job.t0) / step_days);
    }
    if (total_planned == 0) total_planned = 1;

    std::atomic<long> total_steps_done{0};
    std::mutex progress_mutex;
    auto start_time = std::chrono::steady_clock::now();

    std::cout << "[Batch Progress] 0.0% | Elapsed: 0s | ETA: ---s" << std::flush;

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)jobs.size(); ++i) {
        const auto& job = jobs[i];
        std::string asteroid_id = job.id;
        Eigen::VectorXd current_y = job.y0;
        double current_t = job.t0;

        std::string output_file = vm.count("output") ? vm["output"].as<std::string>() : (asteroid_id + ".csv");
        if (jobs.size() > 1 && vm.count("output")) output_file = vm["output"].as<std::string>() + "_" + asteroid_id + ".csv";

        try {
            // Local providers per asteroid to avoid lock contention
            auto local_de441 = std::make_shared<ephemeris::DE441Provider>(bsp_path);
            auto local_ephem = std::make_shared<ephemeris::PlanetaryEphemeris>(local_de441);

            std::shared_ptr<Integrator> integrator;
            if (integrator_name == "AAS") integrator = std::make_shared<AASIntegrator>(tol);
            else if (integrator_name == "RKF78") integrator = std::make_shared<RKF78Integrator>(0.1, tol);
            else if (integrator_name == "GL8") integrator = std::make_shared<GaussIntegrator>(0.01, tol); // Smaller initial step for stability
            else if (integrator_name == "GRKN") integrator = std::make_shared<GRKNIntegrator>(tol, 0.001); // Even smaller for GRKN
            else if (integrator_name == "SABA4") integrator = std::make_shared<SABA4Integrator>(0.01, tol);
            else if (integrator_name == "SABA2") integrator = std::make_shared<SABA4Integrator>(0.01, tol); 
            else if (integrator_name == "IAS15") integrator = std::make_shared<RadauIntegrator>(0.01, tol);
            else continue;

            PropagatorSettings local_settings = settings;
            // Avoid self-perturbation if target is an asteroid in the default list
            try {
                if (!asteroid_id.empty() && std::isdigit(asteroid_id[0])) {
                    int ast_num = std::stoi(asteroid_id);
                    local_settings.exclude_asteroids_list.push_back(ast_num);
                }
            } catch (...) {}

            Propagator propagator(integrator, local_ephem, local_settings);
            
            auto compute_h = [&](const Eigen::VectorXd& state, double mjd) {
                Eigen::Vector3d r = state.head<3>(), v = state.segment<3>(3);
                time::EpochTDB t = time::EpochTDB::from_mjd(mjd);
                auto sun_bary = local_de441->getPosition(ephemeris::CelestialBody::SUN, t).to_eigen_si() * 1e-3 * KM_TO_AU;
                double energy = 0.5 * v.squaredNorm() - GMS / (r - sun_bary).norm();
                if (settings.include_planets) {
                    static const ephemeris::CelestialBody bodies[] = { ephemeris::CelestialBody::MERCURY, ephemeris::CelestialBody::VENUS, ephemeris::CelestialBody::EARTH, ephemeris::CelestialBody::MARS, ephemeris::CelestialBody::JUPITER, ephemeris::CelestialBody::SATURN, ephemeris::CelestialBody::URANUS, ephemeris::CelestialBody::NEPTUNE, ephemeris::CelestialBody::MOON };
                    for (auto b : bodies) {
                        auto p_bary = local_de441->getPosition(b, t).to_eigen_si() * 1e-3 * KM_TO_AU;
                        energy -= ephemeris::PlanetaryEphemeris::planet_gm(b) / (r - p_bary).norm();
                    }
                }
                return energy;
            };

            std::ofstream out(output_file);
            out << std::fixed << std::setprecision(15) << "t_mjd,x_au,y_au,z_au,vx_auday,vy_auday,vz_auday,energy_rel,det_stm\n";
            double h0 = compute_h(current_y, current_t);
            out << current_t << "," << current_y[0] << "," << current_y[1] << "," << current_y[2] << "," << current_y[3] << "," << current_y[4] << "," << current_y[5] << ",0.0,1.0\n";

            // Sliced integration for progress reporting
            double t_target = current_t + step_days;
            while (t_target < tf_mjd + 1e-10) {
                double target = std::min(t_target, tf_mjd);
                current_y = propagator.integrate_raw_au(current_y, current_t, target);
                current_t = target;
                
                double h = compute_h(current_y, current_t);
                out << current_t << "," << current_y[0] << "," << current_y[1] << "," << current_y[2] << "," 
                    << current_y[3] << "," << current_y[4] << "," << current_y[5] << "," << std::abs((h - h0) / h0) << ",1.0\n";
                
                total_steps_done++;
                
                // Report progress frequently (every step)
                if (total_steps_done % 1 == 0) { // Report every step for better feedback
                    std::lock_guard<std::mutex> lock(progress_mutex);
                    auto now = std::chrono::steady_clock::now();
                    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count() / 1000.0;
                    double percent = 100.0 * total_steps_done / total_planned;
                    if (percent > 100.0) percent = 100.0;
                    double eta = (percent > 0.1) ? (elapsed / percent * (100.0 - percent)) : 0.0;
                    
                    std::cout << "\r[Batch Progress] " << std::fixed << std::setprecision(1) << percent << "% "
                              << "| Elapsed: " << (int)elapsed << "s "
                              << "| ETA: " << (int)eta << "s   " << std::flush;
                }
                
                t_target += step_days;
            }
            #pragma omp critical
            std::cout << "\r[trajectory_export]   -> Asteroid " << asteroid_id << " complete.                                \n";
        } catch (const std::exception& e) {
            #pragma omp critical
            std::cerr << "\n[trajectory_export]   -> Error in " << asteroid_id << ": " << e.what() << "\n";
        }
    }

    std::cout << "\r[Batch Progress] 100.0% | Elapsed: " << (int)(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.0) << "s | ETA: 0s      \n";
    std::cout << "[trajectory_export] All exports complete.\n";
    return 0;
}
