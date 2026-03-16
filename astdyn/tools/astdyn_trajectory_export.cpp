/**
 * @file astdyn_trajectory_export.cpp
 * @brief Command-line tool for exporting asteroid trajectories with high-precision physics.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <filesystem>

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

int main(int argc, char** argv) {
    po::options_description desc("astdyn_trajectory_export: Export trajectory to CSV\nUsage: astdyn_trajectory_export --asteroid <id> [options]");
    desc.add_options()
        ("help,h", "produce help message")
        ("asteroid", po::value<std::string>(), "asteroid ID (number or designation)")
        ("t0", po::value<double>()->default_value(60310.0), "initial epoch (MJD TDB)")
        ("tf", po::value<double>(), "final epoch (MJD TDB)")
        ("step", po::value<double>()->default_value(30.0), "output step (days)")
        ("integrator", po::value<std::string>()->default_value("AAS"), "AAS | RKF78 | GL8")
        ("tolerance", po::value<double>()->default_value(1e-4), "tolerance (for RKF78/GL8) or precision (for AAS)")
        ("forces", po::value<std::string>(), "full | twobody (shortcut to toggle all)")
        ("output", po::value<std::string>()->default_value("trajectory.csv"), "output CSV file")
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

    std::string asteroid_id = vm["asteroid"].as<std::string>();
    double t0_mjd = vm["t0"].as<double>();
    double tf_mjd = vm["tf"].as<double>();
    double step_days = vm["step"].as<double>();
    std::string output_file = vm["output"].as<std::string>();
    std::string integrator_name = vm["integrator"].as<std::string>();
    double tol = vm["tolerance"].as<double>();

    // --- 1. Fetch Initial State from Horizons ---
    HorizonsClient horizons;
    time::EpochTDB t0 = time::EpochTDB::from_mjd(t0_mjd);
    
    std::cout << "[trajectory_export] Fetching initial state for " << asteroid_id << " at MJD " << t0_mjd << " (SSB) from Horizons...\n";
    auto state_res = horizons.query_vectors(asteroid_id, t0, "@0"); // @0 = SSB
    if (!state_res) {
        std::cerr << "Error: Failed to fetch state from Horizons: " << error_to_string(state_res.error()) << "\n";
        return 1;
    }
    
    // Convert to AU and AU/day (Horizons returns km and km/s in GCRF/SSB)
    Eigen::VectorXd y0_au = state_res->to_eigen_au_aud();

    // --- 2. Setup Ephemeris and Propagator ---
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
    } catch (const std::exception& e) {
        std::cerr << "Error: Ephemeris setup failed: " << e.what() << "\n";
        return 1;
    }
    auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>(de441);

    PropagatorSettings settings;
    if (vm.count("forces") && vm["forces"].as<std::string>() == "twobody") {
        settings.include_planets = false;
        settings.include_moon = false;
        settings.include_relativity = false;
        settings.include_asteroids = false;
        settings.include_sun_j2 = false;
        settings.include_earth_j2 = false;
    } else {
        settings.include_planets = true;
        settings.include_moon = true;
        settings.include_relativity = vm["relativity"].as<bool>();
        settings.include_asteroids = vm["asteroids"].as<bool>();
        settings.include_sun_j2 = vm["sun-j2"].as<bool>();
        settings.include_earth_j2 = vm["earth-j2"].as<bool>();
        
        if (settings.include_asteroids) {
            std::string set = vm["asteroid-set"].as<std::string>();
            if (set == "17") {
                settings.use_default_asteroid_set = true;
                settings.use_default_30_set = false;
            } else if (set == "30") {
                settings.use_default_asteroid_set = false;
                settings.use_default_30_set = true;
            }
            
            // Default asteroid ephemeris
            settings.asteroid_ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/sb441-n16.bsp";
            if (!std::filesystem::exists(settings.asteroid_ephemeris_file)) {
                settings.asteroid_ephemeris_file = "sb441-n16.bsp";
            }
        }
        
        settings.baricentric_integration = true;
    }

    std::shared_ptr<Integrator> integrator;
    if (integrator_name == "AAS") {
        auto aas = std::make_shared<AASIntegrator>(tol);
        integrator = aas;
    } else if (integrator_name == "RKF78") {
        integrator = std::make_shared<RKF78Integrator>(0.1, tol);
    } else {
        std::cerr << "Unsupported integrator: " << integrator_name << "\n";
        return 1;
    }

    Propagator propagator(integrator, ephem, settings);

    // --- 3. Energy Calculation Helper ---
    auto compute_h = [&](const Eigen::VectorXd& state, double mjd) -> double {
        Eigen::Vector3d r = state.head<3>();
        Eigen::Vector3d v = state.segment<3>(3);
        time::EpochTDB t = time::EpochTDB::from_mjd(mjd);
        
        // Sun-relative for GM_sun/r part
        auto sun_bary = de441->getPosition(ephemeris::CelestialBody::SUN, t).to_eigen_si() * 1e-3 * KM_TO_AU;
        Eigen::Vector3d r_sun_rel = r - sun_bary;
        
        double energy = 0.5 * v.squaredNorm() - GMS / r_sun_rel.norm();
        
        if (settings.include_planets) {
            std::vector<ephemeris::CelestialBody> bodies = {
                ephemeris::CelestialBody::MERCURY, ephemeris::CelestialBody::VENUS,
                ephemeris::CelestialBody::EARTH, ephemeris::CelestialBody::MARS,
                ephemeris::CelestialBody::JUPITER, ephemeris::CelestialBody::SATURN,
                ephemeris::CelestialBody::URANUS, ephemeris::CelestialBody::NEPTUNE,
                ephemeris::CelestialBody::MOON
            };
            for (auto body : bodies) {
                auto planet_bary = de441->getPosition(body, t).to_eigen_si() * 1e-3 * KM_TO_AU;
                double gm = 0;
                switch(body) {
                    case ephemeris::CelestialBody::MERCURY: gm = GM_MERCURY_AU; break;
                    case ephemeris::CelestialBody::VENUS: gm = GM_VENUS_AU; break;
                    case ephemeris::CelestialBody::EARTH: gm = GM_EARTH_AU; break;
                    case ephemeris::CelestialBody::MARS: gm = GM_MARS_AU; break;
                    case ephemeris::CelestialBody::JUPITER: gm = GM_JUPITER_AU; break;
                    case ephemeris::CelestialBody::SATURN: gm = GM_SATURN_AU; break;
                    case ephemeris::CelestialBody::URANUS: gm = GM_URANUS_AU; break;
                    case ephemeris::CelestialBody::NEPTUNE: gm = GM_NEPTUNE_AU; break;
                    case ephemeris::CelestialBody::MOON: gm = GM_MOON_AU; break;
                    default: break;
                }
                energy -= gm / (r - planet_bary).norm();
            }
        }
        return energy;
    };

    // --- 4. Propagation and Export Loop ---
    std::ofstream out(output_file);
    out << "t_mjd,x_au,y_au,z_au,vx_auday,vy_auday,vz_auday,energy_rel,det_stm\n";
    out << std::fixed << std::setprecision(15);

    double h0 = compute_h(y0_au, t0_mjd);
    Eigen::VectorXd current_y = y0_au;
    double current_t = t0_mjd;

    // Output initial state
    out << current_t << "," << current_y[0] << "," << current_y[1] << "," << current_y[2] << ","
        << current_y[3] << "," << current_y[4] << "," << current_y[5] << "," << 0.0 << "," << 1.0 << "\n";

    while (current_t < tf_mjd) {
        double next_t = std::min(current_t + step_days, tf_mjd);
        try {
            // Print progress
            std::cout << "\r[trajectory_export] Propagating... MJD " << std::fixed << std::setprecision(2) << current_t << " / " << tf_mjd << std::flush;

            // Propagate in AU/day frame (Baricentric)
            current_y = propagator.integrate_raw_au(current_y, current_t, next_t);
            current_t = next_t;
            
            double h = compute_h(current_y, current_t);
            double energy_rel = std::abs((h - h0) / h0);
            
            if (std::isnan(energy_rel) || energy_rel > 1.0) {
                 std::cerr << "Error: Energy drift exceeded limit or NaN at t=" << current_t << " (H0=" << h0 << ", H=" << h << ")\n";
                 // Write last corrupted line but don't proceed
                 out << current_t << "," << current_y[0] << "," << current_y[1] << "," << current_y[2] << ","
                     << current_y[3] << "," << current_y[4] << "," << current_y[5] << "," << energy_rel << "," << 1.0 << "\n";
                 return 1;
            }

            out << current_t << "," << current_y[0] << "," << current_y[1] << "," << current_y[2] << ","
                << current_y[3] << "," << current_y[4] << "," << current_y[5] << "," << energy_rel << "," << 1.0 << "\n";
        } catch (const std::exception& e) {
            std::cerr << "Integration error: " << e.what() << "\n";
            return 1;
        }
    }
    std::cout << "\n";

    std::cout << "[trajectory_export] Export complete: " << output_file << "\n";
    return 0;
}
