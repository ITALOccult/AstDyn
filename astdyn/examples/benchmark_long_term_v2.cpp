/**
 * @file benchmark_long_term_v2.cpp
 * @brief AAS vs SABA4 vs RKF78 — 9-body heliocentric mutual integration (8 planets + Asteroid).
 * Includes Full Force Model: Mutual Gravity, Sun J2, Sun-Asteroid 1PN.
 */
#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/propagation/saba4_integrator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <filesystem>
#include <memory>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::constants;

// --- Data Structures ---
struct AsteroidSpec {
    std::string name;
    std::string family;
    Eigen::VectorXd y0_hel;
    double T_yr;
    int lyap_N;
};

// Planet GMs (excluding Sun)
static const std::vector<double> PLANET_GMS = {
    GM_MERCURY_AU, GM_VENUS_AU, GM_EARTH_MOON_AU, GM_MARS_AU,
    GM_JUPITER_AU, GM_SATURN_AU, GM_URANUS_AU, GM_NEPTUNE_AU,
    0.0 // 8: Asteroid
};
static const std::vector<std::string> PLANET_CODES = { "199", "299", "3", "499", "599", "699", "799", "899" };

// --- Force Model (Heliocentric Mutual) ---
DerivativeFunction make_heliocentric_dynamics() {
    return [](double /*t*/, const Eigen::VectorXd& y) {
        const int n = 9; // 8 planets + 1 asteroid
        Eigen::VectorXd dy = Eigen::VectorXd::Zero(54);
        
        for (int i = 0; i < n; ++i) {
            Eigen::Vector3d ri = y.segment<3>(i * 6);
            Eigen::Vector3d vi = y.segment<3>(i * 6 + 3);
            dy.segment<3>(i * 6) = vi;
            
            // 1. Primary Gravity (Sun at origin)
            double r_mag = ri.norm();
            double r3 = r_mag * r_mag * r_mag;
            dy.segment<3>(i * 6 + 3) = -GMS * ri / r3;
            
            // 2. Solar J2
            double r2 = r_mag * r_mag;
            double z = ri.z();
            double j2_factor = -1.5 * GMS * SUN_J2 * std::pow(R_SUN_AU, 2) / (r2 * r2 * r_mag);
            dy.segment<3>(i * 6 + 3) += j2_factor * Eigen::Vector3d(
                ri.x() * (1.0 - 5.0 * z * z / r2),
                ri.y() * (1.0 - 5.0 * z * z / r2),
                ri.z() * (3.0 - 5.0 * z * z / r2)
            );
            
            // 3. Relativistic Correction (1PN, Sun-centric, on all for consistency)
            double c2 = constants::SPEED_OF_LIGHT_AU_PER_DAY * constants::SPEED_OF_LIGHT_AU_PER_DAY;
            double v2 = vi.squaredNorm();
            double r_dot_v = ri.dot(vi);
            double gr_factor = GMS / (c2 * r3);
            dy.segment<3>(i * 6 + 3) += gr_factor * ((4.0 * GMS / r_mag - v2) * ri + 4.0 * r_dot_v * vi);
            
                // 4. Planetary Perturbations (Mutual + Indirect)
                for (int j = 0; j < 8; ++j) {
                    if (i == j) continue;
                    Eigen::Vector3d rj = y.segment<3>(j * 6);
                    Eigen::Vector3d dr = ri - rj;
                    double dr_mag = dr.norm();
                    double rj_mag = rj.norm();
                    if (dr_mag > 1e-12 && rj_mag > 1e-12) {
                        dy.segment<3>(i * 6 + 3) += PLANET_GMS[j] * (-dr / (dr_mag*dr_mag*dr_mag) - rj / (rj_mag*rj_mag*rj_mag));
                    }
                }
        }
        return dy;
    };
}

// --- Energy Diagnostic (Heliocentric with Indirect Term) ---
double energy_relative_asteroid(const Eigen::VectorXd& y) {
    // y is 54-dim heliocentric. Asteroid is at index 8 (segment 48)
    Eigen::Vector3d r_ast = y.segment<3>(48);
    Eigen::Vector3d v_ast = y.segment<3>(51);
    double r = r_ast.norm();
    double v2 = v_ast.squaredNorm();
    double E = 0.5 * v2 - GMS / r;
    
    for (int j = 0; j < 8; ++j) {
        Eigen::Vector3d rj = y.segment<3>(j * 6);
        double rp = (r_ast - rj).norm();
        double rs = rj.norm();
        E -= PLANET_GMS[j] / rp;
        E += PLANET_GMS[j] / rs; // Indirect term (consistent with force model)
    }
    return E;
}

// --- Benchmark Suite ---
std::vector<AsteroidSpec> load_bench_states_v2(const std::string& path) {
    std::vector<AsteroidSpec> asteroids;
    std::ifstream file(path); std::string line; std::getline(file, line);
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line); std::string name, num, x, y, z, vx, vy, vz;
        std::getline(ss, name, ','); std::getline(ss, num, ',');
        std::getline(ss, x, ','); std::getline(ss, y, ','); std::getline(ss, z, ',');
        std::getline(ss, vx, ','); std::getline(ss, vy, ','); std::getline(ss, vz, ',');
        AsteroidSpec a; a.name = name; a.y0_hel.resize(6);
        a.y0_hel << std::stod(x), std::stod(y), std::stod(z), std::stod(vx), std::stod(vy), std::stod(vz);
        if (name == "Apophis" || name == "Icarus" || name == "Phaethon") { a.family = "NEA"; a.T_yr = 50; a.lyap_N = 5; }
        else if (name == "Achilles" || name == "Patroclus" || name == "Hektor") { a.family = "Trojan"; a.T_yr = 1000; a.lyap_N = 100; }
        else if (name == "Hilda" || name == "Thule" || name == "Griqua") { a.family = "Resonant"; a.T_yr = 500; a.lyap_N = 50; }
        else { a.family = "TNO"; a.T_yr = 10000; a.lyap_N = 1000; }
        asteroids.push_back(a);
    }
    return asteroids;
}

int main(int argc, char** argv) {
    std::string ast_filter = "";
    double T_filter = -1.0;
    for(int i=1; i<argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--asteroid" && i+1 < argc) ast_filter = argv[++i];
        if (arg == "--T" && i+1 < argc) T_filter = std::stod(argv[++i]);
    }

    std::cout << "[LTv2] 9-Body Heliocentric Mutual Integration Benchmark" << std::endl;
    
    // 1. Initial State Setup (Heliocentric)
    io::HorizonsClient hzn;
    Eigen::VectorXd y0_planets = Eigen::VectorXd::Zero(48); // 8 planets
    for (size_t i = 0; i < 8; ++i) {
        auto res = hzn.query_vectors(PLANET_CODES[i], time::EpochTDB::from_mjd(60310.0), "500@10"); // Sun-centered
        if (res) {
            y0_planets.segment<6>(i * 6) = res->to_eigen_au_aud();
            std::cout << "  Fetched " << PLANET_CODES[i] << std::endl;
        }
    }

    auto asteroids = load_bench_states_v2("../examples/benchmark_results/initial_states_mjd60310.csv");
    std::ofstream out("../examples/benchmark_results/long_term_tests.csv");
    out << "asteroid,integrator,test,T_yr,value,unit,note\n";
    
    auto f = make_heliocentric_dynamics();
    struct IntSpec { std::string name; double prec; };
    std::vector<IntSpec> specs = { {"AAS", 1e-5}, {"SABA4", 0.0}, {"RKF78", 1e-12} };

    for (const auto& a : asteroids) {
        if (!ast_filter.empty() && a.name != ast_filter) continue;
        double T_run = (T_filter > 0) ? T_filter : a.T_yr;
        std::cout << "Processing " << a.name << " (T=" << T_run << " yr)..." << std::endl;

        // Construct 54-dim system state
        Eigen::VectorXd y_sys = Eigen::VectorXd::Zero(54);
        y_sys.head<48>() = y0_planets;
        y_sys.tail<6>() = a.y0_hel;

        for (const auto& spec : specs) {
            std::cout << "  [" << spec.name << "] " << std::flush;
            
            std::unique_ptr<Integrator> integ;
            if (spec.name == "AAS") {
                // AAS needs 9 GMs for 9 bodies
                std::vector<double> gms_9(9, GMS);
                integ = std::make_unique<AASIntegrator>(spec.prec, gms_9);
            } else if (spec.name == "SABA4") {
                integ = std::make_unique<SABA4Integrator>(0.1, 1e-12, 0.1, 0.1);
            } else {
                integ = std::make_unique<RKF78Integrator>(0.1, spec.prec, 1e-6, 10.0);
            }

            // Forward (+ Energy)
            double E0 = energy_relative_asteroid(y_sys);
            Eigen::VectorXd y_final = integ->integrate(f, y_sys, 0.0, T_run * 365.25);
            double dE_final = std::abs((energy_relative_asteroid(y_final) - E0) / E0);

            // Lyapunov (Shadow Trajectory, 1e-9 shadow)
            Eigen::VectorXd y_shadow = y_sys;
            y_shadow.segment<3>(48) += Eigen::Vector3d(1e-9, 0, 0); // Perturb asteroid
            Eigen::VectorXd y_shadow_final = integ->integrate(f, y_shadow, 0.0, T_run * 365.25);
            double dist = (y_shadow_final.segment<3>(48) - y_final.segment<3>(48)).norm();
            double lyap = std::log(dist / 1e-9) / (T_run);

            // Reversibility (Hamiltonian: reverse v, integrate FORWARD)
            Eigen::VectorXd y_rev = y_final;
            for (int i = 0; i < 9; ++i) y_rev.segment<3>(i * 6 + 3) *= -1.0;
            Eigen::VectorXd y_back = integ->integrate(f, y_rev, 0.0, T_run * 365.25);
            for (int i = 0; i < 9; ++i) y_back.segment<3>(i * 6 + 3) *= -1.0;
            
            double rev_r = (y_back.segment<3>(48) - y_sys.segment<3>(48)).norm() / y_sys.segment<3>(48).norm();
            double rev_v = (y_back.segment<3>(51) - y_sys.segment<3>(51)).norm() / y_sys.segment<3>(51).norm();
            
            std::cout << " dE=" << std::scientific << std::setprecision(2) << dE_final << " rev=" << rev_r << " lyap=" << lyap << std::endl;

            out << a.name << "," << spec.name << ",energy_final," << T_run << "," << dE_final << ",dimensionless,\n";
            out << a.name << "," << spec.name << ",reversibility_r," << T_run << "," << rev_r << ",dimensionless,\n";
            out << a.name << "," << spec.name << ",reversibility_v," << T_run << "," << rev_v << ",dimensionless,\n";
            out << a.name << "," << spec.name << ",lyapunov," << T_run << "," << lyap << ",1/yr,\n";
            out.flush();
        }
    }
    return 0;
}
