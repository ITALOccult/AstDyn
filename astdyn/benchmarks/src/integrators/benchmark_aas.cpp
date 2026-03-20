/**
 * @file benchmark_aas.cpp
 * @brief Numerical Validation Suite for AAS Integrator [M3]
 */

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <filesystem>

using namespace astdyn;
using namespace astdyn::physics;
using namespace astdyn::core;

namespace fs = std::filesystem;

// Ceres Initial State (JPL Horizons, 2024-01-01)
CartesianStateTyped<GCRF> get_ceres_initial() {
    auto epoch = time::EpochTDB::from_mjd(60310.0);
    auto pos = math::Vector3<GCRF, Distance>::from_si(
        -1.777322920405232 * (constants::AU * 1000.0),
         2.108027731776856 * (constants::AU * 1000.0),
         1.611110014769046 * (constants::AU * 1000.0)
    );
    auto vel = math::Vector3<GCRF, Velocity>::from_si(
        -8.841551069485747e-3 * (constants::AU * 1000.0 / 86400.0),
        -6.012176882200350e-3 * (constants::AU * 1000.0 / 86400.0),
        -3.882415174092120e-3 * (constants::AU * 1000.0 / 86400.0)
    );
    return CartesianStateTyped<GCRF>(epoch, pos, vel, GravitationalParameter::from_si(constants::GM_SUN * 1e9));
}

// 1566 Icarus (e=0.827)
CartesianStateTyped<GCRF> get_icarus_initial() {
    auto epoch = time::EpochTDB::from_mjd(60479.0);
    auto icarus = KeplerianStateTyped<GCRF>::from_traditional(
        epoch, 1.078, 0.827, 22.8, 88.0, 31.3, 0.0
    );
    return propagation::keplerian_to_cartesian<GCRF>(icarus);
}

void run_convergence_test() {
    std::cout << "Starting Convergence Test..." << std::endl;
    std::ofstream csv("benchmarks/convergence_test.csv");
    csv << "epsilon,relative_pos_error,nfe\n";

    auto initial = get_ceres_initial();
    auto target_epoch = time::EpochTDB::from_mjd(initial.epoch.mjd() + 1682.0); // 1 period

    // Reference: Analytical Two-Body Solution (Exact for this scenario)
    auto initial_kep = propagation::cartesian_to_keplerian(initial);
    auto ref_state_kep = propagation::TwoBodyPropagator::propagate(initial_kep, target_epoch);
    auto ref_state = propagation::keplerian_to_cartesian(ref_state_kep);

    auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>();
    propagation::PropagatorSettings settings;
    settings.include_planets = false;
    settings.central_body_gm = initial.gm.to_au3_d2(); // Ensure consistency

    std::vector<double> epsilons = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
    for (double eps : epsilons) {
        auto integrator = std::make_shared<propagation::AASIntegrator>(eps, constants::GMS);
        propagation::Propagator prop(integrator, ephem, settings);
        auto final_state = prop.propagate_cartesian(initial, target_epoch);
        
        double error = (final_state.position.to_eigen_si() - ref_state.position.to_eigen_si()).norm();
        double rel_error = error / ref_state.position.norm().to_m();
        csv << std::scientific << eps << "," << rel_error << "," << prop.statistics().num_function_evals << "\n";
        std::cout << "  eps=" << eps << " error=" << error << " m" << std::endl;
    }
    csv.close();
}

void run_work_precision_test() {
    std::cout << "Starting Work-Precision Test..." << std::endl;
    std::ofstream csv("benchmarks/work_precision.csv");
    csv << "integrator_name,nfe,final_error,eccentricity\n";

    auto run_comparison = [&](const std::string& name, CartesianStateTyped<GCRF> initial, double eccentricity) {
        auto target_epoch = time::EpochTDB::from_mjd(initial.epoch.mjd() + 365.25); // 1 year
        auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>();
        propagation::PropagatorSettings settings;
        settings.include_planets = false;

        // Reference
        auto integrator_ref = std::make_shared<propagation::AASIntegrator>(1e-12, constants::GMS);
        propagation::Propagator prop_ref(integrator_ref, ephem, settings);
        auto ref_state = prop_ref.propagate_cartesian(initial, target_epoch);

        // AAS
        std::vector<double> epsilons = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
        for (double eps : epsilons) {
            auto integrator = std::make_shared<propagation::AASIntegrator>(eps, constants::GMS);
            propagation::Propagator prop(integrator, ephem, settings);
            auto final_state = prop.propagate_cartesian(initial, target_epoch);
            double error = (final_state.position.to_eigen_si() - ref_state.position.to_eigen_si()).norm();
            csv << "AAS," << prop.statistics().num_function_evals << "," << error << "," << eccentricity << "\n";
        }

        // RK4 as baseline
        std::vector<double> steps = {1.0, 0.5, 0.1, 0.05, 0.01, 0.005};
        for (double h : steps) {
            AstDynConfig cfg;
            cfg.integrator_type = IntegratorType::RK4;
            cfg.initial_step_size = h;
            cfg.propagator_settings.include_planets = false;
            auto engine = std::make_unique<AstDynEngine>(cfg);
            engine->set_initial_orbit(propagation::cartesian_to_keplerian(initial));
            auto final_state = engine->propagate_to(target_epoch);
            auto final_cart = propagation::keplerian_to_cartesian(final_state);
            double error = (final_cart.position.to_eigen_si() - ref_state.position.to_eigen_si()).norm();
            
            // We need NFE for RK4. RK4 uses 4 evaluations per step.
            long nfe = static_cast<long>((365.25 / h) * 4);
            csv << "RK4," << nfe << "," << error << "," << eccentricity << "\n";
        }
    };

    run_comparison("Ceres", get_ceres_initial(), 0.07);
    run_comparison("Icarus", get_icarus_initial(), 0.82);
    csv.close();
}

void run_energy_conservation_test() {
    std::cout << "Starting Energy Conservation Test..." << std::endl;
    std::ofstream csv("benchmarks/energy_conservation.csv");
    csv << "time_years,delta_H,delta_H_shadow\n";

    auto initial = get_ceres_initial();
    double eps = 1e-5;
    auto integrator = std::make_shared<propagation::AASIntegrator>(eps, constants::GMS);
    auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>();
    propagation::PropagatorSettings settings;
    settings.include_planets = false;
    propagation::Propagator prop(integrator, ephem, settings);

    // Initial energy and shadow Hamiltonian
    Eigen::VectorXd q = initial.position.to_eigen_si();
    Eigen::VectorXd p = initial.velocity.to_eigen_si();
    
    // We need to access some private/protected stats or redo calculation
    // Since we are in a tool, we might need to recreate the physics.
    auto calc_energy = [&](const Eigen::VectorXd& rq, const Eigen::VectorXd& vp) {
        double r = rq.norm();
        double mu = initial.gm.to_m3_s2();
        return 0.5 * vp.squaredNorm() - mu / r;
    };

    double E0 = calc_energy(q, p);
    
    // Sequential propagation is MUCH faster than re-integrating from 0 to T each time
    double total_duration = 36525.0; // 100 years
    double step_days = 365.25; 
    
    CartesianStateTyped<GCRF> current = initial;
    for (double t_days = step_days; t_days <= total_duration; t_days += step_days) {
        auto target = time::EpochTDB::from_mjd(initial.epoch.mjd() + t_days);
        current = prop.propagate_cartesian(current, target);
        
        double Ef = calc_energy(current.position.to_eigen_si(), current.velocity.to_eigen_si());
        double dH = std::abs((Ef - E0) / E0);
        double dH_shadow = prop.statistics().shadow_hamiltonian_drift;
        
        csv << t_days / 365.25 << "," << dH << "," << dH_shadow << "\n";
    }
    csv.close();
}

void run_stm_validation() {
    std::cout << "Starting STM Validation..." << std::endl;
    std::ofstream csv("benchmarks/stm_validation.csv");
    csv << "time,det_error,diff_with_FD\n";

    auto initial_obj = get_icarus_initial(); // Use Icarus as NEA
    double eps = 1e-7;
    
    auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>();
    auto integrator = std::make_shared<propagation::AASIntegrator>(eps, constants::GMS);
    propagation::PropagatorSettings settings;
    settings.include_planets = false;
    propagation::Propagator prop(integrator, ephem, settings);

    // Initial state with identity STM
    Eigen::VectorXd y0(42);
    y0.segment<3>(0) = initial_obj.position.to_eigen_si();
    y0.segment<3>(3) = initial_obj.velocity.to_eigen_si();
    for(int i=0; i<6; ++i) {
        for(int j=0; j<6; ++j) y0[6 + i*6 + j] = (i==j ? 1.0 : 0.0);
    }

    auto f = [&](double t, const Eigen::VectorXd& y) -> Eigen::VectorXd {
        // Simple 2-body derivatives for stability
        Eigen::VectorXd dydt(6);
        Eigen::Vector3d r = y.head<3>();
        double r3 = r.norm() * r.norm() * r.norm();
        dydt.head<3>() = y.segment<3>(3);
        dydt.tail<3>() = -initial_obj.gm.to_m3_s2() * r / r3;
        return dydt;
    };

    double duration = 365.25 / 12.0; // 1 month is enough to check validity
    double dt = 1.0; 

    for (double t = 0; t <= duration; t += dt) {
        auto yf = integrator->integrate(f, y0, 0, t * 86400.0);
        
        Eigen::MatrixXd Phi(6, 6);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) Phi(i, j) = yf[6 + i * 6 + j];
        }

        double det_err = std::abs(Phi.determinant() - 1.0);
        
        // Finite Difference comparison (Simple)
        double h_fd = 100.0; // 100 meters
        Eigen::VectorXd y_plus = y0.head(6); y_plus[0] += h_fd;
        auto yf_plus = integrator->integrate(f, y_plus, 0, t * 86400.0);
        double fd_val = (yf_plus[0] - yf[0]) / h_fd;
        double stm_val = Phi(0, 0);
        double fd_diff = std::abs(stm_val - fd_val);

        csv << t << "," << det_err << "," << fd_diff << "\n";
    }
    csv.close();
}

int main() {
    fs::create_directories("benchmarks");
    
    try {
        run_convergence_test();
        run_work_precision_test();
        run_energy_conservation_test();
        run_stm_validation();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\nBenchmarks completed successfully. CSV files saved in /benchmarks\n";
    return 0;
}
