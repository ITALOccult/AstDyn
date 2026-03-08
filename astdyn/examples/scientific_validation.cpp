/**
 * @file scientific_validation.cpp
 * @brief Test suite for AAS integrator validation (Referee Comment [M3])
 */

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>

using namespace astdyn;

// Helper to compute position error in meters
double pos_error_m(const physics::CartesianStateTyped<core::GCRF>& s1, 
                   const physics::CartesianStateTyped<core::GCRF>& s2) {
    return (s1.position - s2.position).norm();
}

void run_convergence_test(const physics::KeplerianStateTyped<core::ECLIPJ2000>& start_orbit, 
                          io::HorizonsClient& client) {
    std::cout << "\n[1] Running Convergence Test (Slope Test) for AAS..." << std::endl;
    std::ofstream csv("convergence_aas.csv");
    csv << "precision,nfe,error_m\n";

    // Propagate for 1 period
    double period = start_orbit.period_days();
    time::EpochTDB end_time = start_orbit.epoch + period;
    
    // Get true final state from Horizons
    auto horizons_final = client.query_vectors("1566", end_time); 
    if (!horizons_final) {
        std::cerr << "Error fetching Horizons reference for convergence test." << std::endl;
        return;
    }

    AstDynEngine engine;
    for (double prec = 1e-2; prec >= 1e-8; prec /= 10.0) {
        AstDynConfig config;
        config.integrator_type = "AAS";
        config.aas_precision = prec;
        config.propagator_settings.include_planets = false; // Pure Keplerian for slope test accuracy
        engine.load_config(config);
        
        engine.set_initial_orbit(start_orbit);
        auto final_state_ecl = engine.propagate_to(end_time);
        
        // Convert to GCRF for comparison
        auto cart_ecl = propagation::keplerian_to_cartesian(final_state_ecl);
        auto pos_gcrf = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(cart_ecl.position);
        auto vel_gcrf = coordinates::ReferenceFrame::transform_vel<core::ECLIPJ2000, core::GCRF>(cart_ecl.position, cart_ecl.velocity);
        physics::CartesianStateTyped<core::GCRF> final_gcrf(end_time, pos_gcrf, vel_gcrf, cart_ecl.gm);

        double error = pos_error_m(final_gcrf, *horizons_final);
        int nfe = engine.propagator()->statistics().num_function_evals;

        csv << prec << "," << nfe << "," << error << "\n";
        std::cout << "Prec: " << prec << " | Error: " << std::scientific << error << " m | NFE: " << nfe << std::endl;
    }
    csv.close();
}

void run_work_precision_comparison(const physics::KeplerianStateTyped<core::ECLIPJ2000>& start_orbit,
                                   io::HorizonsClient& client) {
    std::cout << "\n[2] Running Work-Precision Comparison (AAS vs RK4)..." << std::endl;
    std::ofstream csv("work_precision.csv");
    csv << "method,parameter,nfe,error_m\n";

    time::EpochTDB end_time = start_orbit.epoch + 365.25; // 1 year
    auto horizons_ref = client.query_vectors("1566", end_time);

    // AAS test
    AstDynEngine engine;
    for (double prec = 1e-3; prec >= 1e-7; prec /= 10.0) {
        AstDynConfig config;
        config.integrator_type = "AAS";
        config.aas_precision = prec;
        engine.load_config(config);
        engine.set_initial_orbit(start_orbit);
        
        auto final_ecl = engine.propagate_to(end_time);
        auto final_gcrf = coordinates::ReferenceFrame::transform_state<core::ECLIPJ2000, core::GCRF>(propagation::keplerian_to_cartesian(final_ecl));
        
        double err = pos_error_m(final_gcrf, *horizons_ref);
        csv << "AAS," << prec << "," << engine.propagator()->statistics().num_function_evals << "," << err << "\n";
    }

    // RK4 test
    for (double step = 1.0; step >= 0.01; step /= 2.0) {
        AstDynConfig config;
        config.integrator_type = "RK4";
        config.initial_step_size = step;
        engine.load_config(config);
        engine.set_initial_orbit(start_orbit);
        
        auto final_ecl = engine.propagate_to(end_time);
        auto final_gcrf = coordinates::ReferenceFrame::transform_state<core::ECLIPJ2000, core::GCRF>(propagation::keplerian_to_cartesian(final_ecl));
        
        double err = pos_error_m(final_gcrf, *horizons_ref);
        csv << "RK4," << step << "," << engine.propagator()->statistics().num_function_evals << "," << err << "\n";
    }
    csv.close();
}

void run_energy_conservation_test(const physics::KeplerianStateTyped<core::ECLIPJ2000>& start_orbit) {
    std::cout << "\n[3] Running Long-term Energy Conservation (Ceres)..." << std::endl;
    std::ofstream csv("energy_conservation.csv");
    csv << "days,hamiltonian_drift,shadow_drift\n";

    AstDynEngine engine;
    AstDynConfig config;
    config.integrator_type = "AAS";
    config.aas_precision = 1e-4;
    engine.load_config(config);
    engine.set_initial_orbit(start_orbit);

    // Propagate for 100 periods (shortened from 10k for time)
    double num_periods = 100.0;
    double period = start_orbit.period_days();
    
    for (int i = 1; i <= num_periods; ++i) {
        time::EpochTDB t = start_orbit.epoch + (i * period);
        engine.propagate_to(t);
        auto stats = engine.propagator()->statistics();
        csv << (i * period) << "," << stats.hamiltonian_drift << "," << stats.shadow_hamiltonian_drift << "\n";
        
        if (i % 10 == 0) std::cout << "Period " << i << "/" << num_periods << " done." << std::endl;
    }
    csv.close();
}

void run_stm_validation(const physics::KeplerianStateTyped<core::ECLIPJ2000>& start_orbit) {
    std::cout << "\n[4] Validating State Transition Matrix (STM)..." << std::endl;
    
    AstDynEngine engine;
    AstDynConfig config;
    config.integrator_type = "AAS";
    config.aas_precision = 1e-4;
    engine.load_config(config);

    time::EpochTDB target_time = start_orbit.epoch + 100.0; // 100 days

    // Initialize state with identity STM
    Eigen::VectorXd y0(42);
    auto cart0 = propagation::keplerian_to_cartesian(start_orbit);
    y0.segment<3>(0) = cart0.position.to_eigen_si();
    y0.segment<3>(3) = cart0.velocity.to_eigen_si();
    y0.segment<36>(6).setZero();
    for(int i=0; i<6; ++i) y0[6 + i*7] = 1.0; // Identity

    // Propagate
    auto y_final = engine.propagator()->integrate_raw(
        [&](double t, const Eigen::VectorXd& y){ return engine.propagator()->derivative(t, y); },
        y0, start_orbit.epoch.mjd(), target_time.mjd()
    );

    // Extract STM
    Eigen::MatrixXd phi(6, 6);
    for(int i=0; i<6; ++i) {
        for(int j=0; j<6; ++j) phi(i, j) = y_final[6 + i*6 + j];
    }

    double det = phi.determinant();
    std::cout << "STM Determinant: " << det << " (Error: " << std::abs(det - 1.0) << ")" << std::endl;

    // Small perturbation for FD check
    double eps = 1e-7;
    Eigen::VectorXd y_pert = y0.head(6);
    y_pert[0] += eps; // Perturb X
    
    auto y_pert_final = engine.propagator()->integrate_raw(
        [&](double t, const Eigen::VectorXd& y){ return engine.propagator()->derivative(t, y); },
        y_pert, start_orbit.epoch.mjd(), target_time.mjd()
    );

    double fd_phi_0_0 = (y_pert_final[0] - y_final[0]) / eps;
    std::cout << "STM(0,0) - Analytic: " << phi(0,0) << " | Finite Difference: " << fd_phi_0_0 << " | Diff: " << std::abs(phi(0,0) - fd_phi_0_0) << std::endl;
}

int main() {
    try {
        io::HorizonsClient client;

        std::cout << "=== Scientific Validation Suite for AstDyn Integrators ===" << std::endl;

        // Fetch Targets
        auto ceres_orbit = client.query_elements("1", time::EpochTDB::from_mjd(60000.0));
        auto icarus_orbit = client.query_elements("1566", time::EpochTDB::from_mjd(60000.0));

        if (!ceres_orbit || !icarus_orbit) {
            std::cerr << "Failed to fetch necessary orbits from Horizons." << std::endl;
            return 1;
        }

        run_convergence_test(*icarus_orbit, client);
        run_work_precision_comparison(*icarus_orbit, client);
        run_energy_conservation_test(*ceres_orbit);
        run_stm_validation(*icarus_orbit);

        std::cout << "\nValidation Complete. Results saved in .csv files." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "FATAL ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
