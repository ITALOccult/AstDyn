/**
 * @file test_integrators_fast.cpp
 * @brief Fast test comparing integrators on appropriate problems
 * @date 2025-12-09
 */

#include <iostream>
#include <iomanip>
#include <chrono>
#include <Eigen/Dense>

#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/RadauIntegrator.hpp"

using namespace astdyn::propagation;

// Simple 2-body Kepler problem
Eigen::VectorXd kepler_derivative(double t, const Eigen::VectorXd& y) {
    static constexpr double mu = 2.959122082855911e-4;
    Eigen::Vector3d r = y.head<3>();
    Eigen::Vector3d v = y.tail<3>();
    double r3 = std::pow(r.norm(), 3);
    
    Eigen::VectorXd dydt(6);
    dydt.head<3>() = v;
    dydt.tail<3>() = -mu * r / r3;
    return dydt;
}

int main() {
    std::cout << "========================================\n";
    std::cout << "  Fast Integrator Test\n";
    std::cout << "========================================\n\n";
    
    // Initial conditions: circular orbit at 1 AU
    Eigen::VectorXd y0(6);
    y0 << 1.0, 0.0, 0.0, 0.0, 0.01720209895, 0.0;
    
    double t0 = 0.0;
    double tf = 10.0;  // Just 10 days for fast test
    
    std::cout << "Problem: Kepler orbit (2-body, non-stiff)\n";
    std::cout << "Integration time: " << tf << " days\n\n";
    
    // Test RKF78
    {
        std::cout << "RKF78 (Explicit):\n";
        RKF78Integrator rkf78(0.1, 1e-12);
        
        auto start = std::chrono::high_resolution_clock::now();
        auto y_final = rkf78.integrate(kepler_derivative, y0, t0, tf);
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        const auto& stats = rkf78.statistics();
        
        std::cout << "  Time: " << duration.count() / 1000.0 << " ms\n";
        std::cout << "  Steps: " << stats.num_steps << "\n";
        std::cout << "  Function evals: " << stats.num_function_evals << "\n\n";
    }
    
    // Test Radau15 with relaxed tolerance (faster)
    {
        std::cout << "Radau15 (Implicit, relaxed tolerance):\n";
        RadauIntegrator radau(1.0, 1e-10, 1e-6, 10.0, 3);  // Larger steps, relaxed tol, fewer iterations
        
        auto start = std::chrono::high_resolution_clock::now();
        auto y_final = radau.integrate(kepler_derivative, y0, t0, tf);
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        const auto& stats = radau.statistics();
        
        std::cout << "  Time: " << duration.count() / 1000.0 << " ms\n";
        std::cout << "  Steps: " << stats.num_steps << "\n";
        std::cout << "  Function evals: " << stats.num_function_evals << "\n\n";
    }
    
    std::cout << "========================================\n";
    std::cout << "Note: For non-stiff Kepler problems,\n";
    std::cout << "RKF78 is the better choice (faster).\n";
    std::cout << "Radau15 shines on stiff problems.\n";
    std::cout << "========================================\n";
    
    return 0;
}
