/**
 * @file test_integrators_comparison.cpp
 * @brief Test and compare RKF78, Radau15, and Gauss integrators
 * @author AstDyn Team
 * @date 2025-12-09
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <Eigen/Dense>

#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/RadauIntegrator.hpp"
#include "astdyn/propagation/GaussIntegrator.hpp"

using namespace astdyn::propagation;

// Test problem: Kepler orbit (2-body problem)
class KeplerProblem {
public:
    static constexpr double mu = 2.959122082855911e-4; // GM Sun [AU^3/day^2]
    
    // State: [x, y, z, vx, vy, vz]
    static Eigen::VectorXd derivative(double t, const Eigen::VectorXd& y) {
        Eigen::Vector3d r = y.head<3>();
        Eigen::Vector3d v = y.tail<3>();
        
        double r_norm = r.norm();
        double r3 = r_norm * r_norm * r_norm;
        
        Eigen::VectorXd dydt(6);
        dydt.head<3>() = v;
        dydt.tail<3>() = -mu * r / r3;
        
        return dydt;
    }
    
    // Orbital energy (should be conserved)
    static double energy(const Eigen::VectorXd& y) {
        Eigen::Vector3d r = y.head<3>();
        Eigen::Vector3d v = y.tail<3>();
        
        double r_norm = r.norm();
        double v2 = v.squaredNorm();
        
        return 0.5 * v2 - mu / r_norm;
    }
    
    // Angular momentum (should be conserved)
    static Eigen::Vector3d angular_momentum(const Eigen::VectorXd& y) {
        Eigen::Vector3d r = y.head<3>();
        Eigen::Vector3d v = y.tail<3>();
        
        return r.cross(v);
    }
};

void print_header() {
    std::cout << "========================================\n";
    std::cout << "  Integrator Comparison Test\n";
    std::cout << "========================================\n\n";
}

void test_integrator(const std::string& name,
                    Integrator& integrator,
                    const Eigen::VectorXd& y0,
                    double t0,
                    double tf) {
    std::cout << "Testing: " << name << "\n";
    std::cout << std::string(40, '-') << "\n";
    
    // Initial conditions
    double E0 = KeplerProblem::energy(y0);
    Eigen::Vector3d H0 = KeplerProblem::angular_momentum(y0);
    
    // Integrate
    auto start = std::chrono::high_resolution_clock::now();
    Eigen::VectorXd y_final = integrator.integrate(KeplerProblem::derivative, y0, t0, tf);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // Final conditions
    double E_final = KeplerProblem::energy(y_final);
    Eigen::Vector3d H_final = KeplerProblem::angular_momentum(y_final);
    
    // Conservation errors
    double dE = std::abs(E_final - E0);
    double dH = (H_final - H0).norm();
    
    // Statistics
    const auto& stats = integrator.statistics();
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Time:           " << duration.count() / 1000.0 << " ms\n";
    std::cout << "Steps:          " << stats.num_steps << "\n";
    std::cout << "Function evals: " << stats.num_function_evals << "\n";
    std::cout << "Rejected steps: " << stats.num_rejected_steps << "\n";
    std::cout << "Energy error:   " << dE << " (relative: " << dE/std::abs(E0) << ")\n";
    std::cout << "H error:        " << dH << " (relative: " << dH/H0.norm() << ")\n";
    std::cout << "\n";
}

int main() {
    print_header();
    
    // Initial conditions: circular orbit at 1 AU
    Eigen::VectorXd y0(6);
    y0 << 1.0, 0.0, 0.0,  // position [AU]
          0.0, 0.01720209895, 0.0;  // velocity [AU/day] for circular orbit
    
    double t0 = 0.0;
    double tf = 30.0;  // 30 days (faster test)
    
    std::cout << "Problem: Kepler orbit (2-body)\n";
    std::cout << "Initial position: " << y0.head<3>().transpose() << " AU\n";
    std::cout << "Initial velocity: " << y0.tail<3>().transpose() << " AU/day\n";
    std::cout << "Integration time: " << tf << " days (30 days)\n\n";
    
    // Test RKF78
    {
        RKF78Integrator rkf78(0.1, 1e-12);
        test_integrator("RKF78 (Explicit, Order 7/8)", rkf78, y0, t0, tf);
    }
    
    // Test Radau15
    {
        RadauIntegrator radau(0.1, 1e-13);
        test_integrator("Radau15 (Implicit, Order 15)", radau, y0, t0, tf);
    }
    
    // Test Gauss
    {
        GaussIntegrator gauss(0.1, 1e-12);
        test_integrator("Gauss-Legendre (Symplectic, Order 8)", gauss, y0, t0, tf);
    }
    
    std::cout << "========================================\n";
    std::cout << "Summary:\n";
    std::cout << "- RKF78: Fast, good accuracy for non-stiff problems\n";
    std::cout << "- Radau15: Highest accuracy, best for stiff problems\n";
    std::cout << "- Gauss: Best energy conservation (symplectic)\n";
    std::cout << "========================================\n";
    
    return 0;
}
