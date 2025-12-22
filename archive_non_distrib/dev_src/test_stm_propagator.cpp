/**
 * @file test_stm_propagator.cpp
 * @brief Test State Transition Matrix propagator
 */

#include <iostream>
#include <iomanip>
#include "astdyn/propagation/STMPropagator.hpp"
#include "astdyn/propagation/Integrator.hpp"

using namespace astdyn::propagation;

// Simple 2-body problem
Eigen::Vector<double, 6> kepler_force(double t, const Eigen::Vector<double, 6>& x) {
    constexpr double mu = 2.959122082855911e-4;  // GM Sun
    
    Eigen::Vector3d r = x.head<3>();
    Eigen::Vector3d v = x.tail<3>();
    
    double r3 = std::pow(r.norm(), 3);
    Eigen::Vector3d a = -mu * r / r3;
    
    Eigen::Vector<double, 6> dxdt;
    dxdt.head<3>() = v;
    dxdt.tail<3>() = a;
    
    return dxdt;
}

int main() {
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║          STM Propagator Test                               ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    // Initial state: circular orbit at 1 AU
    Eigen::Vector<double, 6> x0;
    x0 << 1.0, 0.0, 0.0,  // position [AU]
          0.0, 0.01720209895, 0.0;  // velocity [AU/day]
    
    // Create integrator (RKF78)
    auto integrator = std::make_unique<RKF78Integrator>(0.1, 1e-12);
    
    // Create STM propagator
    STMPropagator stm_prop(std::move(integrator), kepler_force);
    
    // Propagate 10 days
    double t0 = 0.0;
    double tf = 10.0;
    
    std::cout << "Propagating from t=" << t0 << " to t=" << tf << " days...\n\n";
    
    auto result = stm_prop.propagate(x0, t0, tf);
    
    // Display results
    std::cout << "Final state:\n";
    std::cout << "  Position: " << result.state.head<3>().transpose() << " AU\n";
    std::cout << "  Velocity: " << result.state.tail<3>().transpose() << " AU/day\n\n";
    
    std::cout << "State Transition Matrix:\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << result.stm << "\n\n";
    
    // Verify STM properties
    std::cout << "STM Verification:\n";
    
    // 1. Determinant should be 1 (volume preservation)
    double det = result.stm.determinant();
    std::cout << "  det(Φ) = " << std::scientific << det << " (should be ≈1)\n";
    
    // 2. Φ Φ⁻¹ should be identity
    auto identity_check = result.stm * result.stm.inverse();
    double identity_error = (identity_check - Eigen::Matrix<double, 6, 6>::Identity()).norm();
    std::cout << "  ||Φ Φ⁻¹ - I|| = " << identity_error << " (should be ≈0)\n\n";
    
    // Statistics
    std::cout << "Integration statistics:\n";
    std::cout << "  Steps: " << result.stats.num_steps << "\n";
    std::cout << "  Function evals: " << result.stats.num_function_evals << "\n";
    std::cout << "  Rejected steps: " << result.stats.num_rejected_steps << "\n";
    
    std::cout << "\n✓ STM Propagator test complete!\n";
    
    return 0;
}
