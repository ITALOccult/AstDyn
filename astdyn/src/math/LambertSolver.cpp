/**
 * @file LambertSolver.cpp
 * @brief Universal variables implementation for Lambert's problem.
 */

#include "astdyn/math/LambertSolver.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::math {

double LambertSolver::C(double z) {
    if (z > 1e-6)
        return (1.0 - std::cos(std::sqrt(z))) / z;
    else if (z < -1e-6)
        return (std::cosh(std::sqrt(-z)) - 1.0) / (-z);
    else
        return 1.0 / 2.0 - z / 24.0 + z * z / 720.0;
}

double LambertSolver::S(double z) {
    if (z > 1e-6) {
        double sz = std::sqrt(z);
        return (sz - std::sin(sz)) / (sz * sz * sz);
    } else if (z < -1e-6) {
        double sz = std::sqrt(-z);
        return (std::sinh(sz) - sz) / (sz * sz * sz);
    } else
        return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0;
}

Eigen::Vector3d LambertSolver::solve(
    const Eigen::Vector3d& r1_vec,
    const Eigen::Vector3d& r2_vec,
    double dt,
    double mu,
    bool retrograde)
{
    double r1 = r1_vec.norm();
    double r2 = r2_vec.norm();
    double cos_theta = r1_vec.dot(r2_vec) / (r1 * r2);

    if (retrograde) {
        if (r1_vec.cross(r2_vec).z() >= 0)
            cos_theta = -cos_theta; 
    }

    double A = std::sqrt(r1 * r2 * (1.0 + cos_theta));
    if (A == 0.0) return Eigen::Vector3d::Zero();

    // Initial guess for z
    double z = 0.0;
    double psi = 0.0;
    double psi_up = 4.0 * M_PI * M_PI;
    double psi_low = -4.0 * M_PI;

    for (int iter = 0; iter < 100; ++iter) {
        double c = C(psi);
        double s = S(psi);
        double y = r1 + r2 + A * (psi * s - 1.0) / std::sqrt(c);

        if (A > 0.0 && y < 0.0) {
            // Newton step or adjustment for negative y
            while (y < 0.0) {
                psi_low += 0.1;
                psi = psi_low;
                c = C(psi);
                s = S(psi);
                y = r1 + r2 + A * (psi * s - 1.0) / std::sqrt(c);
            }
        }

        double chi = std::sqrt(y / c);
        double t = (chi * chi * chi * s + A * std::sqrt(y)) / std::sqrt(mu);

        if (std::abs(t - dt) < 1e-12) break;

        if (t < dt)
            psi_low = psi;
        else
            psi_up = psi;
        
        psi = (psi_up + psi_low) / 2.0;
    }

    // Final result reconstruction
    double c = C(psi);
    double s = S(psi);
    double y = r1 + r2 + A * (psi * s - 1.0) / std::sqrt(c);
    
    double f = 1.0 - y / r1;
    double g = A * std::sqrt(y / mu);
    // double g_dot = 1.0 - y / r2;

    Eigen::Vector3d v1 = (r2_vec - f * r1_vec) / g;
    return v1;
}

} // namespace astdyn::math
