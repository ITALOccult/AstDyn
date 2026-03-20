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

math::Vector3<core::GCRF, physics::Velocity> LambertSolver::solve(
    const math::Vector3<core::GCRF, physics::Distance>& r1_vec,
    const math::Vector3<core::GCRF, physics::Distance>& r2_vec,
    time::TimeDuration dt,
    physics::GravitationalParameter mu,
    bool retrograde)
{
    // Convert to canonical units [AU, days] for internal solver stability
    double r1_au = r1_vec.norm().to_au();
    double r2_au = r2_vec.norm().to_au();
    double mu_au = mu.to_au3_d2();
    double dt_days = dt.to_days();

    Eigen::Vector3d r1v = r1_vec.to_eigen_si() * (1.0 / (constants::AU * 1000.0));
    Eigen::Vector3d r2v = r2_vec.to_eigen_si() * (1.0 / (constants::AU * 1000.0));

    double cos_theta = r1v.dot(r2v) / (r1_au * r2_au);

    if (retrograde) {
        if (r1v.cross(r2v).z() >= 0)
            cos_theta = -cos_theta; 
    }

    double A = std::sqrt(r1_au * r2_au * (1.0 + cos_theta));
    if (A == 0.0) return math::Vector3<core::GCRF, physics::Velocity>();

    // Initial guess for z
    double psi = 0.0;
    double psi_up = 4.0 * constants::PI * constants::PI;
    double psi_low = -4.0 * constants::PI;

    for (int iter = 0; iter < 100; ++iter) {
        double c = C(psi);
        double s = S(psi);
        double y = r1_au + r2_au + A * (psi * s - 1.0) / std::sqrt(c);

        if (A > 0.0 && y < 0.0) {
            // Newton step adjustment for negative y
            while (y < 0.0) {
                psi_low += 0.1;
                psi = psi_low;
                c = C(psi);
                s = S(psi);
                y = r1_au + r2_au + A * (psi * s - 1.0) / std::sqrt(c);
            }
        }

        double chi = std::sqrt(y / c);
        double t = (chi * chi * chi * s + A * std::sqrt(y)) / std::sqrt(mu_au);

        if (std::abs(t - dt_days) < 1e-12) break;

        if (t < dt_days)
            psi_low = psi;
        else
            psi_up = psi;
        
        psi = (psi_up + psi_low) / 2.0;
    }

    // Final result reconstruction
    double c = C(psi);
    double s = S(psi);
    double y = r1_au + r2_au + A * (psi * s - 1.0) / std::sqrt(c);
    
    double f = 1.0 - y / r1_au;
    double g = A * std::sqrt(y / mu_au);

    // v1 [AU/day]
    Eigen::Vector3d v1_aupd = (r2v - f * r1v) / g;
    
    // Convert back to SI [m/s]
    double v_scale = (constants::AU * 1000.0) / constants::SECONDS_PER_DAY;
    return math::Vector3<core::GCRF, physics::Velocity>::from_si(
        v1_aupd.x() * v_scale, 
        v1_aupd.y() * v_scale, 
        v1_aupd.z() * v_scale);
}

} // namespace astdyn::math
