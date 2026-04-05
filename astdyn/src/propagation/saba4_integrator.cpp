/**
 * @file saba4_integrator.cpp
 * @brief SABA4: Symplectic order-4 with adaptive step (RKF78-style)
 */

#include "astdyn/propagation/saba4_integrator.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::propagation {

using namespace astdyn::constants;

const double SABA4Integrator::c_coeffs_[5] = {
    0.1610444315367623,
    0.3389555684632377,
    0.0, // Calculated as 0.5 - (c1+c2) = 0
    0.3389555684632377,
    0.1610444315367623
};

const double SABA4Integrator::d_coeffs_[num_d_] = {
    0.5144243472499723,
    0.5 - 0.5144243472499723,
    0.5 - 0.5144243472499723,
    0.5144243472499723
};

SABA4Integrator::SABA4Integrator(double initial_step, double tolerance,
                                 double min_step, double max_step)
    : h_initial_(initial_step)
    , tolerance_(tolerance)
    , h_min_(min_step)
    , h_max_(max_step) {}

void SABA4Integrator::split_state(const Eigen::VectorXd& y,
                                  Eigen::VectorXd& q, Eigen::VectorXd& p) const {
    const int n_bodies = y.size() / 6;
    q.resize(n_bodies * 3);
    p.resize(n_bodies * 3);
    for (int i = 0; i < n_bodies; ++i) {
        q.segment<3>(i * 3) = y.segment<3>(i * 6);
        p.segment<3>(i * 3) = y.segment<3>(i * 6 + 3);
    }
}

Eigen::VectorXd SABA4Integrator::join_state(const Eigen::VectorXd& q,
                                            const Eigen::VectorXd& p) const {
    const int n_bodies = q.size() / 3;
    Eigen::VectorXd y(n_bodies * 6);
    for (int i = 0; i < n_bodies; ++i) {
        y.segment<3>(i * 6) = q.segment<3>(i * 3);
        y.segment<3>(i * 6 + 3) = p.segment<3>(i * 3);
    }
    return y;
}

Eigen::VectorXd SABA4Integrator::compute_force(const DerivativeFunction& f,
                                               double t,
                                               const Eigen::VectorXd& q,
                                               const Eigen::VectorXd& p) const {
    Eigen::VectorXd y = join_state(q, p);
    Eigen::VectorXd ydot = f(t, y);
    return ydot.tail(q.size());
}

Eigen::VectorXd SABA4Integrator::saba2_step(const DerivativeFunction& f,
                                            const Eigen::VectorXd& y,
                                            double t, double h) {
    Eigen::VectorXd q, p;
    split_state(y, q, p);

    Eigen::VectorXd force = compute_force(f, t, q, p);
    p += 0.5 * h * force;

    q += h * p;

    force = compute_force(f, t + h, q, p);
    p += 0.5 * h * force;

    stats_.num_function_evals += 2;
    return join_state(q, p);
}

Eigen::VectorXd SABA4Integrator::saba4_step(const DerivativeFunction& f,
                                            const Eigen::VectorXd& y,
                                            double t, double h) {
    // SABA4 (Blanes & Moan 2002): D(c1) K(d1) D(c2) K(d2) D(c3) K(d2) D(c2) K(d1) D(c1)
    Eigen::VectorXd q, p; split_state(y, q, p);
    auto kick = [&](double h_k) { 
        p += h_k * compute_force(f, t, q, p); 
        stats_.num_function_evals++;
    };
    auto drift = [&](double h_d) { 
        q += h_d * p; 
        t += h_d; 
    };

    drift(c_coeffs_[0] * h); kick(d_coeffs_[0] * h); 
    drift(c_coeffs_[1] * h); kick(d_coeffs_[1] * h);
    drift(c_coeffs_[2] * h); kick(d_coeffs_[2] * h); // Middle kick (d2) and no drift (c3=0)
    drift(c_coeffs_[3] * h); kick(d_coeffs_[3] * h);
    drift(c_coeffs_[4] * h);

    return join_state(q, p);
}

double SABA4Integrator::adapt_step_size(double h_current,
                                        double error_estimate,
                                        bool last_accepted) {
    double safety = 0.9;
    double order = 4.0;

    double scale = safety * std::pow(tolerance_ / std::max(error_estimate, 1e-30), 1.0 / (order + 1.0));
    scale = std::max(0.2, std::min(5.0, scale));

    double h_new = h_current * scale;
    h_new = std::max(h_min_, std::min(h_max_, h_new));

    return h_new;
}

double SABA4Integrator::compute_energy(const Eigen::VectorXd& y) const {
    // Energy computation only defined for 6D Keplerian state in unified AU/day
    if (y.size() != 6) return 0.0;

    Eigen::VectorXd q = y.head<3>();
    Eigen::VectorXd p = y.tail<3>();

    double r = q.norm();
    if (r < 1e-10) return 0.0;

    double kinetic = 0.5 * p.squaredNorm();
    double mu = constants::GMS;
    double potential = -mu / r;

    return kinetic + potential;
}

Eigen::VectorXd SABA4Integrator::integrate(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, double tf) {
    Integrator::reset_statistics();
    double t = t0, h = ((tf > t0) ? 1.0 : -1.0) * std::abs(h_initial_);
    Eigen::VectorXd y = y0; initial_energy_ = compute_energy(y0);
    while (std::abs(tf - t) > 1e-14 && stats_.num_steps < 5000000) {
        double current_h = (std::abs(tf - t) < std::abs(h)) ? (tf - t) : h;
        Eigen::VectorXd y_saba4 = saba4_step(f, y, t, current_h);
        Eigen::VectorXd y_saba2 = saba2_step(f, y, t, current_h);
        double err = (y_saba4 - y_saba2).norm() / std::max(std::abs(current_h), 1e-20);
        if (y_saba4.allFinite() && (err <= tolerance_ || std::abs(current_h) <= h_min_)) {
            y = y_saba4; t += current_h; stats_.num_steps++;
            h = adapt_step_size(current_h, err, true);
        } else {
            stats_.num_rejected_steps++; h = adapt_step_size(current_h, err, false);
        }
    }
    stats_.final_time = t; return y;
}

void SABA4Integrator::integrate_steps(const DerivativeFunction& f,
                                      const Eigen::VectorXd& y0,
                                      double t0, double tf,
                                      std::vector<double>& t_out,
                                      std::vector<Eigen::VectorXd>& y_out) {
    Integrator::reset_statistics();

    double t = t0;
    Eigen::VectorXd y = y0;
    double direction = (tf > t0) ? 1.0 : -1.0;
    double h = direction * std::abs(h_initial_);
    
    t_out.clear();
    y_out.clear();
    t_out.push_back(t);
    y_out.push_back(y);

    const int max_steps = 5000000;

    while (std::abs(tf - t) > 1e-14 && stats_.num_steps < max_steps) {
        double current_h = (std::abs(tf - t) < std::abs(h)) ? (tf - t) : h;

        Eigen::VectorXd y_saba4 = saba4_step(f, y, t, current_h);
        Eigen::VectorXd y_saba2 = saba2_step(f, y, t, current_h);

        double error = (y_saba4 - y_saba2).norm() / std::max(std::abs(current_h), 1e-20);

        if (y_saba4.allFinite() && (error <= tolerance_ || std::abs(current_h) <= h_min_)) {
            y = y_saba4;
            t += current_h;
            stats_.num_steps++;

            t_out.push_back(t);
            y_out.push_back(y);

            h = adapt_step_size(current_h, error, true);
        } else {
            stats_.num_rejected_steps++;
            h = adapt_step_size(current_h, error, false);
        }
    }

    stats_.final_time = t;
}

std::vector<Eigen::VectorXd> SABA4Integrator::integrate_at(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, const std::vector<double>& t_targets) {
    Integrator::reset_statistics();
    std::vector<Eigen::VectorXd> res; res.reserve(t_targets.size());
    double t = t0, h = ((t_targets.empty() || t_targets[0] >= t0) ? 1.0 : -1.0) * std::abs(h_initial_);
    Eigen::VectorXd y = y0; initial_energy_ = compute_energy(y0);
    for (double tf : t_targets) {
        while (std::abs(tf - t) > 1e-14 && stats_.num_steps < 5000000) {
            double current_h = (std::abs(tf - t) < std::abs(h)) ? (tf - t) : h;
            Eigen::VectorXd y_saba4 = saba4_step(f, y, t, current_h);
            Eigen::VectorXd y_saba2 = saba2_step(f, y, t, current_h);
            double err = (y_saba4 - y_saba2).norm() / std::max(std::abs(current_h), 1e-20);
            if (y_saba4.allFinite() && (err <= tolerance_ || std::abs(current_h) <= h_min_)) {
                y = y_saba4; t += current_h; stats_.num_steps++;
                h = adapt_step_size(current_h, err, true);
            } else {
                stats_.num_rejected_steps++; h = adapt_step_size(current_h, err, false);
            }
        }
        res.push_back(y);
    }
    stats_.final_time = t; return res;
}

} // namespace astdyn::propagation
