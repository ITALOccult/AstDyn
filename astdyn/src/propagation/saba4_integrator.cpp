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

const double SABA4Integrator::c_coeffs_[num_stages_] = {
    0.0792036964311957,
    0.353172906049774,
    -0.0420650803577195,
    0.2193769557534996,
    -0.0420650803577195,
    0.353172906049774,
    0.0792036964311957
};

const double SABA4Integrator::d_coeffs_[num_stages_] = {
    0.209515106613362,
    -0.143851773179818,
    0.434336666566456,
    0.0,
    0.434336666566456,
    -0.143851773179818,
    0.209515106613362
};

SABA4Integrator::SABA4Integrator(double initial_step, double tolerance,
                                 double min_step, double max_step)
    : h_initial_(initial_step)
    , tolerance_(tolerance)
    , h_min_(min_step)
    , h_max_(max_step) {}

void SABA4Integrator::split_state(const Eigen::VectorXd& y,
                                  Eigen::VectorXd& q, Eigen::VectorXd& p) const {
    const int n = y.size() / 2;
    q = y.head(n);
    p = y.tail(n);
}

Eigen::VectorXd SABA4Integrator::join_state(const Eigen::VectorXd& q,
                                            const Eigen::VectorXd& p) const {
    Eigen::VectorXd y(q.size() + p.size());
    y.head(q.size()) = q;
    y.tail(p.size()) = p;
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
    Eigen::VectorXd q, p;
    split_state(y, q, p);

    double t_stage = t;

    for (int i = 0; i < num_stages_; ++i) {
        if (c_coeffs_[i] != 0.0) {
            q += h * c_coeffs_[i] * p;
            t_stage += h * c_coeffs_[i];
        }

        if (d_coeffs_[i] != 0.0) {
            Eigen::VectorXd force = compute_force(f, t_stage, q, p);
            p += h * d_coeffs_[i] * force;
            stats_.num_function_evals++;
        }
    }

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

Eigen::VectorXd SABA4Integrator::integrate(const DerivativeFunction& f,
                                           const Eigen::VectorXd& y0,
                                           double t0, double tf) {
    Integrator::reset_statistics();

    double t = t0;
    Eigen::VectorXd y = y0;
    double h = h_initial_;

    initial_energy_ = compute_energy(y0);

    const int max_steps = 5000000;
    static const int progress_interval = 3000;
    double direction = (tf > t0) ? 1.0 : -1.0;
    h = std::abs(h) * direction;

    while (std::abs(tf - t) > 1e-14) {
        if (stats_.num_steps.load() >= max_steps) {
            std::cerr << "[SABA4] Warning: max steps " << max_steps << " reached.\n";
            break;
        }
        
        if (std::abs(tf - t) < std::abs(h)) h = tf - t;

        Eigen::VectorXd y_saba4 = saba4_step(f, y, t, h);
        Eigen::VectorXd y_saba2 = saba2_step(f, y, t, h);

        double error = (y_saba4 - y_saba2).norm() / std::max(h, 1e-20);

        if (error <= tolerance_) {
            y = y_saba4;
            t += h;
            stats_.num_steps++;

            double current_energy = compute_energy(y);
            stats_.hamiltonian_drift = std::abs(initial_energy_) > 1e-30
                ? std::abs(current_energy - initial_energy_) / std::abs(initial_energy_)
                : 0.0;

            h = adapt_step_size(h, error, true);
        } else {
            stats_.num_rejected_steps++;
            h = adapt_step_size(h, error, false);

            if (h < h_min_) {
                h = h_min_;
                y = y_saba4;
                t += h;
                stats_.num_steps++;
            }
        }
    }

    stats_.final_time = t;
    return y;
}

void SABA4Integrator::integrate_steps(const DerivativeFunction& f,
                                      const Eigen::VectorXd& y0,
                                      double t0, double tf,
                                      std::vector<double>& t_out,
                                      std::vector<Eigen::VectorXd>& y_out) {
    Integrator::reset_statistics();

    double t = t0;
    Eigen::VectorXd y = y0;
    double h = h_initial_;

    t_out.clear();
    y_out.clear();
    t_out.push_back(t);
    y_out.push_back(y);

    initial_energy_ = compute_energy(y0);

    const int max_steps = 5000000;
    while (t < tf) {
        if (stats_.num_steps.load() >= max_steps) break;
        if (t + h > tf) h = tf - t;

        Eigen::VectorXd y_saba4 = saba4_step(f, y, t, h);
        Eigen::VectorXd y_saba2 = saba2_step(f, y, t, h);

        double error = (y_saba4 - y_saba2).norm() / std::max(h, 1e-20);

        if (error <= tolerance_ || h <= h_min_) {
            y = y_saba4;
            t += h;
            stats_.num_steps++;

            t_out.push_back(t);
            y_out.push_back(y);

            h = adapt_step_size(h, error, true);
        } else {
            stats_.num_rejected_steps++;
            h = adapt_step_size(h, error, false);
        }
    }

    stats_.final_time = t;
}

std::vector<Eigen::VectorXd> SABA4Integrator::integrate_at(const DerivativeFunction& f,
                                                       const Eigen::VectorXd& y0,
                                                       double t0,
                                                       const std::vector<double>& t_targets) {
    Integrator::reset_statistics();
    std::vector<Eigen::VectorXd> results;
    results.reserve(t_targets.size());
    
    double t = t0;
    Eigen::VectorXd y = y0;
    double h = h_initial_;
    
    initial_energy_ = compute_energy(y0);
    
    for (double tf : t_targets) {
        if (std::abs(tf - t) < 1e-14) {
            results.push_back(y);
            continue;
        }
        
        while (std::abs(tf - t) > 1e-14) {
            if (std::abs(tf - t) < std::abs(h)) h = tf - t;
            
            Eigen::VectorXd y_saba4 = saba4_step(f, y, t, h);
            Eigen::VectorXd y_saba2 = saba2_step(f, y, t, h);
            double error = (y_saba4 - y_saba2).norm() / std::max(h, 1e-20);

            if (error <= tolerance_ || h <= h_min_) {
                y = y_saba4;
                t += h;
                stats_.num_steps++;
                h = adapt_step_size(h, error, true);
            } else {
                stats_.num_rejected_steps++;
                h = adapt_step_size(h, error, false);
            }
        }
        results.push_back(y);
    }
    
    stats_.final_time = t;
    return results;
}

} // namespace astdyn::propagation
