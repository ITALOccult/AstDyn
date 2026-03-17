/**
 * @file Integrator.cpp
 * @brief Implementation of numerical integrators
 */

#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/utils/Atomics.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace astdyn::propagation {

using namespace astdyn::utils;

// ============================================================================
// RK4Integrator Implementation
// ============================================================================

RK4Integrator::RK4Integrator(double step_size) : h_(step_size) {
    if (h_ <= 0.0) {
        throw std::invalid_argument("Step size must be positive");
    }
}

Eigen::VectorXd RK4Integrator::step(const DerivativeFunction& f,
                                     double t,
                                     const Eigen::VectorXd& y,
                                     double h) {
    // Classic RK4: k1, k2, k3, k4
    Eigen::VectorXd k1 = f(t, y);
    Eigen::VectorXd k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
    Eigen::VectorXd k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
    Eigen::VectorXd k4 = f(t + h, y + h * k3);
    
    stats_.num_function_evals += 4;
    
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

Eigen::VectorXd RK4Integrator::integrate(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, double tf) {
    stats_.reset();
    double t = t0; Eigen::VectorXd y = y0;
    double h_dir = ((tf > t0) ? 1.0 : -1.0) * std::abs(h_);
    stats_.min_step_size = stats_.max_step_size = std::abs(h_);
    while (std::abs(tf - t) > 1e-14) {
        double current_h = (std::abs(tf - t) < std::abs(h_dir)) ? (tf - t) : h_dir;
        y = step(f, t, y, current_h);
        t += current_h; stats_.num_steps++;
    }
    stats_.final_time = t; return y;
}

void RK4Integrator::integrate_steps(const DerivativeFunction& f,
                                    const Eigen::VectorXd& y0,
                                    double t0,
                                    double tf,
                                    std::vector<double>& t_out,
                                    std::vector<Eigen::VectorXd>& y_out) {
    stats_.reset();
    
    t_out.clear();
    y_out.clear();
    
    double t = t0;
    Eigen::VectorXd y = y0;
    
    // Store initial state
    t_out.push_back(t);
    y_out.push_back(y);
    
    double direction = (tf > t0) ? 1.0 : -1.0;
    double h = direction * std::abs(h_);
    
    stats_.min_step_size = std::abs(h);
    stats_.max_step_size = std::abs(h);
    
    while (std::abs(tf - t) > 1e-14) {
        if (std::abs(tf - t) < std::abs(h)) {
            h = tf - t;
        }
        
        y = step(f, t, y, h);
        t += h;
        stats_.num_steps++;
        
        t_out.push_back(t);
        y_out.push_back(y);
    }
    
    stats_.final_time = t;
}

// ============================================================================
// RKF78Integrator Implementation
// ============================================================================

// Fehlberg 7(8) coefficients (13 stages)
const double RKF78Integrator::c_[13] = {
    0.0,
    2.0/27.0,
    1.0/9.0,
    1.0/6.0,
    5.0/12.0,
    0.5,
    5.0/6.0,
    1.0/6.0,
    2.0/3.0,
    1.0/3.0,
    1.0,
    0.0,
    1.0
};

const double RKF78Integrator::a_[13][12] = {
    {},
    {2.0/27.0},
    {1.0/36.0, 1.0/12.0},
    {1.0/24.0, 0.0, 1.0/8.0},
    {5.0/12.0, 0.0, -25.0/16.0, 25.0/16.0},
    {1.0/20.0, 0.0, 0.0, 1.0/4.0, 1.0/5.0},
    {-25.0/108.0, 0.0, 0.0, 125.0/108.0, -65.0/27.0, 125.0/54.0},
    {31.0/300.0, 0.0, 0.0, 0.0, 61.0/225.0, -2.0/9.0, 13.0/900.0},
    {2.0, 0.0, 0.0, -53.0/6.0, 704.0/45.0, -107.0/9.0, 67.0/90.0, 3.0},
    {-91.0/108.0, 0.0, 0.0, 23.0/108.0, -976.0/135.0, 311.0/54.0, -19.0/60.0, 17.0/6.0, -1.0/12.0},
    {2383.0/4100.0, 0.0, 0.0, -341.0/164.0, 4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0},
    {3.0/205.0, 0.0, 0.0, 0.0, 0.0, -6.0/41.0, -3.0/205.0, -3.0/41.0, 3.0/41.0, 6.0/41.0, 0.0},
    {-1777.0/4100.0, 0.0, 0.0, -341.0/164.0, 4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0, 0.0, 1.0}
};

const double RKF78Integrator::b7_[13] = {
    41.0/840.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0,
    9.0/35.0, 9.0/280.0, 9.0/280.0, 41.0/840.0, 0.0, 0.0
};

const double RKF78Integrator::b8_[13] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0,
    9.0/35.0, 9.0/280.0, 9.0/280.0, 0.0, 41.0/840.0, 41.0/840.0
};

RKF78Integrator::RKF78Integrator(double initial_step,
                                 double tolerance,
                                 double min_step,
                                 double max_step)
    : h_initial_(initial_step),
      tolerance_(tolerance),
      h_min_(min_step),
      h_max_(max_step) {
    if (h_initial_ <= 0.0) {
        throw std::invalid_argument("Initial step size must be positive");
    }
    if (tolerance_ <= 0.0) {
        throw std::invalid_argument("Tolerance must be positive");
    }
}

void RKF78Integrator::compute_stages(const DerivativeFunction& f, double t, const Eigen::VectorXd& y, double h, std::vector<Eigen::VectorXd>& k) {
    k[0] = f(t, y); stats_.num_function_evals++;
    for (int i = 1; i < 13; ++i) {
        Eigen::VectorXd y_tmp = y;
        for (int j = 0; j < i; ++j) y_tmp += h * a_[i][j] * k[j];
        k[i] = f(t + c_[i] * h, y_tmp); stats_.num_function_evals++;
    }
}

double RKF78Integrator::estimate_error(const Eigen::VectorXd& y, const Eigen::VectorXd& y7, const Eigen::VectorXd& y8, double h) const {
    double err = 0.0;
    for (int i = 0; i < y.size(); ++i) {
        double scale = std::max({std::abs(y[i]), std::abs(y8[i]), 1.0});
        err = std::max(err, std::abs(y8[i] - y7[i]) / scale);
    }
    return err;
}

bool RKF78Integrator::adaptive_step(const DerivativeFunction& f, double& t, Eigen::VectorXd& y, double& h, double t_target) {
    std::vector<Eigen::VectorXd> k(13); compute_stages(f, t, y, h, k);
    Eigen::VectorXd y7 = y, y8 = y;
    for (int i = 0; i < 13; ++i) { y7 += h * b7_[i] * k[i]; y8 += h * b8_[i] * k[i]; }
    double rel_err = estimate_error(y, y7, y8, h);
    double dir = (h >= 0.0) ? 1.0 : -1.0;
    if (rel_err > tolerance_) {
        h = dir * std::max(h_min_, 0.9 * std::abs(h) * std::pow(tolerance_/rel_err, 0.125));
        stats_.num_rejected_steps++; return false;
    }
    y = y8; t += h; stats_.num_steps++;
    atomic_min(stats_.min_step_size, std::abs(h)); atomic_max(stats_.max_step_size, std::abs(h));
    h = dir * std::clamp(0.9 * std::abs(h) * std::pow(tolerance_/std::max(rel_err, 1e-20), 0.125), h_min_, h_max_);
    return true;
}

Eigen::VectorXd RKF78Integrator::integrate(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, double tf) {
    stats_.reset();
    double t = t0, h = ((tf > t0) ? 1.0 : -1.0) * std::abs(h_initial_);
    Eigen::VectorXd y = y0;
    int iteration_count = 0, consecutive_rejections = 0;
    while (std::abs(tf - t) > 1e-14) {
        verify_iteration_limits(++iteration_count, consecutive_rejections, t, h);
        double step_h = (std::abs(tf - t) < std::abs(h)) ? (tf - t) : h;
        if (adaptive_step(f, t, y, step_h, tf)) {
            consecutive_rejections = 0; h = step_h;
        } else {
            consecutive_rejections++; h = step_h;
        }
    }
    stats_.final_time = t; return y;
}

void RKF78Integrator::verify_iteration_limits(int count, int rejections, double t, double h) const {
    if (count > 1000000) throw std::runtime_error("RKF78Integrator: Max iterations exceeded");
    if (rejections >= 1000) throw std::runtime_error("RKF78Integrator: Too many consecutive rejections");
}

void RKF78Integrator::integrate_steps(const DerivativeFunction& f,
                                      const Eigen::VectorXd& y0,
                                      double t0,
                                      double tf,
                                      std::vector<double>& t_out,
                                      std::vector<Eigen::VectorXd>& y_out) {
    stats_.reset();
    stats_.min_step_size = std::abs(h_initial_);
    stats_.max_step_size = std::abs(h_initial_);
    
    t_out.clear();
    y_out.clear();
    
    double t = t0;
    Eigen::VectorXd y = y0;
    
    t_out.push_back(t);
    y_out.push_back(y);
    
    double direction = (tf > t0) ? 1.0 : -1.0;
    double h = direction * std::abs(h_initial_);
    
    // Safety limit to prevent infinite loops
    const int max_iterations = 1000000;
    int iteration_count = 0;
    int consecutive_rejections = 0;
    
    while (std::abs(tf - t) > 1e-14) {
        // Check iteration limit
        if (++iteration_count > max_iterations) {
            throw std::runtime_error(
                "RKF78Integrator: Maximum iterations exceeded in integrate_steps");
        }
        
        if (std::abs(tf - t) < std::abs(h)) {
            h = tf - t;
        }
        
        bool accepted = adaptive_step(f, t, y, h, tf);
        
        if (accepted) {
            t_out.push_back(t);
            y_out.push_back(y);
            consecutive_rejections = 0;
        } else {
            consecutive_rejections++;
            
            if (consecutive_rejections >= 1000) {
                throw std::runtime_error(
                    "RKF78Integrator: Too many consecutive step rejections in integrate_steps");
            }
        }
    }
    
    stats_.final_time = t;
}

std::vector<Eigen::VectorXd> RK4Integrator::integrate_at(const DerivativeFunction& f,
                                                     const Eigen::VectorXd& y0,
                                                     double t0,
                                                     const std::vector<double>& t_targets) {
    stats_.reset();
    std::vector<Eigen::VectorXd> results;
    results.reserve(t_targets.size());
    
    double t = t0;
    Eigen::VectorXd y = y0;
    
    for (double tf : t_targets) {
        if (std::abs(tf - t) < 1e-14) {
            results.push_back(y);
            continue;
        }
        
        double direction = (tf > t) ? 1.0 : -1.0;
        double h = direction * std::abs(h_);
        
        while (std::abs(tf - t) > 1e-14) {
            if (std::abs(tf - t) < std::abs(h)) {
                h = tf - t;
            }
            y = step(f, t, y, h);
            t += h;
            stats_.num_steps++;
        }
        results.push_back(y);
    }
    
    stats_.final_time = t;
    return results;
}

std::vector<Eigen::VectorXd> RKF78Integrator::integrate_at(const DerivativeFunction& f,
                                                       const Eigen::VectorXd& y0,
                                                       double t0,
                                                       const std::vector<double>& t_targets) {
    stats_.reset();
    std::vector<Eigen::VectorXd> results; results.reserve(t_targets.size());
    double t = t0, h = ((t_targets.empty() || t_targets[0] >= t0) ? 1.0 : -1.0) * std::abs(h_initial_);
    Eigen::VectorXd y = y0; int total_iterations = 0;
    for (double tf : t_targets) {
        int rejections = 0;
        while (std::abs(tf - t) > 1e-14) {
            verify_iteration_limits(++total_iterations, rejections, t, h);
            double step_h = (std::abs(tf - t) < std::abs(h)) ? (tf - t) : h;
            if (adaptive_step(f, t, y, step_h, tf)) { rejections = 0; h = step_h; }
            else { rejections++; h = step_h; }
        }
        results.push_back(y);
    }
    stats_.final_time = t;
    return results;
}

} // namespace astdyn::propagation
