/**
 * @file Integrator.cpp
 * @brief Implementation of numerical integrators
 */

#include "astdyn/propagation/Integrator.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace astdyn::propagation {

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

Eigen::VectorXd RK4Integrator::integrate(const DerivativeFunction& f,
                                          const Eigen::VectorXd& y0,
                                          double t0,
                                          double tf) {
    stats_.reset();
    
    double t = t0;
    Eigen::VectorXd y = y0;
    
    double direction = (tf > t0) ? 1.0 : -1.0;
    double h = direction * std::abs(h_);
    
    stats_.min_step_size = std::abs(h);
    stats_.max_step_size = std::abs(h);
    
    while (std::abs(tf - t) > 1e-14) {
        // Don't overshoot final time
        if (std::abs(tf - t) < std::abs(h)) {
            h = tf - t;
        }
        
        y = step(f, t, y, h);
        t += h;
        stats_.num_steps++;
    }
    
    stats_.final_time = t;
    return y;
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

bool RKF78Integrator::adaptive_step(const DerivativeFunction& f,
                                    double& t,
                                    Eigen::VectorXd& y,
                                    double& h,
                                    double t_target) {
    int n = y.size();
    std::vector<Eigen::VectorXd> k(13);
    
    // Determine direction for backward/forward integration
    double direction = (h >= 0.0) ? 1.0 : -1.0;
    
    // Compute stages
    k[0] = f(t, y);
    stats_.num_function_evals++;
    
    for (int i = 1; i < 13; ++i) {
        Eigen::VectorXd y_temp = y;
        for (int j = 0; j < i; ++j) {
            y_temp += h * a_[i][j] * k[j];
        }
        k[i] = f(t + c_[i] * h, y_temp);
        stats_.num_function_evals++;
    }
    
    // 7th order solution
    Eigen::VectorXd y7 = y;
    for (int i = 0; i < 13; ++i) {
        y7 += h * b7_[i] * k[i];
    }
    
    // 8th order solution
    Eigen::VectorXd y8 = y;
    for (int i = 0; i < 13; ++i) {
        y8 += h * b8_[i] * k[i];
    }
    
    // Error estimate using mixed component-wise norm for robust scaling
    double relative_error = 0.0;
    for (int i = 0; i < n; ++i) {
        const double scale_i = std::max({std::abs(y[i]), std::abs(y8[i]), 1.0});
        const double error_i = std::abs(y8[i] - y7[i]) / scale_i;
        relative_error = std::max(relative_error, error_i);
    }
    
    // Step size control
    double safety_factor = 0.9;
    double h_new_abs = std::abs(h);
    
    if (relative_error > tolerance_) {
        // Reject step, reduce step size
        h_new_abs = safety_factor * std::abs(h) * std::pow(tolerance_ / relative_error, 1.0 / 8.0);
        h_new_abs = std::max(h_new_abs, h_min_);
        h = direction * h_new_abs;  // Preserve direction
        stats_.num_rejected_steps++;
        return false;
    } else {
        // Accept step
        y = y8; // Use higher-order solution
        t += h;
        stats_.num_steps++;
        
        // Increase step size for next iteration
        if (relative_error > 0.0) {
            h_new_abs = safety_factor * std::abs(h) * std::pow(tolerance_ / relative_error, 1.0 / 8.0);
        } else {
            h_new_abs = h_max_;
        }
        h_new_abs = std::min(h_new_abs, h_max_);
        h_new_abs = std::max(h_new_abs, h_min_);
        
        // Update statistics
        stats_.min_step_size = std::min(stats_.min_step_size, std::abs(h));
        stats_.max_step_size = std::max(stats_.max_step_size, std::abs(h));
        
        h = direction * h_new_abs;  // Preserve direction
        return true;
    }
}

Eigen::VectorXd RKF78Integrator::integrate(const DerivativeFunction& f,
                                            const Eigen::VectorXd& y0,
                                            double t0,
                                            double tf) {
    stats_.reset();
    stats_.min_step_size = std::abs(h_initial_);
    stats_.max_step_size = std::abs(h_initial_);
    
    double t = t0;
    Eigen::VectorXd y = y0;
    
    double direction = (tf > t0) ? 1.0 : -1.0;
    double h = direction * std::abs(h_initial_);
    
    // Safety limit to prevent infinite loops
    const int max_iterations = 1000000;  // 1 million iterations
    int iteration_count = 0;
    int consecutive_rejections = 0;
    
    while (std::abs(tf - t) > 1e-14) {
        // Check iteration limit
        if (++iteration_count > max_iterations) {
            throw std::runtime_error(
                "RKF78Integrator: Maximum iterations (" + std::to_string(max_iterations) + 
                ") exceeded. Integration from t=" + std::to_string(t0) + 
                " to t=" + std::to_string(tf) + " failed at t=" + std::to_string(t) + 
                ". Steps: " + std::to_string(stats_.num_steps) + 
                ", Rejections: " + std::to_string(stats_.num_rejected_steps) +
                ", Current h: " + std::to_string(h));
        }
        
        // Don't overshoot target
        if (std::abs(tf - t) < std::abs(h)) {
            h = tf - t;
        }
        
        // Attempt step (may be rejected)
        bool accepted = adaptive_step(f, t, y, h, tf);
        
        if (accepted) {
            consecutive_rejections = 0;
        } else {
            consecutive_rejections++;
            
            // Warn if too many consecutive rejections (likely a problem)
            if (consecutive_rejections >= 1000) {
                throw std::runtime_error(
                    "RKF78Integrator: Too many consecutive step rejections (" + 
                    std::to_string(consecutive_rejections) + 
                    "). Integration may be stuck. Current t=" + std::to_string(t) + 
                    ", h=" + std::to_string(h) + ", h_min=" + std::to_string(h_min_));
            }
        }
        
        // If step was rejected, t and y are unchanged, h is reduced
        // Loop will retry with smaller h
    }
    
    stats_.final_time = t;
    return y;
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
    std::vector<Eigen::VectorXd> results;
    results.reserve(t_targets.size());
    
    double t = t0;
    Eigen::VectorXd y = y0;
    
    // Initial guess for step size
    double h = (t_targets.empty() || t_targets[0] >= t0) ? std::abs(h_initial_) : -std::abs(h_initial_);
    
    const int max_iterations = 1000000;
    int total_iterations = 0;

    for (double tf : t_targets) {
        if (std::abs(tf - t) < 1e-14) {
            results.push_back(y);
            continue;
        }
        
        // Adjust direction if necessary
        double direction = (tf > t) ? 1.0 : -1.0;
        if (direction * h < 0) h = -h;

        int consecutive_rejections = 0;
        
        while (std::abs(tf - t) > 1e-14) {
            if (++total_iterations > max_iterations) {
                throw std::runtime_error("RKF78Integrator: Maximum iterations exceeded in integrate_at");
            }
            
            // Don't overshoot target, but save old h to restore it after hitting the target
            double old_h = h;
            bool target_limited = false;
            if (std::abs(tf - t) < std::abs(h)) {
                h = tf - t;
                target_limited = true;
            }
            
            bool accepted = adaptive_step(f, t, y, h, tf);
            
            if (accepted) {
                consecutive_rejections = 0;
                // If we were target-limited, we might want to scale h back up for the next target
                // but adaptive_step already suggests a new h. 
                // However, h was tf-t. So suggested_h might be small.
                // For now, let's just let adaptive_step do its thing.
            } else {
                consecutive_rejections++;
                if (consecutive_rejections >= 1000) {
                    throw std::runtime_error("RKF78Integrator: Too many consecutive rejections in integrate_at");
                }
            }
        }
        results.push_back(y);
    }
    
    stats_.final_time = t;
    results.shrink_to_fit();
    return results;
}

} // namespace astdyn::propagation
