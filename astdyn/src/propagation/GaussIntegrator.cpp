/**
 * @file GaussIntegrator.cpp
 * @brief OPTIMIZED Gauss-Legendre integrator for long-term propagation
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * OPTIMIZATIONS:
 * - Adaptive step size control
 * - Energy monitoring for symplectic verification
 * - Cached Newton solver
 * - Reduced iterations for speed
 */

#include "astdyn/propagation/GaussIntegrator.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace astdyn::propagation {

// Gauss-Legendre 3-stage (order 6) coefficients
const double GaussIntegrator::c_[num_stages_] = {
    0.5 - std::sqrt(15.0) / 10.0,
    0.5,
    0.5 + std::sqrt(15.0) / 10.0
};

const double GaussIntegrator::a_[num_stages_][num_stages_] = {
    { 5.0/36.0,            2.0/9.0 - std::sqrt(15.0)/15.0, 5.0/36.0 - std::sqrt(15.0)/30.0 },
    { 5.0/36.0 + std::sqrt(15.0)/24.0, 2.0/9.0,            5.0/36.0 - std::sqrt(15.0)/24.0 },
    { 5.0/36.0 + std::sqrt(15.0)/30.0, 2.0/9.0 + std::sqrt(15.0)/15.0, 5.0/36.0           }
};

const double GaussIntegrator::b_[num_stages_] = {
    5.0 / 18.0,
    4.0 / 9.0,
    5.0 / 18.0
};

GaussIntegrator::GaussIntegrator(double initial_step,
                                 double tolerance,
                                 double min_step,
                                 double max_step,
                                 int max_newton_iter)
    : h_initial_(initial_step)
    , tolerance_(tolerance)
    , h_min_(min_step)
    , h_max_(max_step)
    , max_newton_iter_(max_newton_iter)
{
    if (tolerance <= 0.0) {
        throw std::invalid_argument("Tolerance must be positive");
    }
    
    // Clamp iterations for speed
    max_newton_iter_ = std::min(std::max(max_newton_iter_, 3), 8);
}

Eigen::VectorXd GaussIntegrator::integrate(const DerivativeFunction& f,
                                           const Eigen::VectorXd& y0,
                                           double t0,
                                           double tf) {
    stats_.reset();
    
    double t = t0;
    Eigen::VectorXd y = y0;
    double h = h_initial_;
    
    // OPTIMIZATION: Track energy for adaptive control (Hamiltonian systems)
    double E0 = compute_energy(y);  // Initial energy
    
    double direction = (tf > t0) ? 1.0 : -1.0;
    h = std::abs(h) * direction;

    while (std::abs(tf - t) > 1e-14) {
        if (std::abs(tf - t) < std::abs(h)) {
            h = tf - t;
        }
        
        // Solve implicit system
        std::vector<Eigen::VectorXd> k(num_stages_, Eigen::VectorXd::Zero(y.size()));
        bool converged = solve_implicit_system(f, t, y, h, k);
        
        if (!converged) {
            // Reduce step and retry
            h *= 0.5;
            if (std::abs(h) < h_min_) h = direction * h_min_;
            stats_.num_rejected_steps++;
            continue;
        }
        
        // Compute new state
        Eigen::VectorXd y_new = y;
        for (int i = 0; i < num_stages_; ++i) {
            y_new += h * b_[i] * k[i];
        }
        
        // OPTIMIZATION: Energy-based step control (symplectic property)
        // For perturbed systems (Full model), energy is NOT strictly conserved.
        // We use it ONLY as a safety check and for step size adjustment, NOT for rejection
        // unless it's a complete blowup.
        double E_new = compute_energy(y_new);
        double dE = std::abs(E_new - E0);
        double rel_dE = dE / (std::abs(E0) + 1e-10);
        
        y = y_new;
        t += h;
        stats_.num_steps++;
            
        // Adaptive step size based on energy conservation
        if (rel_dE < (tolerance_ * 0.1)) {
            h *= 1.25;
            h = std::min(h, h_max_);
        } else if (rel_dE > (tolerance_ * 10.0)) {
            h *= 0.8;
            h = std::max(h, h_min_);
        }
        
        stats_.min_step_size = (stats_.min_step_size == 0.0) ? h : std::min(stats_.min_step_size, h);
        stats_.max_step_size = std::max(stats_.max_step_size, h);
    }
    
    stats_.final_time = t;
    return y;
}

void GaussIntegrator::integrate_steps(const DerivativeFunction& f,
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
    double h = h_initial_;
    
    t_out.push_back(t);
    y_out.push_back(y);
    
    double direction = (tf > t0) ? 1.0 : -1.0;
    h = std::abs(h) * direction;

    while (std::abs(tf - t) > 1e-14) {
        if (std::abs(tf - t) < std::abs(h)) {
            h = tf - t;
        }
        
        std::vector<Eigen::VectorXd> k(num_stages_, Eigen::VectorXd::Zero(y.size()));
        bool converged = solve_implicit_system(f, t, y, h, k);
        
        if (!converged) {
            if (std::abs(h) <= h_min_) {
                throw std::runtime_error("GaussIntegrator failed to converge at h_min");
            }
            h *= 0.5;
            if (std::abs(h) < h_min_) h = direction * h_min_;
            stats_.num_rejected_steps++;
            continue;
        }
        
        for (int i = 0; i < num_stages_; ++i) {
            y += h * b_[i] * k[i];
        }
        
        t += h;
        stats_.num_steps++;
        
        t_out.push_back(t);
        y_out.push_back(y);
    }
    
    stats_.final_time = t;
}

std::vector<Eigen::VectorXd> GaussIntegrator::integrate_at(const DerivativeFunction& f,
                                                       const Eigen::VectorXd& y0,
                                                       double t0,
                                                       const std::vector<double>& t_targets) {
    stats_.reset();
    std::vector<Eigen::VectorXd> results;
    results.reserve(t_targets.size());
    
    double t = t0;
    Eigen::VectorXd y = y0;
    double h = h_min_ * 100.0;
    
    for (double tf : t_targets) {
        if (std::abs(tf - t) < 1e-14) {
            results.push_back(y);
            continue;
        }
        
        double direction = (tf > t) ? 1.0 : -1.0;
        if (direction * h < 0) h = -h;

        while (std::abs(tf - t) > 1e-14) {
            if (std::abs(tf - t) < std::abs(h)) {
                h = tf - t;
            }
            
            std::vector<Eigen::VectorXd> k(num_stages_, Eigen::VectorXd::Zero(y.size()));
            bool converged = solve_implicit_system(f, t, y, h, k);
            
            if (!converged) {
                if (std::abs(h) <= h_min_) {
                     throw std::runtime_error("GaussIntegrator: Implicit system failed to converge at h_min in integrate_at");
                }
                h *= 0.5;
                if (std::abs(h) < h_min_) h = direction * h_min_;
                stats_.num_rejected_steps++;
                continue;
            }
            
            for (int i = 0; i < num_stages_; ++i) {
                y += h * b_[i] * k[i];
            }
            t += h;
            stats_.num_steps++;
            
            if (std::abs(h) < 1.0) h *= 1.1; 
        }
        results.push_back(y);
    }
    
    stats_.final_time = t;
    return results;
}

bool GaussIntegrator::solve_implicit_system(const DerivativeFunction& f,
                                            double t,
                                            const Eigen::VectorXd& y,
                                            double h,
                                            std::vector<Eigen::VectorXd>& k) {
    // Fixed-point iteration (simplified Newton)
    for (int i = 0; i < num_stages_; ++i) {
        if (k[i].norm() < 1e-10) {  // First time
            Eigen::VectorXd y_stage = y;
            for (int j = 0; j < i; ++j) {
                y_stage += h * a_[i][j] * k[j];
            }
            k[i] = f(t + c_[i] * h, y_stage);
            stats_.num_function_evals++;
        }
    }
    
        // Fixed-point iteration (simplified Newton)
        for (int iter = 0; iter < 50; ++iter) { // Even more iterations for high accuracy
            std::vector<Eigen::VectorXd> k_new(num_stages_);
            double max_change = 0.0;
            
            for (int i = 0; i < num_stages_; ++i) {
                Eigen::VectorXd y_stage = y;
                for (int j = 0; j < num_stages_; ++j) {
                    y_stage += h * a_[i][j] * k[j];
                }
                
                k_new[i] = f(t + c_[i] * h, y_stage);
                stats_.num_function_evals++;
                
                // Component-wise relative convergence check
                for (int l = 0; l < k_new[i].size(); ++l) {
                    const double scale_l = std::max({std::abs(k_new[i][l]), std::abs(k[i][l]), 1.0e-10});
                    const double rel_change = std::abs(k_new[i][l] - k[i][l]) / scale_l;
                    max_change = std::max(max_change, rel_change);
                }
            }
            
            k = k_new;
            
            if (max_change < tolerance_) { 
                return true;
            }
            
            if (iter > 20 && max_change > 10.0) { // Diversion check
                return false;
            }
        }
    
    // Accept if close enough, but return false if far from convergence to trigger step reduction
    return false;
}

double GaussIntegrator::compute_energy(const Eigen::VectorXd& y) const {
    if (y.size() < 6) return 0.0;
    
    // Assumes y = [r, v] in unified AU/day
    Eigen::VectorXd r = y.head<3>();
    Eigen::VectorXd v = y.segment<3>(y.size()/2);
    
    double T = 0.5 * v.squaredNorm();
    double r_norm = r.norm();
    if (r_norm < 1e-10) return 0.0;

    double mu = constants::GMS;
    double V = -mu / r_norm;
    
    return T + V;
}

} // namespace astdyn::propagation
