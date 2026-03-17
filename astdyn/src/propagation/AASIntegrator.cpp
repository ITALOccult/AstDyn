/**
 * @file AASIntegrator.cpp
 * @brief AAS Adaptive Symplectic Integrator Implementation (Time-Transformed)
 */

#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/utils/Atomics.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::propagation {

AASIntegrator::AASIntegrator(double precision, double mu, double j2, double r_eq)
    : precision_(precision),
      mu_(mu),
      j2_(j2),
      r_eq_(r_eq) 
{
    // Exact Yoshida-4 coefficients derived from 2*w1 + w0 = 1 and 2*w1^3 + w0^3 = 0
    const double k = std::pow(2.0, 1.0/3.0);
    w1 = 1.0 / (2.0 - k);
    w0 = 1.0 - 2.0 * w1;
    
    // c_i (drift) and d_i (kick)
    d1 = w1;
    d2 = w0;
    c1 = w1 / 2.0;
    c2 = (w1 + w0) / 2.0;
}

void AASIntegrator::set_central_body(double mu, double J2, double R_eq) {
    mu_ = mu;
    j2_ = J2;
    r_eq_ = R_eq;
}

void AASIntegrator::split_state(const Eigen::VectorXd& y,
                                Eigen::VectorXd& q, Eigen::VectorXd& p) const {
    const int n = 3; // base physical dimension
    q = y.head(n);
    p = y.segment(n, n);
}

Eigen::VectorXd AASIntegrator::join_state(const Eigen::VectorXd& q,
                                          const Eigen::VectorXd& p) const {
    Eigen::VectorXd y(6);
    y.head(3) = q;
    y.tail(3) = p;
    return y;
}

double AASIntegrator::compute_force_gradient(const Eigen::VectorXd& q) const {
    const double r = q.norm();
    if (r < 1e-3) return 1.0e10; 
    
    // Scale derived from the spectral radius of the Hessian
    double grad = std::abs(2.0 * mu_ / (r * r * r));
    
    if (j2_ != 0.0) {
        const double r2 = r * r;
        const double j2_term = std::abs(4.0 * mu_ * j2_ * r_eq_ * r_eq_ / (r2 * r2 * r));
        grad += j2_term;
    }
    
    return grad;
}

Eigen::Matrix3d AASIntegrator::compute_hessian(const Eigen::VectorXd& q) const {
    const double r = q.norm();
    const double r2 = r * r;
    const double r3 = r2 * r;
    const double r5 = r3 * r2;
    
    // 1. Point mass Hessian: H_pm = -mu/r^3 I + 3mu/r^5 (q ⊗ q)
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d H = (-mu_ / r3) * I + (3.0 * mu_ / r5) * (q * q.transpose());
    
    // 2. J2 Hessian Contribution
    if (j2_ != 0.0) {
        const double req2 = r_eq_ * r_eq_;
        const double z = q[2];
        const double z2 = z * z;
        const double r7 = r5 * r2;
        const double r9 = r7 * r2;
        
        // Potential V_J2 = C (r^-3 - 3z^2 r^-5)
        const double C = 0.5 * mu_ * j2_ * req2;
        
        Eigen::Matrix3d H_j2;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double delta_ij = (i == j) ? 1.0 : 0.0;
                double delta_iz = (i == 2) ? 1.0 : 0.0;
                double delta_jz = (j == 2) ? 1.0 : 0.0;
                
                double term1 = (15.0 * z2 / r7 - 3.0 / r5) * delta_ij;
                double term2 = q[i] * (30.0 * z * delta_jz / r7 - 105.0 * z2 * q[j] / r9 + 15.0 * q[j] / r7);
                double term3 = -6.0 * delta_iz * (delta_jz / r5 - 5.0 * z * q[j] / r7);
                
                H_j2(i, j) = C * (term1 + term2 + term3);
            }
        }
        H += H_j2;
    }
    
    return H;
}

/**
 * @brief Hamiltonian (Total Energy) calculation for drift monitoring
 */
double AASIntegrator::compute_total_energy(const Eigen::VectorXd& q, const Eigen::VectorXd& p) const {
    double r = q.norm();
    double kinetic = 0.5 * p.squaredNorm();
    double potential = -mu_ / r;
    if (j2_ != 0.0) {
        double r2 = r * r;
        double z = q[2];
        potential += (mu_ * j2_ * r_eq_ * r_eq_ / (2.0 * r2 * r)) * (3.0 * (z * z / r2) - 1.0);
    }
    return kinetic + potential;
}

/**
 * @brief Compute Shadow Hamiltonian (Modified Hamiltonian) for Yoshida-4
 * The dominance of O(h^4) terms provides a non-secular energy error metric.
 */
double AASIntegrator::compute_modified_hamiltonian(const Eigen::VectorXd& q,
                                                   const Eigen::VectorXd& p,
                                                   double h) const {
    double H0 = compute_total_energy(q, p);
    
    // Shadow Hamiltonian Lead Terms (C1 and C2 for Yoshida-4)
    // Formula derived from Poisson Algebra expansion of H_shadow
    // H_shadow = H + h^4 * [ sigma_V * {V,{V,{V,{V,T}}}} + sigma_T * {T,{T,{T,{T,V}}}} ]
    // where {T,{T,{T,{T,V}}}} relates to the Hessian and velocities.
    
    // exact coefficient for Yoshida-4 (Shadow Hamiltonian lead term)
    // Yoshida 1990, eq. (3.7)
    const double sigma_4 = (2.0 * std::pow(w1, 4) + std::pow(w0, 4)) / 24.0;
    
    const double r = q.norm();
    const double r3 = r * r * r;
    const double V_pp = mu_ / r3; // V'' for point mass
    
    // Poisson bracket terms (Order h^4)
    const double gradV = mu_ / (r * r);
    const double term_V = (gradV * gradV) * V_pp;
    
    const Eigen::Matrix3d H_hess = compute_hessian(q);
    const double pHp = p.dot(H_hess * p);
    const double term_T = pHp * V_pp; 
    
    const double correction = std::pow(h, 4) * sigma_4 * (term_V + term_T);
    
    return H0 + correction;
}

void AASIntegrator::update_phi_kick(const DerivativeFunction&, double, double step, const Eigen::VectorXd& q, const Eigen::VectorXd&, Eigen::MatrixXd& phi) {
    if (phi.size() == 0) return;
    
    // Gravitational force model: f only depends on position q.
    // JV (Jacobian w.r.t velocity) is zero.
    // JF (Jacobian w.r.t position) is the Hessian of the potential.
    Eigen::Matrix3d JF = compute_hessian(q);
    
    // Symplectic update of the Momentum-Position and Momentum-Momentum STM blocks
    phi.block<3, 6>(3, 0) += step * JF * phi.block<3, 6>(0, 0);
}

void AASIntegrator::symplectic_step(const DerivativeFunction& f, double t, Eigen::VectorXd& q, Eigen::VectorXd& p, Eigen::MatrixXd& phi, double ds) {
    auto drift = [&](double h) { q += p * h; if (phi.size() > 0) phi.block<3, 6>(0, 0) += h * phi.block<3, 6>(3, 0); };
    auto kick = [&](double h) { Eigen::VectorXd force = f(t, join_state(q, p)).tail(3); p += force * h; update_phi_kick(f, t, h, q, p, phi); stats_.num_function_evals++; };
    drift(c1*ds); kick(d1*ds); drift(c2*ds); kick(d2*ds); drift(c2*ds); kick(d1*ds); drift(c1*ds);
}

double AASIntegrator::estimate_step_size(const Eigen::VectorXd& q, const Eigen::VectorXd& p, double target_dt) const {
    double r = q.norm();
    if (!std::isfinite(r)) return 1e-8; // Default safe step on failure
    
    // Scale precision to be more practical for O(4) symplectic integration
    // 1e-6 typically would result in extremely small steps (0.01 days), 
    // we want 1e-6 to map to something like 0.1-0.5 days for Ceres.
    double scaled_precision = precision_ * 50.0;
    double g_n = scaled_precision / std::sqrt(std::max(compute_force_gradient(q), 1e-20));
    double g_avg = g_n;
    
    // Symmetrization loop (1 iteration is sufficient for most cases)
    // Point on trajectory shifted by Sundman step
    Eigen::VectorXd q_next = q + p * g_avg;
    double g_p = scaled_precision / std::sqrt(std::max(compute_force_gradient(q_next), 1e-20));
    
    // Use harmonic mean for time-reversibility preservation
    g_avg = 2.0 * g_n * g_p / (g_n + g_p);
    
    // Apply physical constraints
    // 1. Limit step to 1/20th of local orbital period
    double T_local = 2.0 * constants::PI * std::sqrt(r * r * r / mu_);
    g_avg = std::min(g_avg, T_local / 20.0);
    
    // 2. Global bounds
    g_avg = std::clamp(g_avg, 1e-8, 100.0);
    
    // Perihelion protection (r < 2 Solar Radii)
    if (r < 0.0093) {
        g_avg = std::min(g_avg, 1.0 / 1440.0); // 1 minute limit (hardcoded for stability)
    }
    
    return std::min(g_avg, target_dt);
}

Eigen::VectorXd AASIntegrator::integrate(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, double tf) {
    stats_.reset();
    double t = t0, dir = (tf > t0) ? 1.0 : -1.0;
    Eigen::VectorXd q, p; split_state(y0, q, p);
    Eigen::MatrixXd phi; if (y0.size() == 42) { phi.resize(6, 6); for (int i=0; i<36; ++i) phi(i/6, i%6) = y0[6+i]; }
    double H_i = compute_total_energy(q, p);
    while (std::abs(tf - t) > 1e-14 && stats_.num_steps < 1000000) {
        double dt = estimate_step_size(q, p, std::abs(tf - t));
        symplectic_step(f, t, q, p, phi, dt * dir);
        t += dt * dir; stats_.num_steps++;
    }
    stats_.hamiltonian_drift = std::abs((compute_total_energy(q, p) - H_i) / H_i);
    stats_.final_time = t;
    return (phi.size() == 36) ? finalize_state_phi(y0, q, p, phi) : join_state(q, p);
}

Eigen::VectorXd AASIntegrator::finalize_state_phi(const Eigen::VectorXd& y0, const Eigen::VectorXd& q, const Eigen::VectorXd& p, const Eigen::MatrixXd& phi) {
    Eigen::VectorXd r(42);
    r << q, p, Eigen::Map<const Eigen::VectorXd>(phi.data(), 36);
    return r;
}

void AASIntegrator::integrate_steps(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, double tf, std::vector<double>& t_out, std::vector<Eigen::VectorXd>& y_out) {
    stats_.reset(); t_out.clear(); y_out.clear();
    double t = t0, dir = (tf > t0) ? 1.0 : -1.0;
    Eigen::VectorXd q, p; split_state(y0, q, p);
    Eigen::MatrixXd phi; if (y0.size() == 42) { phi.resize(6, 6); for (int i=0; i<36; ++i) phi(i/6, i%6) = y0[6+i]; }
    t_out.push_back(t); y_out.push_back(y0);
    while (std::abs(tf - t) > 1e-14 && stats_.num_steps < 1000000) {
        double dt = estimate_step_size(q, p, std::abs(tf - t));
        symplectic_step(f, t, q, p, phi, dt * dir); t += dt * dir; stats_.num_steps++;
        t_out.push_back(t); 
        if (phi.size() == 36) { Eigen::VectorXd r(42); r << q, p, Eigen::Map<Eigen::VectorXd>(phi.data(), 36); y_out.push_back(r); }
        else y_out.push_back(join_state(q, p));
    }
}

std::vector<Eigen::VectorXd> AASIntegrator::integrate_at(const DerivativeFunction& f, const Eigen::VectorXd& y0, double t0, const std::vector<double>& t_targets) {
    stats_.reset();
    std::vector<Eigen::VectorXd> res; res.reserve(t_targets.size());
    double t = t0; Eigen::VectorXd q, p; split_state(y0, q, p);
    Eigen::MatrixXd phi; if (y0.size() == 42) { phi.resize(6, 6); for (int i=0; i<36; ++i) phi(i/6, i%6) = y0[6+i]; }
    for (double target : t_targets) {
        double dir = (target > t) ? 1.0 : -1.0;
        while (std::abs(target - t) > 1e-14 && stats_.num_steps < 1000000) {
            double dt = estimate_step_size(q, p, std::abs(target - t));
            symplectic_step(f, t, q, p, phi, dt * dir); t += dt * dir; stats_.num_steps++;
        }
        res.push_back((phi.size() == 36) ? finalize_state_phi(y0, q, p, phi) : join_state(q, p));
    }
    return res;
}

} // namespace astdyn::propagation
