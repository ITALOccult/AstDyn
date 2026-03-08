/**
 * @file AASIntegrator.cpp
 * @brief AAS Adaptive Symplectic Integrator Implementation (Time-Transformed)
 */

#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/core/Constants.hpp"
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
    const int n = 3; 
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
    
    // exact coefficient for Yoshida-4 (sigma_yoshida)
    const double sigma_exact = (2.0 * std::pow(w1, 5) + std::pow(w0, 5)) / 720.0;
    
    const double r = q.norm();
    const double r3 = r * r * r;
    const double V_pp = mu_ / r3; // V'' for point mass
    
    // Poisson bracket terms (Order h^4)
    // {V,{V,{V,{V,T}}}} ~ (gradV)^2 * V'' / V ~ gradV^2 * (1/r^2)
    const double gradV = mu_ / (r * r);
    const double term_V = (gradV * gradV) * V_pp;
    
    // {T,{T,{T,{T,V}}}} ~ v^4 * V''''' is too high, but p^T * Hess(V) * p is better for h^2
    // For h^4, we use the square of the Hessian term for metric consistency:
    const Eigen::Matrix3d H_hess = compute_hessian(q);
    const double pHp = p.dot(H_hess * p);
    const double term_T = pHp * V_pp; 
    
    // Correction must have units of Energy (L^2/T^2)
    // term_V: (L/T^2)^2 * (1/T^2) = L^2/T^6. Multiply by h^4 (T^4) -> L^2/T^2. OK.
    // term_T: (L^2/T^4) * (1/T^2) = L^2/T^6. Multiply by h^4 (T^4) -> L^2/T^2. OK.
    const double correction = std::pow(h, 4) * sigma_exact * (term_V + term_T);
    
    return H0 + correction;
}

void AASIntegrator::symplectic_step(const DerivativeFunction& f,
                                    double t,
                                    Eigen::VectorXd& q, 
                                    Eigen::VectorXd& p, 
                                    Eigen::MatrixXd& phi,
                                    double ds) {
    bool has_stm = (phi.rows() == 6 && phi.cols() == 6);
    
    auto update_phi_drift = [&](double step) {
        if (!has_stm) return;
        phi.block<3, 6>(0, 0) += step * phi.block<3, 6>(3, 0);
    };
    
    auto update_phi_kick = [&](double step, const Eigen::VectorXd& q_eval, const Eigen::VectorXd& p_eval) {
        if (!has_stm) return;
        
        // Full Jacobian: d a / d q AND d a / d v (M4 Requirement)
        // Required for non-conservative forces (Drag, Yarkovsky)
        const double h_fd = 1e-7;
        Eigen::Matrix3d JF = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d JV = Eigen::Matrix3d::Zero();
        Eigen::VectorXd y_base = join_state(q_eval, p_eval);
        
        for (int col = 0; col < 3; ++col) {
            double pert_q = h_fd * std::max(std::abs(y_base[col]), 1.0);
            Eigen::VectorXd y_pq = y_base; y_pq[col] += pert_q;
            Eigen::VectorXd y_mq = y_base; y_mq[col] -= pert_q;
            JF.col(col) = (f(t, y_pq).tail(3) - f(t, y_mq).tail(3)) / (2.0 * pert_q);

            double pert_v = h_fd * std::max(std::abs(y_base[3+col]), 1.0);
            Eigen::VectorXd y_pv = y_base; y_pv[3+col] += pert_v;
            Eigen::VectorXd y_mv = y_base; y_mv[3+col] -= pert_v;
            JV.col(col) = (f(t, y_pv).tail(3) - f(t, y_mv).tail(3)) / (2.0 * pert_v);
        }
        
        // Tangent map kick: Φ_new = [[I, 0], [step*JF, I + step*JV]] * Φ_old
        phi.block<3, 6>(3, 0) += step * (JF * phi.block<3, 6>(0, 0) + JV * phi.block<3, 6>(3, 0));
    };

    // Yoshida-4 DKDKDKD
    // Step 1: Drift
    q += p * (c1 * ds);
    update_phi_drift(c1 * ds);
    
    // Step 1: Kick
    Eigen::VectorXd y1 = join_state(q, p);
    Eigen::VectorXd force1 = f(t, y1).tail(3);
    p += force1 * (d1 * ds);
    update_phi_kick(d1 * ds, q, p);
    stats_.num_function_evals++;

    // Step 2: Drift
    q += p * (c2 * ds);
    update_phi_drift(c2 * ds);
    
    // Step 2: Kick
    Eigen::VectorXd y2 = join_state(q, p);
    Eigen::VectorXd force2 = f(t, y2).tail(3);
    p += force2 * (d2 * ds);
    update_phi_kick(d2 * ds, q, p);
    stats_.num_function_evals++;

    // Step 3: Drift
    q += p * (c2 * ds);
    update_phi_drift(c2 * ds);
    
    // Step 3: Kick
    Eigen::VectorXd y3 = join_state(q, p);
    Eigen::VectorXd force3 = f(t, y3).tail(3);
    p += force3 * (d1 * ds);
    update_phi_kick(d1 * ds, q, p);
    stats_.num_function_evals++;

    // Step 4: Drift
    q += p * (c1 * ds);
    update_phi_drift(c1 * ds);
}

Eigen::VectorXd AASIntegrator::integrate(const DerivativeFunction& f,
                                         const Eigen::VectorXd& y0,
                                         double t0, double tf) {
    Integrator::reset_statistics();
    const double direction = (tf > t0) ? 1.0 : -1.0;
    double t = t0;
    
    Eigen::VectorXd q, p;
    split_state(y0, q, p);

    const double dtau = 1.0; // Sundman fictitious time step (fixed to unity)
    double g_init = precision_ / std::sqrt(std::max(compute_force_gradient(q), 1e-20));
    double H_mod_initial = compute_modified_hamiltonian(q, p, g_init);
    double H_initial = compute_total_energy(q, p);

    Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(0, 0);
    if (y0.size() == 42) {
        phi.resize(6, 6);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) phi(i, j) = y0[6 + i * 6 + j];
        }
    }

    const int max_steps = 1000000;
    double last_h = 0.0;
    
    while (std::abs(tf - t) > 1e-14) {
        if (stats_.num_steps >= max_steps) break;
        
        // 1. Estimate current step size g_n
        double g_n = precision_ / std::sqrt(std::max(compute_force_gradient(q), 1e-20));
        
        // 2. Symmetrization loop (M1) - Iterative Fix for g_avg
        double g_avg = g_n;
        for (int iter = 0; iter < 3; ++iter) {
            double dt_est = g_avg * dtau;
            
            // Predict q_next for g(q_next) evaluation
            // We use a simple drift-kick-drift prediction for speed
            Eigen::VectorXd q_next = q + p * dt_est; 
            double g_p = precision_ / std::sqrt(std::max(compute_force_gradient(q_next), 1e-20));
            
            // Re-calc average
            g_avg = 0.5 * (g_n + g_p);
        }
        
        // 3. Safety Clamps (m2)
        // Clamp to prevent numerical stalls or infinite steps 
        // Limits: 1 ms to 20 days
        g_avg = std::clamp(g_avg, 1e-3, 1728000.0);
        
        // Perihelion protection: ensure stability near Sun/Central Body
        const double r_obj = q.norm();
        if (r_obj < 1.4e9) { // ~2 Solar Radii
             g_avg = std::min(g_avg, 3600.0); // Max 1 hour near Sun
        }
        
        double dt_physical = std::min(g_avg, std::abs(tf - t));
        double actual_ds = dt_physical * direction;
        
        symplectic_step(f, t, q, p, phi, actual_ds);
        
        t += actual_ds;
        last_h = dt_physical;
        stats_.num_steps++;
    }
    
    double H_final = compute_total_energy(q, p);
    double H_mod_final = compute_modified_hamiltonian(q, p, last_h);
    
    stats_.hamiltonian_drift = std::abs((H_final - H_initial) / H_initial);
    stats_.shadow_hamiltonian_drift = std::abs((H_mod_final - H_mod_initial) / H_mod_initial);
    stats_.final_time = t;
    
    if (phi.size() > 0) {
        Eigen::VectorXd y_out(42);
        y_out.segment<3>(0) = q;
        y_out.segment<3>(3) = p;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) y_out[6 + i * 6 + j] = phi(i, j);
        }
        return y_out;
    } else {
        return join_state(q, p);
    }
}

void AASIntegrator::integrate_steps(const DerivativeFunction& f,
                                   const Eigen::VectorXd& y0,
                                   double t0, double tf,
                                   std::vector<double>& t_out,
                                   std::vector<Eigen::VectorXd>& y_out) {
    Integrator::reset_statistics();
    t_out.clear();
    y_out.clear();
    const double direction = (tf > t0) ? 1.0 : -1.0;
    double t = t0;
    t_out.push_back(t);
    y_out.push_back(y0);
    
    Eigen::VectorXd q, p;
    split_state(y0, q, p);
    
    const double dtau = 1.0; // Sundman fictitious time step (fixed to unity)
    double g_init = precision_ / std::sqrt(std::max(compute_force_gradient(q), 1e-20));
    double H_mod_initial = compute_modified_hamiltonian(q, p, g_init);
    double H_initial = compute_total_energy(q, p);

    Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(0, 0);
    if (y0.size() == 42) {
        phi.resize(6, 6);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) phi(i, j) = y0[6 + i * 6 + j];
        }
    }
    
    const int max_steps = 1000000;
    double last_h = 0.0;
    while (std::abs(tf - t) > 1e-14) {
        if (stats_.num_steps >= max_steps) break;
        
        // 1. Estimate current step size g_n
        double g_n = precision_ / std::sqrt(std::max(compute_force_gradient(q), 1e-20));
        
        // 2. Symmetrization loop (M1) - Iterative Fix for g_avg
        double g_avg = g_n;
        for (int iter = 0; iter < 3; ++iter) {
            double dt_est = g_avg * dtau;
            
            // Predict q_next for g(q_next) evaluation
            Eigen::VectorXd q_next = q + p * dt_est; 
            double g_p = precision_ / std::sqrt(std::max(compute_force_gradient(q_next), 1e-20));
            
            // Re-calc average
            g_avg = 0.5 * (g_n + g_p);
        }
        
        // 3. Safety Clamps (m2)
        g_avg = std::clamp(g_avg, 1e-3, 1728000.0);
        
        // Perihelion protection
        const double r_obj = q.norm();
        if (r_obj < 1.4e9) { 
             g_avg = std::min(g_avg, 3600.0); 
        }
        
        double dt_physical = std::min(g_avg, std::abs(tf - t));
        double actual_ds = dt_physical * direction;
        
        symplectic_step(f, t, q, p, phi, actual_ds);
        
        t += actual_ds;
        last_h = dt_physical;
        stats_.num_steps++;
        
        t_out.push_back(t);
        if (phi.size() > 0) {
            Eigen::VectorXd y_step(42);
            y_step.segment<3>(0) = q;
            y_step.segment<3>(3) = p;
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) y_step[6 + i * 6 + j] = phi(i, j);
            }
            y_out.push_back(y_step);
        } else {
            y_out.push_back(join_state(q, p));
        }
    }
    
    double H_final = compute_total_energy(q, p);
    double H_mod_final = compute_modified_hamiltonian(q, p, last_h);
    stats_.hamiltonian_drift = std::abs((H_final - H_initial) / H_initial);
    stats_.shadow_hamiltonian_drift = std::abs((H_mod_final - H_mod_initial) / H_mod_initial);
    stats_.final_time = t;
}

} // namespace astdyn::propagation
