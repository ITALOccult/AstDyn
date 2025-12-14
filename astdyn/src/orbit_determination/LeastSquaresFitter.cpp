/**
 * @file LeastSquaresFitter.cpp
 * @brief Implementation of least squares orbit fitter
 */

#include "astdyn/orbit_determination/LeastSquaresFitter.hpp"
#include <cmath>
#include <algorithm>

namespace astdyn::orbit_determination {

LeastSquaresFitter::LeastSquaresFitter() {}

Eigen::MatrixXd LeastSquaresFitter::build_design_matrix(
    const std::vector<ObservationResidual>& residuals,
    const Eigen::Vector<double, 6>& state,
    double epoch_mjd,
    STMFunction stm_func
) {
    // Design matrix A: each observation contributes 2 rows (RA, Dec)
    int n_obs = residuals.size();
    Eigen::MatrixXd A(2 * n_obs, 6);
    
    // For each observation, compute: A_i = (∂ρ/∂x) × Φ(t_i, t_0)
    // where ∂ρ/∂x is partial of RA/Dec w.r.t. state at obs time
    // and Φ is STM from epoch to obs time
    
    // State is already in Equatorial Frame (Standard for Propagator)
    // We compute partial derivatives of RA/Dec w.r.t Equatorial State.
    
    constexpr double rad_to_arcsec = 206265.04606;

    for (int i = 0; i < n_obs; ++i) {
        double t_obs = residuals[i].epoch_mjd;
        
        // Get STM from epoch to observation time
        auto [state_at_obs, stm] = stm_func(state, epoch_mjd, t_obs);
        
        // Calculate Topocentric position if possible? 
        // For partial derivatives, using Heliocentric/Barycentric State at obs time 
        // is the standard approximation for the design matrix (H matrix).
        // The residuals themselves (O-C) MUST use topocentric, but H can be geocentric/heliocentric approximate.
        
        // state_at_obs is [x, y, z, vx, vy, vz] in Equatorial Frame
        Eigen::Vector3d r = state_at_obs.head<3>();
        
        double x = r(0), y = r(1), z = r(2);
        double r_norm = r.norm();
        double rho_xy = std::sqrt(x*x + y*y);
        
        // ∂RA/∂r (Equatorial)
        Eigen::Vector3d dRA_dr; 
        if (rho_xy > 1e-10) {
            dRA_dr(0) = -y / (rho_xy * rho_xy);
            dRA_dr(1) =  x / (rho_xy * rho_xy);
            dRA_dr(2) =  0.0;
        } else {
            dRA_dr.setZero();
        }
        
        // ∂Dec/∂r (Equatorial)
        Eigen::Vector3d dDec_dr;
        if (r_norm > 1e-10) {
            dDec_dr(0) = -x * z / (r_norm * r_norm * rho_xy);
            dDec_dr(1) = -y * z / (r_norm * r_norm * rho_xy);
            dDec_dr(2) =  rho_xy / (r_norm * r_norm);
        } else {
            dDec_dr.setZero();
        }
        
        // Convert radian gradients to arcsec gradients
        dRA_dr *= rad_to_arcsec;
        dDec_dr *= rad_to_arcsec;
        
        // Design matrix row H_i = [dRA/dr, 0_v]
        Eigen::Matrix<double, 2, 6> dRho_dx_state;
        dRho_dx_state.row(0) << dRA_dr.transpose(), 0, 0, 0;   // ∂RA/∂x
        dRho_dx_state.row(1) << dDec_dr.transpose(), 0, 0, 0;  // ∂Dec/∂x
        
        // Map to initial epoch: A_i = H_i * STM
        Eigen::Matrix<double, 2, 6> A_i = dRho_dx_state * stm;
        
        A.row(2*i) = A_i.row(0);
        A.row(2*i+1) = A_i.row(1);
    }
    
    return A;
}

Eigen::Vector<double, 6> LeastSquaresFitter::solve_normal_equations(
    const Eigen::MatrixXd& A,
    const Eigen::VectorXd& residuals,
    const Eigen::VectorXd& weights,
    Eigen::Matrix<double, 6, 6>& covariance
) {
    // Weight matrix W = diag(weights)
    Eigen::MatrixXd W = weights.asDiagonal();
    
    // Normal equations: (A^T W A) δx = A^T W Δρ
    Eigen::Matrix<double, 6, 6> N = A.transpose() * W * A;
    Eigen::Vector<double, 6> b = A.transpose() * W * residuals;
    
    // Solve using LU decomposition
    Eigen::Vector<double, 6> dx = N.ldlt().solve(b);
    
    // Covariance: (A^T W A)^{-1}
    covariance = N.inverse();
    
    return dx;
}

int LeastSquaresFitter::reject_outliers(std::vector<ObservationResidual>& residuals) {
    if (!outlier_rejection_) return 0;
    
    // Compute RMS
    double sum_sq = 0.0;
    int count = 0;
    
    for (const auto& res : residuals) {
        if (!res.rejected) {
            sum_sq += res.ra_residual_arcsec * res.ra_residual_arcsec;
            sum_sq += res.dec_residual_arcsec * res.dec_residual_arcsec;
            count += 2;
        }
    }
    
    double rms = std::sqrt(sum_sq / count);
    
    // Safety: Do not reject outliers if the fit is still very poor (RMS > 100 arcsec)
    // This allows the fitter to converge from a bad guess without losing data.
    if (rms > 100.0) {
        return 0; 
    }
    
    double threshold = outlier_threshold_ * rms;
    
    // Safety: Ensure threshold is not too small (e.g. < 0.1 arcsec)
    threshold = std::max(threshold, 0.1);
    
    // Reject outliers
    int num_rejected = 0;
    for (auto& res : residuals) {
        if (!res.rejected) {
            double res_norm = std::sqrt(
                res.ra_residual_arcsec * res.ra_residual_arcsec +
                res.dec_residual_arcsec * res.dec_residual_arcsec
            );
            
            if (res_norm > threshold) {
                res.rejected = true;
                num_rejected++;
            }
        }
    }
    
    return num_rejected;
}

void LeastSquaresFitter::compute_statistics(
    const std::vector<ObservationResidual>& residuals,
    FitResult& result
) {
    double sum_ra_sq = 0.0;
    double sum_dec_sq = 0.0;
    int count = 0;
    
    for (const auto& res : residuals) {
        if (!res.rejected) {
            sum_ra_sq += res.ra_residual_arcsec * res.ra_residual_arcsec;
            sum_dec_sq += res.dec_residual_arcsec * res.dec_residual_arcsec;
            count++;
        }
    }
    
    result.num_observations = residuals.size();
    result.num_rejected = std::count_if(residuals.begin(), residuals.end(),
                                        [](const auto& r) { return r.rejected; });
    
    if (count > 0) {
        result.rms_ra_arcsec = std::sqrt(sum_ra_sq / count);
        result.rms_dec_arcsec = std::sqrt(sum_dec_sq / count);
        result.rms_total_arcsec = std::sqrt((sum_ra_sq + sum_dec_sq) / (2 * count));
    } else {
        result.rms_ra_arcsec = 0.0;
        result.rms_dec_arcsec = 0.0;
        result.rms_total_arcsec = 0.0;
    }
    
    result.chi_squared = (sum_ra_sq + sum_dec_sq) / (2 * count - 6);  // DOF = 2*n - 6
}

FitResult LeastSquaresFitter::fit(
    const Eigen::Vector<double, 6>& initial_state,
    double epoch_mjd,
    ResidualFunction residual_func,
    STMFunction stm_func
) {
    FitResult result;
    result.state = initial_state;
    result.converged = false;
    result.num_iterations = 0;
    
    // Helper to calc RMS
    auto calc_rms = [](const std::vector<ObservationResidual>& res) -> double {
        double sum = 0.0;
        int count = 0;
        for (const auto& r : res) {
            if (!r.rejected) {
                sum += r.ra_residual_arcsec * r.ra_residual_arcsec + 
                       r.dec_residual_arcsec * r.dec_residual_arcsec;
                count += 2;
            }
        }
        return (count > 0) ? std::sqrt(sum / count) : 0.0;
    };

    // Initial residuals
    auto current_residuals = residual_func(result.state, epoch_mjd);
    double current_rms = calc_rms(current_residuals);
    
    for (int iter = 0; iter < max_iterations_; ++iter) {
        result.num_iterations = iter + 1;
        
        // Outlier Rejection (on accepted state)
        if (iter > 0) {
            reject_outliers(current_residuals);
            // Re-calc RMS after rejection might change "current_rms", but let's keep it consistent
            current_rms = calc_rms(current_residuals); 
        }
        
        // Build Design Matrix (A) and Normal Equations
        auto A = build_design_matrix(current_residuals, result.state, epoch_mjd, stm_func);
        
        // Pack residuals
        int n_obs = current_residuals.size();
        Eigen::VectorXd res_vec(2 * n_obs);
        Eigen::VectorXd weights(2 * n_obs);
        
        for (int i = 0; i < n_obs; ++i) {
            if (!current_residuals[i].rejected) {
                // IMPORTANT: Normal Equation is A^T * W * (y - y_model)
                // residuals vector here contains (Observed - Computed) = y - f(x)
                res_vec(2*i) = current_residuals[i].ra_residual_arcsec;
                res_vec(2*i+1) = current_residuals[i].dec_residual_arcsec;
                weights(2*i) = current_residuals[i].weight_ra;
                weights(2*i+1) = current_residuals[i].weight_dec;
            } else {
                res_vec(2*i) = 0; res_vec(2*i+1) = 0;
                weights(2*i) = 0; weights(2*i+1) = 0;
            }
        }
        
        // Solve for full step dx
        Eigen::Vector<double, 6> dx_full = solve_normal_equations(A, res_vec, weights, result.covariance);
        
        // Safety: Limit step size to avoid divergence (Trust Region)
        // Especially for position components [0,1,2] (AU)
        double pos_step_norm = dx_full.head<3>().norm();
        constexpr double MAX_STEP_AU = 0.1; // Max 0.1 AU correction per step
        
        if (pos_step_norm > MAX_STEP_AU) {
            double scale_factor = MAX_STEP_AU / pos_step_norm;
            dx_full *= scale_factor;
            // std::cout << "Step limited: " << pos_step_norm << " AU -> " << MAX_STEP_AU << " AU\n";
        }
        
        // --- LINE SEARCH (Step Halving) ---
        double scale = 1.0;
        bool step_accepted = false;
        Eigen::Vector<double, 6> candidate_state;
        std::vector<ObservationResidual> candidate_residuals;
        double candidate_rms = 0.0;

        // Try reducing step size: 1.0, 0.5, 0.25 ... 
        for (int k = 0; k < 8; ++k) { 
            candidate_state = result.state + dx_full * scale;
            
            // Compute residuals at candidate
            candidate_residuals = residual_func(candidate_state, epoch_mjd);
            
            // Apply SAME rejection mask as baseline to compare apples to apples?
            // Or re-eval rejection? Usually re-eval happens next iter.
            // Let's assume rejection mask is fixed for the step calculation.
            for(size_t i=0; i<candidate_residuals.size(); ++i) {
                candidate_residuals[i].rejected = current_residuals[i].rejected;
            }
            
            candidate_rms = calc_rms(candidate_residuals);
            
            if (candidate_rms < current_rms) {
                // Improvement found!
                step_accepted = true;
                break;
            }
            
            // If RMS increased (or is NaN), reduce step
            scale *= 0.5;
        }
        
        // If even small step fails, check if we are already precise enough
        if (!step_accepted) {
            if (dx_full.norm() < tolerance_) {
                // Actually converged locally
                result.converged = true;
                break;
            } else {
                // Stuck in local minimum or divergence
                // Accept the smallest step anyway? Or stop? 
                // Let's accept smallest step to see if it escapes next time, or break.
                // Usually break.
                // std::cerr << "Line search failed. Stopping.\n";
                break;
            }
        }
        
        // Apply Step
        result.state = candidate_state;
        current_residuals = candidate_residuals;
        current_rms = candidate_rms;
        
        // Check Convergence (on the accepted scaled step)
        if ((dx_full * scale).norm() < tolerance_) {
            result.converged = true;
            break;
        }
    }
    
    // Store final results
    result.residuals = current_residuals;
    compute_statistics(current_residuals, result);
    
    return result;
}

} // namespace astdyn::orbit_determination
