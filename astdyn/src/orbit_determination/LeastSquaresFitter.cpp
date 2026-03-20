/**
 * @file LeastSquaresFitter.cpp
 * @brief Implementation of least squares orbit fitter
 */

#include "astdyn/orbit_determination/LeastSquaresFitter.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <cmath>
#include <algorithm>

namespace astdyn::orbit_determination {

using namespace astdyn::constants;

LeastSquaresFitter::LeastSquaresFitter() {}

Eigen::MatrixXd LeastSquaresFitter::build_design_matrix(
    const std::vector<ObservationResidual>& residuals,
    const physics::CartesianStateTyped<core::GCRF>& state,
    time::EpochTDB epoch,
    STMFunction stm_func
) {
    int n_obs = residuals.size();
    Eigen::MatrixXd A(2 * n_obs, 6);
    
    for (int i = 0; i < n_obs; ++i) {
        const auto& res = residuals[i];
        time::EpochTDB t_obs = astdyn::time::to_tdb(res.time);
        
        // Get STM from epoch to observation time
        auto [state_at_obs, stm] = stm_func(state, epoch, t_obs);
        
        // MUST use topocentric vector for RA/Dec partials!
        // We can reconstruct it from range and computed angles in the residual
        double cra = std::cos(res.computed_ra.to_rad());
        double sra = std::sin(res.computed_ra.to_rad());
        double cdc = std::cos(res.computed_dec.to_rad());
        double sdc = std::sin(res.computed_dec.to_rad());
        double rho = res.range.to_m();
        
        double x_topo = rho * cdc * cra;
        double y_topo = rho * cdc * sra;
        double z_topo = rho * sdc;
        
        double rho_xy2 = x_topo * x_topo + y_topo * y_topo;
        double rho2 = rho * rho;
        
        // ∂RA/∂r (rad/m) projected: ∂(RA*cos(dec))/∂r
        Eigen::Vector3d dRA_drho; 
        if (rho_xy2 > 1e-6) { 
            // The residual is projected: res_ra = delta_ra * cos(dec)
            // So the design matrix must also be projected to match.
            dRA_drho(0) = (-y_topo / rho_xy2) * cdc;
            dRA_drho(1) = ( x_topo / rho_xy2) * cdc;
            dRA_drho(2) =  0.0;
        } else {
            dRA_drho.setZero();
        }
        
        // ∂Dec/∂r (rad/m)
        Eigen::Vector3d dDec_drho;
        double rho_xy = std::sqrt(rho_xy2);
        if (rho2 > 1e-6 && rho_xy > 1e-6) {
            dDec_drho(0) = -x_topo * z_topo / (rho2 * rho_xy);
            dDec_drho(1) = -y_topo * z_topo / (rho2 * rho_xy);
            dDec_drho(2) =  rho_xy / rho2;
        } else {
            dDec_drho.setZero();
        }
        
        // Design matrix row H_i = [dRA/drho, 0_v]
        Eigen::Matrix<double, 2, 6> dRho_dx_state;
        dRho_dx_state.row(0) << dRA_drho.transpose(), 0, 0, 0;
        dRho_dx_state.row(1) << dDec_drho.transpose(), 0, 0, 0;

        // STM scaling: The STM (phi) is integrated in internal units [AU, AU/day].
        // To build the design matrix in SI [meters, meters/second], we scale the 
        // time units in the cross-blocks (day -> seconds):
        astdyn::Matrix6d phi_si = stm;
        phi_si.block<3, 3>(0, 3) *= constants::SECONDS_PER_DAY;
        phi_si.block<3, 3>(3, 0) /= constants::SECONDS_PER_DAY;
        
        // Map to initial epoch: A_i = H_i * STM_si
        Eigen::Matrix<double, 2, 6> A_i = dRho_dx_state * phi_si;
        
        A.row(2*i) = A_i.row(0);
        A.row(2*i+1) = A_i.row(1);
    }
    
    return A;
}

Eigen::Vector<double, 6> LeastSquaresFitter::solve_normal_equations(
    const Eigen::MatrixXd& A,
    const Eigen::VectorXd& residuals,
    const Eigen::VectorXd& weights,
    astdyn::Matrix6d& covariance
) {
    Eigen::MatrixXd W = weights.asDiagonal();
    astdyn::Matrix6d N = A.transpose() * W * A;
    Eigen::Vector<double, 6> b = A.transpose() * W * residuals;
    
    Eigen::Vector<double, 6> dx = N.ldlt().solve(b);
    covariance = N.inverse();
    return dx;
}

int LeastSquaresFitter::reject_outliers(std::vector<ObservationResidual>& residuals) {
    if (!outlier_rejection_) return 0;
    
    // Compute RMS
    double sum_sq = 0.0;
    int count = 0;
    
    for (const auto& res : residuals) {
        if (!res.outlier) {
            sum_sq += res.residual_ra.to_arcsec() * res.residual_ra.to_arcsec();
            sum_sq += res.residual_dec.to_arcsec() * res.residual_dec.to_arcsec();
            count += 2;
        }
    }
    
    if (count == 0) return 0;
    
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
        if (!res.outlier) {
            double res_norm = std::sqrt(
                res.residual_ra.to_arcsec() * res.residual_ra.to_arcsec() +
                res.residual_dec.to_arcsec() * res.residual_dec.to_arcsec()
            );
            
            if (res_norm > threshold) {
                res.outlier = true;
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
    
    result.chi_squared = 0.0;
    for (const auto& res : residuals) {
        if (!res.outlier) {
            sum_ra_sq += res.residual_ra.to_arcsec() * res.residual_ra.to_arcsec();
            sum_dec_sq += res.residual_dec.to_arcsec() * res.residual_dec.to_arcsec();
            // Chi-squared is sum of normalized residuals squared
            result.chi_squared += res.chi_squared;
            count++;
        }
    }
    
    result.num_observations = residuals.size();
    result.num_rejected = std::count_if(residuals.begin(), residuals.end(),
                                        [](const auto& r) { return r.outlier; });
    
    if (count > 0) {
        result.rms_ra = astrometry::Angle::from_arcsec(std::sqrt(sum_ra_sq / count));
        result.rms_dec = astrometry::Angle::from_arcsec(std::sqrt(sum_dec_sq / count));
        result.rms_total = astrometry::Angle::from_arcsec(std::sqrt((sum_ra_sq + sum_dec_sq) / (2 * count)));
        
        // Reduced chi-squared
        if (2 * count > 6) {
           result.chi_squared /= (2 * count - 6);
        }
    } else {
        result.rms_ra = astrometry::Angle::from_arcsec(0.0);
        result.rms_dec = astrometry::Angle::from_arcsec(0.0);
        result.rms_total = astrometry::Angle::from_arcsec(0.0);
        result.chi_squared = 0.0;
    }
}

FitResult LeastSquaresFitter::fit(
    const physics::CartesianStateTyped<core::GCRF>& initial_state,
    time::EpochTDB epoch,
    ResidualFunction residual_func,
    STMFunction stm_func
) {
    FitResult result;
    result.state = initial_state;
    result.converged = false;
    result.num_iterations = 0;
    
    auto calc_rms = [](const std::vector<ObservationResidual>& res) -> double {
        double sum = 0.0;
        int count = 0;
        for (const auto& r : res) {
            if (!r.outlier) {
                sum += std::pow(r.residual_ra.to_arcsec(), 2) + 
                       std::pow(r.residual_dec.to_arcsec(), 2);
                count += 2;
            }
        }
        return (count > 0) ? std::sqrt(sum / count) : 0.0;
    };

    auto current_residuals = residual_func(result.state, epoch);
    double current_rms = calc_rms(current_residuals);
    
    for (int iter = 0; iter < max_iterations_; ++iter) {
        result.num_iterations = iter + 1;
        
        if (iter > 0) {
            reject_outliers(current_residuals);
            current_rms = calc_rms(current_residuals); 
        }
        
        auto A = build_design_matrix(current_residuals, result.state, epoch, stm_func);
        
        int n_obs = current_residuals.size();
        Eigen::VectorXd res_vec(2 * n_obs);
        Eigen::VectorXd weights(2 * n_obs);
        
        for (int i = 0; i < n_obs; ++i) {
            if (!current_residuals[i].outlier) {
                // IMPORTANT: A is in [rad/m], so res_vec MUST be in [rad]
                res_vec(2*i) = current_residuals[i].residual_ra.to_rad();
                res_vec(2*i+1) = current_residuals[i].residual_dec.to_rad();
                weights(2*i) = current_residuals[i].weight_ra; // weight is in 1/rad^2 ? check Residuals.hpp
                weights(2*i+1) = current_residuals[i].weight_dec;
            } else {
                res_vec(2*i) = 0; res_vec(2*i+1) = 0;
                weights(2*i) = 0; weights(2*i+1) = 0;
            }
        }
        
        // Solve for dx in SI (m, m/s)
        Eigen::Vector<double, 6> dx_full = solve_normal_equations(A, res_vec, weights, result.covariance);
        
        // --- LINE SEARCH ---
        double scale = 1.0;
        bool step_accepted = false;
        physics::CartesianStateTyped<core::GCRF> candidate_state;
        std::vector<ObservationResidual> candidate_residuals;
        double candidate_rms = 0.0;

        for (int k = 0; k < 8; ++k) { 
            Eigen::Vector<double, 6> candidate_vec = result.state.to_eigen_si();
            candidate_vec.head<6>() += dx_full * scale;

            candidate_state = physics::CartesianStateTyped<core::GCRF>::from_si(
                result.state.epoch,
                candidate_vec[0], candidate_vec[1], candidate_vec[2],
                candidate_vec[3], candidate_vec[4], candidate_vec[5],
                result.state.gm.to_m3_s2()
            );
            
            candidate_residuals = residual_func(candidate_state, epoch);
            for(size_t i=0; i<candidate_residuals.size(); ++i) {
                candidate_residuals[i].outlier = current_residuals[i].outlier;
            }
            
            candidate_rms = calc_rms(candidate_residuals);
            if (candidate_rms < current_rms) {
                step_accepted = true;
                break;
            }
            scale *= 0.5;
        }
        
        if (!step_accepted) {
            if (dx_full.norm() < tolerance_) {
                result.converged = true;
                break;
            }
            break;
        }
        
        result.state = candidate_state;
        current_residuals = candidate_residuals;
        current_rms = candidate_rms;
        
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
