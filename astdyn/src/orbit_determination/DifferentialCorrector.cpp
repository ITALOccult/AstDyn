/**
 * @file DifferentialCorrector.cpp
 * @brief Implementation of differential corrections for orbit determination
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#include "astdyn/orbit_determination/DifferentialCorrector.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

namespace astdyn::orbit_determination {

using namespace astdyn::observations;
using namespace astdyn::propagation;

// ============================================================================
// DifferentialCorrector Implementation
// ============================================================================

DifferentialCorrector::DifferentialCorrector(
    std::shared_ptr<ResidualCalculator> residual_calc,
    std::shared_ptr<StateTransitionMatrix> stm_computer)
    : residual_calc_(residual_calc),
      stm_computer_(stm_computer) {
}

DifferentialCorrectorResult DifferentialCorrector::fit(
    const std::vector<OpticalObservation>& observations,
    const CartesianElements& initial_guess,
    const DifferentialCorrectorSettings& settings) {
    
    DifferentialCorrectorResult result;
    result.final_state = initial_guess;
    result.converged = false;
    result.iterations = 0;
    
    CartesianElements current_state = initial_guess;
    
    if (settings.verbose) {
        std::cout << "\n========================================\n";
        std::cout << "Differential Corrections\n";
        std::cout << "========================================\n";
        std::cout << "Observations: " << observations.size() << "\n";
        std::cout << "Max iterations: " << settings.max_iterations << "\n";
        std::cout << "Convergence: " << settings.convergence_tolerance << " AU\n\n";
    }
    
    // Main iteration loop
    for (int iter = 0; iter < settings.max_iterations; ++iter) {
        result.iterations = iter + 1;
        
        // Perform one iteration
        Eigen::VectorXd correction;
        std::vector<ObservationResidual> residuals;
        
        bool iter_success = iteration(observations, current_state, 
                                      correction, residuals);
        
        if (!iter_success) {
            if (settings.verbose) {
                std::cout << "Iteration " << iter + 1 << " FAILED (singular matrix)\n";
            }
            break;
        }
        
        // Reject outliers if requested
        if (settings.reject_outliers && iter > 0) {
            ResidualCalculator::identify_outliers(residuals, settings.outlier_sigma);
        }
        
        // Compute statistics
        auto stats = ResidualCalculator::compute_statistics(residuals, 6);
        result.rms_history.push_back(stats.rms_total);
        result.correction_norm.push_back(correction.norm());
        
        // Callback for progress monitoring
        if (iteration_callback_) {
            iteration_callback_(iter + 1, stats);
        }
        
        if (settings.verbose) {
            std::cout << "Iter " << std::setw(2) << iter + 1 << ": "
                      << "RMS = " << std::fixed << std::setprecision(3) 
                      << stats.rms_total << " arcsec, "
                      << "||Δx|| = " << std::scientific << std::setprecision(2)
                      << correction.norm() << " AU, "
                      << "Outliers = " << stats.num_outliers << "\n";
        }
        
        // Apply correction
        current_state.position += correction.segment<3>(0);
        current_state.velocity += correction.segment<3>(3);
        
        // Check convergence
        if (check_convergence(correction, settings.convergence_tolerance)) {
            result.converged = true;
            result.final_state = current_state;
            result.residuals = residuals;
            result.statistics = stats;
            
            if (settings.verbose) {
                std::cout << "\n✓ CONVERGED after " << iter + 1 << " iterations\n";
            }
            break;
        }
    }
    
    if (!result.converged && settings.verbose) {
        std::cout << "\n✗ NOT CONVERGED after " << result.iterations << " iterations\n";
    }
    
    // Compute final covariance if requested
    if (settings.compute_covariance && result.converged) {
        result.covariance = compute_covariance(
            observations, result.final_state, result.residuals);
        
        result.formal_uncertainties.resize(6);
        for (int i = 0; i < 6; ++i) {
            result.formal_uncertainties[i] = std::sqrt(result.covariance(i, i));
        }
        
        result.correlation = compute_correlation(result.covariance);
    }
    
    if (settings.verbose) {
        std::cout << "========================================\n\n";
    }
    
    return result;
}

bool DifferentialCorrector::iteration(
    const std::vector<OpticalObservation>& observations,
    const CartesianElements& current_state,
    Eigen::VectorXd& correction,
    std::vector<ObservationResidual>& residuals) {
    
    // Compute residuals
    residuals = residual_calc_->compute_residuals(observations, current_state);
    
    if (residuals.empty()) {
        return false;
    }
    
    // Build design matrix
    auto design_result = build_design_matrix(observations, current_state, residuals);
    
    if (design_result.valid_indices.empty()) {
        return false;
    }
    
    // Solve normal equations
    Matrix6d normal_matrix, normal_inv;
    auto correction_opt = solve_normal_equations(
        design_result.A, design_result.b, design_result.weights,
        normal_matrix, normal_inv);
    
    if (!correction_opt) {
        return false;
    }
    
    correction = *correction_opt;
    return true;
}

DifferentialCorrector::DesignMatrixResult 
DifferentialCorrector::build_design_matrix(
    const std::vector<OpticalObservation>& observations,
    const CartesianElements& state,
    const std::vector<ObservationResidual>& residuals) {
    
    DesignMatrixResult result;
    result.valid_indices.reserve(observations.size());
    
    // Count non-outlier observations
    int n_valid = 0;
    for (size_t i = 0; i < residuals.size(); ++i) {
        if (!residuals[i].outlier) {
            n_valid++;
            result.valid_indices.push_back(i);
        }
    }
    
    if (n_valid == 0) {
        return result;
    }
    
    // Allocate matrices
    int n_equations = 2 * n_valid;  // 2 equations per observation (RA, Dec)
    result.A.resize(n_equations, 6);
    result.b.resize(n_equations);
    result.weights.resize(n_equations);
    
    // Fill design matrix row by row
    int row = 0;
    for (size_t idx : result.valid_indices) {
        const auto& obs = observations[idx];
        const auto& res = residuals[idx];
        
        // Get observer position for this observation
        auto obs_pos_opt = residual_calc_->get_observer_position(obs);
        if (!obs_pos_opt) {
            // Skip if observer position unavailable
            continue;
        }
        Vector3d observer_pos = *obs_pos_opt;
        
        // Convert observation time to TDB (STM needs TDB)
        double obs_mjd_tdb = ResidualCalculator::utc_to_tdb(obs.mjd_utc);
        
        // Compute STM and observation partials
        auto partials = stm_computer_->compute_with_partials(
            state, obs_mjd_tdb, observer_pos);
        
        // Design matrix row: ∂(RA,Dec)/∂x₀ = ∂(RA,Dec)/∂x * Φ(t,t₀)
        Eigen::Matrix<double, 2, 6> A_obs = partials.partial_radec * partials.phi;
        
        // Fill A matrix (2 rows per observation)
        result.A.row(row) = A_obs.row(0);      // RA equation
        result.A.row(row + 1) = A_obs.row(1);  // Dec equation
        
        // Residual vector (O-C)
        result.b[row] = res.residual_ra;
        result.b[row + 1] = res.residual_dec;
        
        // Weights (inverse variance)
        result.weights[row] = 1.0 / (obs.sigma_ra * obs.sigma_ra);
        result.weights[row + 1] = 1.0 / (obs.sigma_dec * obs.sigma_dec);
        
        row += 2;
    }
    
    return result;
}

std::optional<Eigen::VectorXd> DifferentialCorrector::solve_normal_equations(
    const Eigen::MatrixXd& A,
    const Eigen::VectorXd& b,
    const Eigen::VectorXd& W,
    Matrix6d& normal_matrix,
    Matrix6d& normal_inv) {
    
    // Weight matrix (diagonal)
    Eigen::MatrixXd W_diag = W.asDiagonal();
    
    // Normal matrix: N = AᵀWA
    normal_matrix = A.transpose() * W_diag * A;
    
    // Right-hand side: AᵀWb
    Eigen::VectorXd rhs = A.transpose() * W_diag * b;
    
    // Solve: N * Δx = rhs
    Eigen::LDLT<Matrix6d> ldlt(normal_matrix);
    
    if (ldlt.info() != Eigen::Success) {
        return std::nullopt;
    }
    
    Eigen::VectorXd correction = ldlt.solve(rhs);
    
    if (ldlt.info() != Eigen::Success) {
        return std::nullopt;
    }
    
    // Compute inverse for covariance
    normal_inv = ldlt.solve(Matrix6d::Identity());
    
    return correction;
}

Matrix6d DifferentialCorrector::compute_covariance(
    const std::vector<OpticalObservation>& observations,
    const CartesianElements& final_state,
    const std::vector<ObservationResidual>& residuals) {
    
    // Rebuild design matrix for final state
    auto design_result = build_design_matrix(observations, final_state, residuals);
    
    // Compute normal matrix
    Matrix6d normal_matrix, normal_inv;
    Eigen::VectorXd dummy_correction;
    
    auto solution = solve_normal_equations(
        design_result.A, design_result.b, design_result.weights,
        normal_matrix, normal_inv);
    
    if (!solution) {
        return Matrix6d::Zero();
    }
    
    // Covariance = σ² * (AᵀWA)⁻¹
    // where σ² is the variance of unit weight
    
    auto stats = ResidualCalculator::compute_statistics(residuals, 6);
    double sigma_0_squared = stats.chi_squared / stats.degrees_of_freedom;
    
    Matrix6d covariance = sigma_0_squared * normal_inv;
    
    return covariance;
}

Matrix6d DifferentialCorrector::compute_correlation(
    const Matrix6d& covariance) const {
    
    Matrix6d correlation = Matrix6d::Zero();
    
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            double denom = std::sqrt(covariance(i, i) * covariance(j, j));
            if (denom > 0.0) {
                correlation(i, j) = covariance(i, j) / denom;
            }
        }
    }
    
    return correlation;
}

bool DifferentialCorrector::check_convergence(
    const Eigen::VectorXd& correction,
    double tolerance) const {
    
    // Check if position correction is below tolerance
    double pos_correction = correction.segment<3>(0).norm();
    return pos_correction < tolerance;
}

// ============================================================================
// Result Output
// ============================================================================

void DifferentialCorrectorResult::print_summary() const {
    std::cout << "\n========================================\n";
    std::cout << "Differential Corrections Result\n";
    std::cout << "========================================\n";
    
    std::cout << "Status: " << (converged ? "✓ CONVERGED" : "✗ NOT CONVERGED") << "\n";
    std::cout << "Iterations: " << iterations << "\n\n";
    
    std::cout << "Final RMS:\n";
    std::cout << "  RA:    " << std::fixed << std::setprecision(3) 
              << statistics.rms_ra << " arcsec\n";
    std::cout << "  Dec:   " << statistics.rms_dec << " arcsec\n";
    std::cout << "  Total: " << statistics.rms_total << " arcsec\n\n";
    
    std::cout << "Observations:\n";
    std::cout << "  Total:    " << statistics.num_observations << "\n";
    std::cout << "  Outliers: " << statistics.num_outliers << "\n";
    std::cout << "  Used:     " << (statistics.num_observations - statistics.num_outliers) << "\n\n";
    
    std::cout << "Chi-squared:\n";
    std::cout << "  χ²:     " << std::fixed << std::setprecision(2) 
              << statistics.chi_squared << "\n";
    std::cout << "  χ²/dof: " << statistics.reduced_chi_squared << "\n\n";
    
    if (formal_uncertainties.size() == 6) {
        std::cout << "Formal Uncertainties (1σ):\n";
        std::cout << "  Δx:  " << std::scientific << std::setprecision(2) 
                  << formal_uncertainties[0] << " AU\n";
        std::cout << "  Δy:  " << formal_uncertainties[1] << " AU\n";
        std::cout << "  Δz:  " << formal_uncertainties[2] << " AU\n";
        std::cout << "  Δvx: " << formal_uncertainties[3] << " AU/day\n";
        std::cout << "  Δvy: " << formal_uncertainties[4] << " AU/day\n";
        std::cout << "  Δvz: " << formal_uncertainties[5] << " AU/day\n";
    }
    
    std::cout << "========================================\n\n";
}

} // namespace astdyn::orbit_determination
