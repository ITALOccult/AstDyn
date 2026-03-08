/**
 * @file DifferentialCorrector.hpp
 * @brief Differential corrections for orbit determination
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP
#define ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <memory>
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace astdyn::orbit_determination {

/**
 * @brief Settings for differential corrections
 */
struct DifferentialCorrectorSettings {
    int max_iterations = 20;             ///< Maximum iterations
    double convergence_tolerance = 1e-6; ///< Convergence threshold [AU]
    double outlier_sigma = 3.0;          ///< Target Sigma threshold (default)
    
    // Carpentry Settings (Iterative Rejection)
    double outlier_max_sigma = 10.0;     ///< Starting loose sigma
    double outlier_min_sigma = 3.0;      ///< Final tight sigma
    
    bool reject_outliers = true;         ///< Automatically reject outliers
    bool compute_covariance = true;      ///< Compute covariance matrix
    bool verbose = false;                ///< Print iteration details
};

/**
 * @brief Result of differential corrections
 */
template <typename Frame>
struct DifferentialCorrectorResult {
    physics::CartesianStateTyped<Frame> final_state;
    bool converged;
    int iterations;
    
    // Quality metrics
    ResidualStatistics statistics;
    std::vector<ObservationResidual> residuals;
    
    // Uncertainty
    astdyn::Matrix6d covariance;                 ///< Covariance matrix [AU², AU²/day², ...]
    Eigen::VectorXd formal_uncertainties; ///< sqrt(diag(covariance))
    
    // Correlation matrix (normalized covariance)
    astdyn::Matrix6d correlation;
    
    // Normal matrix and its inverse
    astdyn::Matrix6d normal_matrix;              ///< AᵀWA
    astdyn::Matrix6d normal_matrix_inv;          ///< (AᵀWA)⁻¹
    
    // Convergence history
    std::vector<double> rms_history;     ///< RMS at each iteration
    std::vector<double> correction_norm; ///< ||Δx|| at each iteration
    
    /**
     * @brief Get 1-sigma uncertainty for parameter
     */
    double get_uncertainty(int param_index) const {
        return (param_index < 6) ? formal_uncertainties[param_index] : 0.0;
    }
    
    /**
     * @brief Print summary to stdout
     */
    void print_summary() const {
        std::cout << "\n========================================\n";
        std::cout << "Differential Corrections Result\n";
        std::cout << "========================================\n";
        std::cout << "Status: " << (converged ? "✓ CONVERGED" : "✗ NOT CONVERGED") << "\n";
        std::cout << "Iterations: " << iterations << "\n\n";
        std::cout << "Final RMS Total: " << std::fixed << std::setprecision(3) << statistics.rms_total << " arcsec\n";
        std::cout << "Observations Used: " << (statistics.num_observations - statistics.num_outliers) << "\n";
        std::cout << "χ²/dof: " << statistics.reduced_chi_squared << "\n\n";
        std::cout << "========================================\n\n";
    }
};

/**
 * @brief Differential corrector for orbit determination
 * 
 * @tparam Frame Reference frame for integration
 */
template <typename Frame>
class DifferentialCorrector {
public:
    DifferentialCorrector(
        std::shared_ptr<ResidualCalculator<Frame>> residual_calc,
        std::shared_ptr<StateTransitionMatrix<Frame>> stm_computer)
        : residual_calc_(residual_calc),
          stm_computer_(stm_computer) {}
    
    DifferentialCorrectorResult<Frame> fit(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const physics::CartesianStateTyped<Frame>& initial_guess,
        const DifferentialCorrectorSettings& settings = {}) 
    {
        DifferentialCorrectorResult<Frame> result;
        result.final_state = initial_guess;
        result.converged = false;
        result.iterations = 0;
        physics::CartesianStateTyped<Frame> current_state = initial_guess;
        
        if (settings.verbose) {
            std::cout << "\n========================================\n";
            std::cout << "Differential Corrections\n";
            std::cout << "========================================\n";
            std::cout << "Observations: " << observations.size() << "\n";
            std::cout << "Max iterations: " << settings.max_iterations << "\n";
            std::cout << "Convergence: " << settings.convergence_tolerance << " AU\n\n";
        }
        
        double current_sigma = std::max(settings.outlier_sigma, settings.outlier_max_sigma);
        
        for (int iter = 0; iter < settings.max_iterations; ++iter) {
            result.iterations = iter + 1;
            Eigen::VectorXd correction;
            std::vector<ObservationResidual> residuals;
            
            bool iter_success = iteration(observations, current_state, correction, residuals);
            if (!iter_success) break;
            
            if (settings.reject_outliers && iter > 0) {
                ResidualCalculator<Frame>::identify_outliers(residuals, current_sigma);
            }
            
            auto stats = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
            result.rms_history.push_back(stats.rms_total);
            result.correction_norm.push_back(correction.norm());
            
            if (iteration_callback_) iteration_callback_(iter + 1, stats);
            
            // Apply correction (AU -> SI)
            current_state = physics::CartesianStateTyped<Frame>::from_si(
                current_state.epoch,
                current_state.position.x_si() + correction[0] * constants::AU * 1000.0,
                current_state.position.y_si() + correction[1] * constants::AU * 1000.0,
                current_state.position.z_si() + correction[2] * constants::AU * 1000.0,
                current_state.velocity.x_si() + correction[3] * constants::AU * 1000.0 / 86400.0,
                current_state.velocity.y_si() + correction[4] * constants::AU * 1000.0 / 86400.0,
                current_state.velocity.z_si() + correction[5] * constants::AU * 1000.0 / 86400.0,
                current_state.gm.to_m3_s2()
            );
            
            if (check_convergence(correction, settings.convergence_tolerance)) {
                 result.converged = true;
                 result.final_state = current_state;
                 result.residuals = residuals;
                 result.statistics = stats;
                 break;
            }
            result.final_state = current_state;
            result.residuals = residuals;
            result.statistics = stats;
        }
        
        if (settings.compute_covariance && result.converged) {
            result.covariance = compute_covariance(observations, result.final_state, result.residuals);
            result.formal_uncertainties.resize(6);
            for (int i = 0; i < 6; ++i) result.formal_uncertainties[i] = std::sqrt(result.covariance(i, i));
            result.correlation = compute_correlation(result.covariance);
        }
        return result;
    }
    
    bool iteration(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const physics::CartesianStateTyped<Frame>& current_state,
        Eigen::VectorXd& correction,
        std::vector<ObservationResidual>& residuals) 
    {
        residuals = residual_calc_->compute_residuals(observations, current_state);
        if (residuals.empty()) return false;
        
        auto design_result = build_design_matrix(observations, current_state, residuals);
        if (design_result.valid_indices.empty()) return false;
        
        astdyn::Matrix6d normal_matrix, normal_inv;
        auto correction_opt = solve_normal_equations(
            design_result.A, design_result.b, design_result.weights,
            normal_matrix, normal_inv);
        
        if (!correction_opt) return false;
        correction = *correction_opt;
        return true;
    }
    
    astdyn::Matrix6d compute_covariance(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const physics::CartesianStateTyped<Frame>& final_state,
        const std::vector<ObservationResidual>& residuals) 
    {
        auto design_result = build_design_matrix(observations, final_state, residuals);
        astdyn::Matrix6d normal_matrix, normal_inv;
        auto solution = solve_normal_equations(
            design_result.A, design_result.b, design_result.weights,
            normal_matrix, normal_inv);
        
        if (!solution) return astdyn::Matrix6d::Zero();
        auto stats = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
        double sigma_0_squared = (stats.degrees_of_freedom > 0) ? stats.chi_squared / stats.degrees_of_freedom : 1.0;
        return sigma_0_squared * normal_inv;
    }
    
    using IterationCallback = std::function<void(int, const ResidualStatistics&)>;
    void set_iteration_callback(IterationCallback callback) {
        iteration_callback_ = callback;
    }

private:
    struct DesignMatrixResult {
        Eigen::MatrixXd A;
        Eigen::VectorXd b;
        Eigen::VectorXd weights;
        std::vector<size_t> valid_indices;
    };
    
    DesignMatrixResult build_design_matrix(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const physics::CartesianStateTyped<Frame>& state,
        const std::vector<ObservationResidual>& residuals) 
    {
        DesignMatrixResult result;
        result.valid_indices.reserve(observations.size());
        int n_valid = 0;
        for (size_t i = 0; i < residuals.size(); ++i) {
            if (!residuals[i].outlier) {
                n_valid++;
                result.valid_indices.push_back(i);
            }
        }
        if (n_valid == 0) return result;
        
        int n_equations = 2 * n_valid;
        result.A.resize(n_equations, 6);
        result.b.resize(n_equations);
        result.weights.resize(n_equations);
        
        int row = 0;
        for (size_t idx : result.valid_indices) {
            const auto& obs = observations[idx];
            const auto& res = residuals[idx];
            auto obs_pos_opt = residual_calc_->get_observer_position(obs);
            if (!obs_pos_opt) continue;
            
            auto observer_pos = *obs_pos_opt;
            auto obs_time_tdb = astdyn::time::to_tdb(obs.time);
            auto partials = stm_computer_->compute_with_partials(state, obs_time_tdb, observer_pos);
            
            Eigen::Matrix<double, 2, 6> A_obs = partials.partial_radec * partials.phi;
            result.A.row(row) = A_obs.row(0);
            result.A.row(row + 1) = A_obs.row(1);
            result.b[row] = res.residual_ra;
            result.b[row + 1] = res.residual_dec;
            result.weights[row] = 1.0 / (obs.sigma_ra * obs.sigma_ra);
            result.weights[row + 1] = 1.0 / (obs.sigma_dec * obs.sigma_dec);
            row += 2;
        }
        return result;
    }
    
    std::optional<Eigen::VectorXd> solve_normal_equations(
        const Eigen::MatrixXd& A,
        const Eigen::VectorXd& b,
        const Eigen::VectorXd& W,
        astdyn::Matrix6d& normal_matrix,
        astdyn::Matrix6d& normal_inv) 
    {
        Eigen::MatrixXd W_diag = W.asDiagonal();
        normal_matrix = A.transpose() * W_diag * A;
        Eigen::VectorXd rhs = A.transpose() * W_diag * b;
        
        Eigen::LDLT<astdyn::Matrix6d> ldlt(normal_matrix);
        if (ldlt.info() != Eigen::Success) return std::nullopt;
        
        Eigen::VectorXd correction = ldlt.solve(rhs);
        if (ldlt.info() != Eigen::Success) return std::nullopt;
        
        normal_inv = ldlt.solve(astdyn::Matrix6d::Identity());
        return correction;
    }
    
    astdyn::Matrix6d compute_correlation(const astdyn::Matrix6d& covariance) const {
        astdyn::Matrix6d correlation = astdyn::Matrix6d::Zero();
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                double denom = std::sqrt(covariance(i, i) * covariance(j, j));
                if (denom > 0.0) correlation(i, j) = covariance(i, j) / denom;
            }
        }
        return correlation;
    }
    
    bool check_convergence(const Eigen::VectorXd& correction, double tolerance) const {
        return correction.head<3>().norm() < tolerance;
    }

private:
    std::shared_ptr<ResidualCalculator<Frame>> residual_calc_;
    std::shared_ptr<StateTransitionMatrix<Frame>> stm_computer_;
    IterationCallback iteration_callback_;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP
