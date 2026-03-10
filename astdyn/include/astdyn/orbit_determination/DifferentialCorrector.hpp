/**
 * @file DifferentialCorrector.hpp
 * @brief Differential corrections for orbit determination
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP
#define ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/orbit_determination/ODSmartPolicy.hpp"
#include <memory>
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace astdyn::orbit_determination {

/**
 * @brief Settings for differential corrections
 */
struct DifferentialCorrectorSettings {
    int max_iterations = 20;             ///< Maximum iterations
    physics::Distance convergence_tolerance = physics::Distance::from_au(1e-6); ///< Convergence threshold
    double outlier_sigma = 3.0;          ///< Target Sigma threshold (default)
    
    // Carpentry Settings (Iterative Rejection)
    double outlier_max_sigma = 10.0;     ///< Starting loose sigma
    double outlier_min_sigma = 3.0;      ///< Final tight sigma
    
    bool reject_outliers = true;         ///< Automatically reject outliers
    bool compute_covariance = true;      ///< Compute covariance matrix
    bool verbose = false;                ///< Print iteration details
    bool use_line_search = true;         ///< Backtracking line search to prevent divergence
    double line_search_min_alpha = 1e-4; ///< Minimum step fraction before giving up
    double rms_tolerance_arcsec  = 0.001; ///< Absolute ΔRMS convergence threshold [arcsec]

    // Rule 2 — Energy barrier (OrbFit-style non-physical solution rejection)
    bool  check_energy_barrier    = true;  ///< Reject solution if SMA drifts > threshold
    double energy_barrier_fraction = 0.5;  ///< Max allowed |Δa|/a₀ before rejection (0.5 = 50%)
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

    // Rule 2 — Energy barrier
    std::string rejection_reason;        ///< Set if convergence blocked by energy barrier
    
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
        std::cout << "Final RMS Total: " << std::fixed << std::setprecision(3) << statistics.rms_total.to_arcsec() << " arcsec\n";
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
        
        // Ensure observations are sorted by time (structural requirement for sequential propagation)
        std::vector<observations::OpticalObservation> sorted_obs = observations;
        std::sort(sorted_obs.begin(), sorted_obs.end(), [](const auto& a, const auto& b) {
            return a.time < b.time;
        });
        
        physics::CartesianStateTyped<Frame> current_state = initial_guess;

        // Rule 2 — Energy barrier: record initial semi-major axis (AU)
        // Only meaningful for GCRF frames; for other frames compute_sma_au still
        // works because the speed/distance magnitudes are identical.
        double a0_au = 0.0;
        if (settings.check_energy_barrier) {
            a0_au = ODPolicyEngine::compute_sma_au(
                physics::CartesianStateTyped<core::GCRF>::from_si(
                    initial_guess.epoch,
                    initial_guess.position.x_si(), initial_guess.position.y_si(), initial_guess.position.z_si(),
                    initial_guess.velocity.x_si(), initial_guess.velocity.y_si(), initial_guess.velocity.z_si(),
                    initial_guess.gm.to_m3_s2()));
        }

        if (settings.verbose) {
            std::cout << "\n========================================\n";
            std::cout << "Differential Corrections\n";
            std::cout << "========================================\n";
            std::cout << "Observations: " << observations.size() << "\n";
            std::cout << "Max iterations: " << settings.max_iterations << "\n";
            std::cout << "Convergence: " << settings.convergence_tolerance.to_au() << " AU\n\n";
        }
        
        double current_sigma = std::max(settings.outlier_sigma, settings.outlier_max_sigma);
        double prev_iter_rms = std::numeric_limits<double>::max();

        for (int iter = 0; iter < settings.max_iterations; ++iter) {
            result.iterations = iter + 1;
            Eigen::VectorXd correction;
            std::vector<ObservationResidual> residuals;

            bool iter_success = iteration(sorted_obs, current_state, correction, residuals);
            if (!iter_success) break;
            
            if (settings.reject_outliers) {
                if (iter == 0) {
                    // Rule 3 — After iteration 1, if RMS > 3" reject the single worst outlier.
                    // This prevents one bad observation from biasing the normal equations.
                    auto stats0 = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
                    if (stats0.rms_total.to_arcsec() > 3.0) {
                        size_t worst_idx = 0;
                        double worst_chi2 = -1.0;
                        for (size_t i = 0; i < residuals.size(); ++i) {
                            if (!residuals[i].outlier && residuals[i].chi_squared > worst_chi2) {
                                worst_chi2 = residuals[i].chi_squared;
                                worst_idx = i;
                            }
                        }
                        if (worst_chi2 > 0.0) {
                            residuals[worst_idx].outlier = true;
                            if (settings.verbose) {
                                std::cout << "  [Rule3] RMS=" << stats0.rms_total.to_arcsec()
                                          << "\" > 3\": rejected worst obs (chi²="
                                          << worst_chi2 << ")\n";
                            }
                        }
                    }
                } else {
                    // Carpentry: progressively tighten sigma from max to min
                    double t = static_cast<double>(iter) / std::max(settings.max_iterations - 1, 1);
                    current_sigma = settings.outlier_max_sigma + t * (settings.outlier_min_sigma - settings.outlier_max_sigma);
                    ResidualCalculator<Frame>::identify_outliers(residuals, current_sigma);
                }
            }
            
            auto stats = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
            result.rms_history.push_back(stats.rms_total.to_arcsec());
            result.correction_norm.push_back(correction.norm());

            double cur_rms = stats.rms_total.to_arcsec();

            if (settings.verbose) {
                std::cout << "  - Iteration " << std::setw(2) << (iter + 1)
                          << ": RMS = " << std::setw(8) << std::fixed << std::setprecision(3) << cur_rms << "\""
                          << " (used " << (stats.num_observations - stats.num_outliers) << " obs)" << std::endl;
                std::cout.flush();
            }

            // RMS stagnation: converge when |ΔRMS| < rms_tolerance_arcsec (absolute)
            if (iter > 0 && std::isfinite(cur_rms) && std::isfinite(prev_iter_rms) &&
                std::abs(prev_iter_rms - cur_rms) < settings.rms_tolerance_arcsec) {
                result.converged = true;
                result.final_state = current_state;
                result.residuals = residuals;
                result.statistics = stats;
                if (settings.verbose) {
                    std::cout << "  → Converged (ΔRMS=" << std::fixed << std::setprecision(4)
                              << std::abs(prev_iter_rms - cur_rms)
                              << "\" < " << settings.rms_tolerance_arcsec << "\")." << std::endl;
                }
                break;
            }
            prev_iter_rms = cur_rms;

            if (iteration_callback_) iteration_callback_(iter + 1, stats);
            
            // Apply correction: design matrix is now in [rad/m], so correction is in (m, m/s).
            // This matches our SI-based state storage.
            auto make_trial_state = [&](double alpha) {
                return physics::CartesianStateTyped<Frame>(
                    current_state.epoch,
                    astdyn::math::Vector3<Frame, physics::Distance>::from_si(
                        current_state.position.x_si() + alpha * correction[0],
                        current_state.position.y_si() + alpha * correction[1],
                        current_state.position.z_si() + alpha * correction[2]),
                    astdyn::math::Vector3<Frame, physics::Velocity>::from_si(
                        current_state.velocity.x_si() + alpha * correction[3],
                        current_state.velocity.y_si() + alpha * correction[4],
                        current_state.velocity.z_si() + alpha * correction[5]),
                    current_state.gm
                );
            };

            if (settings.use_line_search) {
                double prev_rms = stats.rms_total.to_arcsec();
                double alpha = 1.0;
                bool step_accepted = false;
                while (alpha >= settings.line_search_min_alpha) {
                    auto trial = make_trial_state(alpha);
                    auto trial_res = residual_calc_->compute_residuals(sorted_obs, trial);
                    auto trial_stats = ResidualCalculator<Frame>::compute_statistics(trial_res, 6);
                    double trial_rms = trial_stats.rms_total.to_arcsec();
                    if (std::isfinite(trial_rms) && trial_rms < prev_rms) {
                        current_state = trial;
                        if (settings.verbose && alpha < 1.0) {
                            std::cout << "    (line search: alpha=" << std::setprecision(4) << alpha
                                      << ", rms " << prev_rms << "\" → " << trial_rms << "\")" << std::endl;
                        }
                        step_accepted = true;
                        break;
                    }
                    alpha *= 0.5;
                }
                if (!step_accepted) {
                    // Line search couldn't improve RMS — check if correction is negligible
                    // Tolerance is in AU, so convert to meters
                    if (correction.head<3>().norm() < 10.0 * settings.convergence_tolerance.to_m()) {
                        result.converged = true;
                        result.final_state = current_state;
                        result.residuals = residuals;
                        result.statistics = stats;
                        if (settings.verbose) std::cout << "  → Converged (stagnated)." << std::endl;
                    } else {
                        if (settings.verbose) std::cout << "  → Line search failed, stopping." << std::endl;
                    }
                    break;
                }
            } else {
                current_state = make_trial_state(1.0);
            }

            // Rule 2 — Energy barrier: reject if SMA drifts > 50% from initial value
            if (settings.check_energy_barrier && a0_au > 0.0) {
                double a_now = ODPolicyEngine::compute_sma_au(
                    physics::CartesianStateTyped<core::GCRF>::from_si(
                        current_state.epoch,
                        current_state.position.x_si(), current_state.position.y_si(), current_state.position.z_si(),
                        current_state.velocity.x_si(), current_state.velocity.y_si(), current_state.velocity.z_si(),
                        current_state.gm.to_m3_s2()));
                double drift = std::abs(a_now - a0_au) / a0_au;
                if (drift > settings.energy_barrier_fraction) {
                    result.converged = false;
                    result.rejection_reason = "energy_barrier: a drifted " +
                        std::to_string(static_cast<int>(drift * 100)) + "% (>" +
                        std::to_string(static_cast<int>(settings.energy_barrier_fraction * 100)) + "%)";
                    if (settings.verbose) {
                        std::cout << "  → [EnergyBarrier] Non-physical solution rejected: "
                                  << result.rejection_reason << "\n";
                    }
                    break;
                }
            }

            // Convergence check: correction is in meters
            if (correction.head<3>().norm() < settings.convergence_tolerance.to_m()) {
                 result.converged = true;
                 result.final_state = current_state;
                 result.residuals = residuals;
                 result.statistics = stats;
                 if (settings.verbose) std::cout << "  → Converged!" << std::endl;
                 break;
            }
            result.final_state = current_state;
            result.residuals = residuals;
            result.statistics = stats;
        }
        
        if (settings.compute_covariance && result.converged) {
            result.covariance = compute_covariance(sorted_obs, result.final_state, result.residuals);
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
    
    std::shared_ptr<ResidualCalculator<Frame>> get_residual_calculator() const {
        return residual_calc_;
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
        physics::CartesianStateTyped<Frame> running_state = state;
        astdyn::Matrix6d cumulative_phi = astdyn::Matrix6d::Identity();

        for (size_t idx : result.valid_indices) {
            const auto& obs = observations[idx];
            const auto& res = residuals[idx];
            auto obs_pos_opt = residual_calc_->get_observer_position(obs);
            if (!obs_pos_opt) continue;
            
            auto observer_pos = *obs_pos_opt;
            auto obs_time_tdb = astdyn::time::to_tdb(obs.time);
            
            // Sequential propagation of state and STM (STM comes in AU/Normal units)
            auto partials = stm_computer_->compute_with_partials(running_state, obs_time_tdb, observer_pos);
            cumulative_phi = partials.phi * cumulative_phi;
            running_state = partials.final_state;

            // Transform cumulative_phi (AU-based) to SI-based (meters, m/s)
            astdyn::Matrix6d phi_si = cumulative_phi;
            // dr_si = au_m * dr_au
            // dv_si = v_scale * dv_au
            // d(dr_si)/d(dv0_si) = phi_12 * (dx_si/dx_au) * (dv0_au_d/dv0_si)
            
            const double v_scale_phi = physics::Velocity::from_au_d(1.0).to_ms(); // (AU*1000)/86400
            const double au_m = physics::Distance::from_au(1.0).to_m();
            
            // d(pos_si)/d(vel0_si) = d(pos_au*au_m) / d(vel0_ms) 
            // = au_m * d(pos_au)/d(vel0_aud) * d(vel0_aud)/d(vel0_ms)
            // = au_m * phi_12 * (1.0 / v_scale_phi)
            // = phi_12 * (au_m / v_scale_phi) = phi_12 * 86400
            
            phi_si.block<3, 3>(0, 3) *= (au_m / v_scale_phi);
            phi_si.block<3, 3>(3, 0) *= (v_scale_phi / au_m);

            Eigen::Matrix<double, 2, 6> A_obs = partials.partial_radec * phi_si;
            result.A.row(row) = A_obs.row(0);
            result.A.row(row + 1) = A_obs.row(1);
            result.b[row] = res.residual_ra.to_rad();
            result.b[row + 1] = res.residual_dec.to_rad();
            // Use the weights pre-computed in ResidualCalculator, which already account
            // for the cos(dec) projection of the RA sigma (projected sigma).
            result.weights[row]     = res.weight_ra;
            result.weights[row + 1] = res.weight_dec;
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
