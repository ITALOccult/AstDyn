/**
 * @file DifferentialCorrector.hpp
 * @brief Differential corrector for orbit determination
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP
#define ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include <memory>
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace astdyn::orbit_determination {

struct DifferentialCorrectorSettings {
    int max_iterations = 20;
    physics::Distance convergence_tolerance = physics::Distance::from_au(1e-6);
    double outlier_sigma = 3.0;
    double outlier_max_sigma = 10.0;
    double outlier_min_sigma = 3.0;
    bool reject_outliers = true;
    bool compute_covariance = true;
    bool verbose = false;
    bool use_line_search = true;
    double line_search_min_alpha = 1e-4;
    double rms_tolerance_arcsec = 0.001;
    bool check_energy_barrier = true;
    double energy_barrier_fraction = 0.5;
};

template <typename Frame>
struct DifferentialCorrectorResult {
    physics::CartesianStateTyped<Frame> final_state;
    std::vector<ObservationResidual> residuals;
    ResidualStatistics statistics;
    bool converged = false;
    int iterations = 0;
    std::vector<double> rms_history;
    std::vector<double> correction_norm;
    astdyn::Matrix6d covariance = astdyn::Matrix6d::Zero();
    std::string rejection_reason;

    void print() const {
        std::cout << "\n========================================\n";
        std::cout << "OD Result Summary\n";
        std::cout << "========================================\n";
        std::cout << "Status: " << (converged ? "✓ CONVERGED" : "✗ NOT CONVERGED") << "\n";
        std::cout << "Iterations: " << iterations << "\n";
        std::cout << "Final RMS: " << statistics.rms_total.to_arcsec() << " arcsec\n";
        std::cout << "========================================\n\n";
    }
};

template <typename Frame>
class DifferentialCorrector {
public:
    DifferentialCorrector(
        std::shared_ptr<ResidualCalculator<Frame>> rc,
        std::shared_ptr<StateTransitionMatrix<Frame>> stm)
        : residual_calc_(rc), stm_computer_(stm) {}

    std::shared_ptr<ResidualCalculator<Frame>> get_residual_calculator() const { return residual_calc_; }

    DifferentialCorrectorResult<Frame> fit(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const physics::CartesianStateTyped<Frame>& initial_guess,
        const DifferentialCorrectorSettings& settings = {}) 
    {
        DifferentialCorrectorResult<Frame> result;
        result.final_state = initial_guess;
        
        std::vector<observations::OpticalObservation> sorted_obs = observations;
        std::sort(sorted_obs.begin(), sorted_obs.end(), [](const auto& a, const auto& b) {
            return a.time < b.time;
        });

        physics::CartesianStateTyped<Frame> current_state = initial_guess;
        double prev_iter_rms = 1e18;
        double current_sigma = std::max(settings.outlier_sigma, settings.outlier_max_sigma);

        for (int iter = 0; iter < settings.max_iterations; ++iter) {
            result.iterations = iter + 1;
            
            // 1. Compute residuals
            std::vector<ObservationResidual> residuals = residual_calc_->compute_residuals(sorted_obs, current_state);
            
            // 2. Identify outliers (Carpentry)
            if (settings.reject_outliers) {
                auto stats = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
                if (iter > 0 && stats.rms_total.to_arcsec() < 1000.0) {
                    double t = (double)iter / (double)settings.max_iterations;
                    current_sigma = settings.outlier_max_sigma + t * (settings.outlier_min_sigma - settings.outlier_max_sigma);
                    ResidualCalculator<Frame>::identify_outliers(residuals, current_sigma);
                }
            }

            // 3. Build Normal Equations
            auto dm = build_design_matrix(sorted_obs, current_state, residuals);
            if (dm.valid_indices.empty()) break;

            Eigen::MatrixXd AtW = dm.A.transpose() * dm.weights.asDiagonal();
            Eigen::MatrixXd NormalMat = AtW * dm.A;
            Eigen::VectorXd RHS = AtW * dm.b;
            Eigen::VectorXd correction = NormalMat.ldlt().solve(RHS);

            // 4. Line Search
            auto stats_before = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
            double cur_rms = stats_before.rms_total.to_arcsec();
            
            if (settings.use_line_search) {
                double alpha = 1.0;
                bool success = false;
                if (settings.verbose) std::cout << "    [LSQ] Line search (cur_rms=" << cur_rms << "\")" << std::endl;
                while (alpha >= settings.line_search_min_alpha) {
                    auto trial = apply_correction(current_state, correction, alpha);
                    auto trial_res = residual_calc_->compute_residuals(sorted_obs, trial);
                    for (size_t i=0; i<residuals.size(); ++i) trial_res[i].outlier = residuals[i].outlier;
                    auto trial_stats = ResidualCalculator<Frame>::compute_statistics(trial_res, 6);
                    double trial_rms = trial_stats.rms_total.to_arcsec();
                    
                    if (settings.verbose) std::cout << "      alpha=" << alpha << " -> trial_rms=" << trial_rms << "\"" << std::endl;

                    if (trial_stats.rms_total.to_arcsec() < cur_rms) {
                        current_state = trial;
                        success = true;
                        break;
                    }
                    alpha *= 0.5;
                }
                if (!success) {
                    if (settings.verbose) std::cout << "    [LSQ] Line search failed (reached min alpha)." << std::endl;
                    break;
                }
            } else {
                current_state = apply_correction(current_state, correction, 1.0);
            }

            // 5. Convergence check
            auto final_stats = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
            double final_rms = final_stats.rms_total.to_arcsec();
            if (settings.verbose) {
                std::cout << "  - Iteration " << (iter+1) << ": RMS " << final_rms << "\"" 
                          << " (correction norm: " << correction.norm() << ")" << std::endl;
            }
            if (iter > 0 && std::abs(prev_iter_rms - final_rms) < settings.rms_tolerance_arcsec) {
                result.converged = true;
                if (settings.verbose) std::cout << "    [LSQ] Converged by ΔRMS." << std::endl;
                break;
            }
            if (correction.norm() < 1e-12) {
                result.converged = true;
                if (settings.verbose) std::cout << "    [LSQ] Converged by correction norm." << std::endl;
                break;
            }
            prev_iter_rms = final_rms;
            result.statistics = final_stats;
            result.residuals = residuals;
        }

        result.final_state = current_state;
        
        if (settings.compute_covariance) {
            result.covariance = compute_covariance_internal(sorted_obs, current_state, result.residuals);
        }

        return result;
    }

private:
    physics::CartesianStateTyped<Frame> apply_correction(const physics::CartesianStateTyped<Frame>& s, const Eigen::VectorXd& c, double a) {
        Eigen::VectorXd y = s.to_eigen_au_aud();
        y += a * c;
        return physics::CartesianStateTyped<Frame>::from_au_aud(s.epoch, y, s.gm);
    }

    struct DesignMatrix {
        Eigen::MatrixXd A;
        Eigen::VectorXd b;
        Eigen::VectorXd weights;
        std::vector<size_t> valid_indices;
    };

    DesignMatrix build_design_matrix(const std::vector<observations::OpticalObservation>& obs, const physics::CartesianStateTyped<Frame>& s, const std::vector<ObservationResidual>& res) {
        DesignMatrix dm;
        for (size_t i=0; i<res.size(); ++i) if (!res[i].outlier) dm.valid_indices.push_back(i);
        if (dm.valid_indices.empty()) return dm;

        int n = dm.valid_indices.size();
        dm.A.resize(2*n, 6); dm.b.resize(2*n); dm.weights.resize(2*n);

        std::vector<time::EpochTDB> times;
        std::vector<math::Vector3<core::GCRF, physics::Distance>> positions;
        for (size_t idx : dm.valid_indices) {
            times.push_back(astdyn::time::to_tdb(obs[idx].time));
            auto p = residual_calc_->get_observer_position(obs[idx]);
            positions.push_back(p ? *p : math::Vector3<core::GCRF, physics::Distance>());
        }

        auto batch = stm_computer_->compute_batch(s, times, positions);
        for (int i=0; i<n; ++i) {
            dm.A.template block<2,6>(2*i, 0) = batch[i].partial_radec * batch[i].phi;
            dm.b[2*i] = res[dm.valid_indices[i]].residual_ra.to_rad();
            dm.b[2*i+1] = res[dm.valid_indices[i]].residual_dec.to_rad();
            dm.weights[2*i] = res[dm.valid_indices[i]].weight_ra;
            dm.weights[2*i+1] = res[dm.valid_indices[i]].weight_dec;
        }
        return dm;
    }

    astdyn::Matrix6d compute_covariance_internal(const std::vector<observations::OpticalObservation>& obs, const physics::CartesianStateTyped<Frame>& s, const std::vector<ObservationResidual>& res) {
        auto dm = build_design_matrix(obs, s, res);
        if (dm.valid_indices.empty()) return astdyn::Matrix6d::Zero();
        Eigen::MatrixXd AtW = dm.A.transpose() * dm.weights.asDiagonal();
        Eigen::MatrixXd NormalMat = AtW * dm.A;
        return NormalMat.inverse();
    }

    std::shared_ptr<ResidualCalculator<Frame>> residual_calc_;
    std::shared_ptr<StateTransitionMatrix<Frame>> stm_computer_;
};

} // namespace astdyn::orbit_determination

#endif
