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

    DifferentialCorrectorResult<Frame> fit(const std::vector<astdyn::observations::OpticalObservation>& observations, 
                                           const physics::CartesianStateTyped<Frame>& initial_state, 
                                           const DifferentialCorrectorSettings& settings = {}) {
        DifferentialCorrectorResult<Frame> res;
        auto sorted_obs = prepare_observations(observations);
        physics::CartesianStateTyped<Frame> current_state = initial_state;
        double current_sigma = std::max(settings.outlier_sigma, settings.outlier_max_sigma);
        double prev_rms = 1e18;

        for (int i = 0; i < settings.max_iterations; ++i) {
            res.iterations = i + 1;
            auto residuals = residual_calc_->compute_residuals(sorted_obs, current_state);
            handle_carpentry_sigma(i, settings, current_sigma, residuals);
            
            auto dm = build_design_matrix(sorted_obs, current_state, residuals);
            if (dm.valid_indices.empty()) { res.rejection_reason = "No valid observations remaining"; break; }
            
            Eigen::VectorXd corr = solve_normal_equations(dm);
            double cur_rms = ResidualCalculator<Frame>::compute_statistics(residuals, 6).rms_total.to_arcsec();
            
            if (settings.use_line_search) {
                if (!perform_line_search(sorted_obs, corr, settings, current_state, residuals, cur_rms)) {
                    res.rejection_reason = "Line search failed to find better solution"; break;
                }
            } else {
                current_state = apply_correction(current_state, corr, 1.0);
            }
            
            res.statistics = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
            res.residuals = residuals;
            res.rms_history.push_back(res.statistics.rms_total.to_arcsec());
            res.correction_norm.push_back(corr.norm());

            if (check_convergence(corr.norm(), res.statistics.rms_total.to_arcsec(), prev_rms, settings)) {
                res.converged = true; break;
            }
            prev_rms = res.statistics.rms_total.to_arcsec();
        }
        res.final_state = current_state;
        if (settings.compute_covariance) res.covariance = compute_covariance_internal(sorted_obs, current_state, res.residuals);
        return res;
    }

private:
    struct DesignMatrix {
        Eigen::MatrixXd A; Eigen::VectorXd b; Eigen::VectorXd weights; std::vector<size_t> valid_indices;
    };

    std::vector<observations::OpticalObservation> prepare_observations(const std::vector<observations::OpticalObservation>& obs) {
        auto sorted = obs;
        std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) { return a.time < b.time; });
        return sorted;
    }

    bool check_convergence(double corr_norm, double cur_rms, double prev_rms, const DifferentialCorrectorSettings& settings) {
        if (corr_norm < 1e-12) return true;
        if (std::abs(prev_rms - cur_rms) < settings.rms_tolerance_arcsec) return true;
        return false;
    }

    void handle_carpentry_sigma(int iter, const DifferentialCorrectorSettings& settings, double& current_sigma, std::vector<ObservationResidual>& residuals) {
        if (!settings.reject_outliers) return;
        auto stats = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
        if (iter > 0 && stats.rms_total.to_arcsec() < 1000.0) {
            double t = (double)iter / (double)settings.max_iterations;
            current_sigma = settings.outlier_max_sigma + t * (settings.outlier_min_sigma - settings.outlier_max_sigma);
            ResidualCalculator<Frame>::identify_outliers(residuals, current_sigma);
        }
    }

    Eigen::VectorXd solve_normal_equations(const DesignMatrix& dm) {
        Eigen::MatrixXd AtW = dm.A.transpose() * dm.weights.asDiagonal();
        return (AtW * dm.A).ldlt().solve(AtW * dm.b);
    }

    bool perform_line_search(const std::vector<observations::OpticalObservation>& obs, const Eigen::VectorXd& correction, const DifferentialCorrectorSettings& settings, physics::CartesianStateTyped<Frame>& current_state, std::vector<ObservationResidual>& residuals, double& cur_rms) {
        double alpha = 1.0;
        while (alpha >= settings.line_search_min_alpha) {
            auto trial = apply_correction(current_state, correction, alpha);
            auto trial_res = residual_calc_->compute_residuals(obs, trial);
            for (size_t i=0; i<residuals.size(); ++i) trial_res[i].outlier = residuals[i].outlier;
            auto t_stats = ResidualCalculator<Frame>::compute_statistics(trial_res, 6);
            if (t_stats.rms_total.to_arcsec() < cur_rms) { current_state = trial; residuals = trial_res; cur_rms = t_stats.rms_total.to_arcsec(); return true; }
            alpha *= 0.5;
        }
        return false;
    }

    physics::CartesianStateTyped<Frame> apply_correction(const physics::CartesianStateTyped<Frame>& s, const Eigen::VectorXd& c, double a) {
        Eigen::VectorXd y = s.to_eigen_au_aud(); y += a * c;
        return physics::CartesianStateTyped<Frame>::from_au_aud(s.epoch, y, s.gm);
    }

    DesignMatrix build_design_matrix(const std::vector<observations::OpticalObservation>& obs, const physics::CartesianStateTyped<Frame>& s, const std::vector<ObservationResidual>& res) {
        DesignMatrix dm;
        for (size_t i=0; i<res.size(); ++i) if (!res[i].outlier) dm.valid_indices.push_back(i);
        if (dm.valid_indices.empty()) return dm;

        int n = dm.valid_indices.size();
        dm.A.resize(2*n, 6); dm.b.resize(2*n); dm.weights.resize(2*n);
        std::vector<time::EpochTDB> times; std::vector<math::Vector3<core::GCRF, physics::Distance>> positions;
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
        return (AtW * dm.A).inverse();
    }

    std::shared_ptr<ResidualCalculator<Frame>> residual_calc_;
    std::shared_ptr<StateTransitionMatrix<Frame>> stm_computer_;
};

    std::shared_ptr<ResidualCalculator<Frame>> residual_calc_;
    std::shared_ptr<StateTransitionMatrix<Frame>> stm_computer_;
};

} // namespace astdyn::orbit_determination

#endif
