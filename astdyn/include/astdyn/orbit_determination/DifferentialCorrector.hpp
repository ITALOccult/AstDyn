/**
 * @file DifferentialCorrector.hpp
 * @brief Differential corrector for orbit determination
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP
#define ASTDYN_ORBIT_DETERMINATION_DIFFERENTIAL_CORRECTOR_HPP

#include "astdyn/core/Types.hpp"
#include <chrono>
#include <iostream>
#include "astdyn/core/physics_state.hpp"
#include "astdyn/core/IOCConfig.hpp"
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
    int outlier_ramp_steps = 3;              ///< Iterazioni per passare da
                                             ///< outlier_max_sigma a outlier_min_sigma.
                                             ///< Slegato da max_iterations, che e' un
                                             ///< limite di sicurezza e non deve
                                             ///< influenzare quali osservazioni si scartano.
    int max_line_search_steps = 5;           ///< Tetto ai dimezzamenti di alpha: ogni
                                             ///< tentativo costa una propagazione completa,
                                             ///< e da 1 a 1e-4 sarebbero 14 tentativi.
    double rms_tolerance_arcsec = 0.001;
    bool check_energy_barrier = true;
    double energy_barrier_fraction = 0.5;

    static DifferentialCorrectorSettings from_config(const core::IOCConfig& cfg) {
        DifferentialCorrectorSettings s;
        s.max_iterations       = cfg.get<int>   ("diffcorr.max_iter",            s.max_iterations);
        s.outlier_sigma        = cfg.get<double>("diffcorr.outlier_sigma",        s.outlier_sigma);
        s.outlier_max_sigma    = cfg.get<double>("diffcorr.outlier_max_sigma",    s.outlier_max_sigma);
        s.outlier_min_sigma    = cfg.get<double>("diffcorr.outlier_min_sigma",    s.outlier_min_sigma);
        s.outlier_ramp_steps   = cfg.get<int>   ("diffcorr.outlier_ramp_steps",   s.outlier_ramp_steps);
        s.rms_tolerance_arcsec = cfg.get<double>("diffcorr.rms_tolerance_arcsec",s.rms_tolerance_arcsec);
        s.reject_outliers      = cfg.get<bool>  ("diffcorr.reject_outliers",      s.reject_outliers);
        s.compute_covariance   = cfg.get<bool>  ("diffcorr.compute_covariance",   s.compute_covariance);
        s.use_line_search      = cfg.get<bool>  ("diffcorr.use_line_search",      s.use_line_search);
        s.check_energy_barrier = cfg.get<bool>  ("diffcorr.check_energy_barrier", s.check_energy_barrier);
        s.energy_barrier_fraction = cfg.get<double>("diffcorr.energy_barrier_fraction", s.energy_barrier_fraction);
        s.verbose              = cfg.get<bool>  ("diffcorr.verbose",              s.verbose);
        return s;
    }
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
            // Le statistiche si calcolano QUI, prima del line search: i residui
            // sono gia' disponibili, e un'uscita anticipata lasciava altrimenti
            // la struttura vuota. Il chiamante leggeva zero osservazioni e RMS
            // zero da un fit che si dichiarava riuscito, e sostituiva comunque
            // l'orbita.
            res.statistics = ResidualCalculator<Frame>::compute_statistics(residuals, 6);
            res.residuals = residuals;
            // Il line search deve valutare la stessa quantita' che le equazioni
            // normali minimizzano: il chi quadro, non l'RMS in arcsec.
            //
            // Su (1216) Askania, con osservazioni dal 1936, i due criteri
            // divergono: una lastra con residuo 30" e sigma 20" contribuisce 900
            // alla somma dei quadrati e 2.25 al chi quadro. La correzione ottima
            // nel senso pesato faceva crescere l'RMS aritmetico e veniva
            // rifiutata a ogni alpha, cosi' che il fit non partiva affatto.
            double cur_rms = res.statistics.weighted_rms > 0.0
                           ? res.statistics.weighted_rms
                           : res.statistics.rms_total.to_arcsec();

            if (settings.use_line_search) {
                if (!perform_line_search(sorted_obs, corr, settings, current_state, residuals, cur_rms)) {
                    if (i == 0) {
                        // Alla PRIMA iterazione nessun passo e' ancora stato
                        // compiuto: se gia' la correzione iniziale non migliora
                        // nulla, il fit non e' partito. Non e' un minimo
                        // raggiunto, e l'orbita di partenza va conservata.
                        res.rejection_reason = "il line search non trova alcun passo utile alla "
                                               "prima iterazione: il fit non e' partito";
                        res.converged = false;
                    } else {
                        // Dopo almeno un passo, "nessun miglioramento possibile"
                        // e' la definizione di minimo raggiunto.
                        res.rejection_reason = "minimo raggiunto: nessuna correzione migliora la soluzione";
                        res.converged = true;
                    }
                    break;
                }
            } else {
                current_state = apply_correction(current_state, corr, 1.0);
            }
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

    std::shared_ptr<ResidualCalculator<Frame>> get_residual_calculator() const { return residual_calc_; }

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
            // La soglia scende di un passo fisso a ogni iterazione, da
            // outlier_max_sigma a outlier_min_sigma in outlier_ramp_steps passi.
            // NON dipende da max_iterations: quello e' un limite di sicurezza, e
            // un fit ben condizionato che converge in due iterazioni deve
            // scartare gli outlier come uno che ne impiega venti.
            const int passi = std::max(1, settings.outlier_ramp_steps);
            const double t = std::min(1.0, (double)iter / (double)passi);
            current_sigma = settings.outlier_max_sigma
                          + t * (settings.outlier_min_sigma - settings.outlier_max_sigma);
            ResidualCalculator<Frame>::identify_outliers(residuals, current_sigma);
        }
    }

    Eigen::VectorXd solve_normal_equations(const DesignMatrix& dm) {
        Eigen::MatrixXd AtW = dm.A.transpose() * dm.weights.asDiagonal();
        return (AtW * dm.A).ldlt().solve(AtW * dm.b);
    }

    bool perform_line_search(const std::vector<observations::OpticalObservation>& obs, const Eigen::VectorXd& correction, const DifferentialCorrectorSettings& settings, physics::CartesianStateTyped<Frame>& current_state, std::vector<ObservationResidual>& residuals, double& cur_rms) {
        double alpha = 1.0;
        int tentativi = 0;
        while (alpha >= settings.line_search_min_alpha && tentativi < settings.max_line_search_steps) {
            ++tentativi;
            auto trial = apply_correction(current_state, correction, alpha);
            auto trial_res = residual_calc_->compute_residuals(obs, trial);
            for (size_t i=0; i<residuals.size(); ++i) trial_res[i].outlier = residuals[i].outlier;
            auto t_stats = ResidualCalculator<Frame>::compute_statistics(trial_res, 6);
            // Stessa metrica del chiamante: chi quadro dove disponibile.
            const double t_val = t_stats.weighted_rms > 0.0
                               ? t_stats.weighted_rms
                               : t_stats.rms_total.to_arcsec();
            if (t_val < cur_rms) { current_state = trial; residuals = trial_res; cur_rms = t_val; return true; }
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

private:
    std::shared_ptr<ResidualCalculator<Frame>> residual_calc_;
    std::shared_ptr<StateTransitionMatrix<Frame>> stm_computer_;
};

} // namespace astdyn::orbit_determination

#endif
