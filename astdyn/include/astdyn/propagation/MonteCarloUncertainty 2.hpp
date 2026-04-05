/**
 * @file MonteCarloUncertainty.hpp
 * @brief Monte Carlo uncertainty propagation via multivariate normal sampling.
 *
 * Samples from N(mean, P) using Cholesky decomposition, propagates the
 * ensemble in parallel (OpenMP), and computes summary statistics.
 */

#ifndef ASTDYN_PROPAGATION_MONTE_CARLO_UNCERTAINTY_HPP
#define ASTDYN_PROPAGATION_MONTE_CARLO_UNCERTAINTY_HPP

#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/core/physics_state_au.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/core/Types.hpp"
#include <cassert>
#include <memory>
#include <random>
#include <vector>

namespace astdyn::propagation {

class MonteCarloUncertainty {
public:
    struct Sample {
        physics::KeplerianStateTyped<core::ECLIPJ2000> elements;
        double weight = 1.0;
    };

    struct Statistics {
        physics::CartesianStateTyped<core::ECLIPJ2000> mean;
        physics::CovarianceMatrixAU<core::ECLIPJ2000> covariance;
        double position_rms_km;  // from CovarianceMatrixAU::sigma_pos_km()
        double velocity_rms_ms;  // from CovarianceMatrixAU::sigma_vel_ms()
    };

    [[nodiscard]] static MonteCarloUncertainty make(std::shared_ptr<Propagator> propagator);

    // Samples from N(nominal, P) using Cholesky; each sample converted to Keplerian
    [[nodiscard]] std::vector<Sample> generate_samples(
        const physics::CartesianStateAU<core::ECLIPJ2000>& nominal,
        const physics::CovarianceMatrixAU<core::ECLIPJ2000>& covariance,
        size_t n_samples = 10000,
        uint64_t seed = 42);

    // Propagates ensemble in parallel if OpenMP is available.
    // NOTE: requires the Propagator to be thread-safe for OpenMP execution.
    [[nodiscard]] std::vector<physics::CartesianStateTyped<core::ECLIPJ2000>> propagate_ensemble(
        const std::vector<Sample>& samples,
        time::EpochTDB target_epoch,
        Propagator& propagator);

    [[nodiscard]] Statistics compute_statistics(
        const std::vector<physics::CartesianStateTyped<core::ECLIPJ2000>>& ensemble);

private:
    explicit MonteCarloUncertainty(std::shared_ptr<Propagator> propagator);

    [[nodiscard]] static Matrix6d compute_cholesky_factor(
        const physics::CovarianceMatrixAU<core::ECLIPJ2000>& covariance);

    [[nodiscard]] static Eigen::VectorXd draw_normal_vector(
        std::mt19937_64& rng,
        std::normal_distribution<double>& standard_normal);

    [[nodiscard]] static Sample perturb_and_convert(
        const physics::CartesianStateAU<core::ECLIPJ2000>& nominal,
        const Matrix6d& cholesky_factor,
        const Eigen::VectorXd& z);

    [[nodiscard]] static physics::CartesianStateTyped<core::ECLIPJ2000> compute_ensemble_mean(
        const std::vector<physics::CartesianStateTyped<core::ECLIPJ2000>>& ensemble);

    [[nodiscard]] static physics::CovarianceMatrixAU<core::ECLIPJ2000> compute_ensemble_covariance(
        const std::vector<physics::CartesianStateTyped<core::ECLIPJ2000>>& ensemble,
        const physics::CartesianStateAU<core::ECLIPJ2000>& mean_au);

    std::shared_ptr<Propagator> propagator_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_MONTE_CARLO_UNCERTAINTY_HPP
