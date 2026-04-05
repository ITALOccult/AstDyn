#include "astdyn/propagation/MonteCarloUncertainty.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include <Eigen/Dense>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace astdyn::propagation {

MonteCarloUncertainty MonteCarloUncertainty::make(std::shared_ptr<Propagator> propagator) {
    return MonteCarloUncertainty(std::move(propagator));
}

MonteCarloUncertainty::MonteCarloUncertainty(std::shared_ptr<Propagator> propagator)
    : propagator_(std::move(propagator)) {}

Matrix6d MonteCarloUncertainty::compute_cholesky_factor(
    const physics::CovarianceMatrixAU<core::ECLIPJ2000>& covariance)
{
    if (covariance.matrix().determinant() <= 0.0)
        throw std::invalid_argument("covariance is not positive semi-definite");

    Eigen::LLT<Matrix6d> llt(covariance.matrix());
    if (llt.info() != Eigen::Success)
        throw std::invalid_argument("covariance is not positive semi-definite");

    return llt.matrixL();
}

Eigen::VectorXd MonteCarloUncertainty::draw_normal_vector(
    std::mt19937_64& rng,
    std::normal_distribution<double>& standard_normal)
{
    Eigen::VectorXd z(6);
    for (int k = 0; k < 6; ++k)
        z[k] = standard_normal(rng);
    return z;
}

MonteCarloUncertainty::Sample MonteCarloUncertainty::perturb_and_convert(
    const physics::CartesianStateAU<core::ECLIPJ2000>& nominal,
    const Matrix6d& cholesky_factor,
    const Eigen::VectorXd& z)
{
    const Eigen::VectorXd perturbed_vec = nominal.to_eigen() + cholesky_factor * z;
    const auto perturbed_au = physics::CartesianStateAU<core::ECLIPJ2000>::from_eigen(
        perturbed_vec, nominal.epoch, nominal.gm);
    return Sample{
        .elements = cartesian_to_keplerian<core::ECLIPJ2000>(perturbed_au.to_si()),
        .weight   = 1.0
    };
}

std::vector<MonteCarloUncertainty::Sample> MonteCarloUncertainty::generate_samples(
    const physics::CartesianStateAU<core::ECLIPJ2000>& nominal,
    const physics::CovarianceMatrixAU<core::ECLIPJ2000>& covariance,
    size_t n_samples,
    uint64_t seed)
{
    assert(n_samples > 0);

    const Matrix6d cholesky_factor = compute_cholesky_factor(covariance);
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> standard_normal(0.0, 1.0);

    std::vector<Sample> samples;
    samples.reserve(n_samples);

    for (size_t i = 0; i < n_samples; ++i) {
        const Eigen::VectorXd z = draw_normal_vector(rng, standard_normal);
        samples.push_back(perturb_and_convert(nominal, cholesky_factor, z));
    }

    return samples;
}

std::vector<physics::CartesianStateTyped<core::ECLIPJ2000>>
MonteCarloUncertainty::propagate_ensemble(
    const std::vector<Sample>& samples,
    time::EpochTDB target_epoch,
    Propagator& propagator)
{
    std::vector<physics::CartesianStateTyped<core::ECLIPJ2000>> results(samples.size());

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < static_cast<int>(samples.size()); ++i) {
        const auto cart_si = keplerian_to_cartesian<core::ECLIPJ2000>(samples[i].elements);
        results[i] = propagator.propagate_cartesian(cart_si, target_epoch);
    }

    return results;
}

physics::CartesianStateTyped<core::ECLIPJ2000>
MonteCarloUncertainty::compute_ensemble_mean(
    const std::vector<physics::CartesianStateTyped<core::ECLIPJ2000>>& ensemble)
{
    Eigen::VectorXd sum_au = Eigen::VectorXd::Zero(6);
    const double n_samples = static_cast<double>(ensemble.size());

    for (const auto& state : ensemble)
        sum_au += physics::CartesianStateAU<core::ECLIPJ2000>::from_si(state).to_eigen();

    sum_au /= n_samples;
    const auto mean_au = physics::CartesianStateAU<core::ECLIPJ2000>::from_eigen(
        sum_au, ensemble.front().epoch, ensemble.front().gm);
    return mean_au.to_si();
}

physics::CovarianceMatrixAU<core::ECLIPJ2000>
MonteCarloUncertainty::compute_ensemble_covariance(
    const std::vector<physics::CartesianStateTyped<core::ECLIPJ2000>>& ensemble,
    const physics::CartesianStateAU<core::ECLIPJ2000>& mean_au)
{
    Matrix6d covariance_au2 = Matrix6d::Zero();
    const Eigen::VectorXd mean_vec = mean_au.to_eigen();
    const double n_effective = static_cast<double>(ensemble.size()) - 1.0;

    for (const auto& state : ensemble) {
        const Eigen::VectorXd d =
            physics::CartesianStateAU<core::ECLIPJ2000>::from_si(state).to_eigen() - mean_vec;
        covariance_au2 += d * d.transpose();
    }

    covariance_au2 /= n_effective;
    return physics::CovarianceMatrixAU<core::ECLIPJ2000>(covariance_au2);
}

MonteCarloUncertainty::Statistics MonteCarloUncertainty::compute_statistics(
    const std::vector<physics::CartesianStateTyped<core::ECLIPJ2000>>& ensemble)
{
    if (ensemble.size() < 2)
        throw std::invalid_argument("ensemble must contain at least 2 states");

    const auto mean_si = compute_ensemble_mean(ensemble);
    const auto mean_au = physics::CartesianStateAU<core::ECLIPJ2000>::from_si(mean_si);
    const auto covariance_au = compute_ensemble_covariance(ensemble, mean_au);

    return Statistics{
        .mean            = mean_si,
        .covariance      = covariance_au,
        .position_rms_km = covariance_au.sigma_pos_km(),
        .velocity_rms_ms = covariance_au.sigma_vel_ms()
    };
}

} // namespace astdyn::propagation
