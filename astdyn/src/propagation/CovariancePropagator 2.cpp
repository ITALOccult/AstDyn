#include "astdyn/propagation/CovariancePropagator.hpp"

namespace astdyn::propagation {

CovariancePropagator CovariancePropagator::make(std::shared_ptr<Propagator> propagator) {
    return CovariancePropagator(std::move(propagator));
}

CovariancePropagator::CovariancePropagator(std::shared_ptr<Propagator> propagator)
    : propagator_(std::move(propagator)) {}

Eigen::VectorXd CovariancePropagator::perturbed_finite_difference(
    const physics::CartesianStateAU<core::ECLIPJ2000>& base_au,
    int component,
    double delta,
    time::EpochTDB target_epoch)
{
    Eigen::VectorXd y0 = base_au.to_eigen();
    Eigen::VectorXd y_plus  = y0;
    Eigen::VectorXd y_minus = y0;
    y_plus[component]  += delta;
    y_minus[component] -= delta;

    const double t0 = base_au.epoch.mjd();
    const double tf = target_epoch.mjd();
    const auto yf_plus  = propagator_->integrate_raw_au(y_plus,  t0, tf);
    const auto yf_minus = propagator_->integrate_raw_au(y_minus, t0, tf);

    return (yf_plus - yf_minus) / (2.0 * delta);
}

Matrix6d CovariancePropagator::compute_stm(
    const physics::CartesianStateAU<core::ECLIPJ2000>& initial_au,
    time::EpochTDB target_epoch)
{
    Matrix6d stm;

    for (int j = 0; j < 3; ++j)
        stm.col(j) = perturbed_finite_difference(initial_au, j, DELTA_POS_AU, target_epoch);

    for (int j = 3; j < 6; ++j)
        stm.col(j) = perturbed_finite_difference(initial_au, j, DELTA_VEL_AUD, target_epoch);

    return stm;
}

CovariancePropagator::Result CovariancePropagator::propagate_with_covariance(
    const physics::CartesianStateTyped<core::ECLIPJ2000>& initial_state,
    const physics::CovarianceMatrixAU<core::ECLIPJ2000>& initial_covariance,
    time::EpochTDB target_epoch,
    const physics::CovarianceMatrixAU<core::ECLIPJ2000>& process_noise)
{
    const auto initial_au = physics::CartesianStateAU<core::ECLIPJ2000>::from_si(initial_state);
    const Matrix6d stm = compute_stm(initial_au, target_epoch);

    const Eigen::VectorXd yf_au = propagator_->integrate_raw_au(
        initial_au.to_eigen(), initial_state.epoch.mjd(), target_epoch.mjd());
    const auto final_au = physics::CartesianStateAU<core::ECLIPJ2000>::from_eigen(
        yf_au, target_epoch, initial_au.gm);

    const Matrix6d covariance_au2 = stm * initial_covariance.matrix() * stm.transpose()
                                  + process_noise.matrix();

    return Result{
        .mean_state = final_au.to_si(),
        .covariance = physics::CovarianceMatrixAU<core::ECLIPJ2000>(covariance_au2)
    };
}

} // namespace astdyn::propagation
