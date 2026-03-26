/**
 * @file CovariancePropagator.hpp
 * @brief Propagates 6×6 covariance alongside mean state via numerical STM.
 *
 * Equation:  P_f = STM * P_0 * STM^T + Q
 * STM is computed numerically with symmetric ±delta perturbations.
 */

#ifndef ASTDYN_PROPAGATION_COVARIANCE_PROPAGATOR_HPP
#define ASTDYN_PROPAGATION_COVARIANCE_PROPAGATOR_HPP

#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/core/physics_state_au.hpp"
#include "astdyn/core/Types.hpp"
#include <memory>

namespace astdyn::propagation {

// Position and velocity deltas for numerical STM finite differences
constexpr double DELTA_POS_AU  = 1e-6;  // [AU]
constexpr double DELTA_VEL_AUD = 1e-9;  // [AU/day]

class CovariancePropagator {
public:
    struct Result {
        // Mean state in SI units (public API)
        physics::CartesianStateTyped<core::ECLIPJ2000> mean_state;
        // Covariance in AU: pos-block [AU²], vel-block [(AU/day)²]
        physics::CovarianceMatrixAU<core::ECLIPJ2000> covariance;
    };

    [[nodiscard]] static CovariancePropagator make(std::shared_ptr<Propagator> propagator);

    [[nodiscard]] Result propagate_with_covariance(
        const physics::CartesianStateTyped<core::ECLIPJ2000>& initial_state,
        const physics::CovarianceMatrixAU<core::ECLIPJ2000>& initial_covariance,
        time::EpochTDB target_epoch,
        const physics::CovarianceMatrixAU<core::ECLIPJ2000>& process_noise);

private:
    explicit CovariancePropagator(std::shared_ptr<Propagator> propagator);

    // Assemble STM column-by-column from symmetric perturbations
    [[nodiscard]] Matrix6d compute_stm(
        const physics::CartesianStateAU<core::ECLIPJ2000>& initial_au,
        time::EpochTDB target_epoch);

    // Central-difference column: (y_f(+delta) - y_f(-delta)) / (2*delta)
    [[nodiscard]] Eigen::VectorXd perturbed_finite_difference(
        const physics::CartesianStateAU<core::ECLIPJ2000>& base_au,
        int component,
        double delta,
        time::EpochTDB target_epoch);

    std::shared_ptr<Propagator> propagator_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_COVARIANCE_PROPAGATOR_HPP
