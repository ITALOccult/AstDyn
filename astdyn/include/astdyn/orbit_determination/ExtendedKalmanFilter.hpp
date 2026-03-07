/**
 * @file ExtendedKalmanFilter.hpp
 * @brief Recursive orbit estimation using EKF.
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_EKF_HPP
#define ASTDYN_ORBIT_DETERMINATION_EKF_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include <Eigen/Dense>

namespace astdyn::orbit_determination {

/** @brief Results of an EKF update. */
struct EKFResult {
    physics::CartesianStateTyped<core::GCRF> state; ///< Updated state
    astdyn::Matrix6d covariance;                    ///< Updated covariance
    Eigen::Vector2d innovation;                     ///< Residual (O-C)
};

/**
 * @brief Extended Kalman Filter (EKF) for orbit estimation.
 * 
 * Recursive estimator that updates the state and covariance for each new measurement.
 */
class ExtendedKalmanFilter {
public:
    struct Settings {
        astdyn::Matrix6d process_noise;
        double default_ra_sigma;
        double default_dec_sigma;

        Settings() : 
            process_noise(astdyn::Matrix6d::Identity() * 1e-18),
            default_ra_sigma(0.5 * astdyn::constants::ARCSEC_TO_RAD),
            default_dec_sigma(0.5 * astdyn::constants::ARCSEC_TO_RAD) {}
    };

    explicit ExtendedKalmanFilter(
        std::shared_ptr<propagation::Propagator> propagator,
        const Settings& settings = Settings());

    /**
     * @brief Perform EKF prediction and measurement update.
     * 
     * @param prev_state Previous state at t_k-1
     * @param prev_cov Previous covariance at t_k-1
     * @param obs New observation at t_k
     * @return Updated state and covariance.
     */
    EKFResult update(
        const physics::CartesianStateTyped<core::GCRF>& prev_state,
        const astdyn::Matrix6d& prev_cov,
        const observations::OpticalObservation& obs);

private:
    std::shared_ptr<propagation::Propagator> propagator_;
    StateTransitionMatrix stm_engine_;
    Settings settings_;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_EKF_HPP
