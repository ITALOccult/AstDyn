/**
 * @file ExtendedKalmanFilter.cpp
 * @brief Recursive orbit estimation using EKF.
 */

#include "astdyn/orbit_determination/ExtendedKalmanFilter.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <Eigen/Dense>
#include <cmath>

namespace astdyn::orbit_determination {

ExtendedKalmanFilter::ExtendedKalmanFilter(
    std::shared_ptr<propagation::Propagator> propagator,
    const Settings& settings)
    : propagator_(propagator),
      stm_engine_(propagator),
      settings_(settings) {
    residual_calc_ = std::make_unique<ResidualCalculator<core::GCRF>>(propagator->get_ephemeris());
}

EKFResult ExtendedKalmanFilter::update(
    const physics::CartesianStateTyped<core::GCRF>& prev_state,
    const astdyn::Matrix6d& prev_cov,
    const observations::OpticalObservation& obs)
{
    EKFResult result;
    time::EpochTDB t_obs = time::to_tdb(obs.time);

    // 1. Prediction (Propagate state and covariance)
    auto stm_res = stm_engine_.compute(prev_state, t_obs);
    physics::CartesianStateTyped<core::GCRF> x_pred = stm_res.final_state;
    // Map covariance using ΦPΦᵀ + Q
    Matrix6d P_pred = stm_res.phi * prev_cov * stm_res.phi.transpose() + settings_.process_noise;

    // 2. Innovation and Partials using the ResidualCalculator
    // This handles light-time, aberration, and GCRF frames correctly
    auto res_res = residual_calc_->compute_residuals({obs}, x_pred);
    if (res_res.empty()) return result; // Should not happen with valid input
    
    const auto& res = res_res[0];
    
    // Get Partials (H) at predicted state
    auto obs_pos_opt = residual_calc_->get_observer_position(obs);
    if (!obs_pos_opt) return result;
    
    // We need the partials at t_obs. Since stm_engine_.compute() already integrated to t_obs,
    // we use a zero-duration call at t_obs to get local partials.
    auto partials_res = stm_engine_.compute_with_partials(x_pred, t_obs, *obs_pos_opt);
    Eigen::Matrix<double, 2, 6> H = partials_res.partial_radec;

    // 3. Innovation Vector y = z - h(x)
    // Residuals from ResidualCalculator are already O-C and projected: (O-C)*cos(dec)
    Eigen::Vector2d y;
    y << res.residual_ra.to_rad(), res.residual_dec.to_rad();
    result.innovation = y;

    // 4. Kalman Gain
    Eigen::Matrix2d R;
    R << std::pow(obs.sigma_ra.to_rad(), 2), 0.0,
         0.0, std::pow(obs.sigma_dec.to_rad(), 2);
    
    Eigen::Matrix2d S = H * P_pred * H.transpose() + R;
    Eigen::Matrix<double, 6, 2> K = P_pred * H.transpose() * S.inverse();

    // 5. Update state and covariance
    astdyn::Vector6d dx = K * y;
    
    // Apply correction in SI meters/mps
    result.state = physics::CartesianStateTyped<core::GCRF>::from_si(
        t_obs, 
        x_pred.position.x_si() + dx(0),
        x_pred.position.y_si() + dx(1),
        x_pred.position.z_si() + dx(2),
        x_pred.velocity.x_si() + dx(3),
        x_pred.velocity.y_si() + dx(4),
        x_pred.velocity.z_si() + dx(5),
        x_pred.gm.to_m3_s2()
    );

    result.covariance = (Matrix6d::Identity() - K * H) * P_pred;
    return result;
}

} // namespace astdyn::orbit_determination
