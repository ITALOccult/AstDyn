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
    residual_calc_ = std::make_unique<ResidualCalculator<core::GCRF>>(propagator->get_ephemeris(), propagator);
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
    
    // Normalization: Φ is computed in AU/days units by the integrator for stability.
    // Convert to SI (m, m/s) to match the covariance and observation partials.
    astdyn::Matrix6d phi_si = stm_res.phi;
    const double au_m = physics::Distance::from_au(1.0).to_m();
    const double aud_ms = physics::Velocity::from_au_d(1.0).to_ms();
    
    phi_si.block<3, 3>(0, 3) *= (au_m / aud_ms);
    phi_si.block<3, 3>(3, 0) *= (aud_ms / au_m);
    
    // Map covariance using Φ_si * P * Φ_siᵀ + Q
    astdyn::Matrix6d P_pred = phi_si * prev_cov * phi_si.transpose() + settings_.process_noise;

    // 2. Innovation and Partials using the ResidualCalculator
    auto res_res = residual_calc_->compute_residuals({obs}, x_pred);
    if (res_res.empty()) return result;
    
    const auto& res = res_res[0];
    
    // Get Partials (H) at predicted state
    auto obs_pos_opt = residual_calc_->get_observer_position(obs);
    if (!obs_pos_opt) return result;
    
    // compute_with_partials returns H in [rad/m, rad/(m/s)]
    auto partials_res = stm_engine_.compute_with_partials(x_pred, t_obs, *obs_pos_opt);
    Eigen::Matrix<double, 2, 6> H = partials_res.partial_radec;

    // 3. Innovation Vector y = z - h(x)
    Eigen::Vector2d y;
    y << res.residual_ra.to_rad(), res.residual_dec.to_rad();
    result.innovation = y;

    // 4. Kalman Gain
    // R must be consistent with the innovation y: since y[0] = residual_ra = d_ra * cos(dec)
    // (projected), R[0,0] must be (sigma_ra * cos(dec))^2 = 1 / res.weight_ra.
    Eigen::Matrix2d R;
    R << 1.0 / res.weight_ra,  0.0,
         0.0,                   1.0 / res.weight_dec;
    
    Eigen::Matrix2d S = H * P_pred * H.transpose() + R;
    Eigen::Matrix<double, 6, 2> K = P_pred * H.transpose() * S.ldlt().solve(Eigen::Matrix2d::Identity());

    Eigen::Vector<double, 6> dx = (K * y).head<6>();
    Eigen::Vector<double, 6> state_vec = x_pred.to_eigen_si();
    state_vec += dx;

    result.state = physics::CartesianStateTyped<core::GCRF>::from_si(
        t_obs, 
        state_vec[0], state_vec[1], state_vec[2],
        state_vec[3], state_vec[4], state_vec[5],
        x_pred.gm.to_m3_s2()
    );

    // Joseph form for better numerical stability/symmetry
    auto IKH = astdyn::Matrix6d::Identity() - K * H;
    result.covariance = IKH * P_pred * IKH.transpose() + K * R * K.transpose();
    return result;
}

} // namespace astdyn::orbit_determination
