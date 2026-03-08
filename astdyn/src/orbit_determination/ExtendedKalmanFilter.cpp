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
      settings_(settings) {}

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
    Matrix6d P_pred = stm_res.map_covariance(prev_cov) + settings_.process_noise;

    // 2. Observation Partials (H) at observer location
    // Get Earth position (observer)
    auto earth_pos_si = ephemeris::PlanetaryEphemeris::getPosition(ephemeris::CelestialBody::EARTH, t_obs);
    auto R_obs = math::Vector3<core::GCRF, physics::Distance>::from_si(
        earth_pos_si.x_si(),
        earth_pos_si.y_si(),
        earth_pos_si.z_si()
    );

    auto partials_res = stm_engine_.compute_with_partials(prev_state, t_obs, R_obs);
    Eigen::Matrix<double, 2, 6> H = partials_res.partial_radec;

    // 3. Innovation (O - C)
    auto rho_vec = x_pred.position.to_eigen_si() - R_obs.to_eigen_si();
    double rho = rho_vec.norm();
    double ra_pred = std::atan2(rho_vec.y(), rho_vec.x());
    double dec_pred = std::asin(rho_vec.z() / rho);

    double d_ra = obs.ra - ra_pred;
    while (d_ra > constants::PI) d_ra -= constants::TWO_PI;
    while (d_ra < -constants::PI) d_ra += constants::TWO_PI;
    double d_dec = obs.dec - dec_pred;

    Eigen::Vector2d y(d_ra * std::cos(obs.dec), d_dec); // Innovation
    result.innovation = y;

    // 4. Kalman Gain
    Eigen::Matrix2d R;
    R << std::pow(obs.sigma_ra, 2), 0.0,
         0.0, std::pow(obs.sigma_dec, 2);
    
    Eigen::Matrix2d S = H * P_pred * H.transpose() + R;
    Eigen::Matrix<double, 6, 2> K = P_pred * H.transpose() * S.inverse();

    // 5. Update state and covariance
    astdyn::Vector6d dx = K * y;
    astdyn::Vector6d x_vec = x_pred.to_eigen_si() + dx;
    Matrix6d P_new = (Matrix6d::Identity() - K * H) * P_pred;

    result.state = physics::CartesianStateTyped<core::GCRF>::from_si(
        t_obs, 
        x_vec(0), x_vec(1), x_vec(2),
        x_vec(3), x_vec(4), x_vec(5),
        x_pred.gm.to_m3_s2()
    );
    result.covariance = P_new;

    return result;
}

} // namespace astdyn::orbit_determination
