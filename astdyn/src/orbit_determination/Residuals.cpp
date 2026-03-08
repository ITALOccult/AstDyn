/**
 * @file Residuals.cpp
 * @brief Implementation of observation residual calculations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/observations/ObservatoryDatabase.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <iostream>

namespace astdyn::orbit_determination {

using namespace astdyn::observations;
using namespace astdyn::propagation;
using namespace astdyn::constants;

// ============================================================================
// ResidualCalculator Implementation
// ============================================================================

template <typename Frame>
ResidualCalculator<Frame>::ResidualCalculator(
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris,
    std::shared_ptr<astdyn::propagation::Propagator> propagator)
    : ephemeris_(ephemeris),
      propagator_(propagator) {
}

template <typename Frame>
std::vector<ObservationResidual> ResidualCalculator<Frame>::compute_residuals(
    const std::vector<OpticalObservation>& observations,
    const physics::CartesianStateTyped<Frame>& state) const {
    
    std::vector<ObservationResidual> residuals;
    residuals.reserve(observations.size());
    
    for (const auto& obs : observations) {
        time::EpochTDB obs_time_tdb = astdyn::time::to_tdb(obs.time);
        
        auto state_at_obs = state;
        if (propagator_ && std::abs(obs_time_tdb.mjd() - state.epoch.mjd()) > 1e-10) {
            state_at_obs = propagator_->propagate_cartesian(state, obs_time_tdb);
        }
        
        auto residual = compute_residual(obs, state_at_obs);
        if (residual) {
            residuals.push_back(*residual);
        }
    }
    
    return residuals;
}

template <typename Frame>
std::optional<ObservationResidual> ResidualCalculator<Frame>::compute_residual(
    const OpticalObservation& obs,
    const physics::CartesianStateTyped<Frame>& state) const {
    
    ObservationResidual result;
    result.time = obs.time;
    result.observatory_code = obs.observatory_code;
    result.outlier = false;
    
    // 1. Get observer position (Always GCRF for RA/Dec matching)
    auto observer_pos_opt = get_observer_position(obs);
    if (!observer_pos_opt) return std::nullopt;
    math::Vector3<core::GCRF, physics::Distance> observer_pos_gcrf = *observer_pos_opt;
    
    auto observer_vel_opt = get_observer_velocity(obs);
    if (!observer_vel_opt) return std::nullopt;
    math::Vector3<core::GCRF, physics::Velocity> observer_vel_gcrf = *observer_vel_opt;
    
    // 2. Object state in integration frame
    math::Vector3<Frame, physics::Distance> object_pos_frame = state.position;
    math::Vector3<Frame, physics::Velocity> object_vel_frame = state.velocity;
    
    // 3. Transform object to GCRF for light-time and matching
    auto object_pos_gcrf = coordinates::ReferenceFrame::transform_pos<Frame, core::GCRF>(object_pos_frame, state.epoch);
    auto object_vel_gcrf = coordinates::ReferenceFrame::transform_vel<Frame, core::GCRF>(object_pos_frame, object_vel_frame, state.epoch);
    
    // 4. Light-time correction
    if (light_time_correction_) {
        double tau = 0.0;
        constexpr int max_iter = 3;
        constexpr double tau_tol = 1e-10; 
        
        for (int iter = 0; iter < max_iter; ++iter) {
            math::Vector3<core::GCRF, physics::Distance> rho_vec = object_pos_gcrf - observer_pos_gcrf;
            double tau_new_days = rho_vec.norm().to_m() / (constants::C_LIGHT * 86400.0);
            
            if (iter > 0 && std::abs(tau_new_days - tau) < tau_tol) break;
            
            tau = tau_new_days;
            // Linear approximation for light-time (back-propagation)
            auto pos_corr_si = object_pos_gcrf.to_eigen_si() - (object_vel_gcrf.to_eigen_si() * (tau * 86400.0));
            object_pos_gcrf = math::Vector3<core::GCRF, physics::Distance>::from_si(pos_corr_si.x(), pos_corr_si.y(), pos_corr_si.z());
        }
    }
    
    // 5. Compute topocentric vector and direction
    math::Vector3<core::GCRF, physics::Distance> rho_vec = object_pos_gcrf - observer_pos_gcrf;
    double range = rho_vec.norm().to_m();
    auto direction = math::Vector3<core::GCRF, physics::Distance>::from_si(rho_vec.x_si() / range, rho_vec.y_si() / range, rho_vec.z_si() / range);
    
    math::Vector3<core::GCRF, physics::Velocity> rho_dot = object_vel_gcrf - observer_vel_gcrf;
    double range_rate = rho_dot.to_eigen_si().dot(direction.to_eigen_si());
    
    result.range = range;
    result.range_rate = range_rate;
    
    // 6. Convert to RA/Dec (includes Aberration, Light Bending)
    double computed_ra_rad, computed_dec_rad;
    cartesian_to_radec(direction, rho_vec, observer_pos_gcrf, observer_vel_gcrf, computed_ra_rad, computed_dec_rad);
    
    result.computed_ra = computed_ra_rad;
    result.computed_dec = computed_dec_rad;
    
    // 7. Compute residuals O-C (Radians)
    double d_ra = obs.ra - computed_ra_rad;
    while (d_ra > constants::PI) d_ra -= constants::TWO_PI;
    while (d_ra < -constants::PI) d_ra += constants::TWO_PI;
    
    result.residual_ra = d_ra * std::cos(obs.dec);
    result.residual_dec = obs.dec - computed_dec_rad;
    
    result.normalized_ra = result.residual_ra / obs.sigma_ra;
    result.normalized_dec = result.residual_dec / obs.sigma_dec;
    
    double sig_ra = (obs.sigma_ra > 0.0) ? obs.sigma_ra : 1e-5;
    double sig_dec = (obs.sigma_dec > 0.0) ? obs.sigma_dec : 1e-5;
    result.weight_ra = 1.0 / (sig_ra * sig_ra);
    result.weight_dec = 1.0 / (sig_dec * sig_dec);
    
    result.chi_squared = result.normalized_ra * result.normalized_ra +
                        result.normalized_dec * result.normalized_dec;
    
    return result;
}

template <typename Frame>
void ResidualCalculator<Frame>::cartesian_to_radec(
    const math::Vector3<core::GCRF, physics::Distance>& direction_in,
    const math::Vector3<core::GCRF, physics::Distance>& rho_vec,
    const math::Vector3<core::GCRF, physics::Distance>& observer_pos,
    const math::Vector3<core::GCRF, physics::Velocity>& observer_vel,
    double& ra_rad,
    double& dec_rad) const {
    
    auto eval_dir = direction_in.to_eigen_si();
    auto eval_obs = observer_pos.to_eigen_si();
    auto eval_vel = observer_vel.to_eigen_si();

    // Relativistic Light Bending (Sun)
    Eigen::Vector3d p_sun = -eval_obs; 
    double d_sun = p_sun.norm();
    auto u_sun = p_sun.normalized();
    double cos_elongation = u_sun.dot(eval_dir);
    
    if (d_sun > 1e6 && cos_elongation < 0.999) { 
        constexpr double TWO_GM_C2 = 2953.36; 
        auto cross_prod = u_sun.cross(eval_dir);
        if (cross_prod.norm() > 1e-12) { 
            auto delta_dir = (TWO_GM_C2 / d_sun) * (cross_prod.cross(eval_dir)) / (1.0 + cos_elongation);
            eval_dir = (eval_dir + delta_dir).normalized();
        }
    }

    // Stellar Aberration
    if (aberration_correction_) {
        auto beta = eval_vel / (constants::C_LIGHT);
        double beta_sq = beta.squaredNorm();
        if (beta_sq < 1.0) {
            double gamma = 1.0 / std::sqrt(1.0 - beta_sq);
            double beta_dot_u = beta.dot(eval_dir);
            double factor = gamma / (1.0 + gamma);
            auto num = (1.0/gamma) * eval_dir + beta + (factor * beta_dot_u) * beta;
            double den = 1.0 + beta_dot_u;
            eval_dir = (num / den).normalized();
        }
    }
    
    dec_rad = std::asin(eval_dir.z());
    ra_rad = std::atan2(eval_dir.y(), eval_dir.x());
    if (ra_rad < 0.0) ra_rad += constants::TWO_PI;
}

template <typename Frame>
std::optional<math::Vector3<core::GCRF, physics::Distance>> ResidualCalculator<Frame>::get_observer_position(
    const OpticalObservation& obs) const {
    
    time::EpochTDB t_tdb = time::to_tdb(obs.time);
    auto earth_state = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t_tdb);
    auto sun_state = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::SUN, t_tdb);
    
    math::Vector3<core::GCRF, physics::Distance> earth_center_vec = math::Vector3<core::GCRF, physics::Distance>::from_si(
        earth_state.position().x() - sun_state.position().x(),
        earth_state.position().y() - sun_state.position().y(),
        earth_state.position().z() - sun_state.position().z()
    );
    
    const auto& obs_db = observations::ObservatoryDatabase::getInstance();
    auto obs_info_opt = obs_db.getObservatory(obs.observatory_code);
    if (!obs_info_opt) return earth_center_vec;
    
    auto obs_gcrf = obs_info_opt->getPositionGCRF(obs.time);
    return earth_center_vec + math::Vector3<core::GCRF, physics::Distance>::from_si(obs_gcrf.x_si(), obs_gcrf.y_si(), obs_gcrf.z_si());
}

template <typename Frame>
std::optional<math::Vector3<core::GCRF, physics::Velocity>> ResidualCalculator<Frame>::get_observer_velocity(
    const OpticalObservation& obs) const {
    
    time::EpochTDB t_tdb = time::to_tdb(obs.time);
    auto earth_state = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t_tdb);
    
    return math::Vector3<core::GCRF, physics::Velocity>::from_si(
        earth_state.velocity().x(),
        earth_state.velocity().y(),
        earth_state.velocity().z()
    );
}

template <typename Frame>
ResidualStatistics ResidualCalculator<Frame>::compute_statistics(
    const std::vector<ObservationResidual>& residuals,
    int num_parameters) {
    
    ResidualStatistics stats;
    stats.num_observations = residuals.size();
    stats.num_outliers = 0;
    
    double sum_sq_ra = 0.0;
    double sum_sq_dec = 0.0;
    double chi_squared = 0.0;
    int accepted_count = 0;
    
    stats.max_abs_ra = 0.0;
    stats.max_abs_dec = 0.0;
    
    for (const auto& res : residuals) {
        if (res.outlier) {
            stats.num_outliers++;
            continue;
        }
        
        accepted_count++;
        double res_ra_arcsec = res.residual_ra * constants::RAD_TO_ARCSEC;
        double res_dec_arcsec = res.residual_dec * constants::RAD_TO_ARCSEC;
        
        sum_sq_ra += res_ra_arcsec * res_ra_arcsec;
        sum_sq_dec += res_dec_arcsec * res_dec_arcsec;
        chi_squared += res.chi_squared;
        
        stats.max_abs_ra = std::max(stats.max_abs_ra, std::abs(res_ra_arcsec));
        stats.max_abs_dec = std::max(stats.max_abs_dec, std::abs(res_dec_arcsec));
    }
    
    if (accepted_count > 0) {
        stats.rms_ra = std::sqrt(sum_sq_ra / accepted_count);
        stats.rms_dec = std::sqrt(sum_sq_dec / accepted_count);
        stats.rms_total = std::sqrt((sum_sq_ra + sum_sq_dec) / (2.0 * accepted_count));
        
        stats.chi_squared = chi_squared;
        stats.degrees_of_freedom = 2 * accepted_count - num_parameters;
        
        if (stats.degrees_of_freedom > 0) {
            stats.reduced_chi_squared = chi_squared / stats.degrees_of_freedom;
            stats.weighted_rms = std::sqrt(stats.reduced_chi_squared);
        } else {
            stats.reduced_chi_squared = 0.0;
            stats.weighted_rms = 0.0;
        }
    } else {
        stats.rms_ra = 0.0;
        stats.rms_dec = 0.0;
        stats.rms_total = 0.0;
        stats.chi_squared = 0.0;
        stats.degrees_of_freedom = 0;
        stats.reduced_chi_squared = 0.0;
        stats.weighted_rms = 0.0;
    }
    
    return stats;
}

template <typename Frame>
int ResidualCalculator<Frame>::identify_outliers(
    std::vector<ObservationResidual>& residuals,
    double sigma_threshold) {
    
    int outliers_found = 0;
    for (auto& res : residuals) {
        if (!res.outlier) {
            if (res.is_outlier(sigma_threshold)) {
                res.outlier = true;
                outliers_found++;
            }
        }
    }
    return outliers_found;
}

// Explicit instantiations
template class ResidualCalculator<core::GCRF>;
template class ResidualCalculator<core::ECLIPJ2000>;

} // namespace astdyn::orbit_determination
