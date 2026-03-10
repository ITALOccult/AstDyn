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
    
    physics::CartesianStateTyped<Frame> running_state = state;
    
    for (const auto& obs : observations) {
        time::EpochTDB obs_time_tdb = astdyn::time::to_tdb(obs.time);
        
        // Sequential propagation
        if (propagator_) {
            double dt = obs_time_tdb.mjd() - running_state.epoch.mjd();
            if (dt > 1e-12) {
                // Future point: continue from current running_state
                running_state = propagator_->propagate_cartesian(running_state, obs_time_tdb);
            } else if (dt < -1e-12) {
                 // Out of order: restart from the original base state
                 running_state = propagator_->propagate_cartesian(state, obs_time_tdb);
            }
            // If |dt| < 1e-12, keep running_state as is (already at the right epoch)
        }
        
        auto residual = compute_residual(obs, running_state);
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
    time::TimeDuration tau = time::TimeDuration::zero();
    if (light_time_correction_) {
        const int max_iter = 5;
        const double tau_tol_s = 1e-8; // 10ns tolerance
        
        auto base_pos_gcrf = object_pos_gcrf; // Store original heliocentric position
        
        for (int iter = 0; iter < max_iter; ++iter) {
            math::Vector3<core::GCRF, physics::Distance> rho_vec = object_pos_gcrf - observer_pos_gcrf;
            double rho_m = rho_vec.norm().to_m();
            double c_ms = physics::SpeedOfLight::to_ms();
            double tau_new_s = rho_m / c_ms;
            
            if (iter > 0 && std::abs(tau_new_s - tau.to_seconds()) < tau_tol_s) break;
            
            tau = time::TimeDuration::from_seconds(tau_new_s);
            // Back-propagate from base position: r(t - tau) = r(t) - v * tau
            auto pos_corr_si = base_pos_gcrf.to_eigen_si() - (object_vel_gcrf.to_eigen_si() * tau.to_seconds());
            object_pos_gcrf = math::Vector3<core::GCRF, physics::Distance>::from_si(pos_corr_si.x(), pos_corr_si.y(), pos_corr_si.z());
        }
    }
    
    // 5. Compute topocentric vector and direction
    math::Vector3<core::GCRF, physics::Distance> rho_vec = object_pos_gcrf - observer_pos_gcrf;
    double range_m = rho_vec.norm().to_m();
    
    // Safety check for range
    if (range_m < 1.0) range_m = 1.0;

    auto direction = math::Vector3<core::GCRF, physics::Distance>::from_si(
        rho_vec.x_si() / range_m, 
        rho_vec.y_si() / range_m, 
        rho_vec.z_si() / range_m);
    
    math::Vector3<core::GCRF, physics::Velocity> rho_dot = object_vel_gcrf - observer_vel_gcrf;
    double range_rate_ms = rho_dot.to_eigen_si().dot(direction.to_eigen_si());
    
    result.range = physics::Distance::from_si(range_m);
    result.range_rate = physics::Velocity::from_si(range_rate_ms);
    
    // 6. Convert to RA/Dec (includes Aberration, Light Bending)
    astrometry::RightAscension computed_ra;
    astrometry::Declination computed_dec;
    cartesian_to_radec(direction, rho_vec, observer_pos_gcrf, observer_vel_gcrf, computed_ra, computed_dec);
    
    result.computed_ra = computed_ra;
    result.computed_dec = computed_dec;
    
    // 7. Compute residuals O-C
    astrometry::Angle d_ra = obs.ra - computed_ra;
    d_ra = d_ra.wrap_pi();

    // Residual for RA is projected using COMPUTED dec (consistent with the design matrix
    // partial d(RA*cos(dec))/dx which is also evaluated at the computed position).
    const double cos_comp_dec = std::cos(computed_dec.to_rad());
    result.residual_ra = d_ra * cos_comp_dec;
    result.residual_dec = obs.dec - computed_dec;

    double sig_ra_rad = obs.sigma_ra.to_rad();
    double sig_dec_rad = obs.sigma_dec.to_rad();
    if (sig_ra_rad <= 0.0) sig_ra_rad = 1e-10;
    if (sig_dec_rad <= 0.0) sig_dec_rad = 1e-10;

    // Weight for projected RA residual: sigma must also be projected so that
    // chi_sq = residual_ra^2 * weight_ra = (d_ra*cos(dec))^2 / (sigma_ra*cos(dec))^2
    //        = d_ra^2 / sigma_ra^2  (correct unprojected chi-squared contribution).
    const double sig_ra_proj = sig_ra_rad * std::max(cos_comp_dec, 1e-6);
    result.weight_ra  = 1.0 / (sig_ra_proj  * sig_ra_proj);
    result.weight_dec = 1.0 / (sig_dec_rad  * sig_dec_rad);

    // Normalised residuals: divide projected residual by projected sigma → ratio = d_ra / sigma_ra
    result.normalized_ra  = result.residual_ra.to_rad()  / sig_ra_proj;
    result.normalized_dec = result.residual_dec.to_rad() / sig_dec_rad;
    
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
    astrometry::RightAscension& ra,
    astrometry::Declination& dec) const {
    
    auto eval_dir = direction_in.to_eigen_si();
    auto eval_obs = observer_pos.to_eigen_si();
    auto eval_vel = observer_vel.to_eigen_si();

    // Relativistic Light Bending (Sun)
    Eigen::Vector3d p_sun = -eval_obs; 
    double d_sun = p_sun.norm();
    auto u_sun = p_sun.normalized();
    double cos_elongation = u_sun.dot(eval_dir);
    
    if (d_sun > 1e6 && cos_elongation < 0.999) { 
        // BUG-13: Use derived constant for 2GM/c^2 instead of magic number
        const double TWO_GM_C2 = 2.0 * physics::GravitationalParameter::sun().to_m3_s2() / 
                                (physics::SpeedOfLight::to_ms() * physics::SpeedOfLight::to_ms());
        auto cross_prod = u_sun.cross(eval_dir);
        if (cross_prod.norm() > 1e-12) { 
            auto delta_dir = (TWO_GM_C2 / d_sun) * (cross_prod.cross(eval_dir)) / (1.0 + cos_elongation);
            eval_dir = (eval_dir + delta_dir).normalized();
        }
    }

    // Stellar Aberration
    // eval_vel is in m/s; C_LIGHT constant is in km/s, so use SpeedOfLight::to_ms() to stay in SI.
    if (aberration_correction_) {
        auto beta = eval_vel / (physics::SpeedOfLight::to_ms());
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
    
    double dec_rad = std::asin(eval_dir.z());
    double ra_rad = std::atan2(eval_dir.y(), eval_dir.x());
    
    dec = astrometry::Declination(astrometry::Angle::from_rad(dec_rad));
    ra = astrometry::RightAscension(astrometry::Angle::from_rad(ra_rad));
}

template <typename Frame>
std::optional<math::Vector3<core::GCRF, physics::Distance>> ResidualCalculator<Frame>::get_observer_position(
    const OpticalObservation& obs) const {
    
    time::EpochTDB t_tdb = time::to_tdb(obs.time);
    auto earth_state = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t_tdb);
    
    // earth_state is heliocentric GCRF
    math::Vector3<core::GCRF, physics::Distance> earth_center_vec = earth_state.position;
    
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
    
    return earth_state.velocity;
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
    
    stats.max_abs_ra = astrometry::Angle::zero();
    stats.max_abs_dec = astrometry::Angle::zero();
    
    for (const auto& res : residuals) {
        if (res.outlier) {
            stats.num_outliers++;
            continue;
        }
        
        accepted_count++;
        double res_ra_rad = res.residual_ra.to_rad();
        double res_dec_rad = res.residual_dec.to_rad();
        
        sum_sq_ra += res_ra_rad * res_ra_rad;
        sum_sq_dec += res_dec_rad * res_dec_rad;
        chi_squared += res.chi_squared;
        
        stats.max_abs_ra = astrometry::Angle::from_rad(std::max(stats.max_abs_ra.to_rad(), std::abs(res_ra_rad)));
        stats.max_abs_dec = astrometry::Angle::from_rad(std::max(stats.max_abs_dec.to_rad(), std::abs(res_dec_rad)));
    }
    
    if (accepted_count > 0) {
        stats.rms_ra = astrometry::Angle::from_rad(std::sqrt(sum_sq_ra / accepted_count));
        stats.rms_dec = astrometry::Angle::from_rad(std::sqrt(sum_sq_dec / accepted_count));
        stats.rms_total = astrometry::Angle::from_rad(std::sqrt((sum_sq_ra + sum_sq_dec) / (2.0 * accepted_count)));
        
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
        stats.rms_ra = astrometry::Angle::zero();
        stats.rms_dec = astrometry::Angle::zero();
        stats.rms_total = astrometry::Angle::zero();
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
