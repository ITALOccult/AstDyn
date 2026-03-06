/**
 * @file Residuals.cpp
 * @brief Implementation of observation residual calculations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/observations/ObservatoryDatabase.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "src/utils/time_conversions.hpp"
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <iostream>

namespace astdyn::orbit_determination {

using namespace astdyn::observations;
using namespace astdyn::propagation;
using namespace astdyn::constants;

// WGS84 ellipsoid parameters
static constexpr double WGS84_A = 6378.137;        // Semi-major axis [km]

// ============================================================================
// ResidualCalculator Implementation
// ============================================================================

ResidualCalculator::ResidualCalculator(
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris,
    std::shared_ptr<astdyn::propagation::Propagator> propagator)
    : ephemeris_(ephemeris),
      propagator_(propagator) {
}

std::vector<ObservationResidual> ResidualCalculator::compute_residuals(
    const std::vector<OpticalObservation>& observations,
    const physics::CartesianStateTyped<core::GCRF>& state) const {
    
    std::vector<ObservationResidual> residuals;
    residuals.reserve(observations.size());
    
    for (const auto& obs : observations) {
        // Convert observation time from UTC to TDB
        time::EpochTDB obs_time_tdb = utc_to_tdb(obs.time);
        
        // Propagate state to observation epoch if propagator available
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

time::EpochTDB ResidualCalculator::utc_to_tdb(utils::Instant t_utc) {
    // Simplified conversion UTC -> TDB
    // UTC -> TAI (+37s) -> TT (+32.184s) -> TDB (periodic terms)
    
    double mjd_utc = t_utc.mjd.value;
    double delta_at = 37.0; // seconds (2017-2025)
    double tt_offset = 32.184; // seconds
    
    double mjd_tt = mjd_utc + (delta_at + tt_offset) / 86400.0;
    
    // TDB periodic correction
    double jd_tt = mjd_tt + 2400000.5;
    double g = 357.53 + 0.9856003 * (jd_tt - 2451545.0);
    g = std::fmod(g, 360.0) * astdyn::constants::DEG_TO_RAD;
    
    double mjd_tdb = mjd_tt + 0.001658 * std::sin(g + 0.0167 * std::sin(g)) / 86400.0;
    
    return time::EpochTDB::from_mjd(mjd_tdb);
}

std::optional<ObservationResidual> ResidualCalculator::compute_residual(
    const OpticalObservation& obs,
    const physics::CartesianStateTyped<core::GCRF>& state) const {
    
    ObservationResidual result;
    result.time = obs.time;
    result.observatory_code = obs.observatory_code;
    result.outlier = false;
    
    // Get observer position (Heliocentric GCRF Meter)
    auto observer_pos_opt = get_observer_position(obs);
    if (!observer_pos_opt) {
        return std::nullopt;
    }
    types::Vector3<core::GCRF, core::Meter> observer_pos = *observer_pos_opt;
    
    // Get observer velocity (Heliocentric GCRF Meter/Second)
    auto observer_vel_opt = get_observer_velocity(obs);
    if (!observer_vel_opt) {
        return std::nullopt;
    }
    types::Vector3<core::GCRF, core::Meter> observer_vel = *observer_vel_opt;
    
    // Object state
    types::Vector3<core::GCRF, core::Meter> object_pos(state.position.to_eigen_si());
    types::Vector3<core::GCRF, core::Meter> object_vel(state.velocity.to_eigen_si());
    
    // Light-time correction
    if (light_time_correction_) {
        double tau = 0.0; // Light travel time [days]
        constexpr int max_iter = 3;
        constexpr double tau_tol = 1e-10; 
        
        for (int iter = 0; iter < max_iter; ++iter) {
            types::Vector3<core::GCRF, core::Meter> rho_vec = object_pos - observer_pos;
            // rho.norm() is in meters. SPEED_OF_LIGHT is in m/s.
            double tau_new_sec = rho_vec.norm() / constants::C_LIGHT * 1000.0;
            double tau_new_days = tau_new_sec / 86400.0;
            
            if (iter > 0 && std::abs(tau_new_days - tau) < tau_tol) {
                break;
            }
            
            tau = tau_new_days;
            // object_pos ≈ state.position - state.velocity * tau_sec
            object_pos = types::Vector3<core::GCRF, core::Meter>(state.position.to_eigen_si() - (state.velocity.to_eigen_si() * (tau * 86400.0)));
        }
    }
    
    // Compute topocentric vector
    double range, range_rate;
    auto direction = compute_topocentric_direction(object_pos, observer_pos, observer_vel, range, range_rate);
    
    // Re-calculate range rate properly with object velocity
    types::Vector3<core::GCRF, core::Meter> rho_vec = object_pos - observer_pos;
    types::Vector3<core::GCRF, core::Meter> rho_dot = object_vel - observer_vel;
    range_rate = rho_dot.dot(direction);
    
    result.range = range;
    result.range_rate = range_rate;
    
    // Convert to RA/Dec
    double computed_ra_rad, computed_dec_rad;
    cartesian_to_radec(direction, rho_vec, observer_pos, observer_vel, computed_ra_rad, computed_dec_rad);
    
    result.computed_ra = computed_ra_rad;
    result.computed_dec = computed_dec_rad;
    
    // Compute residuals O-C (Radians)
    double d_ra = obs.ra - computed_ra_rad;
    while (d_ra > PI) d_ra -= TWO_PI;
    while (d_ra < -PI) d_ra += TWO_PI;
    
    result.residual_ra = d_ra * std::cos(obs.dec);
    result.residual_dec = obs.dec - computed_dec_rad;
    
    // Normalized residuals
    result.normalized_ra = result.residual_ra / obs.sigma_ra;
    result.normalized_dec = result.residual_dec / obs.sigma_dec;
    
    // Weights
    double sig_ra = (obs.sigma_ra > 0.0) ? obs.sigma_ra : 1e-5;
    double sig_dec = (obs.sigma_dec > 0.0) ? obs.sigma_dec : 1e-5;
    result.weight_ra = 1.0 / (sig_ra * sig_ra);
    result.weight_dec = 1.0 / (sig_dec * sig_dec);
    
    result.chi_squared = result.normalized_ra * result.normalized_ra +
                        result.normalized_dec * result.normalized_dec;
    
    return result;
}

types::Vector3<core::GCRF, core::Meter> ResidualCalculator::compute_topocentric_direction(
    const types::Vector3<core::GCRF, core::Meter>& heliocentric_pos,
    const types::Vector3<core::GCRF, core::Meter>& observer_pos,
    const types::Vector3<core::GCRF, core::Meter>& observer_vel,
    double& range,
    double& range_rate) const {
    
    types::Vector3<core::GCRF, core::Meter> rho = heliocentric_pos - observer_pos;
    range = rho.norm();
    auto direction = rho.normalized();
    range_rate = 0.0; // Placeholder
    return direction;
}

void ResidualCalculator::cartesian_to_radec(
    const types::Vector3<core::GCRF, core::Meter>& direction_in,
    const types::Vector3<core::GCRF, core::Meter>& rho_vec,
    const types::Vector3<core::GCRF, core::Meter>& observer_pos,
    const types::Vector3<core::GCRF, core::Meter>& observer_vel,
    double& ra_rad,
    double& dec_rad) const {
    
    auto direction = direction_in;

    // Relativistic Light Bending (Sun)
    types::Vector3<core::GCRF, core::Meter> p_sun = -observer_pos; 
    double d_sun = p_sun.norm();
    auto u_sun = p_sun.normalized();
    double cos_elongation = u_sun.dot(direction);
    
    if (d_sun > 1e6 && cos_elongation < 0.999) { 
        // 2GM/c^2 in meters ≈ 2953 m
        constexpr double TWO_GM_C2 = 2953.36; 
        auto cross_prod = u_sun.cross(direction);
        if (cross_prod.norm() > 1e-12) { 
            auto delta_dir = (TWO_GM_C2 / d_sun) * (cross_prod.cross(direction)) / (1.0 + cos_elongation);
            direction = (direction + delta_dir).normalized();
        }
    }

    // Stellar Aberration
    if (aberration_correction_) {
        auto beta = observer_vel / constants::C_LIGHT * 1000.0;
        double beta_sq = beta.squaredNorm();
        if (beta_sq < 1.0) {
            double gamma = 1.0 / std::sqrt(1.0 - beta_sq);
            double beta_dot_u = beta.dot(direction);
            double factor = gamma / (1.0 + gamma);
            auto num = (1.0/gamma) * direction + beta + (factor * beta_dot_u) * beta;
            double den = 1.0 + beta_dot_u;
            direction = (num / den).normalized();
        }
    }
    
    dec_rad = std::asin(direction.z);
    ra_rad = std::atan2(direction.y, direction.x);
    if (ra_rad < 0.0) ra_rad += TWO_PI;
}

std::optional<types::Vector3<core::GCRF, core::Meter>> ResidualCalculator::get_observer_position(
    const OpticalObservation& obs) const {
    
    double jd_tdb = utc_to_tdb(obs.time).jd();
    
    auto earth_state = ephemeris::PlanetaryEphemeris::getState(
        ephemeris::CelestialBody::EARTH, utils::Instant::from_tt(utils::ModifiedJulianDate(jd_tdb - 2400000.5)));
    auto sun_state = ephemeris::PlanetaryEphemeris::getState(
        ephemeris::CelestialBody::SUN, utils::Instant::from_tt(utils::ModifiedJulianDate(jd_tdb - 2400000.5)));
    
    types::Vector3<core::GCRF, core::Meter> earth_center_vec(earth_state.position() - sun_state.position());
    
    const auto& obs_db = observations::ObservatoryDatabase::getInstance();
    auto obs_info_opt = obs_db.getObservatory(obs.observatory_code);
    if (!obs_info_opt) return earth_center_vec;
    
    // Simplified: getPositionGCRF handles ITRF -> GCRF
    return earth_center_vec + obs_info_opt->getPositionGCRF(obs.time);
}

std::optional<types::Vector3<core::GCRF, core::Meter>> ResidualCalculator::get_observer_velocity(
    const OpticalObservation& obs) const {
    
    double jd_tdb = utc_to_tdb(obs.time).jd();
    
    auto earth_state = ephemeris::PlanetaryEphemeris::getState(
        ephemeris::CelestialBody::EARTH, utils::Instant::from_tt(utils::ModifiedJulianDate(jd_tdb - 2400000.5)));
    
    // For now, assume topocentric velocity is negligible or handled by ITRF conversion (not yet in Observatory)
    // In a full implementation, we'd add Earth rotation velocity.
    return types::Vector3<core::GCRF, core::Meter>(earth_state.velocity());
}

ResidualStatistics ResidualCalculator::compute_statistics(
    const std::vector<ObservationResidual>& residuals,
    int num_parameters) {
    
    ResidualStatistics stats;
    int n_valid = 0;
    for (const auto& r : residuals) if (!r.outlier) n_valid++;
    
    stats.num_observations = residuals.size();
    stats.num_outliers = stats.num_observations - n_valid;
    stats.degrees_of_freedom = 2 * n_valid - num_parameters;
    
    if (n_valid == 0) return stats;
    
    double sum_ra2 = 0.0, sum_dec2 = 0.0, sum_chi2 = 0.0;
    double max_ra = 0.0, max_dec = 0.0;
    
    for (const auto& r : residuals) {
        if (r.outlier) continue;
        sum_ra2 += r.residual_ra * r.residual_ra;
        sum_dec2 += r.residual_dec * r.residual_dec;
        sum_chi2 += r.chi_squared;
        max_ra = std::max(max_ra, std::abs(r.residual_ra));
        max_dec = std::max(max_dec, std::abs(r.residual_dec));
    }
    
    stats.rms_ra = std::sqrt(sum_ra2 / n_valid) * RAD_TO_ARCSEC;
    stats.rms_dec = std::sqrt(sum_dec2 / n_valid) * RAD_TO_ARCSEC;
    stats.rms_total = std::sqrt((sum_ra2 + sum_dec2) / (2.0 * n_valid)) * RAD_TO_ARCSEC;
    stats.weighted_rms = std::sqrt(sum_chi2 / (2.0 * n_valid));
    stats.chi_squared = sum_chi2;
    stats.reduced_chi_squared = (stats.degrees_of_freedom > 0) ? stats.chi_squared / stats.degrees_of_freedom : 0.0;
    stats.max_abs_ra = max_ra * RAD_TO_ARCSEC;
    stats.max_abs_dec = max_dec * RAD_TO_ARCSEC;
    
    return stats;
}

int ResidualCalculator::identify_outliers(
    std::vector<ObservationResidual>& residuals,
    double sigma_threshold) {
    
    int num_outliers = 0;
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto& r : residuals) {
            if (r.outlier) continue;
            if (r.is_outlier(sigma_threshold)) {
                r.outlier = true;
                changed = true;
                num_outliers++;
            }
        }
    }
    return num_outliers;
}

} // namespace astdyn::orbit_determination
