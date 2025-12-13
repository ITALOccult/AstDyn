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
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

namespace astdyn::orbit_determination {

using namespace astdyn::observations;
using namespace astdyn::propagation;
using namespace astdyn::constants;

// WGS84 ellipsoid parameters
static constexpr double WGS84_A = 6378.137;        // Semi-major axis [km]
// WGS84 flattening (reserved for future use with geodetic coordinates)
// static constexpr double WGS84_F = 1.0 / 298.257223563;

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * @brief Convert UTC to TDB (Barycentric Dynamical Time)
 * 
 * TDB = TT + periodic terms
 * TT = TAI + 32.184 s
 * TAI = UTC + ΔAT (leap seconds)
 * 
 * For dates 2017-2025, ΔAT = 37 seconds
 * Periodic terms ≈ 0.001658 sin(g) + 0.000014 sin(2g) [seconds]
 * where g = 357.53 + 0.9856003 * (JD - 2451545.0) [degrees]
 * 
 * @param mjd_utc Modified Julian Date in UTC
 * @return MJD in TDB time scale
 */
static double utc_to_tdb_internal(double mjd_utc) {
    // Leap seconds TAI-UTC (valid 2017-2025)
    // TODO: Load from leap seconds file for dates outside this range
    double delta_at = 37.0; // seconds
    
    // TT = TAI + 32.184 s
    double tt_offset = 32.184; // seconds
    
    // MJD in TT
    double mjd_tt = mjd_utc + (delta_at + tt_offset) / 86400.0;
    
    // TDB periodic correction (Fairhead & Bretagnon 1990)
    // Simplified formula accurate to ~10 microseconds
    double jd_tt = mjd_tt + 2400000.5;
    double T = (jd_tt - 2451545.0) / 36525.0; // Julian centuries from J2000.0
    
    // Mean anomaly of Sun [degrees]
    double g = 357.53 + 0.9856003 * (jd_tt - 2451545.0);
    g = std::fmod(g, 360.0) * (M_PI / 180.0); // Convert to radians
    
    // TDB correction [seconds]
    double tdb_correction = 0.001658 * std::sin(g) 
                          + 0.000014 * std::sin(2.0 * g);
    
    // MJD in TDB
    double mjd_tdb = mjd_tt + tdb_correction / 86400.0;
    
    return mjd_tdb;
}

/**
 * @brief Compute Greenwich Mean Sidereal Time
 * 
 * Uses IAU 1982 formula (Aoki et al. 1982, A&A 105, 359)
 * Accuracy: ~0.1 seconds for dates 1900-2100
 * 
 * @param mjd_ut1 Modified Julian Date in UT1 time scale
 * @return GMST [radians], normalized to [0, 2π)
 */
static double compute_gmst(double mjd_ut1) {
    // Julian centuries from J2000.0 (UT1)
    double T = (mjd_ut1 - MJD2000) / 36525.0;
    
    // GMST at 0h UT1 (IAU 1982 formula)
    // GMST = 24110.54841 + 8640184.812866 T + 0.093104 T² - 6.2e-6 T³ [seconds]
    double gmst_seconds = 24110.54841 
                        + 8640184.812866 * T
                        + 0.093104 * T * T
                        - 6.2e-6 * T * T * T;
    
    // Fraction of day
    double frac_day = std::fmod(mjd_ut1, 1.0);
    
    // Add Earth rotation for fraction of day (1.00273790935 sidereal/solar day ratio)
    gmst_seconds += frac_day * 86400.0 * 1.00273790935;
    
    // Convert to radians and normalize to [0, 2π)
    double gmst_rad = gmst_seconds * (TWO_PI / 86400.0);
    gmst_rad = std::fmod(gmst_rad, TWO_PI);
    if (gmst_rad < 0.0) gmst_rad += TWO_PI;
    
    return gmst_rad;
}

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
    const CartesianElements& state) const {
    
    std::vector<ObservationResidual> residuals;
    residuals.reserve(observations.size());
    
    for (const auto& obs : observations) {
        // Convert observation time from UTC to TDB
        double obs_mjd_tdb = utc_to_tdb_internal(obs.mjd_utc);
        
        // Propagate state to observation epoch if propagator available
        // State is assumed to be in EQUATORIAL J2000 (ICRF)
        CartesianElements state_at_obs = state;
        if (propagator_ && std::abs(obs_mjd_tdb - state.epoch_mjd_tdb) > 1e-6) {
            // Propagate from reference epoch to observation epoch
            // Propagator works in Equatorial J2000
            state_at_obs = propagator_->propagate_cartesian(state, obs_mjd_tdb);
        }
        
        auto residual = compute_residual(obs, state_at_obs);
        if (residual) {
            residuals.push_back(*residual);
        }
    }
    
    return residuals;
}

// Public static wrapper
double ResidualCalculator::utc_to_tdb(double mjd_utc) {
    return utc_to_tdb_internal(mjd_utc);
}

std::optional<ObservationResidual> ResidualCalculator::compute_residual(
    const OpticalObservation& obs,
    const CartesianElements& state) const {
    
    ObservationResidual result;
    result.mjd_utc = obs.mjd_utc;
    result.observatory_code = obs.observatory_code;
    result.outlier = false;
    
    // Get observer position (Heliocentric EQUATORIAL J2000)
    auto observer_pos_opt = get_observer_position(obs);
    if (!observer_pos_opt) {
        return std::nullopt;
    }
    Vector3d observer_pos = *observer_pos_opt;
    
    // Get observer velocity (Heliocentric EQUATORIAL J2000)
    auto observer_vel_opt = get_observer_velocity(obs);
    if (!observer_vel_opt) {
        return std::nullopt;
    }
    Vector3d observer_vel = *observer_vel_opt;
    
    // Object state is already in EQUATORIAL J2000 (from Propagator)
    Vector3d object_pos = state.position;
    Vector3d object_vel = state.velocity;
    
    // Light-time correction (iterate to find retarded position)
    // The observer sees the object where it was tau = distance/c ago
    if (light_time_correction_) {
        double tau = 0.0; // Light travel time [days]
        constexpr int max_iter = 3;
        constexpr double tau_tol = 1e-10; // ~10 microseconds
        
        for (int iter = 0; iter < max_iter; ++iter) {
            Vector3d rho = object_pos - observer_pos;
            double tau_new = rho.norm() / SPEED_OF_LIGHT_AU_PER_DAY;
            
            // Check convergence
            if (iter > 0 && std::abs(tau_new - tau) < tau_tol) {
                break;
            }
            
            tau = tau_new;
            
            // Propagate state backward by tau to get retarded position
            // Simple approximation: object_pos ≈ state.position - state.velocity * tau
            object_pos = state.position - state.velocity * tau;
        }
    }
    
    // Compute topocentric vector (EQUATORIAL J2000)
    Vector3d rho = object_pos - observer_pos;
    double range = rho.norm();
    Vector3d direction = rho.normalized();
    
    // Compute range rate
    Vector3d rho_dot = object_vel - observer_vel;
    double range_rate = rho_dot.dot(direction);
    
    result.range = range;
    result.range_rate = range_rate;
    
    // Convert to RA/Dec (directly from Equatorial vector)
    double computed_ra, computed_dec;
    cartesian_to_radec(direction, computed_ra, computed_dec);
    
    result.computed_ra = computed_ra;
    result.computed_dec = computed_dec;
    
    // Compute residuals O-C
    result.residual_ra = obs.ra - computed_ra;
    result.residual_dec = obs.dec - computed_dec;
    
    // Normalize RA residual by cos(dec) for spherical geometry
    result.residual_ra *= std::cos(obs.dec);
    
    // Normalized residuals
    result.normalized_ra = result.residual_ra / obs.sigma_ra;
    result.normalized_dec = result.residual_dec / obs.sigma_dec;
    
    // Chi-squared
    result.chi_squared = result.normalized_ra * result.normalized_ra +
                        result.normalized_dec * result.normalized_dec;
    
    return result;
}

Vector3d ResidualCalculator::compute_topocentric_direction(
    const Vector3d& heliocentric_pos,
    const Vector3d& observer_pos,
    const Vector3d& observer_vel,
    double& range,
    double& range_rate) const {
    
    // Topocentric position vector
    Vector3d rho = heliocentric_pos - observer_pos;
    range = rho.norm();
    
    // Unit direction vector
    Vector3d direction = rho / range;
    
    // Range rate computation moved to compute_residual() where object velocity is available
    // Set placeholder value here (will be overwritten by caller)
    range_rate = 0.0;
    
    return direction;
}

void ResidualCalculator::cartesian_to_radec(
    const Vector3d& direction,
    double& ra,
    double& dec) const {
    
    double x = direction[0];
    double y = direction[1];
    double z = direction[2];
    
    // Declination: arcsin(z)
    dec = std::asin(z);
    
    // Right ascension: atan2(y, x)
    ra = std::atan2(y, x);
    
    // Normalize RA to [0, 2π)
    if (ra < 0.0) {
        ra += TWO_PI;
    }
}

std::optional<Vector3d> ResidualCalculator::get_observer_position(
    const OpticalObservation& obs) const {
    
    // Convert observation time from UTC to TDB
    double mjd_tdb = utc_to_tdb_internal(obs.mjd_utc);
    double jd_tdb = mjd_tdb + 2400000.5;
    
    // Get Earth position from ephemeris (PlanetaryEphemeris returns EQUATORIAL J2000)
    auto earth_state = ephemeris::PlanetaryEphemeris::getState(
        ephemeris::CelestialBody::EARTH, jd_tdb);
    
    Vector3d earth_center = earth_state.position();
    
    // Get observatory topocentric position
    const auto& obs_db = observations::ObservatoryDatabase::getInstance();
    auto obs_info_opt = obs_db.getObservatory(obs.observatory_code);
    
    if (!obs_info_opt) {
        // Unknown observatory, use geocenter
        return earth_center;
    }
    
    const auto& obs_info = *obs_info_opt;
    
    // Compute observatory position relative to Earth center
    double rho_cos_phi = obs_info.rho_cos_phi;
    double rho_sin_phi = obs_info.rho_sin_phi;
    double longitude = obs_info.longitude;
    
    // Compute Greenwich Mean Sidereal Time
    // GMST is the angle between the Greenwich meridian and the vernal equinox
    double gmst = compute_gmst(obs.mjd_utc);
    
    // Local sidereal time = GMST + longitude
    double lst = gmst + longitude;
    
    // Observatory position in GEOCENTRIC EQUATORIAL coordinates [Earth radii]
    // The Z axis aligns with Earth's rotation axis
    // The X axis points to the vernal equinox
    double cos_lst = std::cos(lst);
    double sin_lst = std::sin(lst);
    
    Vector3d obs_geocentric_equatorial;
    obs_geocentric_equatorial[0] = rho_cos_phi * cos_lst;
    obs_geocentric_equatorial[1] = rho_cos_phi * sin_lst;
    obs_geocentric_equatorial[2] = rho_sin_phi;
    
    // Convert from Earth radii to AU
    double earth_radius_au = WGS84_A / AU_TO_KM;
    obs_geocentric_equatorial *= earth_radius_au;
    
    // NO ROTATION TO ECLIPTIC NEEDED
    // The system assumes everything is in Equatorial J2000
    
    // Observatory heliocentric position (in EQUATORIAL J2000)
    Vector3d observer_pos = earth_center + obs_geocentric_equatorial;
    
    return observer_pos;
}

std::optional<Vector3d> ResidualCalculator::get_observer_velocity(
    const OpticalObservation& obs) const {
    
    // Convert to TDB
    double mjd_tdb = utc_to_tdb_internal(obs.mjd_utc);
    double jd_tdb = mjd_tdb + 2400000.5;
    
    // Get Earth velocity (EQUATORIAL J2000)
    auto earth_state = ephemeris::PlanetaryEphemeris::getState(
        ephemeris::CelestialBody::EARTH, jd_tdb);
    
    Vector3d earth_vel = earth_state.velocity();
    
    // Get observatory position to compute rotation velocity
    const auto& obs_db = observations::ObservatoryDatabase::getInstance();
    auto obs_info_opt = obs_db.getObservatory(obs.observatory_code);
    
    if (!obs_info_opt) {
        return earth_vel;
    }
    
    // Re-calculate geocentric position (code reuse justified for now)
    const auto& obs_info = *obs_info_opt;
    double rho_cos_phi = obs_info.rho_cos_phi;
    double rho_sin_phi = obs_info.rho_sin_phi;
    double longitude = obs_info.longitude;
    double gmst = compute_gmst(obs.mjd_utc);
    double lst = gmst + longitude;
    
    double cos_lst = std::cos(lst);
    double sin_lst = std::sin(lst);
    
    // Geocentric position in Earth radii
    Vector3d obs_geo_radii;
    obs_geo_radii[0] = rho_cos_phi * cos_lst;
    obs_geo_radii[1] = rho_cos_phi * sin_lst;
    obs_geo_radii[2] = rho_sin_phi;
    
    // Convert to AU
    double earth_radius_au = WGS84_A / AU_TO_KM;
    Vector3d obs_geocentric = obs_geo_radii * earth_radius_au;
    
    // Earth rotation angular velocity [rad/day]
    // ω = 2π/T_sid where T_sid ≈ 0.99726958 solar days
    double omega_earth = TWO_PI / 0.99726958;  // rad/day
    
    // Rotation vector in Equatorial frame is simply along Z
    Vector3d omega_vec(0.0, 0.0, omega_earth);
    
    // Velocity due to rotation: v = ω × r
    Vector3d v_rotation = omega_vec.cross(obs_geocentric);
    
    // Total observer velocity (EQUATORIAL J2000)
    Vector3d observer_vel = earth_vel + v_rotation;
    
    return observer_vel;
}

// ============================================================================
// Statistics
// ============================================================================

ResidualStatistics ResidualCalculator::compute_statistics(
    const std::vector<ObservationResidual>& residuals,
    int num_parameters) {
    
    ResidualStatistics stats;
    
    // Count non-outlier observations
    int n_valid = 0;
    for (const auto& r : residuals) {
        if (!r.outlier) n_valid++;
    }
    
    stats.num_observations = residuals.size();
    stats.num_outliers = stats.num_observations - n_valid;
    stats.degrees_of_freedom = 2 * n_valid - num_parameters;
    
    if (n_valid == 0) {
        stats.rms_ra = stats.rms_dec = stats.rms_total = 0.0;
        stats.weighted_rms = 0.0;
        stats.chi_squared = 0.0;
        stats.reduced_chi_squared = 0.0;
        return stats;
    }
    
    // Compute RMS
    double sum_ra2 = 0.0, sum_dec2 = 0.0;
    double sum_chi2 = 0.0;
    double max_ra = 0.0, max_dec = 0.0;
    
    for (const auto& r : residuals) {
        if (r.outlier) continue;
        
        sum_ra2 += r.residual_ra * r.residual_ra;
        sum_dec2 += r.residual_dec * r.residual_dec;
        sum_chi2 += r.chi_squared;
        
        max_ra = std::max(max_ra, std::abs(r.residual_ra));
        max_dec = std::max(max_dec, std::abs(r.residual_dec));
    }
    
    // RMS in arcseconds
    stats.rms_ra = std::sqrt(sum_ra2 / n_valid) * RAD_TO_ARCSEC;
    stats.rms_dec = std::sqrt(sum_dec2 / n_valid) * RAD_TO_ARCSEC;
    stats.rms_total = std::sqrt((sum_ra2 + sum_dec2) / (2.0 * n_valid)) * RAD_TO_ARCSEC;
    
    // Weighted RMS (dimensionless)
    stats.weighted_rms = std::sqrt(sum_chi2 / (2.0 * n_valid));
    
    // Chi-squared
    stats.chi_squared = sum_chi2;
    stats.reduced_chi_squared = (stats.degrees_of_freedom > 0) ?
        stats.chi_squared / stats.degrees_of_freedom : 0.0;
    
    // Max residuals in arcseconds
    stats.max_abs_ra = max_ra * RAD_TO_ARCSEC;
    stats.max_abs_dec = max_dec * RAD_TO_ARCSEC;
    
    return stats;
}

int ResidualCalculator::identify_outliers(
    std::vector<ObservationResidual>& residuals,
    double sigma_threshold) {
    
    int num_outliers = 0;
    
    // Iterative 3-sigma clipping
    bool changed = true;
    while (changed) {
        changed = false;
        
        // Compute current statistics (excluding already marked outliers)
        (void)compute_statistics(residuals, 6); // Just recompute for next iteration
        
        // Mark new outliers
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
