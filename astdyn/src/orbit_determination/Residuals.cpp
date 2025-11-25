/**
 * @file Residuals.cpp
 * @brief Implementation of observation residual calculations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/observations/ObservatoryDatabase.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

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
static double utc_to_tdb(double mjd_utc) {
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
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris)
    : ephemeris_(ephemeris) {
}

std::vector<ObservationResidual> ResidualCalculator::compute_residuals(
    const std::vector<OpticalObservation>& observations,
    const CartesianElements& state) const {
    
    std::vector<ObservationResidual> residuals;
    residuals.reserve(observations.size());
    
    for (const auto& obs : observations) {
        auto residual = compute_residual(obs, state);
        if (residual) {
            residuals.push_back(*residual);
        }
    }
    
    return residuals;
}

std::optional<ObservationResidual> ResidualCalculator::compute_residual(
    const OpticalObservation& obs,
    const CartesianElements& state) const {
    
    ObservationResidual result;
    result.mjd_utc = obs.mjd_utc;
    result.observatory_code = obs.observatory_code;
    result.outlier = false;
    
    // Get observer position (heliocentric)
    auto observer_pos_opt = get_observer_position(obs);
    if (!observer_pos_opt) {
        return std::nullopt;
    }
    Vector3d observer_pos = *observer_pos_opt;
    
    // Get observer velocity (for aberration)
    auto observer_vel_opt = get_observer_velocity(obs);
    if (!observer_vel_opt) {
        return std::nullopt;
    }
    Vector3d observer_vel = *observer_vel_opt;
    
    // Object position at observation time (heliocentric)
    Vector3d object_pos = state.position;
    
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
            // This is valid for short light-time and small accelerations
            // For better accuracy, use full numerical integration
            object_pos = state.position - state.velocity * tau;
            
            // Note: For asteroids near Earth, tau ~ 10 minutes, velocity correction
            // is ~0.001 AU. For more distant objects, would need full propagation.
        }
    }
    
    // Compute topocentric direction
    double range, range_rate;
    Vector3d direction = compute_topocentric_direction(
        object_pos, observer_pos, observer_vel, range, range_rate);
    
    // Compute range rate properly using object velocity
    Vector3d rho = object_pos - observer_pos;
    Vector3d rho_dot = state.velocity - observer_vel;
    range_rate = rho_dot.dot(rho.normalized());
    
    result.range = range;
    result.range_rate = range_rate;
    
    // Apply aberration correction
    if (aberration_correction_) {
        // Annual aberration: shift direction by observer velocity
        // Δθ ≈ (v/c) for small angles
        Vector3d v_over_c = observer_vel / SPEED_OF_LIGHT_AU_PER_DAY;
        
        // Aberration correction (to first order)
        // See Explanatory Supplement to the Astronomical Almanac, Section 7.2
        Vector3d aberration_correction = v_over_c - direction * direction.dot(v_over_c);
        direction += aberration_correction;
        direction.normalize();
    }
    
    // Convert to RA/Dec
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
    double mjd_tdb = utc_to_tdb(obs.mjd_utc);
    double jd_tdb = mjd_tdb + 2400000.5;
    
    // Get Earth position from ephemeris
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
    // This requires rotation from ITRF to ICRF (Earth rotation)
    // For now, simplified calculation using parallax constants
    
    double rho_cos_phi = obs_info.rho_cos_phi;
    double rho_sin_phi = obs_info.rho_sin_phi;
    double longitude = obs_info.longitude;
    
    // Compute Greenwich Mean Sidereal Time
    // Note: For higher accuracy, should convert UTC to UT1 (requires ΔUT1 from IERS)
    // For now, assume UTC ≈ UT1 (error < 1 second typically)
    double gmst = compute_gmst(obs.mjd_utc);
    
    // Local sidereal time = GMST + longitude
    double lst = gmst + longitude;
    
    // Observatory position in equatorial coordinates [Earth radii]
    double cos_lst = std::cos(lst);
    double sin_lst = std::sin(lst);
    
    Vector3d obs_geocentric;
    obs_geocentric[0] = rho_cos_phi * cos_lst;
    obs_geocentric[1] = rho_cos_phi * sin_lst;
    obs_geocentric[2] = rho_sin_phi;
    
    // Convert from Earth radii to AU
    double earth_radius_au = WGS84_A / AU_TO_KM;
    obs_geocentric *= earth_radius_au;
    
    // Observatory heliocentric position
    Vector3d observer_pos = earth_center + obs_geocentric;
    
    return observer_pos;
}

std::optional<Vector3d> ResidualCalculator::get_observer_velocity(
    const OpticalObservation& obs) const {
    
    // Convert to TDB
    double mjd_tdb = utc_to_tdb(obs.mjd_utc);
    double jd_tdb = mjd_tdb + 2400000.5;
    
    // Get Earth velocity
    auto earth_state = ephemeris::PlanetaryEphemeris::getState(
        ephemeris::CelestialBody::EARTH, jd_tdb);
    
    Vector3d earth_vel = earth_state.velocity();
    
    // Get observatory position to compute rotation velocity
    auto obs_pos_opt = get_observer_position(obs);
    if (!obs_pos_opt) {
        return earth_vel;  // Return Earth velocity only
    }
    
    Vector3d observer_helio = *obs_pos_opt;
    Vector3d earth_center = earth_state.position();
    Vector3d obs_geocentric = observer_helio - earth_center;
    
    // Earth rotation angular velocity [rad/day]
    // ω = 2π/T_sid where T_sid ≈ 0.99726958 solar days (sidereal day)
    double omega_earth = TWO_PI / 0.99726958;  // rad/day
    
    // Rotation axis (Earth's north pole in ecliptic coordinates)
    // For simplicity, assume z-axis (valid for low accuracy)
    // More accurate: account for obliquity and precession
    Vector3d omega_vec(0.0, 0.0, omega_earth);
    
    // Velocity due to Earth rotation: v_rot = ω × r_geocentric
    Vector3d v_rotation = omega_vec.cross(obs_geocentric);
    
    // Total observer velocity
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
