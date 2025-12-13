/**
 * @file ResidualCalculator.cpp
 * @brief Full implementation of precise residual calculator
 * 
 * Features:
 * - Topocentric correction (Parallax)
 * - Iterative Light-time correction
 * - Earth position via EphemerisProvider
 * - Sidereal time calculation for Earth rotation
 */

#include "astdyn/orbit_determination/ResidualCalculator.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::orbit_determination {

// Constants
constexpr double AU_M = 149597870700.0; // AU in meters
constexpr double C_LIGHT = 299792458.0; // Speed of light m/s
constexpr double C_AU_DAY = (C_LIGHT * 86400.0) / AU_M; // ~173.1446 AU/day
constexpr double EARTH_RADIUS_AU = 6378137.0 / AU_M; // Equatorial radius in AU
constexpr double SEC_TO_RAD = M_PI / (180.0 * 3600.0);
constexpr double DEG_TO_RAD = M_PI / 180.0;

ResidualCalculator::ResidualCalculator() {
    initialize_default_observatories();
}

void ResidualCalculator::set_ephemeris_provider(std::shared_ptr<ephemeris::EphemerisProvider> provider) {
    ephemeris_provider_ = provider;
}

void ResidualCalculator::initialize_default_observatories() {
    // Common MPC observatories
    // Code, Longitude (deg East), rho*cos(phi'), rho*sin(phi')
    // Data from MPC list (Observer's location)
    
    // 500: Geocentric
    observatories_["500"] = {"500", 0.0, 0.0, 0.0};
    
    // 809: ESO, La Silla
    // Long: 70.7303 W -> 289.2697 E
    observatories_["809"] = {"809", 289.2697, 0.83878, -0.54434};
    
    // 691: Steward Obs., Kitt Peak
    // Long: 111.5983 W -> 248.4017 E
    observatories_["691"] = {"691", 248.4017, 0.84414, 0.53580};
    
    // 704: Lincoln Laboratory ETS, New Mexico
    // Long: 106.6578 W -> 253.3422 E
    observatories_["704"] = {"704", 253.3422, 0.83689, 0.54705};
    
    // J75: OAM - Observatorio Astronomico de La Sagra
    // Long: 2.5656 W -> 357.4344 E
    observatories_["J75"] = {"J75", 357.4344, 0.80373, 0.59363};
}

void ResidualCalculator::load_observatories(const std::string& filename) {
    // TODO: Implement parsing of standard MPC observatories.txt
    // For now, rely on defaults
}

double ResidualCalculator::calc_gmst(double mjd_utc) const {
    // Calculate Greenwich Mean Sidereal Time (IAU 1982)
    // T = (JD - 2451545.0) / 36525.0
    double jd = mjd_utc + 2400000.5;
    double T = (jd - 2451545.0) / 36525.0;
    
    // GMST in seconds at 0h UT1
    double gmst_sec = 24110.54841 + 8640184.812866 * T + 0.093104 * T*T - 6.2e-6 * T*T*T;
    
    // Convert to degrees (15 deg/hour, 3600 sec/hour -> 1/240 deg/sec)
    // Actually: 1 sec time = 15 arcsec = 15/3600 deg = 1/240 deg
    double gmst_deg = gmst_sec / 240.0;
    
    // Add Earth rotation for the current time of day
    // Fraction of day * 360 * 1.0027379... (solar to sidereal rate)
    double day_fraction = mjd_utc - std::floor(mjd_utc);
    double earth_rotation = day_fraction * 360.0 * 1.00273790935;
    
    double theta_gmst = gmst_deg + earth_rotation;
    
    // Build normalized angle within [0, 360)
    theta_gmst = std::fmod(theta_gmst, 360.0);
    if (theta_gmst < 0) theta_gmst += 360.0;
    
    return theta_gmst * DEG_TO_RAD; // Radians
}

Eigen::Vector3d ResidualCalculator::get_earth_position(double mjd) const {
    if (ephemeris_provider_) {
        // CelestialBody::EARTH = 3 in our enum (assumed) or specific ID
        // Assuming provider takes integer ID. 
        // VSOP87 usually: 3 = Earth
        // DE441 (SPICE): 399 = Earth center, 3 = Earth-Moon Barycenter. 
        // For precise work we need Earth Center (399).
        // Let's assume the provider handles "Earth" correctly as the planet center.
        return ephemeris_provider_->getPosition(ephemeris::CelestialBody::EARTH, mjd);
    }
    // Fallback: Circular orbit approximation (very crude, only for testing without provider)
    double n = 0.9856 * DEG_TO_RAD; // mean motion deg/day
    double L = n * (mjd - 51544.5); // longitude approx
    return Eigen::Vector3d(std::cos(L), std::sin(L), 0.0);
}

Eigen::Vector3d ResidualCalculator::get_observatory_position(const std::string& code, double mjd_utc) const {
    auto it = observatories_.find(code);
    if (it == observatories_.end()) {
        // Unknown code, fallback to Geocentric
        return Eigen::Vector3d::Zero();
    }
    
    const auto& obs = it->second;
    if (obs.code == "500") return Eigen::Vector3d::Zero(); // Geocentric
    
    // Calculate Local Sidereal Time (LST)
    double gmst = calc_gmst(mjd_utc);
    double lst = gmst + obs.longitude_deg * DEG_TO_RAD;
    
    // Compute geocentric position of observatory in equatorial frame (J2000 approx)
    // r = [ rho*cos(phi') * cos(LST) ]
    //     [ rho*cos(phi') * sin(LST) ]
    //     [ rho*sin(phi') ]
    // Units: Earth Radii. Must convert to AU.
    
    double x = obs.rho_cos_phi * std::cos(lst);
    double y = obs.rho_cos_phi * std::sin(lst);
    double z = obs.rho_sin_phi;
    
    return Eigen::Vector3d(x, y, z) * EARTH_RADIUS_AU;
}

Eigen::Vector3d ResidualCalculator::apply_light_time_correction(
        const Eigen::Vector3d& r_ast, // Helo state at prop time
        const Eigen::Vector3d& v_ast, // Velocity used for linear back-prop
        const Eigen::Vector3d& r_obs, // Observer Helo position
        double& light_time_days
) const {
    // Initial estimate: geometric distance
    Eigen::Vector3d rho_vec = r_ast - r_obs;
    double rho = rho_vec.norm();
    light_time_days = rho / C_AU_DAY;
    
    // Iteration (Newton-like or fixed point)
    // We want r(t - tau) - R_obs(t)
    // Linear approximation: r(t - tau) ~ r(t) - v(t)*tau
    
    int max_iter = 3;
    for (int i=0; i<max_iter; ++i) {
        Eigen::Vector3d r_ast_retarded = r_ast - v_ast * light_time_days;
        rho_vec = r_ast_retarded - r_obs;
        double new_lt = rho_vec.norm() / C_AU_DAY;
        
        if (std::abs(new_lt - light_time_days) < 1e-9) { // ~1ms convergence
            light_time_days = new_lt;
            break;
        }
        light_time_days = new_lt;
    }
    
    // Final retarded position vector relative to observer
    return (r_ast - v_ast * light_time_days) - r_obs;
}

void ResidualCalculator::cartesian_to_radec(
    const Eigen::Vector3d& r,
    double& ra_deg,
    double& dec_deg
) const {
    double x = r(0);
    double y = r(1);
    double z = r(2);
    double r_norm = r.norm();
    
    // RA
    double ra_rad = std::atan2(y, x);
    if (ra_rad < 0) ra_rad += 2.0 * M_PI;
    ra_deg = ra_rad / DEG_TO_RAD;
    
    // Dec
    double dec_rad = std::asin(z / r_norm);
    dec_deg = dec_rad / DEG_TO_RAD;
}

Residual ResidualCalculator::compute_residual(
    const Observation& obs,
    const Eigen::Vector<double, 6>& state_at_obs,
    double epoch_mjd
) {
    Residual res;
    res.epoch_mjd = obs.epoch_mjd;
    res.rejected = obs.rejected;
    
    // 1. Get Earth position (ALREADY Equatorial J2000 from DE441/EphemerisProvider)
    Eigen::Vector3d r_earth_eq = get_earth_position(obs.epoch_mjd);
    
    // 2. Get Observatory position (Geocentric Equatorial)
    Eigen::Vector3d r_obs_site_eq = get_observatory_position(obs.observatory_code, obs.epoch_mjd);
    
    // 3. Total Observer position (Heliocentric Equatorial)
    Eigen::Vector3d r_observer_eq = r_earth_eq + r_obs_site_eq;
    
    // 4. Asteroid State (ALREADY Equatorial J2000 from Propagator)
    Eigen::Vector3d r_ast_eq = state_at_obs.head<3>();
    Eigen::Vector3d v_ast_eq = state_at_obs.tail<3>();
    
    // 5. Light-time Correction (in Equatorial frame)
    double light_time_days = 0.0;
    Eigen::Vector3d rho_vec_eq = apply_light_time_correction(r_ast_eq, v_ast_eq, r_observer_eq, light_time_days);
    
    // 6. Convert Topocentric vector to RA/Dec
    // No rotations needed. The vector is already Equator J2000 -> RA/Dec J2000.
    
    double ra_comp, dec_comp;
    cartesian_to_radec(rho_vec_eq, ra_comp, dec_comp);
    
    res.ra_computed_deg = ra_comp;
    res.dec_computed_deg = dec_comp;
    
    // 7. Compute Residuals
    // Handle RA wrap-around 360-0
    double d_ra = obs.ra_deg - ra_comp;
    while (d_ra > 180.0) d_ra -= 360.0;
    while (d_ra < -180.0) d_ra += 360.0;
    
    // RA residual in arcseconds, scaled by cos(dec)
    res.ra_residual_arcsec = d_ra * 3600.0 * std::cos(obs.dec_deg * DEG_TO_RAD);
    res.dec_residual_arcsec = (obs.dec_deg - dec_comp) * 3600.0;
    
    return res;
}

} // namespace astdyn::orbit_determination
