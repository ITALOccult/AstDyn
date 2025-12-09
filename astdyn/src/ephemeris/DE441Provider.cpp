/**
 * @file DE441Provider.cpp
 * @brief Implementation of JPL DE441 ephemeris provider
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * COMPILATION:
 * Requires CSPICE library. Link with: -lcspice
 * 
 * SETUP:
 * 1. Download CSPICE: https://naif.jpl.nasa.gov/naif/toolkit_C.html
 * 2. Place de441.bsp in a known location (e.g., ~/data/spice/)
 * 3. Set environment variable: export SPICE_DATA=~/data/spice
 */

#include "astdyn/ephemeris/DE441Provider.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::ephemeris {

// NAIF ID mapping
// Reference: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
int DE441Provider::bodyToNAIFId(CelestialBody body) const {
    switch (body) {
        case CelestialBody::SUN:     return 10;
        case CelestialBody::MERCURY: return 199;  // Mercury barycenter
        case CelestialBody::VENUS:   return 299;  // Venus barycenter
        case CelestialBody::EARTH:   return 399;  // Earth barycenter
        case CelestialBody::MARS:    return 499;  // Mars barycenter
        case CelestialBody::JUPITER: return 599;  // Jupiter barycenter
        case CelestialBody::SATURN:  return 699;  // Saturn barycenter
        case CelestialBody::URANUS:  return 799;  // Uranus barycenter
        case CelestialBody::NEPTUNE: return 899;  // Neptune barycenter
        case CelestialBody::PLUTO:   return 999;  // Pluto barycenter
        case CelestialBody::MOON:    return 301;  // Moon
        default:
            throw std::invalid_argument("Unknown celestial body");
    }
}

double DE441Provider::jdToET(double jd_tdb) const {
    // ET (Ephemeris Time) in SPICE is seconds past J2000
    // JD 2451545.0 = J2000.0 epoch
    constexpr double J2000_JD = 2451545.0;
    constexpr double SECONDS_PER_DAY = 86400.0;
    
    return (jd_tdb - J2000_JD) * SECONDS_PER_DAY;
}

Eigen::Vector3d DE441Provider::equatorialToEcliptic(const Eigen::Vector3d& vec) const {
    // Obliquity of ecliptic at J2000.0
    constexpr double epsilon = 23.4392911 * M_PI / 180.0;
    const double c = std::cos(epsilon);
    const double s = std::sin(epsilon);
    
    // Rotation matrix: Rx(epsilon)
    Eigen::Vector3d result;
    result.x() = vec.x();
    result.y() = c * vec.y() + s * vec.z();
    result.z() = -s * vec.y() + c * vec.z();
    
    return result;
}

DE441Provider::DE441Provider(const std::string& bsp_file)
    : bsp_file_(bsp_file)
{
    try {
        // Load SPK kernel
        furnsh_c(bsp_file.c_str());
        loaded_ = true;
        
        std::cout << "✓ DE441 loaded: " << bsp_file << std::endl;
        
    } catch (const std::exception& e) {
        loaded_ = false;
        throw std::runtime_error("Failed to load DE441: " + std::string(e.what()));
    }
}

Eigen::Vector3d DE441Provider::getPosition(CelestialBody body, double jd_tdb) {
    if (!loaded_) {
        throw std::runtime_error("DE441 not loaded");
    }
    
    // Convert to NAIF ID and ET
    int naif_id = bodyToNAIFId(body);
    double et = jdToET(jd_tdb);
    
    // Get state from SPICE
    double state[6];
    double lt;  // Light time (not used)
    
    // Target body relative to Solar System Barycenter
    char target[32];
    snprintf(target, sizeof(target), "%d", naif_id);
    
    spkezr_c(target, et, "J2000", "NONE", "0", state, &lt);
    
    // Convert from km to AU
    constexpr double KM_PER_AU = 149597870.691;
    Eigen::Vector3d pos_eq(state[0] / KM_PER_AU,
                           state[1] / KM_PER_AU,
                           state[2] / KM_PER_AU);
    
    // Convert from J2000 equatorial to ecliptic
    return equatorialToEcliptic(pos_eq);
}

Eigen::Vector3d DE441Provider::getVelocity(CelestialBody body, double jd_tdb) {
    if (!loaded_) {
        throw std::runtime_error("DE441 not loaded");
    }
    
    int naif_id = bodyToNAIFId(body);
    double et = jdToET(jd_tdb);
    
    double state[6];
    double lt;
    
    char target[32];
    snprintf(target, sizeof(target), "%d", naif_id);
    
    spkezr_c(target, et, "J2000", "NONE", "0", state, &lt);
    
    // Convert from km/s to AU/day
    constexpr double KM_PER_AU = 149597870.691;
    constexpr double SECONDS_PER_DAY = 86400.0;
    const double conversion = SECONDS_PER_DAY / KM_PER_AU;
    
    Eigen::Vector3d vel_eq(state[3] * conversion,
                           state[4] * conversion,
                           state[5] * conversion);
    
    return equatorialToEcliptic(vel_eq);
}

} // namespace astdyn::ephemeris
