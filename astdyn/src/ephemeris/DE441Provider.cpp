/**
 * @file DE441Provider.cpp
 * @brief Implementation of JPL DE441 ephemeris provider (Native Reader)
 * @author AstDyn Team
 * @date 2025-12-09
 */

#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/io/SPKReader.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::ephemeris {

// NAIF ID mapping
int DE441Provider::bodyToNAIFId(CelestialBody body) const {
    switch (body) {
        case CelestialBody::SUN:     return 10;
        case CelestialBody::MERCURY: return 1; // Mercury Barycenter (199 usually not in DE441 basics)
        case CelestialBody::VENUS:   return 2; // Venus Barycenter
        case CelestialBody::EARTH:   return 399; // Earth chained via 3
        case CelestialBody::MARS:    return 4; // Mars Barycenter
        case CelestialBody::JUPITER: return 5; // Jupiter Barycenter
        case CelestialBody::SATURN:  return 6; // Saturn Barycenter
        case CelestialBody::URANUS:  return 7; // Uranus Barycenter
        case CelestialBody::NEPTUNE: return 8; // Neptune Barycenter
        case CelestialBody::PLUTO:   return 9; // Pluto Barycenter
        case CelestialBody::MOON:    return 301;
        default:
            throw std::invalid_argument("Unknown celestial body ID: " + std::to_string((int)body));
    }
}

double DE441Provider::jdToET(double jd_tdb) const {
    constexpr double J2000_JD = 2451545.0;
    constexpr double SECONDS_PER_DAY = 86400.0;
    
    // Heuristic: If jd is small (e.g. < 2M), assume Modified Julian Day (MJD)
    // and convert to Full JD. AstDyn Propagator typically uses MJD.
    double full_jd = jd_tdb;
    if (full_jd < 2000000.0) {
        full_jd += 2400000.5;
    }
    
    return (full_jd - J2000_JD) * SECONDS_PER_DAY;
}

Eigen::Vector3d DE441Provider::equatorialToEcliptic(const Eigen::Vector3d& vec) const {
    constexpr double epsilon = 23.4392911 * M_PI / 180.0;
    const double c = std::cos(epsilon);
    const double s = std::sin(epsilon);
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
        reader_ = std::make_unique<io::SPKReader>(bsp_file);
        loaded_ = true;
        std::cout << "✓ DE441 loaded (Native Reader): " << bsp_file << std::endl;
    } catch (const std::exception& e) {
        loaded_ = false;
        throw std::runtime_error("Failed to load DE441: " + std::string(e.what()));
    }
}

DE441Provider::~DE441Provider() = default;

Eigen::Vector3d DE441Provider::getPosition(CelestialBody body, double jd_tdb) {
    if (!loaded_) throw std::runtime_error("DE441 not loaded");
    
    int target = bodyToNAIFId(body);
    // std::cout << "DEBUG: DE441 GetPosition Body=" << (int)body << " NAIF=" << target << " JD=" << jd_tdb << "\n";
    
    double et = jdToET(jd_tdb);
    
    // Logic for Element Chaining (Basic High Precision for Planets)
    // Most planets (Mercury-Pluto) in DE4xx are defined as Barycenters (1, 2, ... 9) relative to SSB (0).
    // Earth (399) is defined relative to EMB (3). EMB (3) is relative to SSB (0).
    // Sun (10) is relative to SSB (0).
    
    Eigen::VectorXd state(6);
    
    try {
        // Handle Earth Special Case (399 -> 3 -> 0)
        if (target == 399) { // Earth
            Eigen::VectorXd s_earth_emb = reader_->getState(399, et); // Earth wrt EMB
            Eigen::VectorXd s_emb_ssb = reader_->getState(3, et);     // EMB wrt SSB
            state = s_earth_emb + s_emb_ssb;
        }
        else if (target == 301) { // Moon (301 -> 3 -> 0)
             Eigen::VectorXd s_moon_emb = reader_->getState(301, et);
             Eigen::VectorXd s_emb_ssb = reader_->getState(3, et);
             state = s_moon_emb + s_emb_ssb;
        }
        else if (target < 10 && target > 0) {
             // Barycenters 1..9 are usually relative to SSB (0) directly in DE4xx?
             // Actually, usually headers say: 1, 2, 3 (EMB), 4, 5, 6... relative to 0.
             // But we used 199, 299... in bodyToNAIFId for Barycenters in CSPICE code?
             // Let's check typical DE4xx contents: 
             // 1 (Mercury Barycenter), 2 (Venus B), 3 (Earth-Moon B), 4 (Mars B)...
             // 199 (Merc), 299 (Venus), 399 (Earth) ... are typically NOT in DE4xx basic, or are relative to their barycenters.
             // For simplicity, we assume we want BARYCENTRIC position of the planet system for outer planets, 
             // but for precision we want the BODY center.
             // DE4xx usually provides:
             // 1,2,3,4,5,6,7,8,9 relative to 0.
             // And sometimes centers relative to barycenters if integrated.
             // Let's map CelestiaBody::JUPITER to 5 (Jupiter Barycenter) for now as default "Position".
             
             // Mapping Logic Update for Native Reader:
             // MERCURY -> 1 (Barycenter)
             // VENUS -> 2
             // EARTH -> handled above (399)
             // MARS -> 4
             // JUPITER -> 5
             // ...
             
             // Remap ID locally
             // Remap ID locally if we get high IDs but need barycenters
             int seg_id = target;
             if (target == 199) seg_id = 1;
             else if (target == 299) seg_id = 2;
             else if (target == 499) seg_id = 4;
             else if (target == 599) seg_id = 5;
             else if (target == 699) seg_id = 6;
             else if (target == 799) seg_id = 7;
             else if (target == 899) seg_id = 8;
             else if (target == 999) seg_id = 9;
             
             state = reader_->getState(seg_id, et);
             
             // If we really want body center (e.g. 599) and only 5 is available, we use 5 (Barycenter).
             // Error is small for distance planets.
             // For Earth, we did the precise chain.
        } 
        else {
             // Direct lookup (Sun, etc)
             state = reader_->getState(target, et);
        }
    } catch (const std::exception& e) {
        // Fallback or rethrow
        throw std::runtime_error("Native SPK Error: " + std::string(e.what()));
    }
    
    // Convert km to AU
    constexpr double KM_PER_AU = 149597870.691;
    Eigen::Vector3d pos_eq(state[0] / KM_PER_AU,
                           state[1] / KM_PER_AU,
                           state[2] / KM_PER_AU);
                           
    return pos_eq; // Return Equatorial J2000 directly
}

Eigen::Vector3d DE441Provider::getVelocity(CelestialBody body, double jd_tdb) {
    if (!loaded_) throw std::runtime_error("DE441 not loaded");
    // (Similar logic to getPosition but extracting indices 3,4,5)
    // Since getState reads both, we can't optimize much without refactoring.
    // Just call the same chaining logic (lazy implementation) or refactor helper.
    // ... For brevity, duplicated logic for now or implement getState helper.
    
    // Let's implement full getState helper to avoid duplication
    // But overriding here due to interface.
    // Reuse logic:
    // ... (same ID logic)
    int target = bodyToNAIFId(body);
    double et = jdToET(jd_tdb);
    Eigen::VectorXd state(6);
    
    if (target == 399) { // Earth
        state = reader_->getState(399, et) + reader_->getState(3, et);
    } else if (target == 301) {
        state = reader_->getState(301, et) + reader_->getState(3, et);
    } else {
         int seg_id = target;
         if (target == 199) seg_id = 1;
         else if (target == 299) seg_id = 2;
         else if (target == 499) seg_id = 4;
         else if (target == 599) seg_id = 5;
         else if (target == 699) seg_id = 6;
         else if (target == 799) seg_id = 7;
         else if (target == 899) seg_id = 8;
         else if (target == 999) seg_id = 9;
         state = reader_->getState(seg_id, et);
    }

    constexpr double KM_PER_AU = 149597870.691;
    constexpr double SECONDS_PER_DAY = 86400.0;
    const double conversion = SECONDS_PER_DAY / KM_PER_AU;
    
    Eigen::Vector3d vel_eq(state[3] * conversion,
                           state[4] * conversion,
                           state[5] * conversion);
                           
    return vel_eq; // Return Equatorial J2000 directly
}

} // namespace astdyn::ephemeris
