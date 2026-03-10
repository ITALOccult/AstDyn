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
#include "astdyn/core/Constants.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"

namespace astdyn::ephemeris {

using Vector3d = Eigen::Vector3d;

// NAIF ID mapping
int DE441Provider::bodyToNAIFId(CelestialBody body) const {
    switch (body) {
        case CelestialBody::SUN:     return 10;
        case CelestialBody::MERCURY: return 1; 
        case CelestialBody::VENUS:   return 2; 
        case CelestialBody::EARTH:   return 399; 
        case CelestialBody::MARS:    return 4; 
        case CelestialBody::JUPITER: return 5; 
        case CelestialBody::SATURN:  return 6; 
        case CelestialBody::URANUS:  return 7; 
        case CelestialBody::NEPTUNE: return 8; 
        case CelestialBody::PLUTO:   return 9; 
        case CelestialBody::MOON:    return 301;
        default:
            throw std::invalid_argument("Unknown celestial body ID: " + std::to_string((int)body));
    }
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

Eigen::VectorXd DE441Provider::readState(CelestialBody body, time::EpochTDB t) const {
    if (!loaded_) throw std::runtime_error("DE441 not loaded");
    
    int target = bodyToNAIFId(body);
    // Convert JD (TDB) to ET
    double et = (t.jd() - constants::JD2000) * 86400.0;
    
    Eigen::VectorXd state(6);
    try {
        if (target == 399) { // Earth (399 -> 3 -> 0)
            try {
                Eigen::VectorXd s_earth_emb = reader_->getState(399, et);
                Eigen::VectorXd s_emb_ssb = reader_->getState(3, et);
                state = s_earth_emb + s_emb_ssb;
            } catch (...) {
                // Fallback: If 399 is not in file, use 3 (EMB) as Earth approximation
                state = reader_->getState(3, et);
            }
        }
        else if (target == 301) { // Moon (301 -> 3 -> 0)
             try {
                 Eigen::VectorXd s_moon_emb = reader_->getState(301, et);
                 Eigen::VectorXd s_emb_ssb = reader_->getState(3, et);
                 state = s_moon_emb + s_emb_ssb;
             } catch (...) {
                 state = reader_->getState(3, et);
             }
        }
        else {
             int seg_id = target;
             // Remap body codes to barycenters if needed (Simplified mapping for DE441)
             if (target < 10 && target > 0) {
                 // Already barycenter IDs 1-9
             } else {
                 if (target == 199) seg_id = 1;
                 else if (target == 299) seg_id = 2;
                 else if (target == 499) seg_id = 4;
                 else if (target == 599) seg_id = 5;
                 else if (target == 699) seg_id = 6;
                 else if (target == 799) seg_id = 7;
                 else if (target == 899) seg_id = 8;
                 else if (target == 999) seg_id = 9;
             }
             state = reader_->getState(seg_id, et);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Native SPK Error: " + std::string(e.what()));
    }
    return state;
}

math::Vector3<core::GCRF, physics::Distance> DE441Provider::getPosition(CelestialBody body, time::EpochTDB t) {
    Eigen::VectorXd state = readState(body, t);
    return math::Vector3<core::GCRF, physics::Distance>::from_si(state[0] * 1000.0, 
                                                                state[1] * 1000.0, 
                                                                state[2] * 1000.0);
}

math::Vector3<core::GCRF, physics::Velocity> DE441Provider::getVelocity(CelestialBody body, time::EpochTDB t) {
    Eigen::VectorXd state = readState(body, t);
    return math::Vector3<core::GCRF, physics::Velocity>::from_si(state[3] * 1000.0, 
                                                                state[4] * 1000.0, 
                                                                state[5] * 1000.0);
}

physics::CartesianStateTyped<core::GCRF> DE441Provider::getState(CelestialBody body, time::EpochTDB t) {
    Eigen::VectorXd raw = readState(body, t);
    // SPK is km and km/s -> convert to SI (meters, m/s)
    auto p = math::Vector3<core::GCRF, physics::Distance>::from_si(raw[0] * 1000.0, raw[1] * 1000.0, raw[2] * 1000.0);
    auto v = math::Vector3<core::GCRF, physics::Velocity>::from_si(raw[3] * 1000.0, raw[4] * 1000.0, raw[5] * 1000.0);
    return physics::CartesianStateTyped<core::GCRF>::from_si(t, p, v);
}

} // namespace astdyn::ephemeris
