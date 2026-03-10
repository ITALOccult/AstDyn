/**
 * @file EphemerisProvider.hpp
 * @brief Abstract interface for planetary ephemeris providers
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * This interface allows switching between different ephemeris sources:
 * - VSOP87 (built-in, fast, ~20 arcsec accuracy)
 * - JPL DE441 (optional, ultra-precise, ~cm accuracy)
 */

#ifndef ASTDYN_EPHEMERIS_PROVIDER_HPP
#define ASTDYN_EPHEMERIS_PROVIDER_HPP

#include "astdyn/ephemeris/CelestialBody.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "src/core/frame_tags.hpp"
#include "src/core/units.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/core/physics_state.hpp"
#include <Eigen/Dense>
#include <string>

namespace astdyn::ephemeris {

/**
 * @brief Abstract base class for ephemeris providers
 */
class EphemerisProvider {
public:
    virtual ~EphemerisProvider() = default;
    
    /**
     * @brief Get position of celestial body
     * 
     * @param body Celestial body
     * @param t Time
     * @return Position vector in J2000 equatorial frame (ICRF)
     */
    virtual math::Vector3<core::GCRF, physics::Distance> getPosition(CelestialBody body, time::EpochTDB t) = 0;
    
    /**
     * @brief Get velocity of celestial body
     * 
     * @param body Celestial body
     * @param t Time
     * @return Velocity vector in J2000 equatorial frame (ICRF)
     */
    virtual math::Vector3<core::GCRF, physics::Velocity> getVelocity(CelestialBody body, time::EpochTDB t) = 0;
    
    /**
     * @brief Get full state (position + velocity)
     * 
     * @param body Celestial body
     * @param t Time
     * @return State in J2000 equatorial frame (GCRF)
     */
    virtual physics::CartesianStateTyped<core::GCRF> getState(CelestialBody body, time::EpochTDB t) {
        return physics::CartesianStateTyped<core::GCRF>::from_si(t, getPosition(body, t), getVelocity(body, t));
    }
    
    /**
     * @brief Get provider name
     */
    virtual std::string getName() const = 0;
    
    /**
     * @brief Get accuracy estimate (arcsec)
     */
    virtual double getAccuracy() const = 0;
    
    /**
     * @brief Check if provider is available/loaded
     */
    virtual bool isAvailable() const = 0;
};

} // namespace astdyn::ephemeris

#endif // ASTDYN_EPHEMERIS_PROVIDER_HPP
