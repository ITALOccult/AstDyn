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

#include <Eigen/Dense>
#include <string>

namespace astdyn::ephemeris {

/**
 * @brief Celestial body identifiers
 */
enum class CelestialBody {
    SUN = 10,
    MERCURY = 1,
    VENUS = 2,
    EARTH = 3,
    MARS = 4,
    JUPITER = 5,
    SATURN = 6,
    URANUS = 7,
    NEPTUNE = 8,
    PLUTO = 9,
    MOON = 301
};

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
     * @param jd_tdb Julian Date (TDB time scale)
     * @return Position vector [AU] in J2000 ecliptic frame
     */
    virtual Eigen::Vector3d getPosition(CelestialBody body, double jd_tdb) = 0;
    
    /**
     * @brief Get velocity of celestial body
     * 
     * @param body Celestial body
     * @param jd_tdb Julian Date (TDB time scale)
     * @return Velocity vector [AU/day] in J2000 ecliptic frame
     */
    virtual Eigen::Vector3d getVelocity(CelestialBody body, double jd_tdb) = 0;
    
    /**
     * @brief Get full state (position + velocity)
     * 
     * @param body Celestial body
     * @param jd_tdb Julian Date (TDB time scale)
     * @return State vector [pos (AU), vel (AU/day)] in J2000 ecliptic frame
     */
    virtual Eigen::Matrix<double, 6, 1> getState(CelestialBody body, double jd_tdb) {
        Eigen::Matrix<double, 6, 1> state;
        state.head<3>() = getPosition(body, jd_tdb);
        state.tail<3>() = getVelocity(body, jd_tdb);
        return state;
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
