/**
 * @file EphemerisInterface.hpp
 * @brief Abstract interface for planetary ephemeris providers
 * 
 * Allows switching between different ephemeris sources:
 * - VSOP87 (analytical, built-in)
 * - JPL DE405/DE441 (numerical, high accuracy)
 * - Custom implementations
 */

#ifndef ORBFIT_EPHEMERIS_INTERFACE_HPP
#define ORBFIT_EPHEMERIS_INTERFACE_HPP

#include <orbfit/coordinates/CartesianState.hpp>
#include <Eigen/Dense>
#include <string>
#include <memory>

namespace orbfit {
namespace ephemeris {

/**
 * @brief Celestial body identifiers
 */
enum class CelestialBody {
    SUN = 0,
    MERCURY = 1,
    VENUS = 2,
    EARTH = 3,
    MARS = 4,
    JUPITER = 5,
    SATURN = 6,
    URANUS = 7,
    NEPTUNE = 8,
    PLUTO = 9,
    MOON = 10,
    EARTH_MOON_BARYCENTER = 3  // Same as Earth for some contexts
};

/**
 * @brief Ephemeris source identifier
 */
enum class EphemerisSource {
    VSOP87,      ///< Analytical VSOP87 theory (built-in, 1-20 arcsec)
    JPL_DE405,   ///< JPL DE405 (1600-2200, 1 km accuracy)
    JPL_DE441,   ///< JPL DE441 (1550-2650, cm-level accuracy)
    CUSTOM       ///< User-provided implementation
};

/**
 * @brief Abstract interface for planetary ephemerides
 */
class IEphemerisProvider {
public:
    virtual ~IEphemerisProvider() = default;
    
    /**
     * @brief Get heliocentric position of a celestial body
     * @param body Celestial body
     * @param jd_tdb Julian Date in TDB time scale
     * @return Position vector [AU] in J2000 ecliptic frame
     */
    virtual Eigen::Vector3d getPosition(CelestialBody body, double jd_tdb) = 0;
    
    /**
     * @brief Get heliocentric velocity of a celestial body
     * @param body Celestial body
     * @param jd_tdb Julian Date in TDB time scale
     * @return Velocity vector [AU/day] in J2000 ecliptic frame
     */
    virtual Eigen::Vector3d getVelocity(CelestialBody body, double jd_tdb) = 0;
    
    /**
     * @brief Get heliocentric state (position + velocity)
     * @param body Celestial body
     * @param jd_tdb Julian Date in TDB time scale
     * @return CartesianState in J2000 ecliptic frame
     */
    virtual coordinates::CartesianState getState(CelestialBody body, double jd_tdb) = 0;
    
    /**
     * @brief Get Sun position relative to solar system barycenter
     * @param jd_tdb Julian Date in TDB time scale
     * @return Barycentric position [AU] in J2000 ecliptic frame
     */
    virtual Eigen::Vector3d getSunBarycentricPosition(double jd_tdb) = 0;
    
    /**
     * @brief Get ephemeris source identifier
     */
    virtual EphemerisSource getSource() const = 0;
    
    /**
     * @brief Get valid time range for this ephemeris
     * @return Pair of (min_jd, max_jd)
     */
    virtual std::pair<double, double> getValidRange() const = 0;
    
    /**
     * @brief Check if ephemeris is valid for given time
     */
    virtual bool isValidForTime(double jd_tdb) const {
        auto range = getValidRange();
        return jd_tdb >= range.first && jd_tdb <= range.second;
    }
};

/**
 * @brief Factory for creating ephemeris providers
 */
class EphemerisFactory {
public:
    /**
     * @brief Create ephemeris provider from source type
     * @param source Desired ephemeris source
     * @param file_path Path to ephemeris file (for JPL DE)
     * @return Unique pointer to ephemeris provider
     */
    static std::unique_ptr<IEphemerisProvider> create(
        EphemerisSource source,
        const std::string& file_path = "");
    
    /**
     * @brief Get default ephemeris provider (VSOP87)
     */
    static std::unique_ptr<IEphemerisProvider> createDefault();
};

} // namespace ephemeris
} // namespace orbfit

#endif // ORBFIT_EPHEMERIS_INTERFACE_HPP
