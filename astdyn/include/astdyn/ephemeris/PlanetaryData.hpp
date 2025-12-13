/**
 * @file PlanetaryData.hpp
 * @brief Planetary masses, radii, and physical constants
 * @date 2025-11-23
 */

#ifndef ASTDYN_PLANETARY_DATA_HPP
#define ASTDYN_PLANETARY_DATA_HPP

#include <cmath>
#include <array>
#include <vector>
#include <string>
#include <map>
#include "astdyn/ephemeris/CelestialBody.hpp"

namespace astdyn::ephemeris {

/**
 * @brief Physical and orbital data for a celestial body
 */
struct PlanetaryBody {
    std::string name;
    double gm;              ///< Gravitational parameter GM [km³/s²]
    double radius;          ///< Mean radius [km]
    double mass;            ///< Mass [kg]
    double semi_major_axis; ///< Mean orbital semi-major axis [AU]
    double period;          ///< Orbital period [days]
};

/**
 * @brief Planetary data repository
 */
class PlanetaryData {
public:
    // ========================================================================
    // Gravitational Parameters (GM) [km³/s²]
    // ========================================================================
    
    static constexpr double GM_SUN     = 1.32712440041279419e11; ///< Sun
    static constexpr double GM_MERCURY = 2.2031868551e4;         ///< Mercury
    static constexpr double GM_VENUS   = 3.24858592e5;           ///< Venus
    static constexpr double GM_EARTH   = 3.98600435507e5;        ///< Earth
    static constexpr double GM_MARS    = 4.282837362e4;          ///< Mars
    static constexpr double GM_JUPITER = 1.26712764100000e8;     ///< Jupiter
    static constexpr double GM_SATURN  = 3.79405852000000e7;     ///< Saturn
    static constexpr double GM_URANUS  = 5.79454900700000e6;     ///< Uranus
    static constexpr double GM_NEPTUNE = 6.83653406400000e6;     ///< Neptune
    static constexpr double GM_PLUTO   = 9.755e2;                ///< Pluto
    static constexpr double GM_MOON    = 4.90280007e3;           ///< Moon
    
    // ========================================================================
    // Mean Radii [km]
    // ========================================================================
    
    static constexpr double RADIUS_SUN     = 695700.0;   ///< Sun
    static constexpr double RADIUS_MERCURY = 2439.7;     ///< Mercury
    static constexpr double RADIUS_VENUS   = 6051.8;     ///< Venus
    static constexpr double RADIUS_EARTH   = 6371.0084;  ///< Earth
    static constexpr double RADIUS_MARS    = 3389.5;     ///< Mars
    static constexpr double RADIUS_JUPITER = 69911.0;    ///< Jupiter
    static constexpr double RADIUS_SATURN  = 58232.0;    ///< Saturn
    static constexpr double RADIUS_URANUS  = 25362.0;    ///< Uranus
    static constexpr double RADIUS_NEPTUNE = 24622.0;    ///< Neptune
    static constexpr double RADIUS_PLUTO   = 1188.3;     ///< Pluto
    static constexpr double RADIUS_MOON    = 1737.4;     ///< Moon
    
    // ========================================================================
    // Masses [kg]
    // ========================================================================
    
    static constexpr double MASS_SUN     = 1.98847e30;   ///< Sun
    static constexpr double MASS_MERCURY = 3.30110e23;   ///< Mercury
    static constexpr double MASS_VENUS   = 4.86732e24;   ///< Venus
    static constexpr double MASS_EARTH   = 5.97217e24;   ///< Earth
    static constexpr double MASS_MARS    = 6.41693e23;   ///< Mars
    static constexpr double MASS_JUPITER = 1.89813e27;   ///< Jupiter
    static constexpr double MASS_SATURN  = 5.68319e26;   ///< Saturn
    static constexpr double MASS_URANUS  = 8.68103e25;   ///< Uranus
    static constexpr double MASS_NEPTUNE = 1.02410e26;   ///< Neptune
    static constexpr double MASS_PLUTO   = 1.46000e22;   ///< Pluto
    static constexpr double MASS_MOON    = 7.34630e22;   ///< Moon
    
    // ========================================================================
    // Helper Functions
    // ========================================================================
    
    static double getGM(CelestialBody body) {
        switch (body) {
            case CelestialBody::SUN: return GM_SUN;
            case CelestialBody::MERCURY: return GM_MERCURY;
            case CelestialBody::VENUS: return GM_VENUS;
            case CelestialBody::EARTH: return GM_EARTH;
            case CelestialBody::MARS: return GM_MARS;
            case CelestialBody::JUPITER: return GM_JUPITER;
            case CelestialBody::SATURN: return GM_SATURN;
            case CelestialBody::URANUS: return GM_URANUS;
            case CelestialBody::NEPTUNE: return GM_NEPTUNE;
            case CelestialBody::PLUTO: return GM_PLUTO;
            case CelestialBody::MOON: return GM_MOON;
            default: return 0.0;
        }
    }
    
    static double getRadius(CelestialBody body) {
        switch (body) {
            case CelestialBody::SUN: return RADIUS_SUN;
            case CelestialBody::MERCURY: return RADIUS_MERCURY;
            case CelestialBody::VENUS: return RADIUS_VENUS;
            case CelestialBody::EARTH: return RADIUS_EARTH;
            case CelestialBody::MARS: return RADIUS_MARS;
            case CelestialBody::JUPITER: return RADIUS_JUPITER;
            case CelestialBody::SATURN: return RADIUS_SATURN;
            case CelestialBody::URANUS: return RADIUS_URANUS;
            case CelestialBody::NEPTUNE: return RADIUS_NEPTUNE;
            case CelestialBody::PLUTO: return RADIUS_PLUTO;
            case CelestialBody::MOON: return RADIUS_MOON;
            default: return 0.0;
        }
    }
    
    static std::string getName(CelestialBody body) {
        switch (body) {
            case CelestialBody::SUN: return "Sun";
            case CelestialBody::MERCURY: return "Mercury";
            case CelestialBody::VENUS: return "Venus";
            case CelestialBody::EARTH: return "Earth";
            case CelestialBody::MARS: return "Mars";
            case CelestialBody::JUPITER: return "Jupiter";
            case CelestialBody::SATURN: return "Saturn";
            case CelestialBody::URANUS: return "Uranus";
            case CelestialBody::NEPTUNE: return "Neptune";
            case CelestialBody::PLUTO: return "Pluto";
            case CelestialBody::MOON: return "Moon";
            default: return "Unknown";
        }
    }

    /**
     * @brief Get mass of celestial body [kg]
     */
    static double getMass(CelestialBody body);

    /**
     * @brief Get full planetary data structure
     */
    static PlanetaryBody getBodyData(CelestialBody body);
};

} // namespace astdyn::ephemeris

#endif // ASTDYN_PLANETARY_DATA_HPP
