/**
 * @file Constants.hpp
 * @brief Physical and astronomical constants for orbital mechanics
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Constants converted from Fortran fund_const module.
 * Values from IAU 2015 and JPL DE431 ephemeris.
 */

#ifndef ASTDYN_CORE_CONSTANTS_HPP
#define ASTDYN_CORE_CONSTANTS_HPP

#include <cmath>

namespace astdyn {
namespace constants {

// ============================================================================
// Mathematical Constants
// ============================================================================

/// Pi
constexpr double PI = 3.14159265358979323846;

/// Pi/2
constexpr double HALF_PI = PI / 2.0;

/// 2*Pi
constexpr double TWO_PI = 2.0 * PI;

/// Degrees to radians conversion factor
constexpr double DEG_TO_RAD = PI / 180.0;

/// Radians to degrees conversion factor
constexpr double RAD_TO_DEG = 180.0 / PI;

/// Arcseconds to radians
constexpr double ARCSEC_TO_RAD = PI / (180.0 * 3600.0);

/// Radians to arcseconds
constexpr double RAD_TO_ARCSEC = (180.0 * 3600.0) / PI;

// ============================================================================
// Time Constants
// ============================================================================

/// Julian Date of J2000.0 epoch
constexpr double JD2000 = 2451545.0;

/// Modified Julian Date of J2000.0 epoch
constexpr double MJD2000 = 51544.5;

/// Seconds per day
constexpr double SECONDS_PER_DAY = 86400.0;

/// Alias for SECONDS_PER_DAY
constexpr double DAY = SECONDS_PER_DAY;

/// Days per Julian year
constexpr double DAYS_PER_YEAR = 365.25;

/// Seconds per Julian year
constexpr double YEAR = DAYS_PER_YEAR * SECONDS_PER_DAY;

/// Days per Julian century
constexpr double DAYS_PER_CENTURY = 36525.0;

// ============================================================================
// Fundamental Astronomical Constants (IAU 2015)
// ============================================================================

/// Speed of light in vacuum [km/s]
constexpr double C_LIGHT = 299792.458;

/// Astronomical Unit [km] (IAU 2012/2015)
constexpr double AU = 149597870.700;

/// Gaussian gravitational constant [AU^(3/2) / (M_sun^(1/2) * day)]
constexpr double K_GAUSS = 0.01720209895;

/// Solar radius [km] (IAU 2015)
constexpr double R_SUN = 695700.0;

/// Earth equatorial radius [km] (WGS84)
constexpr double R_EARTH = 6378.137;

/// Earth flattening (WGS84)
constexpr double EARTH_FLATTENING = 1.0 / 298.257223563;

// ============================================================================
// Gravitational Parameters [km^3/s^2] (from JPL DE431)
// ============================================================================

/// GM Sun
constexpr double GM_SUN = 132712440041.939400;

/// GM Mercury
constexpr double GM_MERCURY = 22031.868551;

/// GM Venus
constexpr double GM_VENUS = 324858.592000;

/// GM Earth+Moon system
constexpr double GM_EARTH_MOON = 403503.235502;

/// GM Earth alone (approximately, from system - Moon)
constexpr double GM_EARTH = 398600.435507;

/// GM Moon
constexpr double GM_MOON = 4902.800066;

/// GM Mars system
constexpr double GM_MARS = 42828.375816;

/// GM Jupiter system
constexpr double GM_JUPITER = 126712764.100000;

/// GM Saturn system
constexpr double GM_SATURN = 37940584.841800;

/// GM Uranus system
constexpr double GM_URANUS = 5794556.400000;

/// GM Neptune system
constexpr double GM_NEPTUNE = 6836527.100580;

/// GM Pluto system
constexpr double GM_PLUTO = 975.500000;

// ============================================================================
// Heliocentric Gravitational Constant
// ============================================================================

/// GMS = k^2 [AU^3/day^2]
constexpr double GMS = K_GAUSS * K_GAUSS;

/// GMS in SI units [km^3/s^2]
constexpr double GMS_SI = GM_SUN;

// ============================================================================
// Gravitational Parameters in AU³/day² (for orbit propagation)
// ============================================================================

/// Conversion factor: km³/s² to AU³/day²
constexpr double GM_KM3S2_TO_AU3DAY2 = (SECONDS_PER_DAY * SECONDS_PER_DAY) / (AU * AU * AU);

/// Speed of light in AU/day
constexpr double SPEED_OF_LIGHT_AU_PER_DAY = C_LIGHT * SECONDS_PER_DAY / AU;

// ============================================================================
// Mass Ratios (Planet/Sun)
// ============================================================================

/// Reciprocal of Sun/Mercury mass ratio
constexpr double MERCURY_MASS_RATIO = 1.0 / 6023600.0;

/// Reciprocal of Sun/Venus mass ratio
constexpr double VENUS_MASS_RATIO = 1.0 / 408523.71;

/// Reciprocal of Sun/Earth+Moon mass ratio
constexpr double EARTH_MOON_MASS_RATIO = 1.0 / 328900.56;

/// Reciprocal of Sun/Mars mass ratio
constexpr double MARS_MASS_RATIO = 1.0 / 3098708.0;

/// Reciprocal of Sun/Jupiter mass ratio
constexpr double JUPITER_MASS_RATIO = 1.0 / 1047.3486;

/// Reciprocal of Sun/Saturn mass ratio
constexpr double SATURN_MASS_RATIO = 1.0 / 3497.898;

/// Reciprocal of Sun/Uranus mass ratio
constexpr double URANUS_MASS_RATIO = 1.0 / 22902.98;

/// Reciprocal of Sun/Neptune mass ratio
constexpr double NEPTUNE_MASS_RATIO = 1.0 / 19412.24;

/// Reciprocal of Sun/Pluto mass ratio
constexpr double PLUTO_MASS_RATIO = 1.0 / 1.35e8;

// ============================================================================
// Planetary Data (Mean Orbital Elements - Approximate)
// ============================================================================

namespace planets {
    /// Planet ID codes (compatible with JPL HORIZONS)
    enum class PlanetID {
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
        SUN = 0
    };
} // namespace planets

// ============================================================================
// Relativistic Constants
// ============================================================================

/// Schwarzschild radius of the Sun [km]
constexpr double SCHWARZSCHILD_SUN = 2.0 * GMS_SI / (C_LIGHT * C_LIGHT);

/// 1/c^2 in AU and day units [day^2/AU^2]
constexpr double INV_C2_AU = SECONDS_PER_DAY * SECONDS_PER_DAY / 
                              (C_LIGHT * C_LIGHT * AU * AU);

// ============================================================================
// Conversion Factors
// ============================================================================

/// AU to km
constexpr double AU_TO_KM = AU;

/// km to AU
constexpr double KM_TO_AU = 1.0 / AU;

/// AU/day to km/s
constexpr double AU_PER_DAY_TO_KM_PER_S = AU / SECONDS_PER_DAY;

/// km/s to AU/day
constexpr double KM_PER_S_TO_AU_PER_DAY = SECONDS_PER_DAY / AU;

/// Solar mass in kg
constexpr double SOLAR_MASS_KG = 1.98847e30;

// ============================================================================
// Numerical Tolerances
// ============================================================================

/// Machine epsilon for double precision
constexpr double EPSILON = 1.0e-15;

/// Small value for comparisons
constexpr double SMALL = 1.0e-12;

/// Large value
constexpr double LARGE = 1.0e30;

/// Default convergence tolerance
constexpr double DEFAULT_TOLERANCE = 1.0e-14;

} // namespace constants
} // namespace astdyn

#endif // ASTDYN_CORE_CONSTANTS_HPP
