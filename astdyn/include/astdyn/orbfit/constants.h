#pragma once
#include <cmath>

namespace orbfit {

// Gravitational constant * solar mass (AU^3 / day^2)
constexpr double GM_SUN    = 2.959122082855911e-4;

// Planetary GM values (AU^3 / day^2) — DE431 values
constexpr double GM_MERCURY = 4.912547451450812e-11;
constexpr double GM_VENUS   = 7.243452486162703e-10;
constexpr double GM_EARTH   = 8.887692445125634e-10;
constexpr double GM_MOON    = 1.093189450742374e-11;
constexpr double GM_MARS    = 9.549535105779258e-11;
constexpr double GM_JUPITER = 2.825345909524226e-7;
constexpr double GM_SATURN  = 8.459715185680659e-8;
constexpr double GM_URANUS  = 1.292024916781969e-8;
constexpr double GM_NEPTUNE = 1.524358900784276e-8;

// Speed of light (AU/day)
constexpr double CLIGHT     = 173.14463267424034;

// Astronomical unit (km)
constexpr double AU_KM      = 1.495978707e8;

// Earth radius (AU)
constexpr double EARTH_RADIUS = 4.263523e-5;

// Conversion factors
constexpr double DEG2RAD    = M_PI / 180.0;
constexpr double RAD2DEG    = 180.0 / M_PI;
constexpr double ARCSEC2RAD = DEG2RAD / 3600.0;
constexpr double RAD2ARCSEC = 3600.0 * RAD2DEG;

// J2000.0 epoch (Julian Date)
constexpr double JD_J2000   = 2451545.0;

// Days per Julian year
constexpr double DAY_PER_YEAR = 365.25;

} // namespace orbfit
