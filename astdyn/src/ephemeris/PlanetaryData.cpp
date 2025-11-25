/**
 * @file PlanetaryData.cpp
 * @brief Implementation of planetary data functions
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#include "astdyn/ephemeris/PlanetaryData.hpp"
#include "astdyn/core/Constants.hpp"

namespace astdyn {
namespace ephemeris {

PlanetaryBody PlanetaryData::getBodyData(CelestialBody body) {
    PlanetaryBody data;
    data.name = getName(body);
    data.gm = getGM(body);
    data.radius = getRadius(body);
    data.mass = getMass(body);
    
    // Mean orbital semi-major axes [AU] and periods [days]
    // Source: JPL Horizons, epoch J2000.0
    switch (body) {
        case CelestialBody::SUN:
            data.semi_major_axis = 0.0;
            data.period = 0.0;
            break;
        case CelestialBody::MERCURY:
            data.semi_major_axis = 0.38709927;
            data.period = 87.9691;
            break;
        case CelestialBody::VENUS:
            data.semi_major_axis = 0.72333566;
            data.period = 224.701;
            break;
        case CelestialBody::EARTH:
            data.semi_major_axis = 1.00000261;
            data.period = 365.256;
            break;
        case CelestialBody::MARS:
            data.semi_major_axis = 1.52371034;
            data.period = 686.980;
            break;
        case CelestialBody::JUPITER:
            data.semi_major_axis = 5.20288700;
            data.period = 4332.59;
            break;
        case CelestialBody::SATURN:
            data.semi_major_axis = 9.53667594;
            data.period = 10759.22;
            break;
        case CelestialBody::URANUS:
            data.semi_major_axis = 19.18916464;
            data.period = 30685.4;
            break;
        case CelestialBody::NEPTUNE:
            data.semi_major_axis = 30.06992276;
            data.period = 60189.0;
            break;
        case CelestialBody::PLUTO:
            data.semi_major_axis = 39.48211675;
            data.period = 90560.0;
            break;
        case CelestialBody::MOON:
            data.semi_major_axis = 0.00256955529; // ~384,400 km in AU
            data.period = 27.321661; // Sidereal month
            break;
        default:
            data.semi_major_axis = 0.0;
            data.period = 0.0;
            break;
    }
    
    return data;
}

} // namespace ephemeris
} // namespace astdyn
