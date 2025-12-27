#include <iostream>
#include <iomanip>
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/ephemeris/CelestialBody.hpp"
#include "astdyn/time/TimeScale.hpp"

using namespace astdyn;
using namespace astdyn::ephemeris;

int main() {
    std::string de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    auto provider = std::make_shared<DE441Provider>(de441_path);

    // Date: 2026-01-10 00:00:00 UTC
    // MJD UTC = 61050.0
    double mjd_utc = 61050.0;
    double mjd_tdb = time::utc_to_tdb(mjd_utc);
    double jd_tdb = mjd_tdb + 2400000.5;

    std::cout << "MJD TDB: " << std::fixed << std::setprecision(8) << mjd_tdb << "\n";

    // JPL Horizons Truth for 2026-Jan-10 00:00:00 TDB (from previous check)
    // Earth: X = -3.11196162E-01, Y =  9.08866170E-01, Z =  3.93880344E-01 AU
    // Sun:   X =  6.79092809E-03, Y = -4.39695679E-03, Z = -2.33610058E-03 AU
    
    auto earth_pos = provider->getPosition(CelestialBody::EARTH, mjd_tdb);
    auto sun_pos = provider->getPosition(CelestialBody::SUN, mjd_tdb);

    std::cout << "\n--- Earth Position (Equatorial J2000 wrt SSB) ---\n";
    std::cout << "AstDyn: " << earth_pos.x() << ", " << earth_pos.y() << ", " << earth_pos.z() << " AU\n";
    std::cout << "JPL:    -0.31119616, 0.90886617, 0.39388034 AU\n";

    std::cout << "\n--- Sun Position (Equatorial J2000 wrt SSB) ---\n";
    std::cout << "AstDyn: " << sun_pos.x() << ", " << sun_pos.y() << ", " << sun_pos.z() << " AU\n";
    std::cout << "JPL:     0.00679093, -0.00439696, -0.00233610 AU\n";

    return 0;
}
