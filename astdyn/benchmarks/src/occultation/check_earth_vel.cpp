
#include "astdyn/AstDyn.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    auto ephem = std::make_shared<ephemeris::DE441Provider>(bsp_path);
    ephemeris::PlanetaryEphemeris::setProvider(ephem);

    io::HorizonsClient horizons;
    time::EpochTDB t = time::EpochTDB::from_jd(2461112.615);

    auto h_earth = horizons.query_vectors("399", t, "@ssb");
    auto a_earth = ephem->getState(ephemeris::CelestialBody::EARTH, t);

    std::cout << "--- Earth Velocity Comparison (SSB, GCRF) ---" << std::endl;
    std::cout << "Horizons VX: " << std::fixed << std::setprecision(12) << h_earth->velocity.x_si() << " m/s" << std::endl;
    std::cout << "AstDyn   VX: " << std::fixed << std::setprecision(12) << a_earth.velocity.x_si() << " m/s" << std::endl;
    std::cout << "Delta VX:    " << (a_earth.velocity.x_si() - h_earth->velocity.x_si()) << " m/s" << std::endl;
    std::cout << "Delta VY:    " << (a_earth.velocity.y_si() - h_earth->velocity.y_si()) << " m/s" << std::endl;
    std::cout << "Delta VZ:    " << (a_earth.velocity.z_si() - h_earth->velocity.z_si()) << " m/s" << std::endl;

    return 0;
}
