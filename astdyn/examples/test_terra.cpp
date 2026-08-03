#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/ephemeris/DE441Provider.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
using namespace astdyn;
int main() {
    const char* percorso = std::getenv("ASTDYN_EPHEMERIS_PATH");
    if (!percorso) { std::cout << "serve ASTDYN_EPHEMERIS_PATH\n"; return 1; }
    auto de = std::make_shared<ephemeris::DE441Provider>(percorso);
    ephemeris::PlanetaryEphemeris eph;
    eph.setProvider(de);
    const auto t = time::EpochTDB::from_mjd(61200.0);
    const auto st = eph.getState(ephemeris::CelestialBody::EARTH, t);
    const auto p = st.position.to_eigen_si() / 1.495978707e11;
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "Terra da noi (GCRF equatoriale?), AU:\n";
    std::cout << "  X = " << p.x() << "\n  Y = " << p.y() << "\n  Z = " << p.z() << "\n";
    std::cout << "\nHorizons eclittico:\n";
    std::cout << "  X = -0.213964399546\n  Y = -0.992221951067\n  Z =  0.000059704737\n";
    return 0;
}
