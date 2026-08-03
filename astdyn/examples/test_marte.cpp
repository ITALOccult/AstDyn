#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/ephemeris/DE441Provider.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace astdyn;
int main() {
    const char* p = std::getenv("ASTDYN_EPHEMERIS_PATH");
    if (!p) { std::cout << "serve ASTDYN_EPHEMERIS_PATH\n"; return 1; }
    auto de = std::make_shared<ephemeris::DE441Provider>(p);
    const auto t = time::EpochTDB::from_mjd(61200.0);

    std::cout << std::fixed << std::setprecision(9)
              << "epoca: MJD " << t.mjd() << "  JD " << t.jd() << "\n"
              << "atteso JD 2461200.500000000\n\n";
    auto marte = de->getPosition(ephemeris::CelestialBody::MARS, t).to_eigen_si();
    auto terra = de->getPosition(ephemeris::CelestialBody::EARTH, t).to_eigen_si();
    std::cout << std::setprecision(12)
              << "TERRA (baricentrica, AU)\n"
              << "  nostro   " << terra.x()/1.495978707e11 << " "
              << terra.y()/1.495978707e11 << " " << terra.z()/1.495978707e11 << "\n"
              << "MARTE (baricentrico, AU)\n"
              << "  nostro   " << marte.x()/1.495978707e11 << " "
              << marte.y()/1.495978707e11 << " " << marte.z()/1.495978707e11 << "\n\n";
    Eigen::Vector3d g = marte - terra;

    // tempo-luce: una sola iterazione basta per l'ordine di grandezza
    const double c = 299792458.0;
    for (int i = 0; i < 6; ++i) {
        const double tau = g.norm() / c / 86400.0;
        auto m2 = de->getPosition(ephemeris::CelestialBody::MARS,
                                  time::EpochTDB::from_mjd(61200.0 - tau)).to_eigen_si();
        g = m2 - terra;
        double ra_i = std::atan2(g.y(), g.x()) * 180.0 / M_PI;
        if (ra_i < 0) ra_i += 360.0;
        std::cout << std::setprecision(6) << "  iter " << i
                  << "  tau " << tau*86400.0 << " s   RA " << std::setprecision(8)
                  << ra_i/15.0 << " h\n";
    }
    std::cout << "\n";
    const double r = g.norm();
    double ra = std::atan2(g.y(), g.x()) * 180.0 / M_PI;
    if (ra < 0) ra += 360.0;
    const double dec = std::asin(g.z() / r) * 180.0 / M_PI;

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Marte, MJD 61200, astrometrico geocentrico\n\n";
    std::cout << "  nostro    RA " << ra/15.0 << " h   Dec " << dec << "\n";
    std::cout << "  Horizons  RA " << (2 + 51.0/60 + 58.62/3600) << " h   Dec "
              << (15 + 53.0/60 + 51.6/3600) << "\n\n";
    const double dra = (ra/15.0 - (2 + 51.0/60 + 58.62/3600)) * 15 * 3600
                     * std::cos(dec * M_PI/180.0);
    const double ddec = (dec - (15 + 53.0/60 + 51.6/3600)) * 3600;
    std::cout << "\n  moduli: Terra " << std::setprecision(9)
              << terra.norm()/1.495978707e11 << " AU, Marte "
              << marte.norm()/1.495978707e11 << " AU\n";
    std::cout << "  (Horizons: Terra 1.015 AU dal Sole, Marte 1.577)\n\n";
    std::cout << "  scarto    " << std::setprecision(3) << dra << "\" in RA, "
              << ddec << "\" in Dec\n";
    return 0;
}
