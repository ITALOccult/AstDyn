/**
 * @file diagnostic_xml_tca.cpp
 * @brief Check geometric separation at XML TCA.
 */

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
    
    // XML TCA: 02:46:34 UT -> JD 2461112.615671
    time::EpochTDB t_xml = time::EpochTDB::from_jd(2461112.615671);

    // Star XML (Gaia-like)
    double xml_ra = 204.73217535;
    double xml_dec = -8.8139102;
    astrometry::SkyCoord<core::GCRF> star_xml = astrometry::SkyCoord<core::GCRF>::from_vector(
        math::Vector3<core::GCRF, physics::Distance>::from_si(
            std::cos(xml_dec*constants::DEG_TO_RAD)*std::cos(xml_ra*constants::DEG_TO_RAD),
            std::cos(xml_dec*constants::DEG_TO_RAD)*std::sin(xml_ra*constants::DEG_TO_RAD),
            std::sin(xml_dec*constants::DEG_TO_RAD)
        ));

    // Asteroid Geometric (Horizons)
    auto h_geo = horizons.query_vectors("50936", t_xml, "@399");
    auto sky_geo = astrometry::SkyCoord<core::GCRF>::from_vector(h_geo->position);

    std::cout << "--- Separation at XML TCA (Geometric) ---" << std::endl;
    std::cout << "Asteroid RA:  " << sky_geo.ra().to_deg() << " deg" << std::endl;
    std::cout << "Star RA:      " << star_xml.ra().to_deg() << " deg" << std::endl;
    std::cout << "Separation:   " << sky_geo.separation(star_xml).to_arcsec() << " arcsec" << std::endl;
    
    return 0;
}
