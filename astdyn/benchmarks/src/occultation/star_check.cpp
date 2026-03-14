/**
 * @file star_check.cpp
 * @brief Compare Gaia DR3 star vs XML Tycho-2 star.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::catalog;

int main() {
    std::string config = "{\"catalog_type\":\"online_esa\"}";
    GaiaDR3Catalog::initialize(config);

    int64_t source_id = 3618364055431251072;
    auto star = GaiaDR3Catalog::instance().by_source_id(source_id);

    if (star) {
        std::cout << "Star Found: Gaia DR3 " << source_id << std::endl;
        std::cout << "RA (Gaia J2016):  " << star->ra.to_deg() << " deg" << std::endl;
        std::cout << "Dec (Gaia J2016): " << star->dec.to_deg() << " deg" << std::endl;
        
        // Propagate to 2026-03-13
        time::EpochTDB t_tca = time::EpochTDB::from_jd(2461112.615671);
        auto star_2026 = star->predict_at(t_tca);
        std::cout << "RA (ICRS 2026.2): " << star_2026.ra().to_deg() << " deg" << std::endl;
        std::cout << "Dec (ICRS 2026.2): " << star_2026.dec().to_deg() << " deg" << std::endl;
        
        // XML Reference
        double xml_ra = (13.0 + 38.0/60.0 + 56.576/3600.0) * 15.0;
        double xml_dec = -(8.0 + 48.0/60.0 + 51.12/3600.0);
        std::cout << "RA (XML Ref):     " << xml_ra << " deg" << std::endl;
        std::cout << "Dec (XML Ref):    " << xml_dec << " deg" << std::endl;
        
        double sep = star_2026.separation(astrometry::SkyCoord<core::GCRF>::from_vector(
            math::Vector3<core::GCRF, physics::Distance>::from_si(
                std::cos(xml_dec*astdyn::constants::DEG_TO_RAD)*std::cos(xml_ra*astdyn::constants::DEG_TO_RAD),
                std::cos(xml_dec*astdyn::constants::DEG_TO_RAD)*std::sin(xml_ra*astdyn::constants::DEG_TO_RAD),
                std::sin(xml_dec*astdyn::constants::DEG_TO_RAD)
            ))).to_arcsec();
        std::cout << "Difference:       " << sep << " arcsec" << std::endl;
    } else {
        std::cout << "Star NOT found!" << std::endl;
    }

    GaiaDR3Catalog::shutdown();
    return 0;
}
