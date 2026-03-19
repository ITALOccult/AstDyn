#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>

using namespace astdyn::catalog;

int main() {
    try {
        GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
        auto& catalog = GaiaDR3Catalog::instance();
        
        OrbitQuery oq;
        oq.t_start = astdyn::time::EpochTDB::from_jd(2461170.5);
        oq.t_end = astdyn::time::EpochTDB::from_jd(2461171.5);
        oq.width = astdyn::astrometry::Angle::from_deg(0.1);
        oq.max_magnitude = 15.0;
        oq.step_days = 0.1;
        
        ChebyshevSegment seg;
        seg.t_start = 2461170.5;
        seg.t_end = 2461171.5;
        // Near RA=267, Dec=-28 (Barbara 234 approx position in May 2026)
        // Wait, let's just use RA=0, Dec=0 for simplicity
        seg.ra_coeffs = {0.0};
        seg.dec_coeffs = {0.0};
        oq.segments = {seg};
        
        std::cout << "Starting query_orbit near RA=0, Dec=0...\n";
        auto stars = catalog.query_orbit(oq);
        
        std::cout << "Found " << stars.size() << " stars via query_orbit\n";
        
        GaiaDR3Catalog::shutdown();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
