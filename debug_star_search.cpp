#include "astdyn/AstDyn.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include <iostream>

using namespace astdyn;
using namespace astdyn::catalog;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
        
        io::HorizonsClient horizons;
        time::EpochTDB t = time::EpochTDB::from_jd(2461171.1);
        cout << "Querying JPL elements for Barbara (234)..." << endl;
        auto elements = horizons.query_elements("234", t);
        
        cout << "Fitting Chebyshev segment..." << endl;
        auto segment = catalog::fit_chebyshev(*elements, t, 0.2, cfg);
        
        auto [pos, vel] = segment.evaluate_full(2461171.1);
        cout << "Segment RA=" << std::get<0>(pos) << " deg, Dec=" << std::get<1>(pos) << " deg" << endl;

        Angle radius = Angle::from_arcsec(120.0); 
        auto stars = catalog::find_stars_near_segment(catalog::GaiaDR3Catalog::instance(), segment, radius, 15.0);
        
        cout << "Found " << stars.size() << " stars." << endl;
        for (const auto& star : stars) {
            if (star.source_id == 3737957664602055808) {
                cout << "SUCCESS: Found target star 3737957664602055808!" << endl;
            }
        }
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
    return 0;
}
