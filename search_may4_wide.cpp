#include "astdyn/AstDyn.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
        
        catalog::CorridorQuery q;
        // Barbara path on May 4 
        // 00:00: 197.962948, 13.783696
        // 24:00: 197.700147, 14.159341
        
        q.path.push_back(astrometry::SkyCoord<core::GCRF>::from_deg(197.962948, 13.783696));
        q.path.push_back(astrometry::SkyCoord<core::GCRF>::from_deg(197.700147, 14.159341));
        
        q.width = astrometry::Angle::from_arcsec(300.0); // 5 arcminutes width
        q.max_magnitude = 17.0;
        
        cout << "Querying Gaia corridor for May 4 path (Width=5 arcmin)..." << endl;
        auto stars = catalog::GaiaDR3Catalog::instance().query_corridor(q);
        cout << "Found " << stars.size() << " stars." << endl;
        
        cout << setprecision(10);
        for (const auto& s : stars) {
            cout << " - Source ID: " << s.source_id << " RA=" << s.ra.to_deg() << " Dec=" << s.dec.to_deg() << " G=" << s.g_mag << endl;
        }

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
