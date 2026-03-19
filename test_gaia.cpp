#include "astdyn/AstDyn.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        cout << "Initializing Gaia DR3 Online..." << endl;
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
        
        auto& cat = catalog::GaiaDR3Catalog::instance();
        
        // Query a small area around Barbara's position on May 4
        // RA ~ 197.96, Dec ~ 13.78
        catalog::ConeQuery q;
        q.ra = astrometry::RightAscension::from_deg(197.963);
        q.dec = astrometry::Declination::from_deg(13.784);
        q.radius = astrometry::Angle::from_arcsec(10.0 * 60.0);
        q.max_magnitude = 17.0;
        
        cout << "Querying Gaia..." << endl;
        auto stars = cat.query_cone(q);
        
        cout << "Found " << stars.size() << " stars." << endl;
        cout << setprecision(10);
        for (const auto& s : stars) {
            cout << " - Source ID: " << s.source_id << " RA=" << s.ra.to_deg() << " Dec=" << s.dec.to_deg() << " (G=" << s.g_mag << ")" << endl;
        }

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
