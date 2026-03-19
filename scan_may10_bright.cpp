#include "astdyn/AstDyn.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
        
        catalog::ConeQuery q;
        // Center of May 10 path
        q.ra = astrometry::RightAscension::from_deg(196.705);
        q.dec = astrometry::Declination::from_deg(15.04);
        q.radius = astrometry::Angle::from_deg(0.3); // Covers the day path
        q.max_magnitude = 14.0;
        
        cout << "Querying Gaia for May 10 corridor (G < 14)..." << endl;
        auto stars = catalog::GaiaDR3Catalog::instance().query_cone(q);
        cout << "Found " << stars.size() << " stars." << endl;
        
        cout << setprecision(10);
        for (const auto& s : stars) {
            cout << " - Source ID: " << s.source_id << " RA=" << s.ra.to_deg() << " Dec=" << s.dec.to_deg() << " G=" << s.g_mag << endl;
        }

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
