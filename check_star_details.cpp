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
        q.ra = astrometry::RightAscension::from_deg(197.9243368);
        q.dec = astrometry::Declination::from_deg(13.83229046);
        q.radius = astrometry::Angle::from_arcsec(5.0); // Slightly wider to be sure
        q.max_magnitude = 15.0;
        
        cout << "Querying Gaia star at RA=197.9243368 Dec=13.83229046..." << endl;
        auto stars = catalog::GaiaDR3Catalog::instance().query_cone(q);
        
        if (!stars.empty()) {
            const auto& star = stars[0];
            cout << "Source ID: " << star.source_id << endl;
            cout << "RA  (deg): " << star.ra.to_deg() << endl;
            cout << "Dec (deg): " << star.dec.to_deg() << endl;
            cout << "PM RA (mas/yr): " << star.pm_ra_cosdec.to_mas_yr() << endl;
            cout << "PM Dec (mas/yr): " << star.pm_dec.to_mas_yr() << endl;
            cout << "Parallax (mas): " << star.parallax.to_mas() << endl;
            cout << "G Mag: " << star.g_mag << endl;
            
            // Predict at 2026-05-04 03:33 UTC (MJD 61164.148)
            double jd = 61164.148 + 2400000.5;
            auto coord = star.predict_at(time::EpochTDB::from_jd(jd));
            cout << "Predicted (2026.3): RA=" << coord.ra().to_deg() << " Dec=" << coord.dec().to_deg() << endl;
        } else {
            cout << "Star not found." << endl;
        }

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
