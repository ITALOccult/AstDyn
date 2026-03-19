#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/ClosestApproachFinder.hpp"
#include "astdyn/astrometry/AsteroidChebyshevEphemeris.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace astdyn;
using namespace astdyn::astrometry;

int main() {
    try {
        // 1. Setup Engine & Config with AAS
        AstDynConfig cfg;
        cfg.integrator_type = IntegratorType::AAS;
        cfg.aas_precision = 1e-6; // High precision for AAS
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.verbose = true;

        // 2. Fetch Asteroid Elements (JPL Nominal)
        string id = "234";
        time::EpochTDB t0 = time::EpochTDB::from_jd(2461171.2); // May 11, 2026
        io::HorizonsClient horizons;
        auto elements = horizons.query_elements(id, t0);
        if (!elements) return 1;

        // 3. Define Star (Gaia DR3 3737957664602055808)
        auto star_opt = catalog::GaiaDR3Catalog::instance().by_source_id(3737957664602055808);
        if (!star_opt) return 1;
        const auto& star = *star_opt;

        // 4. Generate Chebyshev Ephemeris with AAS
        time::EpochTDB t_start = time::EpochTDB::from_jd(2461171.1);
        time::EpochTDB t_end = time::EpochTDB::from_jd(2461171.3);
        
        cout << "Generating Chebyshev Ephemeris using AAS integrator..." << endl;
        AsteroidChebyshevEphemeris ephem(*elements, t_start, t_end, cfg);

        // 5. Run ClosestApproachFinder
        ClosestApproachFinder::Config finder_cfg;
        finder_cfg.max_shadow_distance = astrometry::Angle::from_arcsec(120.0);
        
        cout << "Searching for closest approach..." << endl;
        auto results = ClosestApproachFinder::find_in_segment(
            ephem.get_segment(t0), star, t_start, t_end, finder_cfg.max_shadow_distance);

        cout << fixed << setprecision(6);
        cout << "\nResults with AAS Integrator:" << endl;
        cout << "Found " << results.size() << " approach(es)." << endl;
        
        for (const auto& ca : results) {
            cout << "--------------------------------------------" << endl;
            cout << "T_CA (TDB):  " << ca.t_ca.jd() << endl;
            cout << "Separation:  " << ca.separation.to_arcsec() << " arcsec" << endl;
            cout << "Impact:      " << ca.impact_parameter.to_km() << " km" << endl;
            
            // Check if it's an occultation (Earth radius ~6371 km + Barbara radius ~120 km)
            bool is_hit = ca.impact_parameter < core::Distance::from_km(6371.0 + 120.0);
            cout << "Occultation: " << (is_hit ? "YES" : "NO") << endl;
        }

    } catch (const std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
