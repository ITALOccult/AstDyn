con#include "astdyn/AstDyn.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.integrator_type = IntegratorType::RKF78;
        cfg.tolerance = 1e-13;
        cfg.propagator_settings.include_asteroids = true;
        cfg.propagator_settings.use_default_asteroid_set = true;
        
        // 1. Initial State for Barbara on May 3.5
        time::EpochTDB start_epoch = time::EpochTDB::from_jd(2461164.0);
        io::HorizonsClient horizons;
        auto state_res = horizons.query_vectors("234", start_epoch, "@10");
        if (!state_res) return 1;
        
        auto elements = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(state_res->cast_frame<core::ECLIPJ2000>());
        
        // 2. Fit Chebyshev for 1 day (May 3.5 to May 4.5)
        cout << "Fitting Chebyshev for May 4 window..." << endl;
        auto segment = catalog::fit_chebyshev(elements, start_epoch + time::TimeDuration::from_days(0.5), 1.0, cfg);
        
        // 3. Evaluate at May 4.0 00:00 (JD 2461164.5)
        auto [pos, vel] = segment.evaluate_full(2461164.5);
        cout << setprecision(10) << "Barbara at May 4.0: RA=" << std::get<0>(pos) << ", Dec=" << std::get<1>(pos) << endl;
        
        // 4. Query Gaia
        cout << "Initializing Gaia..." << endl;
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
        
        catalog::OrbitQuery oq;
        oq.t_start = start_epoch;
        oq.t_end = start_epoch + time::TimeDuration::from_days(1.0);
        oq.segments = { segment };
        oq.width = astrometry::Angle::from_deg(0.1); // 6 arcmin
        oq.max_magnitude = 16.5;
        oq.step_days = 0.01; // Every 15 minutes
        
        cout << "Querying Gaia with Radius=0.1 deg and Step=0.01 days..." << endl;
        auto stars = catalog::GaiaDR3Catalog::instance().query_orbit(oq);
        cout << "Found " << stars.size() << " stars in corridor." << endl;
        
        for (const auto& s : stars) {
            cout << " - Star: " << s.source_id << " RA=" << s.ra.to_deg() << " Dec=" << s.dec.to_deg() << " G=" << s.g_mag << endl;
        }

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
