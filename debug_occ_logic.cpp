#include "astdyn/AstDyn.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>

using namespace astdyn;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.integrator_type = IntegratorType::RKF78;
        
        cout << "Testing OccultationLogic for Barbara 10 May..." << endl;
        
        AstDynEngine engine;
        engine.set_config(cfg);
        
        io::HorizonsClient horizons;
        time::EpochTDB t0 = time::EpochTDB::from_jd(2461169.5);
        auto elements_res = horizons.query_elements("234", t0);
        if (!elements_res) {
            cerr << "Failed to query Horizons." << endl;
            return 1;
        }
        
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
        
        astrometry::OccultationConfig occ_cfg;
        occ_cfg.max_mag_star = 17.5;
        occ_cfg.max_shadow_distance = physics::Distance::from_km(500000.0);
        
        cout << "Running search..." << endl;
        auto results = astrometry::OccultationLogic::find_occultations(
            "234", *elements_res, t0, t0 + time::TimeDuration::from_days(1.0), 
            occ_cfg, engine);
        
        cout << "Found " << results.size() << " occultations." << endl;
        for (const auto& res : results) {
            cout << " - Event: " << res.star.source_id << " at " << res.params.t_ca.jd() 
                 << " Mag: " << res.star.g_mag << endl;
        }

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
