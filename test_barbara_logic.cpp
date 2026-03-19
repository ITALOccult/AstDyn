#include "astdyn/AstDyn.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace std;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.integrator_type = IntegratorType::RKF78;
        
        AstDynEngine engine;
        engine.set_config(cfg);
        
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
        
        io::HorizonsClient horizons;
        string id = "234";
        time::EpochTDB t_start = time::EpochTDB::from_jd(2461170.5);
        time::EpochTDB t_end = time::EpochTDB::from_jd(2461171.5);
        
        cout << "Searching for Barbara occultations (Nominal JPL) via Core Logic..." << endl;
        
        auto elements = horizons.query_elements(id, t_start);
        auto physical = horizons.query_physical_properties(id);
        
        if (!elements) { cerr << "Failed to get JPL elements" << endl; return 1; }
        
        // Use a WIDE search window and magnitude limit
        auto events = astrometry::OccultationLogic::find_occultations(
            id, *elements, t_start, t_end, 14.0, engine);
            
        cout << "Found " << events.size() << " occultation events." << endl;
        for (const auto& ev : events) {
            cout << " - Star: " << ev.params.star_id 
                 << " Time: " << ev.params.t_ca.mjd() 
                 << " Miss: " << ev.params.impact_parameter.to_km() << " km" << endl;
        }

    } catch (const exception& e) { cerr << "Error: " << e.what() << endl; }
    return 0;
}
