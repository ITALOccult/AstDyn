/**
 * @file occultation_search_demo.cpp
 * @brief Demo for automatic occultation discovery and refinement.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/time/epoch.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::astrometry;

int main() {
    std::cout << "=== ASTDYN DEMO: Automatic Occultation Search ===\n" << std::endl;

    // 1. Setup Engine and Catalog
    AstDynConfig cfg;
    cfg.ephemeris_type = EphemerisType::DE441;
    cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    cfg.verbose = false;

    AstDynEngine engine(cfg);
    
    // Initialize Gaia Catalog
    std::string gaia_config = "{\"catalog_type\":\"online_esa\",\"timeout_seconds\":30}";
    catalog::GaiaDR3Catalog::initialize(gaia_config);

    // 2. Search for Vesta (4)
    std::cout << "--- Searching Occultations for (4) Vesta ---" << std::endl;
    
    // Approximate elements for Vesta at epoch 2026-01-01
    auto vesta_elements = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_jd(2461041.5), // 2026-01-01
        2.3619, 0.0888, 7.14, 103.85, 151.20, 318.52, 
        physics::GravitationalParameter::sun()
    );

    time::EpochTDB start = time::EpochTDB::from_jd(2461122.0); // 2026-03-22 12:00
    time::EpochTDB end   = time::EpochTDB::from_jd(2461122.16); // 2026-03-22 16:00
    
    auto vesta_events = OccultationLogic::find_occultations(
        "4", vesta_elements, start, end, 12.0, engine, OccultationRefinementMode::ChebyshevDaily);

    std::cout << "Found " << vesta_events.size() << " potential occultations for Vesta." << std::endl;
    for (const auto& ev : vesta_events) {
        std::cout << "  - Star: " << ev.params.star_id << " (G=" << ev.params.star_mag << ")"
                  << " Time: " << std::fixed << std::setprecision(5) << ev.params.t_ca.mjd() << " MJD TDB"
                  << " Impact: " << ev.params.impact_parameter.to_km() << " km" << std::endl;
    }

/*
    // 3. Search for (50936) Nireus
...
    }
*/

    catalog::GaiaDR3Catalog::shutdown();
    return 0;
}
