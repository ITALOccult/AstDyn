/**
 * @file gaia_occultation_test.cpp
 * @brief Test for Gaia catalog integration via OrbitQuery.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::catalog;
using namespace astdyn::physics;
using namespace astdyn::astrometry;
using namespace astdyn::ephemeris;

int main() {
    std::string config_path = "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/external/ioc_gaialib/config_online.json";
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    
    std::cout << "--- Gaia/AstDyn beta 0.9 Integration Test ---" << std::endl;
    
    // 1. Initialize Catalog (Online Mode)
    try {
        GaiaDR3Catalog::initialize(config_path);
        std::cout << "✓ Gaia Catalog initialized (online_esa)" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "✗ Catalog init failed: " << e.what() << std::endl;
        return 1;
    }

    // 2. Setup Propagation for Asteroid 79518 (roughly 1 day window)
    time::EpochTDB t_start = time::to_tdb(time::EpochTT::from_mjd(61100.0));
    time::EpochTDB t_end   = time::to_tdb(time::EpochTT::from_mjd(61101.0));
    
    std::cout << "Window: MJD " << std::fixed << std::setprecision(6) << t_start.mjd() << " to " << t_end.mjd() << std::endl;

    // Fetch initial state from Horizons for 79518
    io::HorizonsClient horizons;
    auto orbit_res = horizons.query_vectors("79518", t_start, "@sun");
    if (!orbit_res) {
        std::cerr << "✗ Failed to fetch initial state" << std::endl;
        return 1;
    }
    auto initial_state = *orbit_res;

    // 3. Setup Ephemeris Provider (for Earth position)
    std::unique_ptr<DE441Provider> ephem_provider;
    try {
        ephem_provider = std::make_unique<DE441Provider>(bsp_path);
        std::cout << "✓ DE441 Provider loaded" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "✗ Failed to load DE441: " << e.what() << std::endl;
        return 1;
    }

    // 4. Propagate Orbit
    AstDynConfig cfg;
    cfg.ephemeris_file = bsp_path;
    cfg.propagator_settings.include_planets = true;
    cfg.propagator_settings.include_relativity = true;

    std::cout << "Propagating asteroid and Earth trajectory..." << std::endl;
    
    std::vector<CartesianStateTyped<core::GCRF>> body_states;
    std::vector<CartesianStateTyped<core::GCRF>> earth_states;
    
    double dt_days = 0.1; // 2.4 hours steps for Chebyshev fitting
    for (double jd = t_start.jd(); jd <= t_end.jd() + 0.001; jd += dt_days) {
        time::EpochTDB t = time::EpochTDB::from_jd(jd);
        
        // Body (Heliocentric GCRF)
        auto state = propagation::Propagator::propagate_cartesian(initial_state, t, cfg);
        body_states.push_back(state);
        
        // Earth (Heliocentric GCRF)
        auto earth_pos = ephem_provider->getPosition(CelestialBody::EARTH, t);
        auto earth_vel = ephem_provider->getVelocity(CelestialBody::EARTH, t);
        
        // Use factory method in CartesianStateTyped
        earth_states.push_back(CartesianStateTyped<core::GCRF>::from_si(
            t, earth_pos.x, earth_pos.y, earth_pos.z,
            earth_vel.x, earth_vel.y, earth_vel.z
        ));
    }

    // 5. Query Catalog
    std::cout << "Searching for stars near the orbit corridor (+/- 60 arcsec)..." << std::endl;
    
    try {
        auto candidates = find_occultation_candidates(
            GaiaDR3Catalog::instance(),
            body_states,
            earth_states,
            t_start,
            t_end,
            Angle::from_arcsec(60.0), // search corridor width
            15.5                      // magnitude limit
        );

        std::cout << "✓ Found " << candidates.size() << " occultation candidates." << std::endl;
        
        for (const auto& star : candidates) {
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << "Gaia ID: " << star.source_id << std::endl;
            std::cout << "RA:      " << star.ra.to_hms() << std::endl;
            std::cout << "Dec:     " << star.dec.to_dms() << std::endl;
            std::cout << "G Mag:   " << std::fixed << std::setprecision(2) << star.g_mag << std::endl;
            if (star.has_parallax()) {
                std::cout << "Dist:    " << star.parallax.to_distance().to_au() << " AU" << std::endl;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "✗ Query failed: " << e.what() << std::endl;
    }

    GaiaDR3Catalog::shutdown();
    return 0;
}
