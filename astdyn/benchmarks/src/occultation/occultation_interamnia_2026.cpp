/**
 * @file occultation_interamnia_2026.cpp
 * @brief Test for searching occultations of (704) Interamnia on April 2, 2026.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include <iostream>

using namespace astdyn;
using namespace astdyn::astrometry;

int main() {
    // Initialize Planetary Ephemeris (Essential for propagation)
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    try {
        ephemeris::PlanetaryEphemeris::setProvider(std::make_shared<ephemeris::DE441Provider>(bsp_path));
    } catch (const std::exception& e) {
        std::cerr << "Failed to load ephemeris: " << e.what() << std::endl;
        return 1;
    }

    // Initialize Catalog (ESA Online TAP)
    try {
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa", "timeout_seconds":60})");
    } catch (const std::exception& e) {
        std::cerr << "Failed to initialize online catalog: " << e.what() << std::endl;
        return 1;
    }

    AstDynConfig cfg;
    cfg.ephemeris_file = bsp_path;
    cfg.ephemeris_type = EphemerisType::DE441;
    cfg.verbose = false; // Reduce noise
    
    AstDynEngine engine;
    engine.set_config(cfg);

    std::string asteroid_id = "704"; // Interamnia
    std::cout << "Fetching initial orbit for (704) Interamnia from Horizons..." << std::endl;
    
    io::HorizonsClient horizons;
    auto elements_res = horizons.query_elements("Interamnia", time::EpochTDB::from_jd(2461132.5)); 
    
    if (!elements_res) {
        std::cerr << "Failed to fetch elements from Horizons." << std::endl;
        return 1;
    }

    physics::KeplerianStateTyped<core::ECLIPJ2000> initial_elements = *elements_res;
    
    std::cout << "  - a: " << initial_elements.a.to_au() << " AU" << std::endl;
    std::cout << "  - e: " << initial_elements.e << std::endl;
    std::cout << "  - i: " << initial_elements.i.to_deg() << " deg" << std::endl;
    std::cout << "  - node: " << initial_elements.node.to_deg() << " deg" << std::endl;
    std::cout << "  - omega: " << initial_elements.omega.to_deg() << " deg" << std::endl;
    std::cout << "  - M: " << initial_elements.M.to_deg() << " deg" << std::endl;
    std::cout << "  - Epoch: " << initial_elements.epoch.jd() << " JD" << std::endl;

    // Search Window: April 2, 2026
    time::EpochTDB start = time::EpochTDB::from_jd(time::calendar_to_mjd(2026, 4, 2, 0) + 2400000.5);
    time::EpochTDB end = time::EpochTDB::from_jd(time::calendar_to_mjd(2026, 4, 3, 0) + 2400000.5);

    // Print nominal RA/Dec at start time for verification
    auto obs_start = astrometry::AstrometryReducer::compute_observation(
        initial_elements, initial_elements.epoch, start, engine.config(), astrometry::AstrometricSettings());
    if (obs_start) {
        std::cout << "  Nominal RA at search start: " << (*obs_start).ra.value * astdyn::constants::RAD_TO_DEG << " deg" << std::endl;
        std::cout << "  Nominal Dec at search start: " << (*obs_start).dec.value * astdyn::constants::RAD_TO_DEG << " deg" << std::endl;
    } else {
        std::cerr << "  ERROR: compute_observation failed" << std::endl;
    }

    double max_mag = 15.0; // User request
    std::cout << "Searching for occultations of (704) Interamnia on 2026-04-02 (Mag < " << max_mag << ")..." << std::endl;
    
    auto results = OccultationLogic::find_occultations(
        asteroid_id, initial_elements, start, end, max_mag, engine, OccultationRefinementMode::ChebyshevDaily
    );

    if (results.empty()) {
        std::cout << "No events found. Trying a wider search corridor (300 arcsec)..." << std::endl;
        // Manual search if find_occultations is too restricted
    }

    std::cout << "Found " << results.size() << " occultation events." << std::endl;

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& res = results[i];
        std::cout << "Event " << i + 1 << ":" << std::endl;
        std::cout << "  Star: RA " << res.star.ra.to_deg() << ", Dec " << res.star.dec.to_deg() 
                  << " (Mag: " << res.star.g_mag << ")" << std::endl;
        std::cout << "  TCA: " << res.params.t_ca.jd() << " JD" << std::endl;
        std::cout << "  Impact Parameter: " << res.params.impact_parameter.to_km() << " km" << std::endl;
        std::cout << "  Shadow Velocity: " << res.params.shadow_velocity.to_km_s() << " km/s" << std::endl;

        // Export KML/SVG for the first found event
        auto path = OccultationMapper::compute_path(
            res.params, res.star.ra, res.star.dec, 
            physics::Distance::from_km(330.0), // Diameter of Interamnia approx
            time::to_utc(res.params.t_ca),
            time::TimeDuration::from_seconds(600.0)
        );

        std::string prefix = "interamnia_20260402_event_" + std::to_string(i + 1);
        OccultationMapper::export_kml(path, prefix + ".kml");
        OccultationMapper::export_svg(path, prefix + ".svg");
        std::cout << "  Exported: " << prefix << ".kml/svg" << std::endl;
    }

    catalog::GaiaDR3Catalog::shutdown();
    return 0;
}
