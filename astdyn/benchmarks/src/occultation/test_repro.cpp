/**
 * @file test_repro.cpp
 * @brief Reproduction script for Interamnia comparison.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/time/epoch.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::astrometry;

int main() {
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    try {
        ephemeris::PlanetaryEphemeris::setProvider(std::make_shared<ephemeris::DE441Provider>(bsp_path));
    } catch (...) { return 1; }

    AstDynConfig cfg;
    cfg.ephemeris_file = bsp_path;
    cfg.ephemeris_type = EphemerisType::DE441;
    AstDynEngine engine;
    engine.set_config(cfg);

    // Hardcoded elements from Horizons (which we know work)
    physics::KeplerianStateTyped<core::ECLIPJ2000> elements;
    elements.a = physics::Distance::from_au(3.05657);
    elements.e = 0.154491;
    elements.i = Angle::from_deg(17.3063);
    elements.node = Angle::from_deg(280.292);
    elements.omega = Angle::from_deg(94.0826);
    elements.M = Angle::from_deg(208.582);
    elements.epoch = time::EpochTDB::from_jd(2461132.5);
    elements.gm = physics::GravitationalParameter::from_si(constants::GM_SUN * 1e9);

    time::EpochTDB t_obs = time::EpochTDB::from_jd(2461132.5);

    auto obs = astrometry::AstrometryReducer::compute_observation(
        elements, elements.epoch, t_obs, engine.config(), astrometry::AstrometricSettings());

    if (obs) {
        std::cout << "REPRO - RA:  " << std::fixed << std::setprecision(5) << (*obs).ra.value * constants::RAD_TO_DEG << " deg" << std::endl;
        std::cout << "REPRO - Dec: " << (*obs).dec.value * constants::RAD_TO_DEG << " deg" << std::endl;
    }

    return 0;
}
