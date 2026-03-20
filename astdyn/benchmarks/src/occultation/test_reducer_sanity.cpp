/**
 * @file test_reducer_sanity.cpp
 * @brief Sanity check for AstrometryReducer using verified elements.
 */

#include "astdyn/AstDyn.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    ephemeris::PlanetaryEphemeris::setProvider(std::make_shared<ephemeris::DE441Provider>(bsp_path));

    AstDynConfig cfg;
    cfg.ephemeris_file = bsp_path;
    cfg.ephemeris_type = EphemerisType::DE441;

    // Verified elements for Interamnia
    physics::KeplerianStateTyped<core::ECLIPJ2000> elements;
    elements.a = physics::Distance::from_au(3.05657);
    elements.e = 0.154491;
    elements.i = astrometry::Angle::from_deg(17.3063);
    elements.node = astrometry::Angle::from_deg(280.292);
    elements.omega = astrometry::Angle::from_deg(94.0826);
    elements.M = astrometry::Angle::from_deg(208.582);
    elements.epoch = time::EpochTDB::from_jd(2461132.5);
    elements.gm = physics::GravitationalParameter::sun();

    // Time: Same as epoch (No propagation)
    time::EpochTDB t_obs = elements.epoch;

    auto obs = astrometry::AstrometryReducer::compute_observation(
        elements, elements.epoch, t_obs, cfg, astrometry::AstrometricSettings());

    if (obs) {
        std::cout << "REDUCER TEST (Verified Elements, No Propagation):" << std::endl;
        std::cout << "  RA:  " << std::fixed << std::setprecision(5) << (*obs).ra.value * constants::RAD_TO_DEG << " deg" << std::endl;
        std::cout << "  Dec: " << (*obs).dec.value * constants::RAD_TO_DEG << " deg" << std::endl;
    }

    // Time: Epoch + 1 day
    t_obs = elements.epoch + time::TimeDuration::from_days(1.0);
    auto obs2 = astrometry::AstrometryReducer::compute_observation(
        elements, elements.epoch, t_obs, cfg, astrometry::AstrometricSettings());
    
    if (obs2) {
        std::cout << "REDUCER TEST (Epoch + 1 day):" << std::endl;
        std::cout << "  RA:  " << (*obs2).ra.value * constants::RAD_TO_DEG << " deg" << std::endl;
        std::cout << "  Dec: " << (*obs2).dec.value * constants::RAD_TO_DEG << " deg" << std::endl;
    }

    return 0;
}
