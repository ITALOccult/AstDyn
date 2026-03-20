/**
 * @file test_manual_prop.cpp
 * @brief Manual propagation test to see what happens step by step.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/OccultationXMLIO.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    ephemeris::PlanetaryEphemeris::setProvider(std::make_shared<ephemeris::DE441Provider>(bsp_path));

    AstDynConfig cfg;
    cfg.ephemeris_file = bsp_path;
    cfg.ephemeris_type = EphemerisType::DE441;
    AstDynEngine engine;
    engine.set_config(cfg);

    physics::KeplerianStateTyped<core::ECLIPJ2000> elements;
    elements.a = physics::Distance::from_au(3.05657);
    elements.e = 0.15513;
    elements.i = astrometry::Angle::from_deg(17.3154);
    elements.node = astrometry::Angle::from_deg(280.1663);
    elements.omega = astrometry::Angle::from_deg(94.0826);
    elements.M = astrometry::Angle::from_deg(208.5820);
    elements.epoch = time::EpochTDB::from_jd(2461131.5);
    elements.gm = physics::GravitationalParameter::sun();

    engine.set_initial_orbit(elements);

    std::cout << "MANUAL PROPAGATION TEST:" << std::endl;
    for (double hour = 0; hour <= 24; hour += 4) {
        time::EpochTDB t = elements.epoch + time::TimeDuration::from_hours(hour);
        auto kep = engine.propagate_to(t);
        auto cart = propagation::keplerian_to_cartesian(kep);
        auto pos = cart.position.to_eigen_si();
        
        std::cout << "T + " << std::setw(2) << hour << "h: Helio r = " << (pos.norm() / (constants::AU * 1000.0)) << " AU" 
                  << "  M = " << kep.M.to_deg() << " deg" << std::endl;
    }

    return 0;
}
