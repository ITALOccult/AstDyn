/**
 * @file test_xml_logic.cpp
 * @brief Test loading elements from XML using the logic that worked in repro.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/OccultationXMLIO.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/time/epoch.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::astrometry;
using namespace astdyn::io;

int main() {
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    try {
        ephemeris::PlanetaryEphemeris::setProvider(std::make_shared<ephemeris::DE441Provider>(bsp_path));
    } catch (...) { return 1; }

    AstDynEngine engine;
    AstDynConfig cfg;
    cfg.ephemeris_file = bsp_path;
    cfg.ephemeris_type = EphemerisType::DE441;
    engine.set_config(cfg);

    // 1. Load XML
    auto events = OccultationXMLIO::read_file("/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/interamnia_reference.xml");
    if (events.empty()) return 1;
    const auto& ref = events[0];

    // 2. Map Elements exactly as in the working repro
    physics::KeplerianStateTyped<core::ECLIPJ2000> elements;
    elements.a = physics::Distance::from_au(ref.semi_major_axis_au);
    elements.e = ref.eccentricity;
    elements.i = Angle::from_deg(ref.inclination_deg);
    elements.node = Angle::from_deg(ref.node_deg);
    elements.omega = Angle::from_deg(ref.mean_anomaly_deg); // 94.0826
    elements.M = Angle::from_deg(ref.ra_node_approx);      // 208.5820
    elements.gm = physics::GravitationalParameter::from_si(constants::GM_SUN * 1e9);
    
    // Test Epoch: MJD 61071.0 (midnight)
    elements.epoch = time::EpochTDB::from_jd(2461131.5); 
    time::EpochTDB t_obs = elements.epoch;

    auto obs = astrometry::AstrometryReducer::compute_observation(
        elements, elements.epoch, t_obs, engine.config(), astrometry::AstrometricSettings());

    if (obs) {
        std::cout << "XML DATA TEST AT MJD 61071.0:" << std::endl;
        std::cout << "  RA:  " << std::fixed << std::setprecision(5) << (*obs).ra.value * constants::RAD_TO_DEG << " deg" << std::endl;
        std::cout << "  Dec: " << (*obs).dec.value * constants::RAD_TO_DEG << " deg" << std::endl;
        std::cout << "  Dist: " << (*obs).distance.value / 1000.0 / constants::AU << " AU" << std::endl;
    }

    return 0;
}
