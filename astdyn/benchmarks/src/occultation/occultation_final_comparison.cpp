/**
 * @file occultation_final_comparison.cpp
 * @brief Formal comparison between AstDyn and provided XML data.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/OccultationXMLIO.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
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

    auto events = OccultationXMLIO::read_file("/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/interamnia_reference.xml");
    if (events.empty()) return 1;
    const auto& ref = events[0];

    std::cout << "=========================================================" << std::endl;
    std::cout << " ASTDYN VS OCCULT4 (Reference XML) - (704) Interamnia" << std::endl;
    std::cout << "=========================================================" << std::endl;

    // Orbit from XML
    physics::KeplerianStateTyped<core::ECLIPJ2000> elements;
    elements.a = physics::Distance::from_au(ref.semi_major_axis_au);
    elements.e = ref.eccentricity;
    elements.i = Angle::from_deg(ref.inclination_deg);
    elements.node = Angle::from_deg(ref.node_deg);
    elements.omega = Angle::from_deg(ref.mean_anomaly_deg); // 94.0826
    elements.M = Angle::from_deg(ref.ra_node_approx);      // 208.5820
    elements.gm = physics::GravitationalParameter::from_si(constants::GM_SUN * 1e9);
    elements.epoch = time::EpochTDB::from_jd(2461131.5); // 2026-04-02 00:00

    time::EpochTDB t_obs = time::EpochTDB::from_mjd(ref.mjd);

    auto obs_res = astrometry::AstrometryReducer::compute_observation(
        elements, elements.epoch, t_obs, engine.config(), astrometry::AstrometricSettings());

    if (!obs_res) return 1;
    const auto& obs = *obs_res;

    double star_ra = ref.ra_event_h * 15.0;
    double star_dec = ref.dec_event_deg;
    double ast_ra = obs.ra.value * constants::RAD_TO_DEG;
    double ast_dec = obs.dec.value * constants::RAD_TO_DEG;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Event Time (MJD): " << ref.mjd << std::endl;
    std::cout << "\nCoordinates at TCA:" << std::endl;
    std::cout << "  Source  |     RA (deg)    |     Dec (deg)   |" << std::endl;
    std::cout << "  --------|-----------------|-----------------|" << std::endl;
    std::cout << "  XML Star|    " << std::setw(12) << star_ra << " |    " << std::setw(12) << star_dec << " |" << std::endl;
    std::cout << "  AstDyn  |    " << std::setw(12) << ast_ra << " |    " << std::setw(12) << ast_dec << " |" << std::endl;

    double dra = (ast_ra - star_ra) * 3600.0 * std::cos(star_dec * constants::DEG_TO_RAD);
    double ddec = (ast_dec - star_dec) * 3600.0;

    std::cout << "\nSeparation AstDyn - XML Star:" << std::endl;
    std::cout << "  dRA*cos(dec): " << std::setw(10) << std::setprecision(3) << dra << " arcsec" << std::endl;
    std::cout << "  dDec:         " << std::setw(10) << std::setprecision(3) << ddec << " arcsec" << std::endl;

    // Geometry
    time::EpochTDB t2 = t_obs + time::TimeDuration::from_seconds(1.0);
    auto obs2 = *astrometry::AstrometryReducer::compute_observation(elements, elements.epoch, t2, engine.config(), astrometry::AstrometricSettings());
    Angle dra_dt = Angle::from_rad(obs2.ra.value - obs.ra.value);
    Angle ddec_dt = Angle::from_rad(obs2.dec.value - obs.dec.value);

    auto params = OccultationLogic::compute_parameters(
        RightAscension(Angle::from_deg(star_ra)), Declination(Angle::from_deg(star_dec)),
        RightAscension(Angle::from_rad(obs.ra.value)), Declination(Angle::from_rad(obs.dec.value)),
        physics::Distance::from_m(obs.distance.value),
        dra_dt, ddec_dt, physics::Velocity::zero());

    std::cout << "\nPath Geometry:" << std::endl;
    std::cout << "  Impact Parameter: " << std::setprecision(2) << params.impact_parameter.to_km() << " km" << std::endl;
    std::cout << "  Shadow Velocity:  " << params.shadow_velocity.to_km_s() << " km/s" << std::endl;
    
    if (params.impact_parameter.to_km() < 6378.0 + 330.0/2.0) {
        std::cout << "  Result: HIT (Shadow crosses Earth near reference)" << std::endl;
    } else {
        std::cout << "  Result: MISS (Check elements or earth geometry)" << std::endl;
    }
    
    std::cout << "=========================================================" << std::endl;

    return 0;
}
