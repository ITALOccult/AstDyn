/**
 * @file occultation_comparison_test.cpp
 * @brief Benchmark test to compare AstDyn calculations with user-provided XML (Occult4).
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/OccultationXMLIO.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::astrometry;
using namespace astdyn::io;

int main() {
    // 1. Setup Environment
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    try {
        ephemeris::PlanetaryEphemeris::setProvider(std::make_shared<ephemeris::DE441Provider>(bsp_path));
    } catch (...) {
        std::cerr << "Failed to load ephemeris." << std::endl;
        return 1;
    }

    AstDynConfig cfg;
    cfg.ephemeris_file = bsp_path;
    cfg.ephemeris_type = EphemerisType::DE441;
    AstDynEngine engine;
    engine.set_config(cfg);

    // 2. Load Reference XML
    auto events = OccultationXMLIO::read_file("/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/interamnia_reference.xml");
    if (events.empty()) {
        std::cerr << "Could not load interamnia_reference.xml" << std::endl;
        return 1;
    }

    const auto& ref = events[0];
    std::cout << "--- COMPARING ASTDYN VS XML (REFERENCE ID: " << ref.event_id << ") ---" << std::endl;

    // 3. Setup Orbit from XML
    physics::KeplerianStateTyped<core::ECLIPJ2000> elements;
    elements.a = physics::Distance::from_au(ref.semi_major_axis_au);
    elements.e = ref.eccentricity;
    elements.i = Angle::from_deg(ref.inclination_deg);
    elements.node = Angle::from_deg(ref.node_deg);
    
    // Correction: In XML/Occult4, Orbit[1] is M (208.582) and Orbit[5] is Omega (94.0826)
    // OccultationXMLIO maps: Orbit[1] -> ra_node_approx, Orbit[5] -> mean_anomaly_deg
    elements.M = Angle::from_deg(ref.ra_node_approx);      // 208.5820
    elements.omega = Angle::from_deg(ref.mean_anomaly_deg); // 94.0826
    
    // IMPORTANT: Must set Solar GM for asteroid orbits
    elements.gm = physics::GravitationalParameter::from_si(astdyn::constants::GM_SUN * 1e9); 
    
    elements.epoch = time::EpochTDB::from_jd(time::calendar_to_mjd(ref.epoch_year, ref.epoch_month, ref.epoch_day) + 2400000.5);
    
    std::cout << "Astrometry Reference Time (MJD): " << std::fixed << std::setprecision(5) << ref.mjd << std::endl;

    // 4. Calculate AstDyn Prediction at EPOCH of elements (MJD 61071.0)
    auto obs_t0 = astrometry::AstrometryReducer::compute_observation(
        elements, elements.epoch, elements.epoch, engine.config(), astrometry::AstrometricSettings());
    
    if (obs_t0) {
        std::cout << "\n0. TEST AT EPOCH (MJD 61071.0):" << std::endl;
        std::cout << "   - AstDyn RA:  " << (*obs_t0).ra.value * constants::RAD_TO_DEG << " deg" << std::endl;
        std::cout << "   - AstDyn Dec: " << (*obs_t0).dec.value * constants::RAD_TO_DEG << " deg" << std::endl;
    }

    // 4b. Calculate AstDyn Prediction at EVENT time
    time::EpochTDB event_time = time::EpochTDB::from_mjd(ref.mjd);
    auto obs_res = astrometry::AstrometryReducer::compute_observation(
        elements, elements.epoch, event_time, engine.config(), astrometry::AstrometricSettings());

    if (!obs_res) {
        std::cerr << "compute_observation failed!" << std::endl;
        return 1;
    }
    const auto& obs = *obs_res;

    // Star coordinates from XML
    double star_ra_deg = ref.ra_event_h * 15.0;
    double star_dec_deg = ref.dec_event_deg;
    
    double ast_ra_deg = obs.ra.value * constants::RAD_TO_DEG;
    double ast_dec_deg = obs.dec.value * constants::RAD_TO_DEG;

    std::cout << "\n1. ASTROMETRY COMPARISON (EPOCH OF EVENT):" << std::endl;
    std::cout << "   - XML Star RA:  " << star_ra_deg << " deg" << std::endl;
    std::cout << "   - XML Star Dec: " << star_dec_deg << " deg" << std::endl;
    std::cout << "   - AstDyn RA:    " << ast_ra_deg << " deg" << std::endl;
    std::cout << "   - AstDyn Dec:   " << ast_dec_deg << " deg" << std::endl;
    
    double dra = (ast_ra_deg - star_ra_deg) * 3600.0 * std::cos(star_dec_deg * constants::DEG_TO_RAD);
    double ddec = (ast_dec_deg - star_dec_deg) * 3600.0;
    std::cout << "   - Separation (RA):  " << dra << " arcsec" << std::endl;
    std::cout << "   - Separation (Dec): " << ddec << " arcsec" << std::endl;

    // 5. Calculate Shadow Parameters
    // We need velocities for compute_parameters
    time::EpochTDB t2 = event_time + time::TimeDuration::from_seconds(1.0);
    auto obs2_res = astrometry::AstrometryReducer::compute_observation(
        elements, elements.epoch, t2, engine.config(), astrometry::AstrometricSettings());
    const auto& obs2 = *obs2_res;

    Angle dra_dt = Angle::from_rad(obs2.ra.value - obs.ra.value);
    Angle ddec_dt = Angle::from_rad(obs2.dec.value - obs.dec.value);

    auto params = OccultationLogic::compute_parameters(
        RightAscension(Angle::from_deg(star_ra_deg)), Declination(Angle::from_deg(star_dec_deg)),
        RightAscension(Angle::from_rad(obs.ra.value)), Declination(Angle::from_rad(obs.dec.value)),
        physics::Distance::from_m(obs.distance.value),
        dra_dt, ddec_dt, physics::Velocity::zero());

    std::cout << "\n2. SHADOW PATH PARAMETERS:" << std::endl;
    std::cout << "   - XML Lat/Lon: " << ref.latitude_deg << ", " << ref.longitude_deg << std::endl;
    std::cout << "   - AstDyn Impact Parameter: " << params.impact_parameter.to_km() << " km" << std::endl;
    
    // The impact parameter should be < 6371 (Earth radius) for a hit.
    if (params.impact_parameter.to_km() < 6400.0) {
        std::cout << "   - Prediction: HIT (Potential Occultation Found)" << std::endl;
    } else {
        std::cout << "   - Prediction: MISS (Asteroid misses Earth or is off-track)" << std::endl;
    }

    return 0;
}
