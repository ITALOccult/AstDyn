/**
 * @file debug_vectors.cpp
 * @brief Heavy debug of vectors during propagation and reduction.
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

    physics::KeplerianStateTyped<core::ECLIPJ2000> elements;
    elements.a = physics::Distance::from_au(3.05657);
    elements.e = 0.15513;
    elements.i = astrometry::Angle::from_deg(17.3154);
    elements.node = astrometry::Angle::from_deg(280.1663);
    elements.omega = astrometry::Angle::from_deg(94.0826);
    elements.M = astrometry::Angle::from_deg(208.5820);
    elements.epoch = time::EpochTDB::from_jd(2461131.5); // MJD 61071.0
    elements.gm = physics::GravitationalParameter::sun();

    auto run_debug = [&](const std::string& label, time::EpochTDB t) {
        std::cout << "\nDEBUG [" << label << "] MJD: " << std::fixed << std::setprecision(6) << t.mjd() << std::endl;
        
        // 1. Earth Position (Heliocentric Ecliptic)
        auto de441 = std::make_shared<ephemeris::DE441Provider>(bsp_path);
        auto star_eq = de441->getPosition(ephemeris::CelestialBody::EARTH, t);
        auto sun_eq = de441->getPosition(ephemeris::CelestialBody::SUN, t);
        auto earth_h_eq = star_eq.to_eigen_si() - sun_eq.to_eigen_si();
        auto earth_h_ecl_vec = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(
            math::Vector3<core::GCRF, physics::Distance>::from_si(earth_h_eq.x(), earth_h_eq.y(), earth_h_eq.z())).to_eigen_si();
        
        std::cout << "  Earth Helio Ecl: " << earth_h_ecl_vec.transpose() / (constants::AU*1000.0) << " AU" << std::endl;

        // 2. Asteroid Position (Heliocentric Ecliptic)
        AstDynEngine engine(cfg);
        engine.set_initial_orbit(elements);
        auto ast_kep = engine.propagate_to(t);
        auto ast_cart = propagation::keplerian_to_cartesian(ast_kep);
        auto ast_h_ecl = ast_cart.position.to_eigen_si();

        std::cout << "  Asteroid Helio Ecl: " << ast_h_ecl.transpose() / (constants::AU*1000.0) << " AU" << std::endl;
        
        // 3. Observed Vector
        auto rho_ecl = ast_h_ecl - earth_h_ecl_vec;
        std::cout << "  Rho Ecl: " << rho_ecl.transpose() / (constants::AU*1000.0) << " AU" << std::endl;

        // 4. RA/Dec
        auto obs = astrometry::AstrometryReducer::compute_observation(elements, elements.epoch, t, cfg, astrometry::AstrometricSettings());
        if (obs) {
             std::cout << "  -> Final RA/Dec: " << (*obs).ra.value * constants::RAD_TO_DEG << ", " << (*obs).dec.value * constants::RAD_TO_DEG << std::endl;
        }
    };

    run_debug("EPOCH", elements.epoch);
    run_debug("EVENT", time::EpochTDB::from_mjd(61071.61));

    return 0;
}
