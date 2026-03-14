/**
 * @file diagnostic_50936.cpp
 * @brief Diagnostic to compare AstDyn Apparent position vs Horizons Apparent position.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    auto ephem = std::make_shared<ephemeris::DE441Provider>(bsp_path);
    ephemeris::PlanetaryEphemeris::setProvider(ephem);

    io::HorizonsClient horizons;
    
    // TCA from previous run: 2026-03-13 02:47:43.18 UT
    // JD: 2461112.616472
    time::EpochTDB t_tca = time::EpochTDB::from_jd(2461112.616472);

    std::cout << "--- Diagnostic Asteroid 50936 at 02:47:43.18 UT ---" << std::endl;

    auto horizons_geo = horizons.query_vectors("50936", t_tca, "@399"); // Geocentric Geometric
    auto horizons_helio = horizons.query_vectors("50936", t_tca, "@sun"); // Heliocentric Geometric
    
    // 2. AstDyn Propagated
    time::EpochTDB t_start = time::to_tdb(time::EpochUTC::from_mjd(61112.041667)); // 01:00 UT
    auto start_state = horizons.query_vectors("50936", t_start, "@sun");
    
    propagation::PropagatorSettings settings;
    settings.include_planets = true;
    settings.integrate_in_ecliptic = true;
    propagation::Propagator prop(std::make_shared<propagation::RKF78Integrator>(0.001), 
                               std::make_shared<ephemeris::PlanetaryEphemeris>(), settings);
    
    auto ast_final_helio = prop.propagate_cartesian(*start_state, t_tca);
    
    std::cout << "\n1. Geometric Comparison (at TCA):" << std::endl;
    std::cout << "Horizons Helio: " << horizons_helio->position.to_eigen_si().transpose() << std::endl;
    std::cout << "AstDyn   Helio: " << ast_final_helio.position.to_eigen_si().transpose() << std::endl;
    std::cout << "Diff (Helio):   " << (ast_final_helio.position.to_eigen_si() - horizons_helio->position.to_eigen_si()).norm() << " meters" << std::endl;

    // Geocentric Geometric
    auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t_tca);
    auto sun_ssb = ephem->getState(ephemeris::CelestialBody::SUN, t_tca);
    Eigen::Vector3d p_ast_geo_astdyn = (ast_final_helio.position.to_eigen_si() + sun_ssb.position.to_eigen_si()) - earth_ssb.position.to_eigen_si();
    
    std::cout << "Horizons Geo:   " << horizons_geo->position.to_eigen_si().transpose() << std::endl;
    std::cout << "AstDyn   Geo:   " << p_ast_geo_astdyn.transpose() << std::endl;
    std::cout << "Diff (Geo):     " << (p_ast_geo_astdyn - horizons_geo->position.to_eigen_si()).norm() << " meters" << std::endl;

    return 0;
}
