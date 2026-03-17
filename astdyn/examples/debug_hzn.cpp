/**
 * @file debug_hzn.cpp
 * @brief Comprehensive validation of all AstDyn APIs against JPL Horizons
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    std::string de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    auto provider = std::make_shared<ephemeris::DE441Provider>(de441_path);
    auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>(provider);

    io::HorizonsClient horizons;
    std::string asteroid = "21105";
    time::EpochTDB t_ref = time::EpochTDB::from_jd(2460320.5);

    auto obs_res = horizons.query_observation(asteroid, t_ref);
    if (obs_res) {
        std::cout << "HZN-OBS: RA=" << obs_res->ra.to_deg() << " deg\n";
        std::cout << "HZN-OBS: Dec=" << obs_res->dec.to_deg() << " deg\n";
    }

    auto vec_res = horizons.query_vectors(asteroid, t_ref);
    if (vec_res) {
        std::cout << "HZN-VEC: Pos=" << vec_res->position.to_eigen_si().transpose() << " m\n";
    }

    // Manual O-C with proper Heliocentric Earth
    if (obs_res && vec_res) {
        auto earth_ssb = ephem->getPosition(ephemeris::CelestialBody::EARTH, t_ref);
        auto sun_ssb = ephem->getPosition(ephemeris::CelestialBody::SUN, t_ref);
        
        Eigen::Vector3d earth_helio(earth_ssb.x_si(), earth_ssb.y_si(), earth_ssb.z_si());
        Eigen::Vector3d sun_helio(sun_ssb.x_si(), sun_ssb.y_si(), sun_ssb.z_si());
        
        Eigen::Vector3d earth_to_sun = sun_helio - earth_helio;
                                      
        Eigen::Vector3d rho = vec_res->position.to_eigen_si() - earth_helio;
        double ra_calc = std::atan2(rho.y(), rho.x());
        if (ra_calc < 0) ra_calc += constants::TWO_PI;
        double dec_calc = std::asin(rho.z() / rho.norm());

        std::cout << "ASTDYN-CALC: RA=" << ra_calc * constants::RAD_TO_DEG << " deg\n";
        std::cout << "ASTDYN-CALC: Dec=" << dec_calc * constants::RAD_TO_DEG << " deg\n";
        
        double dRA = (obs_res->ra.to_rad() - ra_calc) * constants::RAD_TO_DEG * 3600.0;
        double dDec = (obs_res->dec.to_rad() - dec_calc) * constants::RAD_TO_DEG * 3600.0;
        
        std::cout << "DIFF: dRA=" << dRA << " arcsec\n";
        std::cout << "DIFF: dDec=" << dDec << " arcsec\n";
    }

    return 0;
}
