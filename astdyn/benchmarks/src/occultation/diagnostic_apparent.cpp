/**
 * @file diagnostic_apparent.cpp
 * @brief Compare AstDyn Apparent position vs Horizons Apparent position.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

// Helper from my test
Eigen::Vector3d apply_aberration_test(const Eigen::Vector3d& rho, const Eigen::Vector3d& v_earth) {
    double r = rho.norm();
    Eigen::Vector3d u = rho / r;
    Eigen::Vector3d v_c = v_earth / (astdyn::constants::C_LIGHT * 1000.0);
    double d_uv = u.dot(v_c);
    double beta_inv = std::sqrt(1.0 - v_c.squaredNorm());
    return (beta_inv * u + (1.0 + d_uv / (1.0 + beta_inv)) * v_c) / (1.0 + d_uv);
}

int main() {
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    auto ephem = std::make_shared<ephemeris::DE441Provider>(bsp_path);
    ephemeris::PlanetaryEphemeris::setProvider(ephem);

    io::HorizonsClient horizons;
    time::EpochTDB t_tca = time::EpochTDB::from_jd(2461112.615); 

    std::cout << "--- Apparent Position Diagnostic: 50936 ---" << std::endl;

    auto h_geo = horizons.query_vectors("50936", t_tca, "@399");
    
    auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t_tca);
    
    double tau = 0.0;
    Eigen::Vector3d p_ast_ssb;
    for(int i=0; i<3; ++i) {
        time::EpochTDB t_emit = time::EpochTDB::from_jd(t_tca.jd() - tau);
        auto h_helio_emit = horizons.query_vectors("50936", t_emit, "@sun");
        auto sun_ssb_emit = ephem->getState(ephemeris::CelestialBody::SUN, t_emit);
        p_ast_ssb = h_helio_emit->position.to_eigen_si() + sun_ssb_emit.position.to_eigen_si();
        tau = (p_ast_ssb - earth_ssb.position.to_eigen_si()).norm() / (astdyn::constants::C_LIGHT * 1000.0 * 86400.0);
    }
    
    Eigen::Vector3d rho_geo = p_ast_ssb - earth_ssb.position.to_eigen_si();
    Eigen::Vector3d u_app = apply_aberration_test(rho_geo, earth_ssb.velocity.to_eigen_si());
    
    double ra_app = std::atan2(u_app.y(), u_app.x()) * constants::RAD_TO_DEG;
    double dec_app = std::asin(u_app.z()) * constants::RAD_TO_DEG;
    if (ra_app < 0) ra_app += 360.0;

    std::cout << "AstDyn Apparent RA:  " << std::fixed << std::setprecision(8) << ra_app << " deg" << std::endl;
    std::cout << "AstDyn Apparent Dec: " << std::fixed << std::setprecision(8) << dec_app << " deg" << std::endl;

    std::cout << "\nCorrections applied by AstDyn:" << std::endl;
    Eigen::Vector3d u_geo = h_geo->position.to_eigen_si().normalized();
    double ra_geo = std::atan2(u_geo.y(), u_geo.x()) * constants::RAD_TO_DEG;
    double dec_geo = std::asin(u_geo.z()) * constants::RAD_TO_DEG;
    if (ra_geo < 0) ra_geo += 360.0;
    
    std::cout << "Geometric RA:       " << ra_geo << " deg" << std::endl;
    std::cout << "Shift:              " << (ra_app - ra_geo)*3600.0 << " arcsec" << std::endl;

    return 0;
}
