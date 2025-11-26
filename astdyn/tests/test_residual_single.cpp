/**
 * @file test_residual_single.cpp
 * @brief Detailed debug of single residual calculation
 */

#include <gtest/gtest.h>
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::observations;
using namespace astdyn::orbit_determination;

TEST(ResidualSingle, DetailedDebug) {
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         Detailed Residual Calculation Debug               ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    // 1. Load ONE real observation from Pompeja
    std::cout << "1. Loading one observation from 203.rwo...\n";
    
    auto obs_list = RWOReader::readFile("tools/203_recent100.rwo");
    ASSERT_FALSE(obs_list.empty()) << "No observations loaded";
    
    const auto& obs = obs_list[0];
    
    std::cout << "   Parsed " << obs_list.size() << " observations total\n";
    std::cout << "   Observation details (first):\n";
    std::cout << "   • Object: " << obs.object_designation << "\n";
    std::cout << "   • MJD UTC: " << std::fixed << std::setprecision(6) << obs.mjd_utc << "\n";
    std::cout << "   • Observatory: " << obs.observatory_code << "\n";
    std::cout << "   • RA observed:  " << std::setprecision(8) << obs.ra << " rad = " 
              << std::setprecision(6) << obs.ra * constants::RAD_TO_DEG << " deg\n";
    std::cout << "   • Dec observed: " << std::setprecision(8) << obs.dec << " rad = " 
              << std::setprecision(6) << obs.dec * constants::RAD_TO_DEG << " deg\n";
    
    // Print also last observation
    if (obs_list.size() > 1) {
        const auto& last_obs = obs_list[obs_list.size() - 1];
        std::cout << "\n   Last observation:\n";
        std::cout << "   • MJD UTC: " << std::fixed << std::setprecision(6) << last_obs.mjd_utc << "\n";
        std::cout << "   • RA:  " << std::setprecision(6) << last_obs.ra * constants::RAD_TO_DEG << " deg\n";
        std::cout << "   • Dec: " << std::setprecision(6) << last_obs.dec * constants::RAD_TO_DEG << " deg\n";
    }
    std::cout << "   • Sigma RA:  " << obs.sigma_ra << " rad = " 
              << obs.sigma_ra * constants::RAD_TO_ARCSEC << " arcsec\n";
    std::cout << "   • Sigma Dec: " << obs.sigma_dec << " rad = " 
              << obs.sigma_dec * constants::RAD_TO_ARCSEC << " arcsec\n";
    
    // 2. Create orbit at observation epoch (from OrbFit)
    std::cout << "\n2. Creating Pompeja orbit (from OrbFit)...\n";
    
    // Parse OrbFit elements from .oel file  
    std::ifstream oel_file("tools/203_astdys.oel");
    ASSERT_TRUE(oel_file.is_open()) << "Cannot open .oel file";
    
    double a = 0, h = 0, k = 0, p = 0, q = 0, lambda = 0, mjd_ref = 0;
    std::string oel_line;
    while (std::getline(oel_file, oel_line)) {
        if (oel_line.empty() || oel_line[0] == '!') continue;
        std::istringstream iss(oel_line);
        std::string key;
        double value;
        if (iss >> key >> value) {
            if (key == "a") a = value;
            else if (key == "h") h = value;
            else if (key == "k") k = value;
            else if (key == "p") p = value;
            else if (key == "q") q = value;
            else if (key == "lambda") lambda = value;
            else if (key == "MJD") mjd_ref = value;
        }
    }
    
    // Convert equinoctial to Keplerian
    double e = std::sqrt(h*h + k*k);
    double i = 2.0 * std::atan(std::sqrt(p*p + q*q));
    double Omega = std::atan2(p, q);
    if (Omega < 0) Omega += 2.0 * constants::PI;
    double omega_plus_Omega = std::atan2(h, k);
    if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * constants::PI;
    double omega = omega_plus_Omega - Omega;
    if (omega < 0) omega += 2.0 * constants::PI;
    double M = lambda - omega_plus_Omega;
    while (M < 0) M += 2.0 * constants::PI;
    
    KeplerianElements kep;
    kep.epoch_mjd_tdb = mjd_ref;
    kep.semi_major_axis = a;
    kep.eccentricity = e;
    kep.inclination = i;
    kep.longitude_ascending_node = Omega;
    kep.argument_perihelion = omega;
    kep.mean_anomaly = M;
    kep.gravitational_parameter = constants::GMS;
    
    std::cout << "   OrbFit orbit at MJD " << kep.epoch_mjd_tdb << ":\n";
    std::cout << "   • a = " << kep.semi_major_axis << " AU\n";
    std::cout << "   • e = " << kep.eccentricity << "\n";
    std::cout << "   • i = " << kep.inclination * constants::RAD_TO_DEG << " deg\n";
    std::cout << "   • Ω = " << kep.longitude_ascending_node * constants::RAD_TO_DEG << " deg\n";
    std::cout << "   • ω = " << kep.argument_perihelion * constants::RAD_TO_DEG << " deg\n";
    std::cout << "   • M = " << kep.mean_anomaly * constants::RAD_TO_DEG << " deg\n";
    
    // 3. Convert to Cartesian
    auto cart = keplerian_to_cartesian(kep);
    
    std::cout << "\n   Cartesian state at MJD " << cart.epoch_mjd_tdb << ":\n";
    std::cout << "   • r = [" << std::setprecision(8) << cart.position.transpose() << "] AU\n";
    std::cout << "   • v = [" << cart.velocity.transpose() << "] AU/day\n";
    std::cout << "   • |r| = " << cart.position.norm() << " AU\n";
    std::cout << "   • |v| = " << cart.velocity.norm() << " AU/day\n";
    
    // 4. Propagate to observation epoch
    std::cout << "\n3. Propagating orbit to observation epoch...\n";
    
    double obs_mjd_tdb = ResidualCalculator::utc_to_tdb(obs.mjd_utc);
    std::cout << "   • Observation MJD TDB: " << obs_mjd_tdb << "\n";
    std::cout << "   • Time difference: " << (obs_mjd_tdb - cart.epoch_mjd_tdb) << " days\n";
    
    auto ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
    auto propagator = std::make_shared<Propagator>(
        std::make_unique<RKF78Integrator>(0.1, 1e-12),
        ephemeris,
        PropagatorSettings());
    
    auto cart_at_obs = propagator->propagate_cartesian(cart, obs_mjd_tdb);
    
    std::cout << "\n   Propagated state at observation epoch:\n";
    std::cout << "   • r = [" << cart_at_obs.position.transpose() << "] AU\n";
    std::cout << "   • v = [" << cart_at_obs.velocity.transpose() << "] AU/day\n";
    std::cout << "   • Δr from ref = " << (cart_at_obs.position - cart.position).norm() 
              << " AU\n";
    
    // 5. Get Earth position
    std::cout << "\n4. Computing observer position...\n";
    
    auto residual_calc = std::make_shared<ResidualCalculator>(ephemeris, propagator);
    
    auto observer_pos_opt = residual_calc->get_observer_position(obs);
    ASSERT_TRUE(observer_pos_opt.has_value()) << "Cannot get observer position";
    auto observer_pos = *observer_pos_opt;
    
    std::cout << "   • Observer (heliocentric): [" << observer_pos.transpose() << "] AU\n";
    std::cout << "   • |r_obs| = " << observer_pos.norm() << " AU\n";
    
    // 6. Compute topocentric position
    std::cout << "\n5. Computing topocentric position...\n";
    
    Vector3d rho = cart_at_obs.position - observer_pos;
    double range = rho.norm();
    Vector3d direction = rho / range;
    
    std::cout << "   • Topocentric vector: [" << rho.transpose() << "] AU\n";
    std::cout << "   • Range: " << range << " AU = " << range * constants::AU_TO_KM << " km\n";
    std::cout << "   • Direction (unit): [" << direction.transpose() << "]\n";
    
    // 7. Convert to RA/Dec
    std::cout << "\n6. Converting to RA/Dec...\n";
    
    double computed_ra = std::atan2(direction[1], direction[0]);
    if (computed_ra < 0) computed_ra += 2.0 * constants::PI;
    double computed_dec = std::asin(direction[2]);
    
    std::cout << "   • Computed RA:  " << computed_ra << " rad = " 
              << computed_ra * constants::RAD_TO_DEG << " deg = "
              << computed_ra * constants::RAD_TO_DEG * 24.0 / 360.0 << " hours\n";
    std::cout << "   • Computed Dec: " << computed_dec << " rad = " 
              << computed_dec * constants::RAD_TO_DEG << " deg\n";
    
    // 8. Compute residuals
    std::cout << "\n7. Computing residuals...\n";
    
    double res_ra = obs.ra - computed_ra;
    double res_dec = obs.dec - computed_dec;
    
    // Normalize RA residual by cos(dec)
    double res_ra_normalized = res_ra * std::cos(obs.dec);
    
    std::cout << "   • Residual RA (raw):  " << res_ra << " rad = " 
              << res_ra * constants::RAD_TO_ARCSEC << " arcsec\n";
    std::cout << "   • Residual RA (×cos dec): " << res_ra_normalized << " rad = " 
              << res_ra_normalized * constants::RAD_TO_ARCSEC << " arcsec\n";
    std::cout << "   • Residual Dec: " << res_dec << " rad = " 
              << res_dec * constants::RAD_TO_ARCSEC << " arcsec\n";
    
    // 9. Use ResidualCalculator
    std::cout << "\n8. Using ResidualCalculator::compute_residual()...\n";
    
    auto residual_opt = residual_calc->compute_residual(obs, cart);
    ASSERT_TRUE(residual_opt.has_value()) << "Residual computation failed";
    
    auto residual = *residual_opt;
    
    std::cout << "   • Computed RA:  " << residual.computed_ra * constants::RAD_TO_DEG << " deg\n";
    std::cout << "   • Computed Dec: " << residual.computed_dec * constants::RAD_TO_DEG << " deg\n";
    std::cout << "   • Residual RA:  " << residual.residual_ra * constants::RAD_TO_ARCSEC << " arcsec\n";
    std::cout << "   • Residual Dec: " << residual.residual_dec * constants::RAD_TO_ARCSEC << " arcsec\n";
    std::cout << "   • Normalized RA:  " << residual.normalized_ra << " sigma\n";
    std::cout << "   • Normalized Dec: " << residual.normalized_dec << " sigma\n";
    std::cout << "   • Chi²: " << residual.chi_squared << "\n";
    std::cout << "   • Range: " << residual.range << " AU\n";
    
    // 10. Analysis
    std::cout << "\n9. Analysis:\n";
    
    double res_total_arcsec = std::sqrt(
        residual.residual_ra * residual.residual_ra + 
        residual.residual_dec * residual.residual_dec
    ) * constants::RAD_TO_ARCSEC;
    
    std::cout << "   • Total residual: " << res_total_arcsec << " arcsec\n";
    
    if (res_total_arcsec > 100.0) {
        std::cout << "   ⚠️  RESIDUAL TOO LARGE! Expected < 100 arcsec for good orbit\n";
        std::cout << "   • Possible issues:\n";
        std::cout << "     - Orbit epoch far from observation\n";
        std::cout << "     - Wrong coordinate frame\n";
        std::cout << "     - Observer position error\n";
        std::cout << "     - Propagation error\n";
    } else if (res_total_arcsec > 10.0) {
        std::cout << "   ⚠️  Residual moderate (> 10 arcsec)\n";
    } else {
        std::cout << "   ✓ Residual reasonable (< 10 arcsec)\n";
    }
    
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                    DEBUG COMPLETE                          ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
