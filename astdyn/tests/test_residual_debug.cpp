/**
 * @file test_residual_debug.cpp
 * @brief Debug test for residual calculation
 */

#include <gtest/gtest.h>
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::observations;
using namespace astdyn::orbit_determination;

TEST(ResidualDebug, SimpleTest) {
    // Create a simple orbit (Pompeja-like)
    KeplerianElements kep;
    kep.epoch_mjd_tdb = 60800.0;  // Mid-2025
    kep.semi_major_axis = 2.737;
    kep.eccentricity = 0.061;
    kep.inclination = 3.15 * constants::DEG_TO_RAD;
    kep.longitude_ascending_node = 193.0 * constants::DEG_TO_RAD;
    kep.argument_perihelion = 292.0 * constants::DEG_TO_RAD;
    kep.mean_anomaly = 230.0 * constants::DEG_TO_RAD;
    kep.gravitational_parameter = constants::GMS;
    
    // Convert to Cartesian
    auto cart = keplerian_to_cartesian(kep);
    
    std::cout << "Initial state at MJD " << cart.epoch_mjd_tdb << ":\n";
    std::cout << "  r = [" << cart.position.transpose() << "] AU\n";
    std::cout << "  v = [" << cart.velocity.transpose() << "] AU/day\n";
    
    // Create propagator
    auto ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
    auto propagator = std::make_shared<Propagator>(
        std::make_unique<RKF78Integrator>(0.1, 1e-12),
        ephemeris,
        PropagatorSettings());
    
    // Create residual calculator WITH propagator
    auto residual_calc = std::make_shared<ResidualCalculator>(ephemeris, propagator);
    
    // Create a synthetic observation at same epoch
    OpticalObservation obs;
    obs.mjd_utc = 60800.0;  // Same as orbit epoch
    obs.object_designation = "203";
    obs.observatory_code = "500";  // Geocenter
    
    // Compute expected RA/Dec from state
    // For now, just set rough values
    obs.ra = 3.0;      // ~172 deg
    obs.dec = -0.05;   // ~-3 deg
    obs.sigma_ra = 1.0 / 206265.0;   // 1 arcsec
    obs.sigma_dec = 1.0 / 206265.0;
    
    // Compute residual
    auto residual_opt = residual_calc->compute_residual(obs, cart);
    
    ASSERT_TRUE(residual_opt.has_value()) << "Residual computation failed";
    
    auto residual = *residual_opt;
    
    std::cout << "\nResidual at same epoch:\n";
    std::cout << "  Computed RA:  " << residual.computed_ra << " rad (" 
              << residual.computed_ra * constants::RAD_TO_DEG << " deg)\n";
    std::cout << "  Computed Dec: " << residual.computed_dec << " rad (" 
              << residual.computed_dec * constants::RAD_TO_DEG << " deg)\n";
    std::cout << "  Residual RA:  " << residual.residual_ra << " rad (" 
              << residual.residual_ra * constants::RAD_TO_ARCSEC << " arcsec)\n";
    std::cout << "  Residual Dec: " << residual.residual_dec << " rad (" 
              << residual.residual_dec * constants::RAD_TO_ARCSEC << " arcsec)\n";
    std::cout << "  Range: " << residual.range << " AU\n";
    
    // Test propagation to different epoch
    double target_mjd = 60850.0;  // 50 days later
    std::cout << "\nPropagating to MJD " << target_mjd << "...\n";
    
    auto cart_propagated = propagator->propagate_cartesian(cart, target_mjd);
    
    std::cout << "  r = [" << cart_propagated.position.transpose() << "] AU\n";
    std::cout << "  v = [" << cart_propagated.velocity.transpose() << "] AU/day\n";
    std::cout << "  Δr = " << (cart_propagated.position - cart.position).norm() << " AU\n";
    
    // Test residual computation with propagation
    obs.mjd_utc = target_mjd;
    auto residual_prop = residual_calc->compute_residual(obs, cart);  // Using original state
    
    if (residual_prop) {
        std::cout << "\nResidual with automatic propagation:\n";
        std::cout << "  Computed RA:  " << residual_prop->computed_ra * constants::RAD_TO_DEG << " deg\n";
        std::cout << "  Computed Dec: " << residual_prop->computed_dec * constants::RAD_TO_DEG << " deg\n";
        std::cout << "  Residual RA:  " << residual_prop->residual_ra * constants::RAD_TO_ARCSEC << " arcsec\n";
        std::cout << "  Residual Dec: " << residual_prop->residual_dec * constants::RAD_TO_ARCSEC << " arcsec\n";
    }
    
    // Verify residuals are reasonable (< 1 degree for same epoch)
    EXPECT_LT(std::abs(residual.residual_ra), 0.02) << "RA residual too large";
    EXPECT_LT(std::abs(residual.residual_dec), 0.02) << "Dec residual too large";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
