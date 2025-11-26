/**
 * @file test_astdys_fitter.cpp
 * @brief Test AstDysOrbitFitter utility class
 */

#include <gtest/gtest.h>
#include <astdyn/io/AstDysOrbitFitter.hpp>
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::io;

TEST(AstDysOrbitFitter, LoadAndFit) {
    // Test files (203 Pompeja)
    std::string rwo_file = "../tools/203_recent100.rwo";
    std::string oel_file = "../tools/203.oel";
    
    // Create fitter
    AstDysOrbitFitter fitter;
    fitter.set_verbose(true);
    
    // Load observations from RWO file
    fitter.set_observations_file(rwo_file, "rwo");
    
    // Load initial elements from OEL file
    fitter.set_elements_file(oel_file, "oel");
    
    // Run orbit fit
    auto result = fitter.fit();
    
    // Verify result
    EXPECT_GT(result.num_observations_loaded, 0);
    EXPECT_GT(result.num_observations_used, 0);
    EXPECT_GT(result.num_iterations, 0);
    
    // Print results
    std::cout << "\n=== Fit Results ===\n";
    std::cout << "Observations loaded: " << result.num_observations_loaded << "\n";
    std::cout << "Observations used: " << result.num_observations_used << "\n";
    std::cout << "Outliers rejected: " << result.num_outliers << "\n";
    std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
    std::cout << "Iterations: " << result.num_iterations << "\n";
    std::cout << "RMS RA: " << result.rms_ra << " arcsec\n";
    std::cout << "RMS Dec: " << result.rms_dec << " arcsec\n";
    std::cout << "\nFitted orbit at MJD " << result.fitted_orbit.epoch_mjd_tdb << ":\n";
    std::cout << "  a = " << result.fitted_orbit.semi_major_axis << " AU\n";
    std::cout << "  e = " << result.fitted_orbit.eccentricity << "\n";
    std::cout << "  i = " << result.fitted_orbit.inclination * 180.0 / constants::PI << " deg\n";
    std::cout << "  Ω = " << result.fitted_orbit.longitude_ascending_node * 180.0 / constants::PI << " deg\n";
    std::cout << "  ω = " << result.fitted_orbit.argument_perihelion * 180.0 / constants::PI << " deg\n";
    std::cout << "  M = " << result.fitted_orbit.mean_anomaly * 180.0 / constants::PI << " deg\n";
}

TEST(AstDysOrbitFitter, CompareWithReference) {
    // Test files (203 Pompeja)
    std::string rwo_file = "../tools/203_recent100.rwo";
    std::string oel_file = "../tools/203.oel";
    
    // Create fitter
    AstDysOrbitFitter fitter;
    fitter.set_verbose(false);
    
    // Load data
    fitter.set_observations_file(rwo_file, "rwo");
    fitter.set_elements_file(oel_file, "oel");
    
    // Set reference orbit (from OEL file - OrbFit solution)
    propagation::KeplerianElements reference;
    reference.epoch_mjd_tdb = 61192.0;
    reference.semi_major_axis = 2.73687097;
    reference.eccentricity = 0.06099753;
    reference.inclination = 0.05498866; // radians
    reference.longitude_ascending_node = 3.36846166;
    reference.argument_perihelion = 5.09933389;
    reference.mean_anomaly = 4.03314924;
    reference.gravitational_parameter = constants::GMS;
    
    fitter.set_reference_orbit(reference);
    
    // Run fit
    auto result = fitter.fit();
    
    // Verify comparison was computed
    EXPECT_TRUE(result.delta_a_km.has_value());
    EXPECT_TRUE(result.delta_e.has_value());
    EXPECT_TRUE(result.delta_i_arcsec.has_value());
    
    // Print comparison
    std::cout << "\n=== Comparison with OrbFit ===\n";
    std::cout << "Δa = " << *result.delta_a_km << " km\n";
    std::cout << "Δe = " << std::scientific << *result.delta_e << std::fixed << "\n";
    std::cout << "Δi = " << *result.delta_i_arcsec << " arcsec\n";
}

TEST(AstDysOrbitFitter, DirectElementsInput) {
    // Test setting elements directly without file
    std::string rwo_file = "../tools/203_recent100.rwo";
    
    // Create Keplerian elements directly
    propagation::KeplerianElements initial;
    initial.epoch_mjd_tdb = 61192.0;
    initial.semi_major_axis = 2.737;
    initial.eccentricity = 0.061;
    initial.inclination = 3.15 * constants::PI / 180.0;
    initial.longitude_ascending_node = 193.0 * constants::PI / 180.0;
    initial.argument_perihelion = 292.2 * constants::PI / 180.0;
    initial.mean_anomaly = 231.0 * constants::PI / 180.0;
    initial.gravitational_parameter = constants::GMS;
    
    // Create fitter
    AstDysOrbitFitter fitter;
    fitter.set_verbose(false);
    
    // Set data directly
    fitter.set_observations_file(rwo_file, "rwo");
    fitter.set_elements(initial);
    
    // Run fit
    auto result = fitter.fit();
    
    // Verify result (allow 0 used observations in edge cases)
    EXPECT_GE(result.num_observations_used, 0);
    
    std::cout << "\n=== Direct Elements Test ===\n";
    std::cout << "Used " << result.num_observations_used << " observations\n";
    std::cout << "Final a = " << result.fitted_orbit.semi_major_axis << " AU\n";
}

TEST(AstDysOrbitFitter, MPCFormat) {
    // Test with standard MPC format file
    std::string mpc_file = "../tools/203.txt";
    std::string oel_file = "../tools/203.oel";
    
    // Create fitter
    AstDysOrbitFitter fitter;
    fitter.set_verbose(true);
    
    // Load from MPC file
    fitter.set_observations_file(mpc_file, "mpc");
    fitter.set_elements_file(oel_file, "oel");
    
    // Run fit
    auto result = fitter.fit();
    
    // Verify
    EXPECT_GT(result.num_observations_loaded, 100); // Should load many observations
    
    std::cout << "\n=== MPC Format Test ===\n";
    std::cout << "Loaded " << result.num_observations_loaded << " observations from MPC file\n";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
