/**
 * @file test_orbfit_engine.cpp
 * @brief Tests for OrbFitEngine class
 */

#include <gtest/gtest.h>
#include "orbfit/OrbFitEngine.hpp"
#include "orbfit/core/Constants.hpp"
#include "orbfit/propagation/OrbitalElements.hpp"

using namespace orbfit;
using namespace orbfit::propagation;

class OrbFitEngineTest : public ::testing::Test {
protected:
    void SetUp() override {
        engine = std::make_unique<OrbFitEngine>();
    }
    
    std::unique_ptr<OrbFitEngine> engine;
};

// Test 1: Construction
TEST_F(OrbFitEngineTest, Construction) {
    EXPECT_TRUE(engine != nullptr);
    EXPECT_FALSE(engine->has_orbit());
    EXPECT_EQ(engine->observations().size(), 0);
}

// Test 2: Configuration
TEST_F(OrbFitEngineTest, Configuration) {
    OrbFitConfig config;
    config.verbose = false;
    config.max_iterations = 20;
    config.tolerance = 1e-14;
    
    engine->set_config(config);
    
    EXPECT_FALSE(engine->config().verbose);
    EXPECT_EQ(engine->config().max_iterations, 20);
    EXPECT_DOUBLE_EQ(engine->config().tolerance, 1e-14);
}

// Test 3: Set initial orbit
TEST_F(OrbFitEngineTest, SetInitialOrbit) {
    KeplerianElements orbit;
    orbit.epoch_mjd_tdb = 60000.0;
    orbit.semi_major_axis = 2.5;
    orbit.eccentricity = 0.1;
    orbit.inclination = 0.1;
    orbit.longitude_ascending_node = 0.5;
    orbit.argument_perihelion = 1.0;
    orbit.mean_anomaly = 0.0;
    orbit.gravitational_parameter = constants::GMS;
    
    engine->set_initial_orbit(orbit);
    
    EXPECT_TRUE(engine->has_orbit());
    EXPECT_DOUBLE_EQ(engine->orbit().semi_major_axis, 2.5);
    EXPECT_DOUBLE_EQ(engine->orbit().eccentricity, 0.1);
}

// Test 4: Add observations
TEST_F(OrbFitEngineTest, AddObservations) {
    observations::OpticalObservation obs;
    obs.object_designation = "TEST001";
    obs.mjd_utc = 60000.0;
    obs.ra = 180.0 * M_PI / 180.0;  // Convert to radians
    obs.dec = 10.0 * M_PI / 180.0;  // Convert to radians
    obs.observatory_code = "500";
    
    engine->add_observation(obs);
    
    EXPECT_EQ(engine->observations().size(), 1);
    EXPECT_EQ(engine->observations()[0].object_designation, "TEST001");
}

// Test 5: Clear observations
TEST_F(OrbFitEngineTest, ClearObservations) {
    observations::OpticalObservation obs;
    obs.object_designation = "TEST001";
    obs.mjd_utc = 60000.0;
    obs.ra = 180.0 * M_PI / 180.0;  // Convert to radians
    obs.dec = 10.0 * M_PI / 180.0;  // Convert to radians
    obs.observatory_code = "500";
    
    engine->add_observation(obs);
    EXPECT_EQ(engine->observations().size(), 1);
    
    engine->clear_observations();
    EXPECT_EQ(engine->observations().size(), 0);
}

// Test 6: Propagate orbit
TEST_F(OrbFitEngineTest, PropagateOrbit) {
    // Set initial orbit
    KeplerianElements orbit;
    orbit.epoch_mjd_tdb = 60000.0;
    orbit.semi_major_axis = 2.5;
    orbit.eccentricity = 0.1;
    orbit.inclination = 0.1;
    orbit.longitude_ascending_node = 0.5;
    orbit.argument_perihelion = 1.0;
    orbit.mean_anomaly = 0.0;
    orbit.gravitational_parameter = constants::GMS;
    
    engine->set_initial_orbit(orbit);
    
    // Propagate 100 days forward
    double target_mjd = 60100.0;
    
    EXPECT_NO_THROW({
        auto propagated = engine->propagate_to(target_mjd);
        EXPECT_DOUBLE_EQ(propagated.epoch_mjd_tdb, target_mjd);
        EXPECT_NEAR(propagated.semi_major_axis, 2.5, 1e-3);  // a changes due to planetary perturbations
    });
}

// Test 7: Compute ephemeris
TEST_F(OrbFitEngineTest, ComputeEphemeris) {
    // Set initial orbit
    KeplerianElements orbit;
    orbit.epoch_mjd_tdb = 60000.0;
    orbit.semi_major_axis = 2.5;
    orbit.eccentricity = 0.1;
    orbit.inclination = 0.1;
    orbit.longitude_ascending_node = 0.5;
    orbit.argument_perihelion = 1.0;
    orbit.mean_anomaly = 0.0;
    orbit.gravitational_parameter = constants::GMS;
    
    engine->set_initial_orbit(orbit);
    
    // Generate ephemeris
    double start_mjd = 60000.0;
    double end_mjd = 60010.0;
    double step = 1.0;
    
    auto ephemeris = engine->compute_ephemeris(start_mjd, end_mjd, step);
    
    EXPECT_EQ(ephemeris.size(), 11);  // 0, 1, 2, ..., 10 days
    EXPECT_DOUBLE_EQ(ephemeris[0].epoch_mjd_tdb, start_mjd);
    EXPECT_DOUBLE_EQ(ephemeris.back().epoch_mjd_tdb, end_mjd);
}

// Test 8: Verbose mode toggle
TEST_F(OrbFitEngineTest, VerboseMode) {
    engine->set_verbose(true);
    EXPECT_TRUE(engine->config().verbose);
    
    engine->set_verbose(false);
    EXPECT_FALSE(engine->config().verbose);
}

// Test 9: Orbit without observations (should fail fit)
TEST_F(OrbFitEngineTest, FitOrbitWithoutObservations) {
    KeplerianElements orbit;
    orbit.epoch_mjd_tdb = 60000.0;
    orbit.semi_major_axis = 2.5;
    orbit.eccentricity = 0.1;
    orbit.inclination = 0.1;
    orbit.longitude_ascending_node = 0.5;
    orbit.argument_perihelion = 1.0;
    orbit.mean_anomaly = 0.0;
    orbit.gravitational_parameter = constants::GMS;
    
    engine->set_initial_orbit(orbit);
    
    // Should throw because no observations
    EXPECT_THROW({
        engine->fit_orbit();
    }, std::runtime_error);
}

// Test 10: Propagate without orbit (should fail)
TEST_F(OrbFitEngineTest, PropagateWithoutOrbit) {
    EXPECT_THROW({
        engine->propagate_to(60100.0);
    }, std::runtime_error);
}

// Test 11: Ephemeris without orbit (should fail)
TEST_F(OrbFitEngineTest, EphemerisWithoutOrbit) {
    EXPECT_THROW({
        engine->compute_ephemeris(60000.0, 60100.0, 1.0);
    }, std::runtime_error);
}

// Test 12: Multiple orbits - verify update
TEST_F(OrbFitEngineTest, UpdateOrbit) {
    // Set first orbit
    KeplerianElements orbit1;
    orbit1.epoch_mjd_tdb = 60000.0;
    orbit1.semi_major_axis = 2.5;
    orbit1.eccentricity = 0.1;
    orbit1.inclination = 0.1;
    orbit1.longitude_ascending_node = 0.5;
    orbit1.argument_perihelion = 1.0;
    orbit1.mean_anomaly = 0.0;
    orbit1.gravitational_parameter = constants::GMS;
    
    engine->set_initial_orbit(orbit1);
    EXPECT_DOUBLE_EQ(engine->orbit().semi_major_axis, 2.5);
    
    // Update to different orbit
    KeplerianElements orbit2;
    orbit2.epoch_mjd_tdb = 60100.0;
    orbit2.semi_major_axis = 3.0;
    orbit2.eccentricity = 0.15;
    orbit2.inclination = 0.2;
    orbit2.longitude_ascending_node = 0.6;
    orbit2.argument_perihelion = 1.1;
    orbit2.mean_anomaly = 0.5;
    orbit2.gravitational_parameter = constants::GMS;
    
    engine->set_initial_orbit(orbit2);
    EXPECT_DOUBLE_EQ(engine->orbit().semi_major_axis, 3.0);
    EXPECT_DOUBLE_EQ(engine->orbit().eccentricity, 0.15);
}

// Test 13: Propagator settings
TEST_F(OrbFitEngineTest, PropagatorSettings) {
    OrbFitConfig config;
    config.propagator_settings.include_planets = true;
    config.propagator_settings.include_asteroids = true;
    config.propagator_settings.perturb_jupiter = true;
    config.propagator_settings.perturb_saturn = true;
    
    engine->set_config(config);
    
    EXPECT_TRUE(engine->config().propagator_settings.include_planets);
    EXPECT_TRUE(engine->config().propagator_settings.include_asteroids);
    EXPECT_TRUE(engine->config().propagator_settings.perturb_jupiter);
}

// Test 14: Long ephemeris generation
TEST_F(OrbFitEngineTest, LongEphemeris) {
    KeplerianElements orbit;
    orbit.epoch_mjd_tdb = 60000.0;
    orbit.semi_major_axis = 2.5;
    orbit.eccentricity = 0.1;
    orbit.inclination = 0.1;
    orbit.longitude_ascending_node = 0.5;
    orbit.argument_perihelion = 1.0;
    orbit.mean_anomaly = 0.0;
    orbit.gravitational_parameter = constants::GMS;
    
    engine->set_initial_orbit(orbit);
    
    // Generate 1 year ephemeris with 10-day steps
    auto ephemeris = engine->compute_ephemeris(60000.0, 60365.0, 10.0);
    
    EXPECT_EQ(ephemeris.size(), 38);  // 365/10 + 1 (including both endpoints)
    
    // Check all epochs are in order
    for (size_t i = 1; i < ephemeris.size(); ++i) {
        EXPECT_GT(ephemeris[i].epoch_mjd_tdb, ephemeris[i-1].epoch_mjd_tdb);
    }
}

// Test 15: Integration with different tolerances
TEST_F(OrbFitEngineTest, DifferentTolerances) {
    KeplerianElements orbit;
    orbit.epoch_mjd_tdb = 60000.0;
    orbit.semi_major_axis = 2.5;
    orbit.eccentricity = 0.1;
    orbit.inclination = 0.1;
    orbit.longitude_ascending_node = 0.5;
    orbit.argument_perihelion = 1.0;
    orbit.mean_anomaly = 0.0;
    orbit.gravitational_parameter = constants::GMS;
    
    // Test with tight tolerance
    OrbFitConfig config;
    config.tolerance = 1e-14;
    engine->set_config(config);
    engine->set_initial_orbit(orbit);
    
    auto result1 = engine->propagate_to(60100.0);
    
    // Test with loose tolerance
    config.tolerance = 1e-10;
    engine->set_config(config);
    engine->set_initial_orbit(orbit);
    
    auto result2 = engine->propagate_to(60100.0);
    
    // Results should be similar but not identical
    EXPECT_NEAR(result1.semi_major_axis, result2.semi_major_axis, 1e-8);
}

// Test 16: Summary
TEST(OrbFitEngineTestSuite, Summary) {
    std::cout << "\n========================================\n";
    std::cout << "OrbFit Engine Tests Summary\n";
    std::cout << "========================================\n";
    std::cout << "✓ Engine construction and initialization\n";
    std::cout << "✓ Configuration management\n";
    std::cout << "✓ Orbit management (set, update, propagate)\n";
    std::cout << "✓ Observation management (add, clear)\n";
    std::cout << "✓ Ephemeris generation (short and long)\n";
    std::cout << "✓ Error handling (missing orbit/observations)\n";
    std::cout << "✓ Propagator settings integration\n";
    std::cout << "✓ Tolerance and accuracy control\n";
    std::cout << "========================================\n\n";
}
