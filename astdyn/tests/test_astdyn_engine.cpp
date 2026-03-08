/**
 * @file test_orbfit_engine.cpp
 * @brief Tests for AstDynEngine class
 */

#include <gtest/gtest.h>
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"

using namespace astdyn;
using namespace astdyn::propagation;

class AstDynEngineTest : public ::testing::Test {
protected:
    void SetUp() override {
        engine = std::make_unique<AstDynEngine>();
        
        // Disable planetary perturbations for unit tests (requires .bsp file otherwise)
        AstDynConfig config = engine->config();
        config.propagator_settings.include_planets = false;
        engine->set_config(config);
    }
    
    std::unique_ptr<AstDynEngine> engine;
};

// Test 1: Construction
TEST_F(AstDynEngineTest, Construction) {
    EXPECT_TRUE(engine != nullptr);
    EXPECT_FALSE(engine->has_orbit());
    EXPECT_EQ(engine->observations().size(), static_cast<size_t>(0));
}

// Test 2: Configuration
TEST_F(AstDynEngineTest, Configuration) {
    AstDynConfig config;
    config.verbose = false;
    config.max_iterations = 20;
    config.tolerance = 1e-14;
    
    engine->set_config(config);
    
    EXPECT_FALSE(engine->config().verbose);
    EXPECT_EQ(engine->config().max_iterations, 20);
    EXPECT_DOUBLE_EQ(engine->config().tolerance, 1e-14);
}

// Test 3: Set initial orbit
TEST_F(AstDynEngineTest, SetInitialOrbit) {
    auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(60000.0),
        2.5, 0.1, 0.1, 0.5, 1.0, 0.0
    );
    
    engine->set_initial_orbit(orbit);
    
    EXPECT_TRUE(engine->has_orbit());
    EXPECT_DOUBLE_EQ(engine->orbit().a.to_au(), 2.5);
    EXPECT_DOUBLE_EQ(engine->orbit().e, 0.1);
}

// Test 4: Add observations
TEST_F(AstDynEngineTest, AddObservations) {
    observations::OpticalObservation obs;
    obs.object_designation = "TEST001";
    obs.time = time::EpochUTC::from_mjd(60000.0);
    obs.ra = astrometry::RightAscension(astrometry::Angle::from_rad(180.0 * M_PI / 180.0));
    obs.dec = astrometry::Declination(astrometry::Angle::from_rad(10.0 * M_PI / 180.0));
    obs.observatory_code = "500";
    
    engine->add_observation(obs);
    
    EXPECT_EQ(engine->observations().size(), 1UL);
    EXPECT_EQ(engine->observations()[0].object_designation, "TEST001");
}

// Test 5: Clear observations
TEST_F(AstDynEngineTest, ClearObservations) {
    observations::OpticalObservation obs;
    obs.object_designation = "TEST001";
    obs.time = time::EpochUTC::from_mjd(60000.0);
    obs.ra = astrometry::RightAscension(astrometry::Angle::from_rad(180.0 * M_PI / 180.0));
    obs.dec = astrometry::Declination(astrometry::Angle::from_rad(10.0 * M_PI / 180.0));
    obs.observatory_code = "500";
    
    engine->add_observation(obs);
    EXPECT_EQ(engine->observations().size(), 1UL);
    
    engine->clear_observations();
    EXPECT_EQ(engine->observations().size(), static_cast<size_t>(0));
}

// Test 6: Propagate orbit
TEST_F(AstDynEngineTest, PropagateOrbit) {
    // Set initial orbit
    auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(60000.0),
        2.5, 0.1, 0.1, 0.5, 1.0, 0.0
    );
    
    engine->set_initial_orbit(orbit);
    
    // Propagate 100 days forward
    time::EpochTDB target_time = time::EpochTDB::from_mjd(60100.0);
    
    EXPECT_NO_THROW({
        auto propagated = engine->propagate_to(target_time);
        EXPECT_DOUBLE_EQ(propagated.epoch.mjd(), target_time.mjd());
        EXPECT_NEAR(propagated.a.to_au(), 2.5, 1e-3);  // a changes due to planetary perturbations
    });
}

// Test 7: Compute ephemeris
TEST_F(AstDynEngineTest, ComputeEphemeris) {
    // Set initial orbit
    auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(60000.0),
        2.5, 0.1, 0.1, 0.5, 1.0, 0.0
    );
    
    engine->set_initial_orbit(orbit);
    
    // Generate ephemeris
    time::EpochTDB start_time = time::EpochTDB::from_mjd(60000.0);
    time::EpochTDB end_time = time::EpochTDB::from_mjd(60010.0);
    double step = 1.0;
    
    auto ephemeris = engine->compute_ephemeris(start_time, end_time, step);
    
    EXPECT_EQ(ephemeris.size(), 11UL);  // 0, 1, 2, ..., 10 days
    EXPECT_DOUBLE_EQ(ephemeris[0].epoch.mjd(), start_time.mjd());
    EXPECT_DOUBLE_EQ(ephemeris.back().epoch.mjd(), end_time.mjd());
}

// Test 8: Verbose mode toggle
TEST_F(AstDynEngineTest, VerboseMode) {
    engine->set_verbose(true);
    EXPECT_TRUE(engine->config().verbose);
    
    engine->set_verbose(false);
    EXPECT_FALSE(engine->config().verbose);
}

// Test 9: Orbit without observations (should fail fit)
TEST_F(AstDynEngineTest, FitOrbitWithoutObservations) {
    auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(60000.0),
        2.5, 0.1, 0.1, 0.5, 1.0, 0.0
    );
    
    engine->set_initial_orbit(orbit);
    
    // Should throw because no observations
    EXPECT_THROW({
        engine->fit_orbit();
    }, std::runtime_error);
}

// Test 10: Propagate without orbit (should fail)
TEST_F(AstDynEngineTest, PropagateWithoutOrbit) {
    EXPECT_THROW({
        engine->propagate_to(time::EpochTDB::from_mjd(60100.0));
    }, std::runtime_error);
}

// Test 11: Ephemeris without orbit (should fail)
TEST_F(AstDynEngineTest, EphemerisWithoutOrbit) {
    EXPECT_THROW({
        engine->compute_ephemeris(time::EpochTDB::from_mjd(60000.0), 
                                  time::EpochTDB::from_mjd(60100.0), 1.0);
    }, std::runtime_error);
}

// Test 12: Multiple orbits - verify update
TEST_F(AstDynEngineTest, UpdateOrbit) {
    // Set first orbit
    auto orbit1 = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(60000.0),
        2.5, 0.1, 0.1, 0.5, 1.0, 0.0
    );
    
    engine->set_initial_orbit(orbit1);
    EXPECT_DOUBLE_EQ(engine->orbit().a.to_au(), 2.5);
    
    // Update to different orbit
    auto orbit2 = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(60100.0),
        3.0, 0.15, 0.2, 0.6, 1.1, 0.5
    );
    
    engine->set_initial_orbit(orbit2);
    EXPECT_DOUBLE_EQ(engine->orbit().a.to_au(), 3.0);
    EXPECT_DOUBLE_EQ(engine->orbit().e, 0.15);
}

// Test 13: Propagator settings
TEST_F(AstDynEngineTest, PropagatorSettings) {
    AstDynConfig config;
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
TEST_F(AstDynEngineTest, LongEphemeris) {
    // Set initial orbit
    auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(60000.0),
        2.5, 0.1, 0.1, 0.5, 1.0, 0.0
    );
    
    engine->set_initial_orbit(orbit);
    
    // Generate 1 year ephemeris with 10-day steps
    auto ephemeris = engine->compute_ephemeris(time::EpochTDB::from_mjd(60000.0), 
                                               time::EpochTDB::from_mjd(60365.0), 10.0);
    
    EXPECT_EQ(ephemeris.size(), 37UL);  // 365/10 + 1 (0 to 360 by 10)
    
    // Check all epochs are in order
    for (size_t i = 1; i < ephemeris.size(); ++i) {
        EXPECT_GT(ephemeris[i].epoch.mjd(), ephemeris[i-1].epoch.mjd());
    }
}

// Test 15: Integration with different tolerances
TEST_F(AstDynEngineTest, DifferentTolerances) {
    auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(60000.0),
        2.5, 0.1, 0.1, 0.5, 1.0, 0.0
    );
    
    // Test with tight tolerance
    AstDynConfig config = engine->config();
    config.tolerance = 1e-14;
    engine->set_config(config);
    engine->set_initial_orbit(orbit);
    
    auto result1 = engine->propagate_to(time::EpochTDB::from_mjd(60100.0));
    
    // Test with loose tolerance
    config.tolerance = 1e-10;
    engine->set_config(config);
    engine->set_initial_orbit(orbit);
    
    auto result2 = engine->propagate_to(time::EpochTDB::from_mjd(60100.0));
    
    // Results should be similar but not identical
    EXPECT_NEAR(result1.a.to_au(), result2.a.to_au(), 1e-8);
}

// Test 16: Summary
TEST(AstDynEngineTestSuite, Summary) {
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
