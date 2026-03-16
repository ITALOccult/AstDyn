/**
 * @file test_asteroid_perturbations.cpp
 * @brief Tests for asteroid perturbation module
 */

#include <gtest/gtest.h>
#include "astdyn/ephemeris/AsteroidPerturbations.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/time/epoch.hpp"
#include <iostream>

using namespace astdyn;
using namespace astdyn::ephemeris;
using namespace astdyn::propagation;

TEST(AsteroidPerturbationsTest, LoadDefault17Set) {
    AsteroidPerturbations ast;
    ast.loadAstDynDefaultSet();
    
    const auto& asteroids = ast.getAsteroids();
    
    // Check if we have 17 asteroids + Pluto (custom numbering 10-26)
    EXPECT_EQ(asteroids.size(), 17);
    
    // Check Pluto
    EXPECT_EQ(asteroids[0].name, "Pluto");
    EXPECT_EQ(asteroids[0].number, 10);
    EXPECT_EQ(asteroids[0].naif_id, 999);
    
    // Check Vesta
    EXPECT_EQ(asteroids[16].name, "Vesta");
    EXPECT_EQ(asteroids[16].number, 26);
    EXPECT_EQ(asteroids[16].naif_id, 2000004);
}

TEST(AsteroidPerturbationsTest, LoadDefault30Set) {
    AsteroidPerturbations ast;
    ast.loadDefault30Asteroids();
    
    const auto& asteroids = ast.getAsteroids();
    EXPECT_EQ(asteroids.size(), 30);
    
    // Check some members
    EXPECT_EQ(asteroids[0].name, "Ceres");
    EXPECT_EQ(asteroids[3].name, "Vesta");
}

TEST(AsteroidPerturbationsTest, PropagatorSettingsDefaults) {
    PropagatorSettings settings;
    
    // Check new high-precision defaults
    EXPECT_TRUE(settings.include_asteroids);
    EXPECT_TRUE(settings.use_default_asteroid_set);
    EXPECT_FALSE(settings.use_default_30_set);
    EXPECT_TRUE(settings.include_sun_j2);
    EXPECT_TRUE(settings.include_earth_j2);
    EXPECT_FALSE(settings.asteroid_ephemeris_file.empty());
}

TEST(AsteroidPerturbationsTest, AccelerationCalculation) {
    AsteroidPerturbations ast;
    ast.loadAstDynDefaultSet();
    
    // Test object at 2 AU
    Eigen::Vector3d pos_au(2.0, 0.0, 0.0);
    double mjd = 60000.0;
    Eigen::Vector3d sun_pos_bary_au(0.0, 0.0, 0.0);
    
    // Compute perturbation (analytical fallback since no SPK loaded here)
    Eigen::Vector3d acc = ast.computePerturbationRaw(pos_au, mjd, sun_pos_bary_au, true);
    
    // Acceleration should be non-zero
    EXPECT_GT(acc.norm(), 0.0);
    
    // Should be very small compared to Sun gravity (~mu/r^2 = 3e-4 / 4 = 7e-5)
    // Massive asteroids produce around 1e-12 to 1e-10 AU/d^2
    EXPECT_LT(acc.norm(), 1e-8);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
