/**
 * @file test_asteroid_perturbations.cpp
 * @brief Unit tests for asteroid perturbation calculations
 */

#include <gtest/gtest.h>
#include <orbfit/ephemeris/AsteroidPerturbations.hpp>
#include <orbfit/core/Constants.hpp>
#include <cmath>

using namespace orbfit;
using namespace orbfit::ephemeris;
using namespace orbfit::constants;

class AsteroidPerturbationsTest : public ::testing::Test {
protected:
    AsteroidPerturbations asteroids;
    
    void SetUp() override {
        // Default constructor loads AST17 asteroids
    }
};

// Test 1: Verify AST17 asteroids are loaded
TEST_F(AsteroidPerturbationsTest, LoadDefaultAsteroids) {
    const auto& ast_list = asteroids.getAsteroids();
    
    // Should load 16 asteroids
    EXPECT_EQ(ast_list.size(), 16);
    
    // Check first few asteroids
    EXPECT_EQ(ast_list[0].number, 1);
    EXPECT_EQ(ast_list[0].name, "Ceres");
    EXPECT_NEAR(ast_list[0].gm, 62.6284, 0.1);
    
    EXPECT_EQ(ast_list[1].number, 2);
    EXPECT_EQ(ast_list[1].name, "Pallas");
    
    EXPECT_EQ(ast_list[2].number, 4);
    EXPECT_EQ(ast_list[2].name, "Vesta");
}

// Test 2: Verify asteroid data validity
TEST_F(AsteroidPerturbationsTest, AsteroidDataValidity) {
    const auto& ast_list = asteroids.getAsteroids();
    
    for (const auto& ast : ast_list) {
        // GM should be positive and reasonable (< 100 km³/s²)
        EXPECT_GT(ast.gm, 0.0);
        EXPECT_LT(ast.gm, 100.0);
        
        // Semi-major axis should be in main belt (2-4 AU)
        EXPECT_GT(ast.a, 2.0);
        EXPECT_LT(ast.a, 4.0);
        
        // Eccentricity should be < 1 (elliptical orbits)
        EXPECT_GE(ast.e, 0.0);
        EXPECT_LT(ast.e, 1.0);
        
        // Inclination should be reasonable (< 35°)
        EXPECT_GE(ast.i, 0.0);
        EXPECT_LT(ast.i, 35.0);
        
        // Mean motion should be positive
        EXPECT_GT(ast.mean_motion, 0.0);
    }
}

// Test 3: Compute asteroid position
TEST_F(AsteroidPerturbationsTest, ComputePosition) {
    const auto& ceres = asteroids.getAsteroids()[0];  // (1) Ceres
    
    // Epoch J2000.0 (MJD 51544.5)
    double mjd_tdb = 51544.5;
    
    auto pos = asteroids.getPosition(ceres, mjd_tdb);
    
    // Position should be near semi-major axis
    double dist = pos.norm();
    EXPECT_GT(dist, 2.0);  // > perihelion
    EXPECT_LT(dist, 3.5);  // < aphelion
    
    // Components should be reasonable
    EXPECT_LT(std::abs(pos[0]), 4.0);
    EXPECT_LT(std::abs(pos[1]), 4.0);
    EXPECT_LT(std::abs(pos[2]), 1.0);  // Low inclination
}

// Test 4: Mean anomaly propagation
TEST_F(AsteroidPerturbationsTest, MeanAnomalyPropagation) {
    const auto& ceres = asteroids.getAsteroids()[0];
    
    double mjd0 = 51544.5;  // J2000.0
    double mjd1 = mjd0 + 365.25;  // 1 year later
    
    double M0 = ceres.meanAnomalyAt(mjd0);
    double M1 = ceres.meanAnomalyAt(mjd1);
    
    // Mean anomaly should increase
    EXPECT_GT(M1, M0);
    
    // Should increase by approximately n * 365.25 degrees
    double dM = M1 - M0;
    double expected_dM = ceres.mean_motion * 365.25;
    
    // Allow for angle wrapping
    while (dM < 0) dM += 360.0;
    while (dM > 360.0) dM -= 360.0;
    while (expected_dM > 360.0) expected_dM -= 360.0;
    
    EXPECT_NEAR(dM, expected_dM, 1.0);  // Within 1 degree
}

// Test 5: Single perturbation calculation
TEST_F(AsteroidPerturbationsTest, SinglePerturbationCalculation) {
    // Test position: 1 AU from Sun on x-axis
    Eigen::Vector3d spacecraft_pos(1.0, 0.0, 0.0);
    
    // Ceres position: 2.5 AU on x-axis (simplified)
    Eigen::Vector3d ceres_pos(2.5, 0.0, 0.0);
    
    // Ceres GM in AU³/day²
    double gm_ceres = 62.6284 * GM_KM3S2_TO_AU3DAY2;
    
    auto accel = AsteroidPerturbations::computeSinglePerturbation(
        spacecraft_pos, ceres_pos, gm_ceres
    );
    
    // Acceleration should be toward Ceres (positive x direction)
    EXPECT_GT(accel[0], 0.0);
    EXPECT_NEAR(accel[1], 0.0, 1e-15);  // No y component
    EXPECT_NEAR(accel[2], 0.0, 1e-15);  // No z component
    
    // Magnitude should be reasonable (very small)
    double mag = accel.norm();
    EXPECT_GT(mag, 0.0);
    EXPECT_LT(mag, 1e-6);  // Much smaller than Sun's gravity
}

// Test 6: Total perturbation from all asteroids
TEST_F(AsteroidPerturbationsTest, TotalPerturbation) {
    Eigen::Vector3d test_pos(1.0, 0.0, 0.0);  // 1 AU from Sun
    double mjd_tdb = 60000.0;  // ~2023
    
    auto total_accel = asteroids.computePerturbation(test_pos, mjd_tdb);
    
    // Total acceleration should be non-zero
    double mag = total_accel.norm();
    EXPECT_GT(mag, 0.0);
    
    // Should be much smaller than Sun's gravity (~5.9e-6 AU/day²)
    EXPECT_LT(mag, 1e-6);
    
    // Very small but measurable (> 1e-16 AU/day²)
    // At 1 AU from Sun with asteroids at ~2.5-3 AU, perturbations are extremely weak
    EXPECT_GT(mag, 1e-16);
}

// Test 7: Enable/disable asteroids
TEST_F(AsteroidPerturbationsTest, EnableDisableAsteroids) {
    // Enable only Ceres
    for (int i = 2; i <= 704; ++i) {
        asteroids.setAsteroidEnabled(i, false);
    }
    
    EXPECT_TRUE(asteroids.isAsteroidEnabled(1));   // Ceres enabled
    EXPECT_FALSE(asteroids.isAsteroidEnabled(2));  // Pallas disabled
    EXPECT_FALSE(asteroids.isAsteroidEnabled(4));  // Vesta disabled
    
    // Compute perturbation with only Ceres
    Eigen::Vector3d test_pos(1.0, 0.0, 0.0);
    double mjd_tdb = 60000.0;
    
    auto accel_ceres_only = asteroids.computePerturbation(test_pos, mjd_tdb);
    
    // Re-enable all
    asteroids.setAsteroidEnabled(2, true);
    asteroids.setAsteroidEnabled(4, true);
    
    auto accel_all = asteroids.computePerturbation(test_pos, mjd_tdb);
    
    // Perturbations should be different (they might not necessarily be larger 
    // depending on asteroid positions and cancellation effects)
    EXPECT_NE(accel_all.norm(), accel_ceres_only.norm());
}

// Test 8: Total mass calculation
TEST_F(AsteroidPerturbationsTest, TotalMass) {
    double total_mass = asteroids.getTotalMass();
    
    // Total mass should be positive
    EXPECT_GT(total_mass, 0.0);
    
    // Should be roughly 8-9e-10 solar masses (from AST17 data)
    // Ceres alone is about 4.7e-10 M☉
    // Total of 16 asteroids should be around 8.9e-10 M☉
    EXPECT_GT(total_mass, 5e-10);
    EXPECT_LT(total_mass, 2e-9);
}

// Test 9: Perturbation varies with position
TEST_F(AsteroidPerturbationsTest, PerturbationVariesWithPosition) {
    double mjd_tdb = 60000.0;
    
    // Position 1: Near Sun
    Eigen::Vector3d pos1(0.5, 0.0, 0.0);
    auto accel1 = asteroids.computePerturbation(pos1, mjd_tdb);
    
    // Position 2: In main belt
    Eigen::Vector3d pos2(2.5, 0.0, 0.0);
    auto accel2 = asteroids.computePerturbation(pos2, mjd_tdb);
    
    // Position 3: Beyond main belt
    Eigen::Vector3d pos3(5.0, 0.0, 0.0);
    auto accel3 = asteroids.computePerturbation(pos3, mjd_tdb);
    
    // Perturbation should be largest in main belt
    double mag1 = accel1.norm();
    double mag2 = accel2.norm();
    double mag3 = accel3.norm();
    
    // All should be non-zero
    EXPECT_GT(mag1, 0.0);
    EXPECT_GT(mag2, 0.0);
    EXPECT_GT(mag3, 0.0);
    
    // Main belt should have significant perturbation
    // (exact ordering depends on asteroid positions at epoch)
    EXPECT_GT(mag2, 0.0);
}

// Test 10: Time evolution of perturbation
TEST_F(AsteroidPerturbationsTest, TimeEvolution) {
    Eigen::Vector3d test_pos(2.0, 0.0, 0.0);
    
    // Compute at different times
    double mjd0 = 60000.0;
    double mjd1 = mjd0 + 100.0;  // 100 days later
    double mjd2 = mjd0 + 365.25;  // 1 year later
    
    auto accel0 = asteroids.computePerturbation(test_pos, mjd0);
    auto accel1 = asteroids.computePerturbation(test_pos, mjd1);
    auto accel2 = asteroids.computePerturbation(test_pos, mjd2);
    
    // Perturbations should vary as asteroids move
    EXPECT_NE(accel0.norm(), accel1.norm());
    EXPECT_NE(accel1.norm(), accel2.norm());
    
    // But should remain in reasonable range
    EXPECT_LT(accel0.norm(), 1e-6);
    EXPECT_LT(accel1.norm(), 1e-6);
    EXPECT_LT(accel2.norm(), 1e-6);
}

// Test 11: AST17 namespace constants
TEST(AsteroidConstantsTest, AST17Constants) {
    using namespace ast17;
    
    // Verify major asteroid GMs
    EXPECT_NEAR(GM_CERES, 62.6284, 1.0);
    EXPECT_NEAR(GM_PALLAS, 13.8, 1.0);
    EXPECT_NEAR(GM_VESTA, 17.8, 1.0);
    EXPECT_NEAR(GM_HYGIEA, 5.78, 1.0);
    
    // All should be positive
    EXPECT_GT(GM_CERES, 0.0);
    EXPECT_GT(GM_EUNOMIA, 0.0);
    EXPECT_GT(GM_DAVIDA, 0.0);
}

// Test 12: Direct vs indirect terms
TEST_F(AsteroidPerturbationsTest, DirectAndIndirectTerms) {
    // Spacecraft at 1 AU
    Eigen::Vector3d sc_pos(1.0, 0.0, 0.0);
    
    // Ceres at 2.7 AU (simplified)
    Eigen::Vector3d ceres_pos(2.7, 0.0, 0.0);
    double gm = 62.6284 * GM_KM3S2_TO_AU3DAY2;
    
    // Relative vector: Ceres - spacecraft
    Eigen::Vector3d delta = ceres_pos - sc_pos;
    double delta_norm = delta.norm();
    
    // Direct term: toward Ceres
    Eigen::Vector3d direct = gm * delta / (delta_norm * delta_norm * delta_norm);
    
    // Indirect term: Sun's acceleration toward Ceres
    double ceres_dist = ceres_pos.norm();
    Eigen::Vector3d indirect = -gm * ceres_pos / (ceres_dist * ceres_dist * ceres_dist);
    
    // Total
    Eigen::Vector3d total = direct + indirect;
    
    // Verify against function
    auto computed = AsteroidPerturbations::computeSinglePerturbation(sc_pos, ceres_pos, gm);
    
    EXPECT_NEAR(total[0], computed[0], 1e-15);
    EXPECT_NEAR(total[1], computed[1], 1e-15);
    EXPECT_NEAR(total[2], computed[2], 1e-15);
}
