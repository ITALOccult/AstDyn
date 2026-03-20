/**
 * @file test_constants.cpp
 * @brief Unit tests for physical constants
 */

#include <gtest/gtest.h>
#include <astdyn/core/Constants.hpp>
#include <cmath>

using namespace astdyn::constants;

class ConstantsTest : public ::testing::Test {
protected:
    // Tolerance for floating point comparisons
    static constexpr double TOLERANCE = 1.0e-10;
};

// ============================================================================
// Mathematical Constants Tests
// ============================================================================

TEST_F(ConstantsTest, MathematicalConstants) {
    // Test Pi
    EXPECT_NEAR(PI, 3.14159265358979323846, TOLERANCE);
    
    // Test Pi derivatives
    EXPECT_NEAR(TWO_PI, 2.0 * PI, TOLERANCE);
    EXPECT_NEAR(HALF_PI, PI / 2.0, TOLERANCE);
    
    // Test degree/radian conversions
    EXPECT_NEAR(DEG_TO_RAD * 180.0, PI, TOLERANCE);
    EXPECT_NEAR(RAD_TO_DEG * PI, 180.0, TOLERANCE);
    
    // Test arcsecond conversions
    EXPECT_NEAR(ARCSEC_TO_RAD * RAD_TO_ARCSEC, 1.0, TOLERANCE);
}

// ============================================================================
// Time Constants Tests
// ============================================================================

TEST_F(ConstantsTest, TimeConstants) {
    // Test J2000 epoch
    EXPECT_DOUBLE_EQ(JD2000, 2451545.0);
    EXPECT_DOUBLE_EQ(MJD2000, 51544.5);
    
    // Test relationship between JD and MJD
    EXPECT_DOUBLE_EQ(JD2000 - MJD2000, 2400000.5);
    
    // Test time conversion factors
    EXPECT_DOUBLE_EQ(SECONDS_PER_DAY, 86400.0);
    EXPECT_DOUBLE_EQ(DAYS_PER_YEAR, 365.25);
    EXPECT_DOUBLE_EQ(DAYS_PER_CENTURY, 36525.0);
}

// ============================================================================
// Physical Constants Tests
// ============================================================================

TEST_F(ConstantsTest, PhysicalConstants) {
    // Test speed of light
    EXPECT_DOUBLE_EQ(C_LIGHT, 299792.458);
    
    // Test astronomical unit
    EXPECT_DOUBLE_EQ(AU, 149597870.700);
    
    // Test Gaussian constant
    EXPECT_NEAR(K_GAUSS, 0.01720209895, TOLERANCE);
    
    // Test GMS relationship: GMS = k^2
    EXPECT_NEAR(GMS, K_GAUSS * K_GAUSS, TOLERANCE);
}

// ============================================================================
// Gravitational Parameters Tests
// ============================================================================

TEST_F(ConstantsTest, GravitationalParameters) {
    // Test that GM values are positive
    EXPECT_GT(GM_SUN, 0.0);
    EXPECT_GT(GM_EARTH, 0.0);
    EXPECT_GT(GM_JUPITER, 0.0);
    
    // Test that Jupiter has much larger GM than Earth
    EXPECT_GT(GM_JUPITER, GM_EARTH * 100.0);
    
    // Test Earth+Moon system consistency
    EXPECT_NEAR(GM_EARTH_MOON, GM_EARTH + GM_MOON, 100.0); // Within 100 km^3/s^2
}

// ============================================================================
// Unit Conversion Tests
// ============================================================================

TEST_F(ConstantsTest, UnitConversions) {
    // Test AU conversions
    EXPECT_DOUBLE_EQ(AU_TO_KM, AU);
    EXPECT_NEAR(KM_TO_AU * AU_TO_KM, 1.0, TOLERANCE);
    
    // Test velocity conversions
    EXPECT_NEAR(AU_PER_DAY_TO_KM_PER_S * KM_PER_S_TO_AU_PER_DAY, 1.0, TOLERANCE);
    
    // Test that AU/day to km/s makes sense
    double velocity_au_per_day = 1.0;  // 1 AU/day
    double velocity_km_per_s = velocity_au_per_day * AU_PER_DAY_TO_KM_PER_S;
    EXPECT_NEAR(velocity_km_per_s, AU / SECONDS_PER_DAY, TOLERANCE);
}

// ============================================================================
// Mass Ratio Tests
// ============================================================================

TEST_F(ConstantsTest, MassRatios) {
    // Test that mass ratios are positive and less than 1
    EXPECT_GT(MERCURY_MASS_RATIO, 0.0);
    EXPECT_LT(MERCURY_MASS_RATIO, 1.0);
    
    EXPECT_GT(JUPITER_MASS_RATIO, 0.0);
    EXPECT_LT(JUPITER_MASS_RATIO, 1.0);
    
    // Test that Jupiter has larger mass ratio than Earth
    EXPECT_GT(JUPITER_MASS_RATIO, EARTH_MOON_MASS_RATIO);
}

// ============================================================================
// Relativistic Constants Tests
// ============================================================================

TEST_F(ConstantsTest, RelativisticConstants) {
    // Test Schwarzschild radius of the Sun is positive and reasonable
    EXPECT_GT(SCHWARZSCHILD_SUN, 0.0);
    EXPECT_LT(SCHWARZSCHILD_SUN, 10.0); // Should be around 2.95 km
    
    // Test inverse c-squared in AU units
    EXPECT_GT(INV_C2_AU, 0.0);
    EXPECT_LT(INV_C2_AU, 1.0e-6); // Very small number
}

// ============================================================================
// Tolerance Tests
// ============================================================================

TEST_F(ConstantsTest, NumericalTolerances) {
    EXPECT_GT(EPSILON, 0.0);
    EXPECT_LT(EPSILON, 1.0e-10);
    
    EXPECT_GT(SMALL, 0.0);
    EXPECT_LT(SMALL, 1.0e-6);
    
    EXPECT_GT(LARGE, 1.0e20);
    
    EXPECT_GT(DEFAULT_TOLERANCE, 0.0);
    EXPECT_LT(DEFAULT_TOLERANCE, 1.0e-6);
}

// Run all tests
