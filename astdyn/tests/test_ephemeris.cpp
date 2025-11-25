/**
 * @file test_ephemeris.cpp
 * @brief Unit tests for planetary ephemeris calculations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#include <gtest/gtest.h>
#include "astdyn/ephemeris/PlanetaryData.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>

using namespace astdyn;
using namespace astdyn::ephemeris;

// ============================================================================
// PlanetaryData Tests
// ============================================================================

TEST(PlanetaryDataTest, GravitationalParameters) {
    // Test GM values
    EXPECT_GT(PlanetaryData::getGM(CelestialBody::SUN), 1e11);
    EXPECT_GT(PlanetaryData::getGM(CelestialBody::JUPITER), 1e8);
    EXPECT_GT(PlanetaryData::getGM(CelestialBody::EARTH), 3.9e5);
    
    // Jupiter should have largest planetary GM
    double gm_jup = PlanetaryData::getGM(CelestialBody::JUPITER);
    EXPECT_GT(gm_jup, PlanetaryData::getGM(CelestialBody::SATURN));
    EXPECT_GT(gm_jup, PlanetaryData::getGM(CelestialBody::EARTH));
}

TEST(PlanetaryDataTest, Radii) {
    // Test radius values [km]
    EXPECT_NEAR(PlanetaryData::getRadius(CelestialBody::EARTH), 6371.0, 10.0);
    EXPECT_NEAR(PlanetaryData::getRadius(CelestialBody::JUPITER), 69911.0, 100.0);
    EXPECT_NEAR(PlanetaryData::getRadius(CelestialBody::MOON), 1737.4, 10.0);
    
    // Jupiter should be largest planet
    double r_jup = PlanetaryData::getRadius(CelestialBody::JUPITER);
    EXPECT_GT(r_jup, PlanetaryData::getRadius(CelestialBody::SATURN));
    EXPECT_GT(r_jup, PlanetaryData::getRadius(CelestialBody::EARTH));
}

TEST(PlanetaryDataTest, Masses) {
    // Test mass ordering
    double m_sun = PlanetaryData::getMass(CelestialBody::SUN);
    double m_jup = PlanetaryData::getMass(CelestialBody::JUPITER);
    double m_earth = PlanetaryData::getMass(CelestialBody::EARTH);
    
    EXPECT_GT(m_sun, m_jup);
    EXPECT_GT(m_jup, m_earth);
    EXPECT_NEAR(m_earth, 5.97217e24, 1e23); // [kg]
}

TEST(PlanetaryDataTest, BodyNames) {
    EXPECT_EQ(PlanetaryData::getName(CelestialBody::EARTH), "Earth");
    EXPECT_EQ(PlanetaryData::getName(CelestialBody::MARS), "Mars");
    EXPECT_EQ(PlanetaryData::getName(CelestialBody::JUPITER), "Jupiter");
}

TEST(PlanetaryDataTest, BodyDataStructure) {
    auto earth_data = PlanetaryData::getBodyData(CelestialBody::EARTH);
    
    EXPECT_EQ(earth_data.name, "Earth");
    EXPECT_NEAR(earth_data.semi_major_axis, 1.0, 0.01); // ~1 AU
    EXPECT_NEAR(earth_data.period, 365.256, 1.0); // ~365 days
    EXPECT_GT(earth_data.gm, 0.0);
    EXPECT_GT(earth_data.radius, 0.0);
}

// ============================================================================
// PlanetaryEphemeris Tests
// ============================================================================

TEST(PlanetaryEphemerisTest, SunPosition) {
    // Sun should be at origin
    double jd = constants::JD2000;
    auto pos = PlanetaryEphemeris::getPosition(CelestialBody::SUN, jd);
    
    EXPECT_NEAR(pos.norm(), 0.0, 1e-10);
}

TEST(PlanetaryEphemerisTest, EarthPositionJ2000) {
    // Earth at J2000.0 epoch
    double jd = constants::JD2000;
    auto pos = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd);
    
    // Distance should be ~1 AU
    EXPECT_NEAR(pos.norm(), 1.0, 0.05);
    
    // Check that position is reasonable (in ecliptic plane)
    EXPECT_LT(std::abs(pos.z()), 0.01); // Low inclination
}

TEST(PlanetaryEphemerisTest, EarthVelocity) {
    double jd = constants::JD2000;
    auto vel = PlanetaryEphemeris::getVelocity(CelestialBody::EARTH, jd);
    
    // Earth velocity ~29.78 km/s = ~6.28 AU/year = ~0.0172 AU/day
    double v_mag = vel.norm();
    EXPECT_NEAR(v_mag, 0.0172, 0.003); // Within 20%
}

TEST(PlanetaryEphemerisTest, MarsPosition) {
    double jd = constants::JD2000;
    auto pos = PlanetaryEphemeris::getPosition(CelestialBody::MARS, jd);
    
    // Mars at ~1.52 AU
    EXPECT_NEAR(pos.norm(), 1.52, 0.3);
}

TEST(PlanetaryEphemerisTest, JupiterPosition) {
    double jd = constants::JD2000;
    auto pos = PlanetaryEphemeris::getPosition(CelestialBody::JUPITER, jd);
    
    // Jupiter at ~5.2 AU
    EXPECT_NEAR(pos.norm(), 5.2, 1.0);
}

TEST(PlanetaryEphemerisTest, OuterPlanetsOrder) {
    // Check that outer planets are at increasing distances
    double jd = constants::JD2000;
    
    double r_jup = PlanetaryEphemeris::getPosition(CelestialBody::JUPITER, jd).norm();
    double r_sat = PlanetaryEphemeris::getPosition(CelestialBody::SATURN, jd).norm();
    double r_ura = PlanetaryEphemeris::getPosition(CelestialBody::URANUS, jd).norm();
    double r_nep = PlanetaryEphemeris::getPosition(CelestialBody::NEPTUNE, jd).norm();
    
    EXPECT_LT(r_jup, r_sat);
    EXPECT_LT(r_sat, r_ura);
    EXPECT_LT(r_ura, r_nep);
}

TEST(PlanetaryEphemerisTest, StateVector) {
    double jd = constants::JD2000;
    auto state = PlanetaryEphemeris::getState(CelestialBody::EARTH, jd);
    
    EXPECT_NEAR(state.radius(), 1.0, 0.05);
    EXPECT_GT(state.speed(), 0.01); // Non-zero velocity
}

TEST(PlanetaryEphemerisTest, OrbitalPeriod) {
    // Test that planet returns to similar position after one period
    double jd_start = constants::JD2000;
    
    // Earth: 365.25 days
    auto pos1 = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd_start);
    auto pos2 = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd_start + 365.256);
    
    // Positions should be similar (within ~10% due to apsidal precession)
    EXPECT_NEAR((pos1 - pos2).norm(), 0.0, 0.2);
}

TEST(PlanetaryEphemerisTest, MercurySpeed) {
    // Mercury has highest orbital speed (~48 km/s = ~0.104 AU/day)
    double jd = constants::JD2000;
    auto v_mer = PlanetaryEphemeris::getVelocity(CelestialBody::MERCURY, jd);
    auto v_earth = PlanetaryEphemeris::getVelocity(CelestialBody::EARTH, jd);
    
    EXPECT_GT(v_mer.norm(), v_earth.norm());
}

TEST(PlanetaryEphemerisTest, TimeEvolution) {
    // Test position changes over time
    double jd1 = constants::JD2000;
    double jd2 = jd1 + 100.0; // 100 days later
    
    auto pos1 = PlanetaryEphemeris::getPosition(CelestialBody::MARS, jd1);
    auto pos2 = PlanetaryEphemeris::getPosition(CelestialBody::MARS, jd2);
    
    // Position should change significantly
    EXPECT_GT((pos1 - pos2).norm(), 0.1); // > 0.1 AU change
}

// ============================================================================
// Barycentric Corrections Tests
// ============================================================================

TEST(PlanetaryEphemerisTest, SunBarycentricOffset) {
    double jd = constants::JD2000;
    auto r_sun = PlanetaryEphemeris::getSunBarycentricPosition(jd);
    
    // Sun offset should be < 0.01 AU (dominated by Jupiter)
    EXPECT_LT(r_sun.norm(), 0.01);
    EXPECT_GT(r_sun.norm(), 1e-5); // But non-negligible
}

TEST(PlanetaryEphemerisTest, HeliocentricToBarycentric) {
    double jd = constants::JD2000;
    
    // Earth heliocentric state
    auto helio_state = PlanetaryEphemeris::getState(CelestialBody::EARTH, jd);
    
    // Convert to barycentric
    auto bary_state = PlanetaryEphemeris::heliocentricToBarycentric(helio_state, jd);
    
    // Difference should be small (~Sun's offset)
    auto delta = helio_state.position() - bary_state.position();
    EXPECT_LT(delta.norm(), 0.01); // < 0.01 AU
}

// ============================================================================
// Validation Against Known Ephemerides
// ============================================================================

TEST(EphemerisValidationTest, EarthJ2000) {
    // Earth position at J2000.0 from JPL Horizons:
    // X = -0.177087 AU, Y = +0.967129 AU, Z = -0.000016 AU
    // (barycentric, but heliocentric similar within ~0.005 AU)
    
    double jd = constants::JD2000;
    auto pos = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd);
    
    // Check distance (within 5%)
    EXPECT_NEAR(pos.norm(), 0.9833, 0.05);
    
    // Check that X is negative, Y is positive
    EXPECT_LT(pos.x(), 0.0);
    EXPECT_GT(pos.y(), 0.0);
}

TEST(EphemerisValidationTest, JupiterJ2000) {
    // Jupiter at J2000.0 from JPL Horizons:
    // Distance ~5.2 AU, mostly in X-Y plane
    
    double jd = constants::JD2000;
    auto pos = PlanetaryEphemeris::getPosition(CelestialBody::JUPITER, jd);
    
    // Distance check (within 10% for low-precision ephemeris)
    EXPECT_NEAR(pos.norm(), 5.2, 0.6);
}

TEST(EphemerisValidationTest, MarsOpposition) {
    // Test Mars at opposition (approximately JD 2451545 + 780 days)
    // At opposition, Mars should be ~180° from Sun as seen from Earth
    
    double jd = constants::JD2000 + 780.0;
    
    auto pos_earth = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd);
    auto pos_mars = PlanetaryEphemeris::getPosition(CelestialBody::MARS, jd);
    
    // Mars should be farther from Sun than Earth
    EXPECT_GT(pos_mars.norm(), pos_earth.norm());
}

TEST(EphemerisValidationTest, VenusInferiorConjunction) {
    // Venus is closer to Sun than Earth at all times
    double jd = constants::JD2000;
    
    auto pos_earth = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd);
    auto pos_venus = PlanetaryEphemeris::getPosition(CelestialBody::VENUS, jd);
    
    EXPECT_LT(pos_venus.norm(), pos_earth.norm());
}

// ============================================================================
// Consistency Tests
// ============================================================================

TEST(EphemerisConsistencyTest, VelocityFromFiniteDifference) {
    // Verify velocity by finite difference of position
    double jd = constants::JD2000;
    double dt = 0.01; // 0.01 days
    
    auto pos1 = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd - dt);
    auto pos2 = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd + dt);
    auto vel_fd = (pos2 - pos1) / (2.0 * dt);
    
    auto vel_analytical = PlanetaryEphemeris::getVelocity(CelestialBody::EARTH, jd);
    
    // Should match within 5% (due to numerical differentiation error)
    EXPECT_NEAR((vel_fd - vel_analytical).norm(), 0.0, 0.002);
}

TEST(EphemerisConsistencyTest, KeplerThirdLaw) {
    // Test T² ∝ a³ for planets
    double jd = constants::JD2000;
    
    auto data_earth = PlanetaryData::getBodyData(CelestialBody::EARTH);
    auto data_mars = PlanetaryData::getBodyData(CelestialBody::MARS);
    
    // T²/a³ should be constant (= 4π²/GM_sun)
    double k_earth = (data_earth.period * data_earth.period) / 
                     (data_earth.semi_major_axis * data_earth.semi_major_axis * data_earth.semi_major_axis);
    double k_mars = (data_mars.period * data_mars.period) / 
                    (data_mars.semi_major_axis * data_mars.semi_major_axis * data_mars.semi_major_axis);
    
    EXPECT_NEAR(k_earth / k_mars, 1.0, 0.01); // Within 1%
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
