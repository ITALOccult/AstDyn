/**
 * @file test_validation.cpp
 * @brief Validation tests against known reference data
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * This test suite validates AstDyn C++ against:
 * 1. Published orbital elements (MPC, JPL)
 * 2. Well-known orbits (ISS, GPS, GEO satellites)
 * 3. IAU standard test cases
 */

#include <gtest/gtest.h>
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/coordinates/KeplerianElements.hpp"
#include "astdyn/coordinates/CometaryElements.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"

using namespace astdyn;
using namespace astdyn::coordinates;
using namespace astdyn::constants;

// ========== Validation Against Known Orbits ==========

TEST(ValidationTest, ISSOrbit) {
    /**
     * ISS orbital parameters (approximate, epoch 2020)
     * Source: https://www.heavens-above.com/
     */
    double a = 6778.14;              // Semi-major axis [km]
    double e = 0.0005730;            // Eccentricity
    double i = 51.6433 * DEG_TO_RAD; // Inclination [rad]
    
    KeplerianElements iss(a, e, i, 0, 0, 0, GM_EARTH);
    
    // Validate orbital period (~92.9 minutes)
    double period_minutes = iss.period() / 60.0;
    EXPECT_NEAR(period_minutes, 92.9, 0.5);
    
    // Validate altitude (~408 km)
    double altitude = a - 6371.0;  // Mean Earth radius
    EXPECT_NEAR(altitude, 408.0, 10.0);
    
    // Validate orbital velocity (~7.66 km/s)
    CartesianState state = iss.to_cartesian();
    EXPECT_NEAR(state.speed(), 7.66, 0.1);
}

TEST(ValidationTest, GPSOrbit) {
    /**
     * GPS satellite orbital parameters
     * Source: GPS Interface Specification IS-GPS-200
     */
    double a = 26559.7;              // Semi-major axis [km]
    double e = 0.01;                 // Max eccentricity (typically < 0.02)
    double i = 55.0 * DEG_TO_RAD;    // Inclination [rad]
    
    KeplerianElements gps(a, e, i, 0, 0, 0, GM_EARTH);
    
    // Validate orbital period (~12 hours = 43200 seconds)
    double period_hours = gps.period() / 3600.0;
    EXPECT_NEAR(period_hours, 11.967, 0.1);  // ~11h 58m
    
    // Validate altitude (~20,200 km)
    double altitude = a - 6371.0;
    EXPECT_NEAR(altitude, 20188.7, 10.0);
}

TEST(ValidationTest, GEOOrbit) {
    /**
     * Geostationary orbit
     * Period must equal Earth's sidereal day (23h 56m 4s = 86164 s)
     */
    double period_geo = 86164.0905;  // Sidereal day [s]
    
    // Compute semi-major axis from Kepler's third law: a³ = μ·T²/(4π²)
    double a_geo = std::cbrt(GM_EARTH * period_geo * period_geo / 
                             (4.0 * PI * PI));
    
    KeplerianElements geo(a_geo, 0.0, 0.0, 0, 0, 0, GM_EARTH);
    
    // Validate computed period matches sidereal day
    EXPECT_NEAR(geo.period(), period_geo, 1.0);
    
    // Validate altitude (~35,786 km)
    double altitude = a_geo - 6371.0;
    EXPECT_NEAR(altitude, 35786.0, 10.0);
    
    // Validate orbital velocity (~3.07 km/s)
    CartesianState state = geo.to_cartesian();
    EXPECT_NEAR(state.speed(), 3.07, 0.1);
}

TEST(ValidationTest, MolniyaOrbit) {
    /**
     * Molniya orbit (Russian communications satellites)
     * Source: Wertz, "Mission Geometry; Orbit and Constellation Design"
     */
    double a = 26554.0;              // Semi-major axis [km]
    double e = 0.722;                // Eccentricity (high!)
    double i = 63.4 * DEG_TO_RAD;    // Critical inclination (frozen orbit)
    
    KeplerianElements molniya(a, e, i, 0, 0, 0, GM_EARTH);
    
    // Validate period (~12 hours)
    double period_hours = molniya.period() / 3600.0;
    EXPECT_NEAR(period_hours, 11.97, 0.1);
    
    // Validate periapsis altitude (actual ~1000 km with this a and e)
    double r_peri = molniya.periapsis_distance();
    double alt_peri = r_peri - 6371.0;
    EXPECT_NEAR(alt_peri, 1000.0, 200.0);
    
    // Validate apoapsis altitude (actual ~39,350 km with this a and e)
    double r_apo = molniya.apoapsis_distance();
    double alt_apo = r_apo - 6371.0;
    EXPECT_NEAR(alt_apo, 39350.0, 500.0);
}

// ========== Comet Validation ==========

TEST(ValidationTest, HalleyComet) {
    /**
     * 1P/Halley orbital elements (epoch 1986 perihelion)
     * Source: JPL Small-Body Database
     */
    double q = 0.5871 * AU;          // Perihelion distance
    double e = 0.967143;             // Eccentricity
    double i = 162.2626 * DEG_TO_RAD;// Inclination (retrograde)
    
    CometaryElements halley(q, e, i, 0, 0, 0, GM_SUN);
    
    // Validate semi-major axis (~17.8 AU)
    double a_au = halley.semi_major_axis() / AU;
    EXPECT_NEAR(a_au, 17.83, 0.5);
    
    // Validate period (~75-76 years)
    double period_years = halley.period() / YEAR;
    EXPECT_NEAR(period_years, 75.3, 1.0);
    
    // Validate aphelion distance (~35.1 AU)
    double Q_au = halley.aphelion_distance() / AU;
    EXPECT_NEAR(Q_au, 35.1, 0.5);
}

TEST(ValidationTest, EnckeComet) {
    /**
     * 2P/Encke - shortest period comet
     * Source: JPL Small-Body Database
     */
    double q = 0.3363 * AU;          // Perihelion
    double e = 0.8489;               // Eccentricity
    double i = 11.7808 * DEG_TO_RAD; // Inclination
    
    CometaryElements encke(q, e, i, 0, 0, 0, GM_SUN);
    
    // Validate period (~3.3 years)
    double period_years = encke.period() / YEAR;
    EXPECT_NEAR(period_years, 3.3, 0.1);
    
    // Validate semi-major axis (~2.2 AU)
    double a_au = encke.semi_major_axis() / AU;
    EXPECT_NEAR(a_au, 2.22, 0.1);
}

// ========== Reference Frame Validation ==========

TEST(ValidationTest, EclipticObliquity) {
    /**
     * Mean obliquity of the ecliptic at J2000.0
     * IAU 2000: ε₀ = 23°26'21.406" = 23.439291°
     */
    Matrix3d R = ReferenceFrame::j2000_to_ecliptic();
    
    // Extract rotation angle from matrix
    // For rotation about X-axis: R(1,1) = cos(ε), R(2,1) = -sin(ε)
    double obliquity = std::atan2(-R(2,1), R(1,1));
    double obliquity_deg = obliquity * RAD_TO_DEG;
    
    EXPECT_NEAR(obliquity_deg, 23.439291, 1e-6);
}

TEST(ValidationTest, J2000ICRSFrameBias) {
    /**
     * Frame bias between J2000 (FK5) and ICRS
     * IERS Conventions 2010: ~0.02 arcsec offset
     */
    Matrix3d bias = ReferenceFrame::j2000_to_icrs();
    
    // Bias should be very close to identity
    Matrix3d identity = Matrix3d::Identity();
    Matrix3d diff = bias - identity;
    
    // Maximum element should be ~1e-7 (corresponds to ~0.02 arcsec)
    double max_diff = diff.cwiseAbs().maxCoeff();
    EXPECT_LT(max_diff, 1e-6);  // < 0.2 arcsec
}

TEST(ValidationTest, EarthRotationRate) {
    /**
     * Earth rotation rate: 7.292115e-5 rad/s
     * GMST advances by one full rotation (2π rad) per sidereal day
     * Sidereal day = 23h 56m 4s = 86164 s
     * Solar day = 24h = 86400 s
     * GMST gain per solar day ≈ 2π × (86400/86164) ≈ 2π × 1.00274
     */
    double mjd1 = MJD2000;
    double mjd2 = MJD2000 + 1.0;  // 1 solar day later
    
    double gmst1 = ReferenceFrame::gmst(mjd1);
    double gmst2 = ReferenceFrame::gmst(mjd2);
    
    // GMST advance per solar day
    double dgmst = gmst2 - gmst1;
    // Normalize to [0, 2π)
    while (dgmst < 0) dgmst += 2.0 * PI;
    while (dgmst > 2.0 * PI) dgmst -= 2.0 * PI;
    
    // GMST should advance slightly MORE than 2π per solar day
    // (sidereal day is shorter than solar day)
    // Expected: 2π × 1.00274 = 6.3 rad, but after wrapping should be ~0.017 rad
    double excess = 2.0 * PI * 1.00274 - 2.0 * PI;
    EXPECT_NEAR(dgmst, excess, 0.01);  // Should be ~0.017 rad
}

// ========== Conservation Laws ==========

TEST(ValidationTest, EnergyConservation) {
    /**
     * Specific orbital energy should be conserved under coordinate transformations
     */
    KeplerianElements kep(7000.0, 0.1, 30.0*DEG_TO_RAD, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    CartesianState cart = kep.to_cartesian();
    
    // Energy from Keplerian: E = -μ/(2a)
    double E_kep = -GM_EARTH / (2.0 * kep.semi_major_axis());
    
    // Energy from Cartesian: E = v²/2 - μ/r
    double E_cart = cart.specific_energy();
    
    // Should match within numerical precision
    EXPECT_NEAR(E_kep, E_cart, 1e-6);
}

TEST(ValidationTest, AngularMomentumConservation) {
    /**
     * Angular momentum magnitude should be conserved
     */
    KeplerianElements kep(7000.0, 0.1, 30.0*DEG_TO_RAD, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    CartesianState cart = kep.to_cartesian();
    
    // Angular momentum from Keplerian: h = √(μ·a·(1-e²))
    double h_kep = std::sqrt(GM_EARTH * kep.semi_major_axis() * 
                            (1.0 - kep.eccentricity() * kep.eccentricity()));
    
    // Angular momentum from Cartesian: h = |r × v|
    double h_cart = cart.angular_momentum().norm();
    
    // Should match within numerical precision
    EXPECT_NEAR(h_kep, h_cart, 1e-3);
}

TEST(ValidationTest, MagnitudePreservation) {
    /**
     * Position and velocity magnitudes should be preserved under frame transformations
     */
    Vector3d pos(7000.0, 3000.0, 1000.0);
    Vector3d vel(2.0, 5.0, -3.0);
    CartesianState state_j2000(pos, vel, GM_EARTH);
    
    double r_mag = pos.norm();
    double v_mag = vel.norm();
    
    // Transform to multiple frames
    CartesianState state_icrs = ReferenceFrame::transform_state(
        state_j2000, FrameType::J2000, FrameType::ICRS);
    CartesianState state_ecl = ReferenceFrame::transform_state(
        state_j2000, FrameType::J2000, FrameType::ECLIPTIC);
    
    // Magnitudes should be preserved (rotations only)
    EXPECT_NEAR(state_icrs.radius(), r_mag, 1e-6);
    EXPECT_NEAR(state_icrs.speed(), v_mag, 1e-6);
    EXPECT_NEAR(state_ecl.radius(), r_mag, 1e-6);
    EXPECT_NEAR(state_ecl.speed(), v_mag, 1e-6);
}

// ========== Edge Cases ==========

TEST(ValidationTest, EscapeVelocity) {
    /**
     * Escape velocity at Earth's surface: v_esc = √(2μ/r) ≈ 11.186 km/s
     */
    double r = 6371.0;  // Earth radius
    double v_esc_theoretical = std::sqrt(2.0 * GM_EARTH / r);
    
    // Create state at escape velocity (parabolic orbit, e = 1)
    Vector3d pos(r, 0, 0);
    Vector3d vel(0, v_esc_theoretical, 0);
    CartesianState state(pos, vel, GM_EARTH);
    
    // Eccentricity should be 1.0 (parabolic)
    EXPECT_NEAR(state.eccentricity(), 1.0, 1e-6);
    
    // Specific energy should be zero
    EXPECT_NEAR(state.specific_energy(), 0.0, 1e-3);
}

TEST(ValidationTest, VisVivaEquation) {
    /**
     * Vis-viva equation: v² = μ(2/r - 1/a)
     * Validate at periapsis and apoapsis
     */
    double a = 10000.0;  // km
    double e = 0.3;
    
    KeplerianElements kep(a, e, 0, 0, 0, 0, GM_EARTH);
    
    // At periapsis (ν = 0)
    kep.set_mean_anomaly(0.0);
    CartesianState state_peri = kep.to_cartesian();
    double r_peri = state_peri.radius();
    double v_peri = state_peri.speed();
    
    // Vis-viva at periapsis
    double v_expected = std::sqrt(GM_EARTH * (2.0/r_peri - 1.0/a));
    EXPECT_NEAR(v_peri, v_expected, 1e-3);
    
    // At apoapsis (ν = π)
    kep.set_mean_anomaly(PI);
    CartesianState state_apo = kep.to_cartesian();
    double r_apo = state_apo.radius();
    double v_apo = state_apo.speed();
    
    // Vis-viva at apoapsis
    v_expected = std::sqrt(GM_EARTH * (2.0/r_apo - 1.0/a));
    EXPECT_NEAR(v_apo, v_expected, 1e-3);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
