/**
 * @file test_coordinates.cpp
 * @brief Unit tests for coordinate systems
 */

#include <gtest/gtest.h>
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/coordinates/KeplerianElements.hpp"
#include "astdyn/coordinates/EquinoctialElements.hpp"
#include "astdyn/coordinates/CometaryElements.hpp"

using namespace astdyn;
using namespace astdyn::coordinates;
using namespace astdyn::constants;

// ========== CartesianState Tests ==========

TEST(CartesianStateTest, Construction) {
    Vector3d pos(7000.0, 0.0, 0.0);
    Vector3d vel(0.0, 7.5, 0.0);
    
    CartesianState state(pos, vel, GM_EARTH);
    
    EXPECT_TRUE(state.position().isApprox(pos));
    EXPECT_TRUE(state.velocity().isApprox(vel));
    EXPECT_NEAR(state.mu(), GM_EARTH, 1e-10);
}

TEST(CartesianStateTest, StateVector) {
    Vector6d state_vec;
    state_vec << 7000.0, 0.0, 0.0, 0.0, 7.5, 0.0;
    
    CartesianState state(state_vec, GM_EARTH);
    Vector6d result = state.state_vector();
    
    EXPECT_TRUE(result.isApprox(state_vec));
}

TEST(CartesianStateTest, OrbitalProperties) {
    // Circular orbit at 7000 km
    Vector3d pos(7000.0, 0.0, 0.0);
    Vector3d vel(0.0, 7.546, 0.0); // Circular velocity
    
    CartesianState state(pos, vel, GM_EARTH);
    
    EXPECT_NEAR(state.radius(), 7000.0, 1e-10);
    EXPECT_NEAR(state.speed(), 7.546, 1e-3);
    
    // Angular momentum should be perpendicular to orbit plane
    Vector3d h = state.angular_momentum();
    EXPECT_NEAR(h.x(), 0.0, 1e-6);
    EXPECT_NEAR(h.y(), 0.0, 1e-6);
    EXPECT_GT(h.z(), 0.0);
}

TEST(CartesianStateTest, CircularOrbit) {
    Vector3d pos(7000.0, 0.0, 0.0);
    double v_circ = std::sqrt(GM_EARTH / 7000.0);
    Vector3d vel(0.0, v_circ, 0.0);
    
    CartesianState state(pos, vel, GM_EARTH);
    
    EXPECT_TRUE(state.is_circular(1e-3));
    EXPECT_TRUE(state.is_elliptic());
    EXPECT_FALSE(state.is_hyperbolic());
    
    double e = state.eccentricity();
    EXPECT_LT(e, 1e-3);
}

TEST(CartesianStateTest, EllipticOrbit) {
    // Molniya-like orbit
    Vector3d pos(26554.0, 0.0, 0.0); // Apogee
    Vector3d vel(0.0, 1.5, 0.0);
    
    CartesianState state(pos, vel, GM_EARTH);
    
    EXPECT_TRUE(state.is_elliptic());
    EXPECT_FALSE(state.is_circular());
    EXPECT_GT(state.eccentricity(), 0.1);
    EXPECT_LT(state.eccentricity(), 1.0);
}

TEST(CartesianStateTest, Inclination) {
    // Inclined orbit (45 degrees)
    double r = 7000.0;
    double v = std::sqrt(GM_EARTH / r);
    double i = 45.0 * DEG_TO_RAD;
    
    Vector3d pos(r, 0.0, 0.0);
    Vector3d vel(0.0, v * std::cos(i), v * std::sin(i));
    
    CartesianState state(pos, vel, GM_EARTH);
    
    EXPECT_NEAR(state.inclination(), i, 1e-6);
}

TEST(CartesianStateTest, UnitConversions) {
    Vector3d pos_km(AU, 0.0, 0.0);
    Vector3d vel_km(0.0, 30.0, 0.0);
    
    CartesianState state_km(pos_km, vel_km, GM_SUN);
    CartesianState state_au = state_km.to_AU_per_day();
    
    EXPECT_NEAR(state_au.position().x(), 1.0, 1e-10); // 1 AU
    
    // Convert back
    CartesianState state_km2 = state_au.to_km_per_s();
    EXPECT_NEAR(state_km2.position().x(), AU, 1.0); // Back to km
}

TEST(CartesianStateTest, RelativeState) {
    Vector3d pos1(7000.0, 0.0, 0.0);
    Vector3d vel1(0.0, 7.5, 0.0);
    CartesianState state1(pos1, vel1, GM_EARTH);
    
    Vector3d pos2(7100.0, 0.0, 0.0);
    Vector3d vel2(0.0, 7.4, 0.0);
    CartesianState state2(pos2, vel2, GM_EARTH);
    
    CartesianState rel = relative_state(state1, state2);
    
    EXPECT_NEAR(rel.position().x(), -100.0, 1e-10);
    EXPECT_NEAR(rel.velocity().y(), 0.1, 1e-10);
    
    double dist = distance(state1, state2);
    EXPECT_NEAR(dist, 100.0, 1e-6);
}

// ========== KeplerianElements Tests ==========

TEST(KeplerianElementsTest, Construction) {
    double a = 7000.0;
    double e = 0.1;
    double i = 45.0 * DEG_TO_RAD;
    double Omega = 30.0 * DEG_TO_RAD;
    double omega = 60.0 * DEG_TO_RAD;
    double M = 90.0 * DEG_TO_RAD;
    
    KeplerianElements kep(a, e, i, Omega, omega, M, GM_EARTH);
    
    EXPECT_NEAR(kep.semi_major_axis(), a, 1e-10);
    EXPECT_NEAR(kep.eccentricity(), e, 1e-10);
    EXPECT_NEAR(kep.inclination(), i, 1e-10);
}

TEST(KeplerianElementsTest, DerivedQuantities) {
    double a = 10000.0;
    double e = 0.2;
    
    KeplerianElements kep(a, e, 0, 0, 0, 0, GM_EARTH);
    
    double q = kep.periapsis_distance();
    EXPECT_NEAR(q, a * (1.0 - e), 1e-6);
    
    double Q = kep.apoapsis_distance();
    EXPECT_NEAR(Q, a * (1.0 + e), 1e-6);
    
    double p = kep.semi_latus_rectum();
    EXPECT_NEAR(p, a * (1.0 - e * e), 1e-6);
    
    double T = kep.period();
    double expected_T = 2.0 * PI * std::sqrt(a * a * a / GM_EARTH);
    EXPECT_NEAR(T, expected_T, 1.0);
}

TEST(KeplerianElementsTest, AnomalyConversions) {
    double e = 0.5;
    double M = PI / 4.0; // 45 degrees
    
    // M → E
    double E = KeplerianElements::mean_to_eccentric_anomaly(M, e);
    
    // E → M (should recover original)
    double M_back = KeplerianElements::eccentric_to_mean_anomaly(E, e);
    EXPECT_NEAR(M_back, M, 1e-10);
    
    // E → ν
    double nu = KeplerianElements::eccentric_to_true_anomaly(E, e);
    
    // ν → E (should recover)
    double E_back = KeplerianElements::true_to_eccentric_anomaly(nu, e);
    EXPECT_NEAR(E_back, E, 1e-10);
    
    // M → ν → M
    double nu2 = KeplerianElements::mean_to_true_anomaly(M, e);
    double M2 = KeplerianElements::true_to_mean_anomaly(nu2, e);
    EXPECT_NEAR(M2, M, 1e-10);
}

TEST(KeplerianElementsTest, CircularOrbit) {
    KeplerianElements kep(7000.0, 0.0, 0.0, 0.0, 0.0, 0.0, GM_EARTH);
    
    EXPECT_TRUE(kep.is_circular());
    EXPECT_TRUE(kep.is_elliptic());
    EXPECT_FALSE(kep.is_hyperbolic());
}

TEST(KeplerianElementsTest, HyperbolicOrbit) {
    KeplerianElements kep(7000.0, 1.5, 0.0, 0.0, 0.0, 0.0, GM_EARTH);
    
    EXPECT_TRUE(kep.is_hyperbolic());
    EXPECT_FALSE(kep.is_elliptic());
    
    double Q = kep.apoapsis_distance();
    EXPECT_TRUE(std::isinf(Q));
}

// ========== Conversion Tests ==========

TEST(ConversionTest, CircularEquatorialOrbit) {
    // Simple circular equatorial orbit
    double a = 7000.0;
    KeplerianElements kep(a, 0.0, 0.0, 0.0, 0.0, 0.0, GM_EARTH);
    
    CartesianState cart = kep.to_cartesian();
    
    // Should be at periapsis (x = a, y = 0, z = 0)
    EXPECT_NEAR(cart.radius(), a, 1e-6);
    EXPECT_NEAR(cart.position().z(), 0.0, 1e-6);
    
    // Velocity should be circular
    double v_circ = std::sqrt(GM_EARTH / a);
    EXPECT_NEAR(cart.speed(), v_circ, 1e-3);
    
    // Convert back
    KeplerianElements kep_back = KeplerianElements::from_cartesian(cart);
    
    EXPECT_NEAR(kep_back.semi_major_axis(), a, 1.0);
    EXPECT_LT(kep_back.eccentricity(), 1e-3); // Nearly circular
}

TEST(ConversionTest, EllipticOrbit) {
    double a = 10000.0;
    double e = 0.3;
    double i = 30.0 * DEG_TO_RAD;
    
    KeplerianElements kep(a, e, i, 0.0, 0.0, 0.0, GM_EARTH);
    
    CartesianState cart = kep.to_cartesian();
    KeplerianElements kep_back = KeplerianElements::from_cartesian(cart);
    
    EXPECT_NEAR(kep_back.semi_major_axis(), a, 1.0);
    EXPECT_NEAR(kep_back.eccentricity(), e, 1e-6);
    EXPECT_NEAR(kep_back.inclination(), i, 1e-6);
}

TEST(ConversionTest, InclinedOrbit) {
    double a = 7000.0;
    double i = 63.4 * DEG_TO_RAD; // Molniya
    
    KeplerianElements kep(a, 0.1, i, 0.0, 0.0, 0.0, GM_EARTH);
    
    CartesianState cart = kep.to_cartesian();
    KeplerianElements kep_back = KeplerianElements::from_cartesian(cart);
    
    EXPECT_NEAR(kep_back.inclination(), i, 1e-6);
}

TEST(ConversionTest, RoundTrip) {
    // Complex orbit with all elements non-zero
    double a = 8000.0;
    double e = 0.25;
    double i = 45.0 * DEG_TO_RAD;
    double Omega = 30.0 * DEG_TO_RAD;
    double omega = 60.0 * DEG_TO_RAD;
    double M = 120.0 * DEG_TO_RAD;
    
    KeplerianElements kep_orig(a, e, i, Omega, omega, M, GM_EARTH);
    
    // Kep → Cart → Kep
    CartesianState cart = kep_orig.to_cartesian();
    KeplerianElements kep_final = KeplerianElements::from_cartesian(cart);
    
    EXPECT_NEAR(kep_final.semi_major_axis(), a, 1.0);
    EXPECT_NEAR(kep_final.eccentricity(), e, 1e-6);
    EXPECT_NEAR(kep_final.inclination(), i, 1e-6);
    EXPECT_NEAR(kep_final.RAAN(), Omega, 1e-6);
    EXPECT_NEAR(kep_final.argument_of_periapsis(), omega, 1e-6);
    // Mean anomaly might differ slightly due to numerical precision
    EXPECT_NEAR(kep_final.mean_anomaly(), M, 1e-3);
}

// ========== EquinoctialElements Tests ==========

TEST(EquinoctialElementsTest, Construction) {
    double a = 10000.0;  // km
    double h = 0.15;     // e·sin(ω+Ω)
    double k = 0.2;      // e·cos(ω+Ω)
    double p = 0.1;      // tan(i/2)·sin(Ω)
    double q = 0.2;      // tan(i/2)·cos(Ω)
    double lambda = PI / 4.0; // mean longitude
    
    EquinoctialElements eq(a, h, k, p, q, lambda, GM_EARTH);
    
    EXPECT_NEAR(eq.semi_major_axis(), a, 1e-10);
    EXPECT_NEAR(eq.h(), h, 1e-10);
    EXPECT_NEAR(eq.k(), k, 1e-10);
    EXPECT_NEAR(eq.p(), p, 1e-10);
    EXPECT_NEAR(eq.q(), q, 1e-10);
    EXPECT_NEAR(eq.mean_longitude(), lambda, 1e-10);
}

TEST(EquinoctialElementsTest, DerivedElements) {
    // Near-circular orbit
    double a = 7000.0;
    double h = 0.01;  // small eccentricity
    double k = 0.01;
    double p = 0.0;   // equatorial
    double q = 0.0;
    double lambda = 0.0;
    
    EquinoctialElements eq(a, h, k, p, q, lambda, GM_EARTH);
    
    // Check derived elements
    double e_expected = std::sqrt(h*h + k*k);
    EXPECT_NEAR(eq.eccentricity(), e_expected, 1e-10);
    
    double i_expected = 2.0 * std::atan(std::sqrt(p*p + q*q));
    EXPECT_NEAR(eq.inclination(), i_expected, 1e-10);
    
    EXPECT_TRUE(eq.is_circular(0.02));
    EXPECT_TRUE(eq.is_equatorial(1e-5));
    EXPECT_TRUE(eq.is_elliptic());
}

TEST(EquinoctialElementsTest, KeplerianConversion) {
    // Create Keplerian elements
    double a = 8000.0;
    double e = 0.3;
    double i = 30.0 * DEG_TO_RAD;
    double Omega = 45.0 * DEG_TO_RAD;
    double omega = 60.0 * DEG_TO_RAD;
    double M = 90.0 * DEG_TO_RAD;
    
    KeplerianElements kep(a, e, i, Omega, omega, M, GM_EARTH);
    
    // Convert to equinoctial
    EquinoctialElements eq = EquinoctialElements::from_keplerian(kep);
    
    // Verify derived elements match
    EXPECT_NEAR(eq.semi_major_axis(), a, 1e-6);
    EXPECT_NEAR(eq.eccentricity(), e, 1e-6);
    EXPECT_NEAR(eq.inclination(), i, 1e-6);
    EXPECT_NEAR(eq.RAAN(), Omega, 1e-6);
    
    // Convert back to Keplerian
    KeplerianElements kep2 = eq.to_keplerian();
    
    EXPECT_NEAR(kep2.semi_major_axis(), a, 1e-6);
    EXPECT_NEAR(kep2.eccentricity(), e, 1e-6);
    EXPECT_NEAR(kep2.inclination(), i, 1e-6);
    EXPECT_NEAR(kep2.RAAN(), Omega, 1e-6);
}

TEST(EquinoctialElementsTest, CartesianRoundTrip) {
    // Create circular equatorial orbit
    Vector3d pos(7000.0, 0.0, 0.0);
    Vector3d vel(0.0, 7.546, 0.0);
    CartesianState cart_orig(pos, vel, GM_EARTH);
    
    // Cart → Eq → Cart
    EquinoctialElements eq = EquinoctialElements::from_cartesian(cart_orig);
    CartesianState cart_final = eq.to_cartesian();
    
    EXPECT_TRUE(cart_final.position().isApprox(pos, 1e-3));
    EXPECT_TRUE(cart_final.velocity().isApprox(vel, 1e-3));
}

TEST(EquinoctialElementsTest, NoSingularityCircular) {
    // Test that equinoctial handles circular orbits smoothly (e=0)
    double a = 7000.0;
    double e = 0.0;
    double i = 30.0 * DEG_TO_RAD;
    double Omega = 45.0 * DEG_TO_RAD;
    double omega = 60.0 * DEG_TO_RAD;  // Should be well-defined even for e=0
    double M = 90.0 * DEG_TO_RAD;
    
    KeplerianElements kep(a, e, i, Omega, omega, M, GM_EARTH);
    EquinoctialElements eq = EquinoctialElements::from_keplerian(kep);
    
    // h and k should both be zero
    EXPECT_NEAR(eq.h(), 0.0, 1e-10);
    EXPECT_NEAR(eq.k(), 0.0, 1e-10);
    EXPECT_NEAR(eq.eccentricity(), 0.0, 1e-10);
    
    // Can still convert back
    KeplerianElements kep2 = eq.to_keplerian();
    EXPECT_NEAR(kep2.eccentricity(), 0.0, 1e-10);
}

TEST(EquinoctialElementsTest, NoSingularityEquatorial) {
    // Test that equinoctial handles equatorial orbits smoothly (i=0)
    double a = 7000.0;
    double e = 0.2;
    double i = 0.0;
    double Omega = 45.0 * DEG_TO_RAD;  // Should be well-defined even for i=0
    double omega = 60.0 * DEG_TO_RAD;
    double M = 90.0 * DEG_TO_RAD;
    
    KeplerianElements kep(a, e, i, Omega, omega, M, GM_EARTH);
    EquinoctialElements eq = EquinoctialElements::from_keplerian(kep);
    
    // p and q should both be zero
    EXPECT_NEAR(eq.p(), 0.0, 1e-10);
    EXPECT_NEAR(eq.q(), 0.0, 1e-10);
    EXPECT_NEAR(eq.inclination(), 0.0, 1e-10);
    
    // Can still convert back
    KeplerianElements kep2 = eq.to_keplerian();
    EXPECT_NEAR(kep2.inclination(), 0.0, 1e-10);
}

// ========== CometaryElements Tests ==========

TEST(CometaryElementsTest, Construction) {
    double q = 0.5 * AU;  // perihelion at 0.5 AU
    double e = 0.8;
    double i = 45.0 * DEG_TO_RAD;
    double Omega = 90.0 * DEG_TO_RAD;
    double omega = 180.0 * DEG_TO_RAD;
    double T = 2451545.0;  // JD
    
    CometaryElements comet(q, e, i, Omega, omega, T);
    
    EXPECT_NEAR(comet.perihelion_distance(), q, 1e-10);
    EXPECT_NEAR(comet.eccentricity(), e, 1e-10);
    EXPECT_NEAR(comet.inclination(), i, 1e-10);
    EXPECT_NEAR(comet.RAAN(), Omega, 1e-10);
    EXPECT_NEAR(comet.argument_of_periapsis(), omega, 1e-10);
    EXPECT_NEAR(comet.time_of_perihelion(), T, 1e-10);
}

TEST(CometaryElementsTest, EllipticOrbit) {
    double q = 0.5 * AU;
    double e = 0.6;  // elliptic
    double i = 30.0 * DEG_TO_RAD;
    double Omega = 0.0;
    double omega = 0.0;
    double T = 0.0;
    
    CometaryElements comet(q, e, i, Omega, omega, T);
    
    EXPECT_TRUE(comet.is_elliptic());
    EXPECT_FALSE(comet.is_parabolic());
    EXPECT_FALSE(comet.is_hyperbolic());
    EXPECT_TRUE(comet.is_bound());
    
    // Check derived quantities
    double a_expected = q / (1.0 - e);
    EXPECT_NEAR(comet.semi_major_axis(), a_expected, 1e-6);
    
    double Q_expected = q * (1.0 + e) / (1.0 - e);
    EXPECT_NEAR(comet.aphelion_distance(), Q_expected, 1e-6);
    
    EXPECT_GT(comet.period(), 0.0);
    EXPECT_LT(comet.period(), std::numeric_limits<double>::infinity());
}

TEST(CometaryElementsTest, ParabolicOrbit) {
    double q = 1.0 * AU;
    double e = 1.0;  // parabolic
    double i = 60.0 * DEG_TO_RAD;
    double Omega = 0.0;
    double omega = 0.0;
    double T = 0.0;
    
    CometaryElements comet(q, e, i, Omega, omega, T);
    
    EXPECT_FALSE(comet.is_elliptic());
    EXPECT_TRUE(comet.is_parabolic(1e-10));
    EXPECT_FALSE(comet.is_hyperbolic());
    EXPECT_FALSE(comet.is_bound());
    EXPECT_TRUE(comet.is_unbound());
    
    // Semi-major axis should be infinite
    EXPECT_EQ(comet.semi_major_axis(), std::numeric_limits<double>::infinity());
    EXPECT_EQ(comet.period(), std::numeric_limits<double>::infinity());
}

TEST(CometaryElementsTest, HyperbolicOrbit) {
    double q = 2.0 * AU;
    double e = 1.5;  // hyperbolic
    double i = 90.0 * DEG_TO_RAD;
    double Omega = 0.0;
    double omega = 0.0;
    double T = 0.0;
    
    CometaryElements comet(q, e, i, Omega, omega, T);
    
    EXPECT_FALSE(comet.is_elliptic());
    EXPECT_FALSE(comet.is_parabolic());
    EXPECT_TRUE(comet.is_hyperbolic());
    EXPECT_FALSE(comet.is_bound());
    EXPECT_TRUE(comet.is_unbound());
    
    // Semi-major axis should be negative
    EXPECT_LT(comet.semi_major_axis(), 0.0);
    EXPECT_EQ(comet.period(), std::numeric_limits<double>::infinity());
    EXPECT_EQ(comet.aphelion_distance(), std::numeric_limits<double>::infinity());
}

TEST(CometaryElementsTest, KeplerianConversion) {
    // Elliptic comet
    double q = 1.2 * AU;
    double e = 0.7;
    double i = 45.0 * DEG_TO_RAD;
    double Omega = 30.0 * DEG_TO_RAD;
    double omega = 60.0 * DEG_TO_RAD;
    double T = 2451545.0;
    
    CometaryElements comet(q, e, i, Omega, omega, T);
    
    // Convert to Keplerian at perihelion
    KeplerianElements kep = comet.to_keplerian(T);
    
    // At perihelion, mean anomaly should be 0
    EXPECT_NEAR(kep.mean_anomaly(), 0.0, 1e-6);
    EXPECT_NEAR(kep.eccentricity(), e, 1e-10);
    EXPECT_NEAR(kep.inclination(), i, 1e-10);
    EXPECT_NEAR(kep.RAAN(), Omega, 1e-10);
    EXPECT_NEAR(kep.argument_of_periapsis(), omega, 1e-10);
    
    // Convert back
    CometaryElements comet2 = CometaryElements::from_keplerian(kep, T);
    EXPECT_NEAR(comet2.perihelion_distance(), q, 1e-6);
    EXPECT_NEAR(comet2.eccentricity(), e, 1e-10);
}

TEST(CometaryElementsTest, CartesianConversion) {
    // Create cometary orbit
    double q = 0.8 * AU;
    double e = 0.5;
    double i = 20.0 * DEG_TO_RAD;
    double Omega = 0.0;
    double omega = 0.0;
    double T = 0.0;
    
    CometaryElements comet(q, e, i, Omega, omega, T);
    
    // Convert to Cartesian at perihelion
    CartesianState cart = comet.to_cartesian(T);
    
    // At perihelion: r = q, v perpendicular
    EXPECT_NEAR(cart.radius(), q, 1e-3);
    
    // Velocity should be √(μ(1+e)/q)
    double v_peri = std::sqrt(GM_SUN * (1.0 + e) / q);
    EXPECT_NEAR(cart.speed(), v_peri, 1e-3);
}

