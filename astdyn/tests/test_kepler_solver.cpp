/**
 * @file test_kepler_solver.cpp
 * @brief Unit tests for Kepler equation solvers (elliptic, hyperbolic, parabolic edge cases)
 */

#include <gtest/gtest.h>
#include "astdyn/math/kepler_solver.hpp"
#include <cmath>

using namespace astdyn::math;

// ============================================================================
// Elliptic solver
// ============================================================================

TEST(KeplerEllipticTest, CircularOrbit) {
    KeplerSolverOptions opts{0.0};
    auto E = solve_kepler_elliptic(Radian(1.2345), opts);
    ASSERT_TRUE(E.has_value());
    // For e=0: E == M
    EXPECT_NEAR(E->value, 1.2345, 1e-10);
}

TEST(KeplerEllipticTest, EarthEccentricity) {
    // Earth: e = 0.0167, M = 1 rad
    KeplerSolverOptions opts{0.0167};
    auto E = solve_kepler_elliptic(Radian(1.0), opts);
    ASSERT_TRUE(E.has_value());
    // Verify: E - e*sin(E) == M
    double residual = E->value - 0.0167 * std::sin(E->value) - 1.0;
    EXPECT_NEAR(residual, 0.0, 1e-11);
}

TEST(KeplerEllipticTest, HighEccentricity) {
    // e = 0.95 (near-parabolic), M = 0.1 rad
    KeplerSolverOptions opts{0.95};
    auto E = solve_kepler_elliptic(Radian(0.1), opts);
    ASSERT_TRUE(E.has_value());
    double residual = E->value - 0.95 * std::sin(E->value) - 0.1;
    EXPECT_NEAR(residual, 0.0, 1e-11);
}

TEST(KeplerEllipticTest, RejectsHyperbolicEccentricity) {
    KeplerSolverOptions opts{1.5};
    auto E = solve_kepler_elliptic(Radian(1.0), opts);
    EXPECT_FALSE(E.has_value());
}

TEST(KeplerEllipticTest, RejectsParabolicEccentricity) {
    KeplerSolverOptions opts{1.0};
    auto E = solve_kepler_elliptic(Radian(1.0), opts);
    EXPECT_FALSE(E.has_value());
}

TEST(KeplerEllipticTest, ZeroMeanAnomaly) {
    KeplerSolverOptions opts{0.5};
    auto E = solve_kepler_elliptic(Radian(0.0), opts);
    ASSERT_TRUE(E.has_value());
    EXPECT_NEAR(E->value, 0.0, 1e-12);
}

TEST(KeplerEllipticTest, FullCircle) {
    // M = 2π should wrap to E ≈ 0
    KeplerSolverOptions opts{0.3};
    auto E = solve_kepler_elliptic(Radian(2.0 * M_PI), opts);
    ASSERT_TRUE(E.has_value());
    double residual = E->value - 0.3 * std::sin(E->value);
    EXPECT_NEAR(std::fmod(residual, 2.0 * M_PI), 0.0, 1e-10);
}

// ============================================================================
// Hyperbolic solver
// ============================================================================

TEST(KeplerHyperbolicTest, BasicConvergence) {
    // e = 2.0, Mh = 1.0 — should converge
    KeplerSolverOptions opts{2.0};
    auto H = solve_kepler_hyperbolic(1.0, opts);
    ASSERT_TRUE(H.has_value());
    double residual = 2.0 * std::sinh(*H) - *H - 1.0;
    EXPECT_NEAR(residual, 0.0, 1e-11);
}

TEST(KeplerHyperbolicTest, SmallMhNearPerihelion) {
    // Mh ≈ 0: H ≈ Mh/(e-1)
    KeplerSolverOptions opts{1.5};
    auto H = solve_kepler_hyperbolic(0.01, opts);
    ASSERT_TRUE(H.has_value());
    double residual = 1.5 * std::sinh(*H) - *H - 0.01;
    EXPECT_NEAR(residual, 0.0, 1e-11);
}

TEST(KeplerHyperbolicTest, LargeMh) {
    // Large Mh tests logarithmic initial guess
    KeplerSolverOptions opts{1.2};
    auto H = solve_kepler_hyperbolic(50.0, opts);
    ASSERT_TRUE(H.has_value());
    double residual = 1.2 * std::sinh(*H) - *H - 50.0;
    EXPECT_NEAR(residual, 0.0, 1e-8);
}

TEST(KeplerHyperbolicTest, NegativeMh) {
    // Negative hyperbolic mean anomaly (approaching from other side)
    KeplerSolverOptions opts{2.0};
    auto H_pos = solve_kepler_hyperbolic( 2.0, opts);
    auto H_neg = solve_kepler_hyperbolic(-2.0, opts);
    ASSERT_TRUE(H_pos.has_value());
    ASSERT_TRUE(H_neg.has_value());
    EXPECT_NEAR(*H_pos, -*H_neg, 1e-11); // Antisymmetry
}

TEST(KeplerHyperbolicTest, RejectsEllipticEccentricity) {
    KeplerSolverOptions opts{0.9};
    auto H = solve_kepler_hyperbolic(1.0, opts);
    EXPECT_FALSE(H.has_value());
}

TEST(KeplerHyperbolicTest, RejectsParabolicEccentricity) {
    KeplerSolverOptions opts{1.0};
    auto H = solve_kepler_hyperbolic(1.0, opts);
    EXPECT_FALSE(H.has_value());
}

TEST(KeplerHyperbolicTest, ZeroMh) {
    // H == 0 when Mh == 0 (at perihelion)
    KeplerSolverOptions opts{3.0};
    auto H = solve_kepler_hyperbolic(0.0, opts);
    ASSERT_TRUE(H.has_value());
    EXPECT_NEAR(*H, 0.0, 1e-12);
}
