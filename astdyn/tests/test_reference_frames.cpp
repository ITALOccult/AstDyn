/**
 * @file test_reference_frames.cpp
 * @brief Unit tests for reference frame transformations
 */

#include <gtest/gtest.h>
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/coordinates/CartesianState.hpp"

using namespace astdyn;
using namespace astdyn::coordinates;
using namespace astdyn::constants;

// ========== Rotation Matrix Tests ==========

TEST(ReferenceFrameTest, RotationMatrices) {
    // Test that rotation matrices are orthogonal
    double angle = PI / 4.0; // 45 degrees
    
    Matrix3d Rx = ReferenceFrame::rotation_x(angle);
    Matrix3d Ry = ReferenceFrame::rotation_y(angle);
    Matrix3d Rz = ReferenceFrame::rotation_z(angle);
    

    
    // Check orthogonality: R * R^T = I
    EXPECT_TRUE((Rx * Rx.transpose()).isApprox(Matrix3d::Identity(), 1e-10));
    EXPECT_TRUE((Ry * Ry.transpose()).isApprox(Matrix3d::Identity(), 1e-10));
    EXPECT_TRUE((Rz * Rz.transpose()).isApprox(Matrix3d::Identity(), 1e-10));
    
    // Check determinant = 1
    EXPECT_NEAR(Rx.determinant(), 1.0, 1e-10);
    EXPECT_NEAR(Ry.determinant(), 1.0, 1e-10);
    EXPECT_NEAR(Rz.determinant(), 1.0, 1e-10);
}

TEST(ReferenceFrameTest, RotationXAxis) {
    // Rotation about X-axis should not change X-component
    Vector3d vec(1.0, 2.0, 3.0);
    double angle = PI / 3.0;
    
    Matrix3d R = ReferenceFrame::rotation_x(angle);
    Vector3d rotated = R * vec;
    
    EXPECT_NEAR(rotated.x(), vec.x(), 1e-10);
}

TEST(ReferenceFrameTest, RotationYAxis) {
    // Rotation about Y-axis should not change Y-component
    Vector3d vec(1.0, 2.0, 3.0);
    double angle = PI / 3.0;
    
    Matrix3d R = ReferenceFrame::rotation_y(angle);
    Vector3d rotated = R * vec;
    
    EXPECT_NEAR(rotated.y(), vec.y(), 1e-10);
}

TEST(ReferenceFrameTest, RotationZAxis) {
    // Rotation about Z-axis should not change Z-component
    Vector3d vec(1.0, 2.0, 3.0);
    double angle = PI / 3.0;
    
    Matrix3d R = ReferenceFrame::rotation_z(angle);
    Vector3d rotated = R * vec;
    
    EXPECT_NEAR(rotated.z(), vec.z(), 1e-10);
}

// ========== Frame Bias Tests (J2000 ↔ ICRS) ==========

TEST(ReferenceFrameTest, J2000ToICRS) {
    Matrix3d bias = ReferenceFrame::j2000_to_icrs();
    
    // Frame bias should be very close to identity (offset ~0.02 arcsec)
    EXPECT_TRUE(bias.isApprox(Matrix3d::Identity(), 1e-5));
    
    // Should be orthogonal
    EXPECT_TRUE((bias * bias.transpose()).isApprox(Matrix3d::Identity(), 1e-10));
}

TEST(ReferenceFrameTest, ICRSToJ2000Inverse) {
    Matrix3d j2000_to_icrs = ReferenceFrame::j2000_to_icrs();
    Matrix3d icrs_to_j2000 = ReferenceFrame::icrs_to_j2000();
    
    // Should be inverses
    EXPECT_TRUE((j2000_to_icrs * icrs_to_j2000).isApprox(Matrix3d::Identity(), 1e-10));
}

TEST(ReferenceFrameTest, J2000ICRSRoundTrip) {
    Vector3d pos_j2000(7000.0, 0.0, 0.0);

    // J2000 → ICRS → J2000 via rotation matrices
    Matrix3d R_fwd = ReferenceFrame::get_transformation(FrameType::J2000, FrameType::ICRS);
    Matrix3d R_inv = ReferenceFrame::get_transformation(FrameType::ICRS, FrameType::J2000);
    Vector3d pos_icrs  = R_fwd * pos_j2000;
    Vector3d pos_final = R_inv * pos_icrs;

    EXPECT_TRUE(pos_final.isApprox(pos_j2000, 1e-6));
}

// ========== Ecliptic Transformations ==========

TEST(ReferenceFrameTest, J2000ToEcliptic) {
    Matrix3d R = ReferenceFrame::j2000_to_ecliptic();
    
    // Should be rotation about X-axis by obliquity
    Matrix3d expected = ReferenceFrame::rotation_x(constants::OBLIQUITY_J2000);
    
    EXPECT_TRUE(R.isApprox(expected, 1e-10));
}

TEST(ReferenceFrameTest, EclipticToJ2000Inverse) {
    Matrix3d to_ecliptic = ReferenceFrame::j2000_to_ecliptic();
    Matrix3d from_ecliptic = ReferenceFrame::ecliptic_to_j2000();
    
    // Should be inverses
    EXPECT_TRUE((to_ecliptic * from_ecliptic).isApprox(Matrix3d::Identity(), 1e-10));
}

TEST(ReferenceFrameTest, EclipticCoordinates) {
    // Point on equator in J2000 should have non-zero Z in ecliptic
    Vector3d pos_j2000(AU, 0.0, 0.0);

    Matrix3d R = ReferenceFrame::get_transformation(FrameType::J2000, FrameType::ECLIPTIC);
    Vector3d pos_ecliptic = R * pos_j2000;

    // X-component should be unchanged (rotation about X)
    EXPECT_NEAR(pos_ecliptic.x(), pos_j2000.x(), 1e-6);

    // Should preserve magnitude
    EXPECT_NEAR(pos_ecliptic.norm(), pos_j2000.norm(), 1e-6);
}

TEST(ReferenceFrameTest, EclipticRoundTrip) {
    Vector3d pos_j2000(1.5e8, 5.0e7, 2.0e7); // Arbitrary position

    Matrix3d R_fwd = ReferenceFrame::get_transformation(FrameType::J2000, FrameType::ECLIPTIC);
    Matrix3d R_inv = ReferenceFrame::get_transformation(FrameType::ECLIPTIC, FrameType::J2000);
    Vector3d pos_ecliptic = R_fwd * pos_j2000;
    Vector3d pos_final    = R_inv * pos_ecliptic;

    EXPECT_TRUE(pos_final.isApprox(pos_j2000, 1e-3));
}

// ========== ITRF Transformations ==========

TEST(ReferenceFrameTest, J2000ToITRFSimple) {
    // At J2000.0 epoch
    
    Matrix3d R = ReferenceFrame::j2000_to_itrf_simple(time::EpochUTC::from_mjd(MJD2000));
    
    // Should be orthogonal
    EXPECT_TRUE((R * R.transpose()).isApprox(Matrix3d::Identity(), 1e-10));
    
    // Determinant should be 1
    EXPECT_NEAR(R.determinant(), 1.0, 1e-10);
}

TEST(ReferenceFrameTest, ITRFToJ2000Inverse) {
    time::EpochUTC mjd_instant = time::EpochUTC::from_mjd(MJD2000 + 1.0); // One day after J2000
    
    Matrix3d to_itrf = ReferenceFrame::j2000_to_itrf_simple(mjd_instant);
    Matrix3d from_itrf = ReferenceFrame::itrf_to_j2000_simple(mjd_instant);
    
    // Should be inverses
    EXPECT_TRUE((to_itrf * from_itrf).isApprox(Matrix3d::Identity(), 1e-10));
}

TEST(ReferenceFrameTest, ITRFRotation) {
    // Position at different times should show Earth rotation
    Vector3d pos_j2000(7000.0, 0.0, 0.0);

    time::EpochTDB mjd1_instant = time::EpochTDB::from_mjd(MJD2000);
    time::EpochTDB mjd2_instant = time::EpochTDB::from_mjd(MJD2000 + 0.25); // 6 hours later

    Matrix3d R1 = ReferenceFrame::get_transformation(FrameType::J2000, FrameType::ITRF, mjd1_instant);
    Matrix3d R2 = ReferenceFrame::get_transformation(FrameType::J2000, FrameType::ITRF, mjd2_instant);
    Vector3d pos_itrf1 = R1 * pos_j2000;
    Vector3d pos_itrf2 = R2 * pos_j2000;

    // Positions should be different (Earth has rotated)
    EXPECT_GT((pos_itrf2 - pos_itrf1).norm(), 1000.0);

    // But magnitude should be preserved
    EXPECT_NEAR(pos_itrf1.norm(), pos_j2000.norm(), 1e-6);
    EXPECT_NEAR(pos_itrf2.norm(), pos_j2000.norm(), 1e-6);
}

// ========== Generic Transformation Tests ==========

TEST(ReferenceFrameTest, IdentityTransformation) {
    // Same frame should return identity
    Matrix3d R = ReferenceFrame::get_transformation(
        FrameType::J2000, FrameType::J2000);
    
    EXPECT_TRUE(R.isApprox(Matrix3d::Identity(), 1e-10));
}

TEST(ReferenceFrameTest, ChainedTransformation) {
    // ICRS → ECLIPTIC should equal ICRS → J2000 → ECLIPTIC
    Vector3d pos_icrs(AU, 0.0, 0.0);

    // Direct transformation
    Matrix3d R_direct = ReferenceFrame::get_transformation(FrameType::ICRS, FrameType::ECLIPTIC);
    Vector3d pos_ecliptic1 = R_direct * pos_icrs;

    // Chained through J2000
    Matrix3d R_to_j2000   = ReferenceFrame::get_transformation(FrameType::ICRS, FrameType::J2000);
    Matrix3d R_to_ecl     = ReferenceFrame::get_transformation(FrameType::J2000, FrameType::ECLIPTIC);
    Vector3d pos_j2000    = R_to_j2000 * pos_icrs;
    Vector3d pos_ecliptic2 = R_to_ecl * pos_j2000;

    EXPECT_TRUE(pos_ecliptic1.isApprox(pos_ecliptic2, 1e-3));
}

// ========== State Transformation Tests ==========

// Helper: apply a rotation matrix to both position and velocity of a CartesianState.
// For rotating frames (ITRF) the Coriolis term is not included here;
// use the typed ReferenceFrame::transform_vel<F,T>() API for that.
static CartesianState apply_rotation(const CartesianState& s, const Matrix3d& R) {
    return CartesianState(R * s.position(), R * s.velocity(), s.mu());
}

TEST(ReferenceFrameTest, StateTransformationPosition) {
    // Create a state in J2000
    Vector3d pos(7000.0, 0.0, 0.0);
    Vector3d vel(0.0, 7.5, 0.0);
    CartesianState state_j2000(pos, vel, GM_EARTH);

    // Transform to ICRS via rotation matrix
    Matrix3d R = ReferenceFrame::get_transformation(FrameType::J2000, FrameType::ICRS);
    CartesianState state_icrs = apply_rotation(state_j2000, R);

    // Position should change slightly (frame bias ~0.02 arcsec)
    EXPECT_TRUE(state_icrs.position().isApprox(pos, 1e-2));

    // Magnitude should be preserved
    EXPECT_NEAR(state_icrs.radius(), state_j2000.radius(), 1e-6);
}

TEST(ReferenceFrameTest, StateTransformationVelocity) {
    Vector3d pos(7000.0, 0.0, 0.0);
    Vector3d vel(0.0, 7.5, 0.0);
    CartesianState state_j2000(pos, vel, GM_EARTH);

    // Transform to Ecliptic via rotation matrix
    Matrix3d R = ReferenceFrame::get_transformation(FrameType::J2000, FrameType::ECLIPTIC);
    CartesianState state_ecliptic = apply_rotation(state_j2000, R);

    // Speed should be preserved (pure rotation, no Coriolis)
    EXPECT_NEAR(state_ecliptic.speed(), state_j2000.speed(), 1e-6);
}

TEST(ReferenceFrameTest, StateTransformationRoundTrip) {
    Vector3d pos(7000.0, 3000.0, 1000.0);
    Vector3d vel(2.0, 5.0, -3.0);
    CartesianState state_orig(pos, vel, GM_EARTH);

    // J2000 → ECLIPTIC → J2000
    Matrix3d R_fwd = ReferenceFrame::get_transformation(FrameType::J2000, FrameType::ECLIPTIC);
    Matrix3d R_inv = ReferenceFrame::get_transformation(FrameType::ECLIPTIC, FrameType::J2000);
    CartesianState state_ecliptic = apply_rotation(state_orig, R_fwd);
    CartesianState state_final    = apply_rotation(state_ecliptic, R_inv);

    EXPECT_TRUE(state_final.position().isApprox(pos, 1e-3));
    EXPECT_TRUE(state_final.velocity().isApprox(vel, 1e-3));
}

TEST(ReferenceFrameTest, ITRFVelocityTransformation) {
    // Test that ITRF velocity differs from inertial (rotation effect).
    // Use the typed API which correctly handles the Coriolis term.
    using namespace astdyn::core;
    using namespace astdyn::physics;
    using namespace astdyn::math;

    time::EpochTDB t = time::EpochTDB::from_mjd(MJD2000);

    // Construct typed state in GCRF (≈J2000)
    auto pos_gcrf = Vector3<GCRF, Distance>::from_si(7000.0e3, 0.0, 0.0); // 7000 km in meters
    auto vel_gcrf = Vector3<GCRF, Velocity>::from_si(0.0, 7500.0, 0.0);   // 7.5 km/s in m/s

    // Transform velocity to ITRF (includes Coriolis)
    auto vel_itrf = ReferenceFrame::transform_vel<GCRF, ITRF>(pos_gcrf, vel_gcrf, t);

    // Speed in ITRF should differ from GCRF speed due to Earth rotation
    double speed_gcrf = vel_gcrf.norm().to_ms();
    double speed_itrf = vel_itrf.norm().to_ms();
    EXPECT_NE(speed_gcrf, speed_itrf);
}

// ========== Utility Function Tests ==========

TEST(ReferenceFrameTest, IsInertial) {
    EXPECT_TRUE(ReferenceFrame::is_inertial(FrameType::J2000));
    EXPECT_TRUE(ReferenceFrame::is_inertial(FrameType::ICRS));
    EXPECT_TRUE(ReferenceFrame::is_inertial(FrameType::ECLIPTIC));
    EXPECT_FALSE(ReferenceFrame::is_inertial(FrameType::ITRF));
}

TEST(ReferenceFrameTest, IsRotating) {
    EXPECT_FALSE(ReferenceFrame::is_rotating(FrameType::J2000));
    EXPECT_FALSE(ReferenceFrame::is_rotating(FrameType::ICRS));
    EXPECT_FALSE(ReferenceFrame::is_rotating(FrameType::ECLIPTIC));
    EXPECT_TRUE(ReferenceFrame::is_rotating(FrameType::ITRF));
}

TEST(ReferenceFrameTest, GMSTCalculation) {
    // GMST at J2000.0 should be approximately 6h 41m 50.5s = 100.4606° = 1.753368 rad
    double gmst0 = ReferenceFrame::gmst(time::EpochUTC::from_mjd(MJD2000));
    
    // Check it's in valid range [0, 2π)
    EXPECT_GE(gmst0, 0.0);
    EXPECT_LT(gmst0, 2.0 * PI);
    
    // GMST should increase with time
    double gmst_later = ReferenceFrame::gmst(time::EpochUTC::from_mjd(MJD2000 + 1.0));
    EXPECT_NE(gmst0, gmst_later);
}

TEST(ReferenceFrameTest, FrameTypeToString) {
    EXPECT_EQ(frame_type_to_string(FrameType::J2000), "J2000");
    EXPECT_EQ(frame_type_to_string(FrameType::ICRS), "ICRS");
    EXPECT_EQ(frame_type_to_string(FrameType::ECLIPTIC), "ECLIPTIC");
    EXPECT_EQ(frame_type_to_string(FrameType::ITRF), "ITRF");
}

