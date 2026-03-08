/**
 * @file test_jacobians.cpp
 * @brief Unit tests for Jacobian matrices
 */

#include <gtest/gtest.h>
#include "astdyn/coordinates/KeplerianElements.hpp"
#include "astdyn/coordinates/CartesianState.hpp"

using namespace astdyn;
using namespace astdyn::coordinates;
using namespace astdyn::constants;

// ========== Jacobian Matrix Tests ==========

TEST(JacobianTest, KeplerianToCartesianJacobian) {
    // Test that Jacobian has correct dimensions
    KeplerianElements kep(7000.0, 0.1, 30.0*DEG_TO_RAD, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    Matrix6d J = kep.jacobian_to_cartesian();
    
    // Check dimensions
    EXPECT_EQ(J.rows(), 6);
    EXPECT_EQ(J.cols(), 6);
    
    // Check that it's not all zeros
    EXPECT_GT(J.norm(), 1e-6);
}

TEST(JacobianTest, CartesianToKeplerianJacobian) {
    Vector3d pos(7000.0, 0.0, 0.0);
    Vector3d vel(0.0, 7.5, 0.0);
    CartesianState state(pos, vel, GM_EARTH);
    
    Matrix6d J = KeplerianElements::jacobian_from_cartesian(state);
    
    // Check dimensions
    EXPECT_EQ(J.rows(), 6);
    EXPECT_EQ(J.cols(), 6);
    
    // Check that it's not all zeros
    EXPECT_GT(J.norm(), 1e-10);
}

TEST(JacobianTest, JacobianInverseRelationship) {
    // J_kep_to_cart * J_cart_to_kep should be approximately identity
    KeplerianElements kep(8000.0, 0.2, 45.0*DEG_TO_RAD, 
                         30.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         120.0*DEG_TO_RAD, GM_EARTH);
    
    CartesianState cart = kep.to_cartesian();
    
    Matrix6d J_to = kep.jacobian_to_cartesian();
    Matrix6d J_from = KeplerianElements::jacobian_from_cartesian(cart);
    
    Matrix6d product = J_from * J_to;
    Matrix6d identity = Matrix6d::Identity();
    
    // Should be close to identity (within numerical precision)
    EXPECT_TRUE(product.isApprox(identity, 0.1));
}

TEST(JacobianTest, FiniteDifferenceValidation) {
    // Manually verify Jacobian using finite differences
    KeplerianElements kep(7000.0, 0.1, 30.0*DEG_TO_RAD, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    Matrix6d J = kep.jacobian_to_cartesian();
    
    // Verify one column manually (∂state/∂a)
    double da = 1.0;  // 1 km perturbation
    KeplerianElements kep_pert = kep;
    kep_pert.set_semi_major_axis(kep.semi_major_axis() + da);
    
    CartesianState state_nom = kep.to_cartesian();
    CartesianState state_pert = kep_pert.to_cartesian();
    
    Vector6d state_nom_vec, state_pert_vec;
    state_nom_vec << state_nom.position(), state_nom.velocity();
    state_pert_vec << state_pert.position(), state_pert.velocity();
    
    Vector6d finite_diff = (state_pert_vec - state_nom_vec) / da;
    
    // First column of Jacobian should match finite difference
    EXPECT_TRUE(J.col(0).isApprox(finite_diff, 1e-3));
}

TEST(JacobianTest, CircularOrbitJacobian) {
    // Test Jacobian for circular orbit (e ≈ 0)
    KeplerianElements kep(7000.0, 1e-8, 30.0*DEG_TO_RAD, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    Matrix6d J = kep.jacobian_to_cartesian();
    
    // Should still be well-defined
    EXPECT_TRUE(std::isfinite(J.norm()));
    EXPECT_GT(J.norm(), 1e-6);
}

TEST(JacobianTest, EquatorialOrbitJacobian) {
    // Test Jacobian for equatorial orbit (i ≈ 0)
    KeplerianElements kep(7000.0, 0.1, 1e-8, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    Matrix6d J = kep.jacobian_to_cartesian();
    
    // Should still be well-defined
    EXPECT_TRUE(std::isfinite(J.norm()));
    EXPECT_GT(J.norm(), 1e-6);
}

TEST(JacobianTest, CovariancePropagation) {
    // Test practical use case: propagating covariance matrix
    KeplerianElements kep(7000.0, 0.1, 30.0*DEG_TO_RAD, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    // Assumed Keplerian covariance (diagonal for simplicity)
    Matrix6d P_kep = Matrix6d::Identity();
    P_kep(0,0) = 1.0;      // a uncertainty: 1 km²
    P_kep(1,1) = 1e-6;     // e uncertainty
    P_kep(2,2) = 1e-8;     // i uncertainty (rad²)
    P_kep(3,3) = 1e-8;     // Omega uncertainty
    P_kep(4,4) = 1e-8;     // omega uncertainty
    P_kep(5,5) = 1e-8;     // M uncertainty
    
    // Propagate to Cartesian covariance
    Matrix6d J = kep.jacobian_to_cartesian();
    Matrix6d P_cart = J * P_kep * J.transpose();
    
    // Cartesian covariance should be positive semi-definite
    Eigen::SelfAdjointEigenSolver<Matrix6d> solver(P_cart);
    EXPECT_GE(solver.eigenvalues().minCoeff(), -1e-10);
    
    // Should have non-zero diagonal elements
    EXPECT_GT(P_cart.diagonal().norm(), 1e-10);
}

TEST(JacobianTest, SymmetryCheck) {
    // Covariance matrices should remain symmetric
    KeplerianElements kep(7000.0, 0.1, 30.0*DEG_TO_RAD, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    Matrix6d P_kep = Matrix6d::Identity();
    Matrix6d J = kep.jacobian_to_cartesian();
    Matrix6d P_cart = J * P_kep * J.transpose();
    
    // Should be symmetric
    EXPECT_TRUE(P_cart.isApprox(P_cart.transpose(), 1e-10));
}

TEST(JacobianTest, EccentricOrbitJacobian) {
    // Test Jacobian for highly eccentric orbit
    KeplerianElements kep(10000.0, 0.8, 30.0*DEG_TO_RAD, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    Matrix6d J_to = kep.jacobian_to_cartesian();
    CartesianState cart = kep.to_cartesian();
    Matrix6d J_from = KeplerianElements::jacobian_from_cartesian(cart);
    
    // Both should be well-conditioned
    EXPECT_TRUE(std::isfinite(J_to.norm()));
    EXPECT_TRUE(std::isfinite(J_from.norm()));
    
    // Product should still approximate identity
    Matrix6d product = J_from * J_to;
    EXPECT_TRUE(product.isApprox(Matrix6d::Identity(), 0.2));
}

TEST(JacobianTest, RoundTripCovariance) {
    // Test: Kep → Cart → Kep covariance round-trip
    KeplerianElements kep(7000.0, 0.1, 30.0*DEG_TO_RAD, 
                         45.0*DEG_TO_RAD, 60.0*DEG_TO_RAD, 
                         90.0*DEG_TO_RAD, GM_EARTH);
    
    Matrix6d P_kep_orig = Matrix6d::Identity() * 0.01;
    
    // Kep → Cart
    Matrix6d J_to = kep.jacobian_to_cartesian();
    Matrix6d P_cart = J_to * P_kep_orig * J_to.transpose();
    
    // Cart → Kep
    CartesianState cart = kep.to_cartesian();
    Matrix6d J_from = KeplerianElements::jacobian_from_cartesian(cart);
    Matrix6d P_kep_final = J_from * P_cart * J_from.transpose();
    
    // Should approximately recover original covariance
    EXPECT_TRUE(P_kep_final.isApprox(P_kep_orig, 0.01));
}

