/**
 * @file test_chebyshev.cpp
 * @brief Unit tests for ChebyshevApproximation
 */

#include <gtest/gtest.h>
#include "chebyshev_approximation.h"
#include <Eigen/Dense>
#include <vector>
#include <cmath>

using namespace ioccultcalc;

class ChebyshevTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup a simple parabolic trajectory
        // x = t^2, y = 2t, z = 0.5t^2 + t
        start_epoch = 61000.0;
        end_epoch = 61010.0;
        num_points = 11;
        
        for (int i = 0; i < num_points; ++i) {
            double t = start_epoch + i;
            double rel_t = i;  // relative time for simplicity
            positions.push_back(Eigen::Vector3d(
                rel_t * rel_t,
                2.0 * rel_t,
                0.5 * rel_t * rel_t + rel_t
            ));
            epochs.push_back(t);
        }
    }

    double start_epoch;
    double end_epoch;
    int num_points;
    std::vector<Eigen::Vector3d> positions;
    std::vector<double> epochs;
};

TEST_F(ChebyshevTest, FitAndEvaluate) {
    // We need at least degree 2 to fit a parabola exactly
    // num_coefficients = degree + 1, so 3 coefficients
    ChebyshevApproximation approx(4); // Use 4 to be safe and test higher order
    
    ASSERT_TRUE(approx.fit(positions, start_epoch, end_epoch));
    
    // Evaluate at midpoints
    for (int i = 0; i < num_points - 1; ++i) {
        double t = start_epoch + i + 0.5;
        double rel_t = i + 0.5;
        
        Eigen::Vector3d expected(
            rel_t * rel_t,
            2.0 * rel_t,
            0.5 * rel_t * rel_t + rel_t
        );
        
        Eigen::Vector3d actual = approx.evaluatePosition(t);
        
        EXPECT_NEAR(expected.x(), actual.x(), 1e-10);
        EXPECT_NEAR(expected.y(), actual.y(), 1e-10);
        EXPECT_NEAR(expected.z(), actual.z(), 1e-10);
    }
}

TEST_F(ChebyshevTest, VelocityEvaluation) {
    ChebyshevApproximation approx(4);
    ASSERT_TRUE(approx.fit(positions, start_epoch, end_epoch));
    
    // x = t^2 -> dx/dt = 2t
    // y = 2t  -> dy/dt = 2
    // z = 0.5t^2 + t -> dz/dt = t + 1
    
    for (int i = 0; i < num_points - 1; ++i) {
        double t = start_epoch + i + 0.5;
        double rel_t = i + 0.5;
        
        Eigen::Vector3d expected_vel(
            2.0 * rel_t,
            2.0,
            rel_t + 1.0
        );
        
        Eigen::Vector3d actual_vel = approx.evaluateVelocity(t);
        
        EXPECT_NEAR(expected_vel.x(), actual_vel.x(), 1e-9);
        EXPECT_NEAR(expected_vel.y(), actual_vel.y(), 1e-9);
        EXPECT_NEAR(expected_vel.z(), actual_vel.z(), 1e-9);
    }
}

TEST_F(ChebyshevTest, Serialization) {
    ChebyshevApproximation approx(4);
    ASSERT_TRUE(approx.fit(positions, start_epoch, end_epoch));
    
    std::string test_file = "test_coeffs.txt";
    ASSERT_TRUE(approx.saveToFile(test_file));
    
    ChebyshevApproximation approx2(4);
    ASSERT_TRUE(approx2.loadFromFile(test_file));
    
    double t = start_epoch + 5.5;
    Eigen::Vector3d p1 = approx.evaluatePosition(t);
    Eigen::Vector3d p2 = approx2.evaluatePosition(t);
    
    EXPECT_DOUBLE_EQ(p1.x(), p2.x());
    EXPECT_DOUBLE_EQ(p1.y(), p2.y());
    EXPECT_DOUBLE_EQ(p1.z(), p2.z());
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
