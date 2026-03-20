/**
 * @file test_types.cpp
 * @brief Unit tests for type definitions and utilities
 */

#include <gtest/gtest.h>
#include <astdyn/core/Types.hpp>
#include <cmath>

using namespace astdyn;

class TypesTest : public ::testing::Test {
protected:
    static constexpr double TOLERANCE = 1.0e-10;
};

// ============================================================================
// Eigen Vector Tests
// ============================================================================

TEST_F(TypesTest, Vector3dCreation) {
    Vector3d v(1.0, 2.0, 3.0);
    EXPECT_DOUBLE_EQ(v(0), 1.0);
    EXPECT_DOUBLE_EQ(v(1), 2.0);
    EXPECT_DOUBLE_EQ(v(2), 3.0);
}

TEST_F(TypesTest, Vector6dCreation) {
    Vector6d state;
    state << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    
    EXPECT_EQ(state.size(), 6);
    EXPECT_DOUBLE_EQ(state(0), 1.0);
    EXPECT_DOUBLE_EQ(state(5), 6.0);
}

TEST_F(TypesTest, Matrix3dCreation) {
    Matrix3d m = Matrix3d::Identity();
    
    EXPECT_DOUBLE_EQ(m(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(m(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(m(2, 2), 1.0);
    EXPECT_DOUBLE_EQ(m(0, 1), 0.0);
}

// ============================================================================
// Special Values Tests
// ============================================================================

TEST_F(TypesTest, NaNDetection) {
    double nan_value = NaN;
    double normal_value = 1.0;
    
    EXPECT_TRUE(isNaN(nan_value));
    EXPECT_FALSE(isNaN(normal_value));
}

TEST_F(TypesTest, InfinityDetection) {
    double inf_value = Infinity;
    double normal_value = 1.0;
    
    EXPECT_FALSE(isFinite(inf_value));
    EXPECT_TRUE(isFinite(normal_value));
}

TEST_F(TypesTest, FiniteDetection) {
    EXPECT_TRUE(isFinite(0.0));
    EXPECT_TRUE(isFinite(1.0));
    EXPECT_TRUE(isFinite(-1.0));
    EXPECT_FALSE(isFinite(Infinity));
    EXPECT_FALSE(isFinite(-Infinity));
    EXPECT_FALSE(isFinite(NaN));
}

// ============================================================================
// Result Type Tests
// ============================================================================

TEST_F(TypesTest, ResultSuccess) {
    auto result = Result<double>::Success(42.0);
    
    EXPECT_TRUE(result.success);
    EXPECT_TRUE(result); // Test implicit bool conversion
    EXPECT_DOUBLE_EQ(result.value, 42.0);
    EXPECT_TRUE(result.error_message.empty());
}

TEST_F(TypesTest, ResultFailure) {
    auto result = Result<double>::Failure("Something went wrong");
    
    EXPECT_FALSE(result.success);
    EXPECT_FALSE(result); // Test implicit bool conversion
    EXPECT_EQ(result.error_message, "Something went wrong");
}

TEST_F(TypesTest, ResultWithComplexType) {
    auto result = Result<Vector3d>::Success(Vector3d(1.0, 2.0, 3.0));
    
    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(result.value(0), 1.0);
    EXPECT_DOUBLE_EQ(result.value(1), 2.0);
    EXPECT_DOUBLE_EQ(result.value(2), 3.0);
}



// Run all tests
// main() is provided by gtest_main
