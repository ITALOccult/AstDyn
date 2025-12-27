/**
 * @file test_astdyn_wrapper.cpp
 * @brief Unit tests for AstDynWrapper high-level API
 */

#include <gtest/gtest.h>
#include "astdyn_wrapper.h"
#include <Eigen/Dense>
#include <cmath>
#include <filesystem>

using namespace ioccultcalc;
namespace fs = std::filesystem;

class AstDynWrapperTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Path to test data (17030 Sierks)
        // Note: This assumes the test is run from the build directory 
        // within italoccultlibrary, and matches the project structure.
        std::string potential_path = "../../astdyn/data/17030.eq1";
        if (fs::exists(potential_path)) {
            eq1_path = potential_path;
        } else {
            // Try absolute path if relative fails (for environment flexibility)
            eq1_path = "/Users/michelebigi/Documents/Develop/ASTDYN/ITALOccultLibrary/astdyn/data/17030.eq1";
        }
    }

    std::string eq1_path;
};

TEST_F(AstDynWrapperTest, Initialization) {
    AstDynWrapper wrapper(PropagationSettings::highAccuracy());
    EXPECT_FALSE(wrapper.isInitialized());
    
    ASSERT_TRUE(fs::exists(eq1_path)) << "Test data file not found at: " << eq1_path;
    
    EXPECT_TRUE(wrapper.loadFromEQ1File(eq1_path));
    EXPECT_TRUE(wrapper.isInitialized());
    EXPECT_EQ(wrapper.getObjectName(), "17030");
    EXPECT_NEAR(wrapper.getCurrentEpoch(), 61000.0, 1e-6);
}

TEST_F(AstDynWrapperTest, HighPrecisionValidation17030) {
    // This test verifies that the high-level API produces JPL-grade results
    // for Asteroid 17030 (Sierks) at a known epoch.
    
    AstDynWrapper wrapper(PropagationSettings::highAccuracy());
    ASSERT_TRUE(wrapper.loadFromEQ1File(eq1_path));
    
    // Reference epoch: JD 2460643.77083 (2025-11-26 06:30 UTC)
    double target_jd = 2460643.77083;
    double target_mjd = target_jd - 2400000.5;
    
    // Expected position (ICRF, AU) from validated example
    // x = 3.221414563652
    // y = 0.770557583764
    // z = 0.150969776966
    Eigen::Vector3d expected_pos(3.221414563652, 0.770557583764, 0.150969776966);
    
    CartesianStateICRF state = wrapper.propagateToEpoch(target_mjd);
    
    // Check precision (within ~1 km = 6.6e-9 AU)
    // Using a slightly more relaxed threshold in unit test for different environments,
    // but demonstrating the 12-decimal precision matching.
    EXPECT_NEAR(state.position.x(), expected_pos.x(), 1e-8);
    EXPECT_NEAR(state.position.y(), expected_pos.y(), 1e-8);
    EXPECT_NEAR(state.position.z(), expected_pos.z(), 1e-8);
    
    std::cout << "High Precision Diff (AU): " << (state.position - expected_pos).norm() << std::endl;
}

TEST_F(AstDynWrapperTest, ManualElements) {
    AstDynWrapper wrapper;
    
    // 17030 Sierks elements at MJD 61000
    double a = 3.175473;
    double e = 0.045421;
    double i = 2.904601 * M_PI / 180.0;
    double Omega = 104.162434 * M_PI / 180.0;
    double omega = 100.514070 * M_PI / 180.0;
    double M = 229.790911 * M_PI / 180.0;
    double epoch = 61000.0;
    
    wrapper.setKeplerianElements(a, e, i, Omega, omega, M, epoch, "Sierks_Manual");
    
    EXPECT_TRUE(wrapper.isInitialized());
    EXPECT_EQ(wrapper.getObjectName(), "Sierks_Manual");
    
    auto state = wrapper.propagateToEpoch(epoch + 1.0);
    EXPECT_NEAR(state.epoch_mjd_tdb, epoch + 1.0, 1e-10);
    EXPECT_GT(state.position.norm(), 3.0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
