#include <gtest/gtest.h>
#include "astdyn/astrometry/OccultationLogic.hpp"
#include <cmath>

using namespace astdyn::astrometry;

/**
 * @brief Test for basic occultation geometry logic
 */
TEST(OccultationTest, BasicLogic) {
    // Example: Asteroid moving across a star
    // Star at RA=0, Dec=0
    auto star_ra = RightAscension(Angle::from_deg(0.0));
    auto star_dec = Declination(Angle::from_deg(0.0));

    // Asteroid approaching from the West, offset slightly North
    // At t=0, it's at RA = -0.01 deg, Dec = 0.001 deg
    auto ast_ra = RightAscension(Angle::from_deg(-0.01));
    auto ast_dec = Declination(Angle::from_deg(0.001));
    
    // Distance = 2 AU
    double dist = 2.0 * 149597870.7 * 1000.0;
    
    // Velocity: moves East at 0.01 deg/hour
    double dra_dt = (0.01 / 3600.0) * (M_PI / 180.0);
    double ddec_dt = 0.0;
    
    auto params = OccultationLogic::compute_parameters(
        star_ra, star_dec,
        ast_ra, ast_dec,
        dist,
        dra_dt, ddec_dt, 0.0
    );

    // Expected:
    // Impact parameter should be dist * tan(0.001 deg)
    double expected_b_km = dist * std::tan(0.001 * M_PI / 180.0) / 1000.0;
    
    EXPECT_NEAR(params.impact_parameter_km, expected_b_km, 1e-3);
    EXPECT_GT(params.shadow_velocity_kms, 0.0);
    EXPECT_NEAR(params.position_angle_deg, 90.0, 1e-1); // Moving strictly East (PA=90)
}
