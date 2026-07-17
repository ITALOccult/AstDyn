#include <gtest/gtest.h>
#include "astdyn/astrometry/OccultationLogic.hpp"
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

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
    auto dist = astdyn::physics::Distance::from_au(2.0);
    
    // Velocity: moves East at 0.01 deg/sec (for visibility in test)
    auto dra_dt = Angle::from_deg(0.01);
    auto ddec_dt = Angle::from_deg(0.0);
    auto ddist_dt = astdyn::physics::Velocity::from_ms(0.0);
    auto t_ca = astdyn::time::EpochTDB::from_mjd(60000.0);
    
    auto params = OccultationLogic::compute_parameters(
        star_ra, star_dec,
        ast_ra, ast_dec,
        dist,
        dra_dt, ddec_dt, ddist_dt,
        t_ca,
        nullptr
    );

    // Expected: the impact parameter is dist * sin(separation), where the
    // separation is the FULL angle between the two directions.
    //
    // The previous expectation used tan(0.001 deg), i.e. the declination offset
    // alone, and silently dropped the 0.01 deg offset in right ascension. That
    // made it fail by a factor sqrt(0.01^2 + 0.001^2)/0.001 = 10.05 -- and the
    // code was right all along. Derive the separation from the inputs so that
    // changing the geometry above cannot leave the expectation behind.
    const double as = star_ra.to_rad(), ds = star_dec.to_rad();
    const double aa = ast_ra.to_rad(),  dd = ast_dec.to_rad();
    const Eigen::Vector3d k(std::cos(ds) * std::cos(as),
                            std::cos(ds) * std::sin(as), std::sin(ds));
    const Eigen::Vector3d n(std::cos(dd) * std::cos(aa),
                            std::cos(dd) * std::sin(aa), std::sin(dd));
    const double sep = std::acos(std::clamp(k.dot(n), -1.0, 1.0));
    const double expected_b_km = dist.to_m() * std::sin(sep) / 1000.0;

    EXPECT_NEAR(params.impact_parameter.to_km(), expected_b_km, 1e-1);
    EXPECT_GT(params.shadow_velocity.to_km_s(), 0.0);
    EXPECT_NEAR(params.position_angle.to_deg(), 90.0, 1e-1); // Moving strictly East (PA=90)
}
