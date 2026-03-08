#include <iostream>
#include <cassert>
#include "astdyn/astrometry/OccultationLogic.hpp"

using namespace astdyn::astrometry;

int main() {
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

    std::cout << "Occultation Parameters:" << std::endl;
    std::cout << "  Impact Parameter: " << params.impact_parameter_km << " km" << std::endl;
    std::cout << "  Shadow Velocity:  " << params.shadow_velocity_kms << " km/s" << std::endl;
    std::cout << "  Position Angle:   " << params.position_angle_deg << " deg" << std::endl;
    std::cout << "  Time to TCA:      " << params.closest_approach_time_offset_sec << " s" << std::endl;
    std::cout << "  RA Velocity (mas/s): " << params.d_ra_cos_dec_mas_sec << std::endl;

    // Expected:
    // Impact parameter should be dist * tan(0.001 deg)
    double expected_b_km = dist * std::tan(0.001 * M_PI / 180.0) / 1000.0;
    std::cout << "  Expected Impact:  " << expected_b_km << " km" << std::endl;

    return 0;
}
