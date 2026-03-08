/**
 * @file OccultationLogic.hpp
 * @brief Logic for computing physical parameters of a stellar occultation.
 */

#ifndef ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP
#define ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP

#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/core/physics_state.hpp"
#include <Eigen/Dense>

namespace astdyn::astrometry {

/**
 * @brief Physical parameters of a stellar occultation.
 */
struct OccultationParameters {
    double impact_parameter_km;    ///< Minimum distance of shadow center from geocenter [km]
    double shadow_velocity_kms;    ///< Velocity of the shadow on the fundamental plane [km/s]
    double position_angle_deg;     ///< Direction of the shadow track (North to East) [deg]
    double d_ra_cos_dec_mas_sec;   ///< Relative velocity in RA*cos(dec) [mas/s]
    double d_dec_mas_sec;          ///< Relative velocity in Declination [mas/s]
    double closest_approach_time_offset_sec; ///< Time offset from input center to TCA [s]
};

/**
 * @brief Computes occultation parameters from apparent sky coordinates.
 */
class OccultationLogic {
public:
    /**
     * @brief Computes physical parameters of the occultation.
     * 
     * @param star_ra Star Right Ascension
     * @param star_dec Star Declination
     * @param ast_ra Asteroid Apparent Right Ascension
     * @param ast_dec Asteroid Apparent Declination
     * @param ast_dist Asteroid Geocentric Distance [m]
     * @param ast_dra_dt Angular RA velocity [rad/s]
     * @param ast_ddec_dt Angular Dec velocity [rad/s]
     * @param ast_ddist_dt Radial velocity [m/s]
     * @return Physical parameters
     */
    static OccultationParameters compute_parameters(
        const RightAscension& star_ra, const Declination& star_dec,
        const RightAscension& ast_ra, const Declination& ast_dec,
        double ast_dist_m,
        double ast_dra_dt_rad_s, double ast_ddec_dt_rad_s,
        double ast_ddist_dt_m_s);

    /**
     * @brief Convenience overload using SkyCoord.
     */
    template <typename Frame>
    static OccultationParameters compute_from_sky(
        const SkyCoord<Frame>& star,
        const SkyCoord<Frame>& asteroid,
        double distance_m,
        double dra_dt_rad_s,
        double ddec_dt_rad_s) 
    {
        return compute_parameters(
            star.ra(), star.dec(),
            asteroid.ra(), asteroid.dec(),
            distance_m,
            dra_dt_rad_s, ddec_dt_rad_s, 0.0
        );
    }
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP
