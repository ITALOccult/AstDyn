/**
 * @file OccultationLogic.hpp
 * @brief Logic for computing physical parameters of a stellar occultation.
 */

#ifndef ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP
#define ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP

#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/catalog/CatalogTypes.hpp"
#include <Eigen/Dense>
#include <vector>
#include <string>

namespace astdyn { class AstDynEngine; }

namespace astdyn::astrometry {

/**
 * @brief Physical parameters of a stellar occultation.
 */
struct OccultationParameters {
    physics::Distance xi_ca;              ///< Signed xi coordinate at TCA
    physics::Distance eta_ca;             ///< Signed eta coordinate at TCA
    physics::Distance impact_parameter;    ///< Minimum distance of shadow center from geocenter
    physics::Velocity shadow_velocity;    ///< Velocity magnitude on fundamental plane
    physics::Velocity dxi_dt;             ///< Velocity component along Xi (East)
    physics::Velocity deta_dt;            ///< Velocity component along Eta (North)
    Angle position_angle;                 ///< Direction of the shadow track (North to East)
    
    // Relative velocities (represented as angular displacement per second)
    Angle d_ra_cos_dec_per_sec;
    Angle d_dec_per_sec;
    
    time::TimeDuration closest_approach_time_offset; ///< Time offset from input center to TCA
    time::TimeDuration time_uncertainty;           ///< 1-sigma uncertainty in time
    physics::Distance cross_track_uncertainty;      ///< 1-sigma uncertainty in cross-track position
    
    // Metadata
    std::string star_id;
    double star_mag;
    // Geometry & Timing
    time::EpochTDB t_ca;  ///< Absolute time of closest approach
    Angle center_lon;     ///< Longitude of the sub-asteroid point at TCA (deg)
    Angle center_lat;     ///< Latitude of the sub-asteroid point at TCA (deg)
    time::TimeDuration max_duration; ///< Maximum occultation duration (sec)
    double total_apparent_rate;    ///< Total angular velocity (arcsec/hr)
    
    // Observability Filters
    bool is_daylight;     ///< True if TCA is during daylight at central point
    double sun_altitude;  ///< Elevation of Sun at central point (deg)
    double moon_altitude; ///< Elevation of Moon at central point (deg)
    double moon_dist;     ///< Angular distance from Moon to Star (deg)
    double mag_drop;      ///< Expected magnitude drop (m_star - m_asteroid)
};

/**
 * @brief Configuration for filtering occultation events.
 */
struct OccultationConfig {
    double min_sun_altitude = -12.0;       ///< Max Sun elevation (civil/nautical twilight)
    double min_object_altitude = 10.0;     ///< Min Object elevation at center
    double min_moon_dist = 5.0;            ///< Min angular distance from Moon (deg)
    double min_mag_drop = 0.05;            ///< Min magnitude drop for detectability
    double max_mag_star = 16.0;            ///< Max star magnitude
    bool filter_daylight = true;           ///< Skip events in daylight
    
    // Refinement settings
    bool use_proper_motion = true;         ///< Apply star proper motion to TCA
    bool use_parallax = true;             ///< Apply star parallax
};

/**
 * @brief Mode for asteroid state calculation.
 */
enum class OccultationRefinementMode {
    ChebyshevDaily,     ///< Use precomputed daily Chebyshev polynomials (Default)
    OrbitFitting,       ///< Download observations and fit orbit
    PropagationOnly     ///< Propagate original elements
};

/**
 * @brief Result for a single body in a multi-body system.
 */
struct BodyOccultation {
    std::string name;
    OccultationParameters params;
    physics::Distance diameter;
};

/**
 * @brief Result of a discovery/refinement run for a multi-body system.
 */
struct OccultationSystemCandidate {
    catalog::Star star;
    std::vector<BodyOccultation> bodies;
};

/**
 * @brief Result of a discovery/refinement run.
 */
struct OccultationCandidate {
    std::string asteroid_id;
    catalog::Star star;
    OccultationParameters params;
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
    /**
     * @brief Computes physical parameters of the occultation.
     * ...
     */
    static OccultationParameters compute_parameters(
        const RightAscension& star_ra, const Declination& star_dec,
        const RightAscension& ast_ra, const Declination& ast_dec,
        const physics::Distance& ast_dist,
        const Angle& ast_dra_dt, const Angle& ast_ddec_dt,
        const physics::Velocity& ast_ddist_dt,
        const time::EpochTDB& t_ca);

    /**
     * @brief Main entry point for searching occultations.
     * 
     * @param asteroid_id       Designation or name (for downloads)
     * @param initial_elements  Seed elements for search
     * @param start             Search window start
     * @param end               Search window end
     * @param max_mag           Magnitude limit for stars
     * @param engine            AstDyn engine (holds config, propagator, etc)
     * @param mode              Refinement mode (Daily Chebyshev, Fit, or Prop)
     * @return List of verified occultations
     */
    static std::vector<OccultationCandidate> find_occultations(
        const std::string& asteroid_id,
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
        time::EpochTDB start,
        time::EpochTDB end,
        const OccultationConfig& config,
        AstDynEngine& engine,
        OccultationRefinementMode mode = OccultationRefinementMode::ChebyshevDaily);

    /**
     * @brief Search for occultations of a whole system (e.g. asteroid + satellites).
     * 
     * @param body_ids Names/IDs of the bodies to search for.
     * @param bsp_path Path to the SPK/BSP file containing the system ephemeris.
     * ...
     */
    static std::vector<OccultationSystemCandidate> find_system_occultations(
        const std::vector<std::string>& body_ids,
        const std::string& bsp_path,
        time::EpochTDB start,
        time::EpochTDB end,
        const OccultationConfig& config,
        AstDynEngine& engine);

    /**
     * @brief Search for occultations for multiple asteroids using pre-calculated polynomials.
     * 
     * @param asteroid_ids List of asteroid designations to search.
     * @param manager      Pre-filled manager containing AsteroidChebyshevEphemeris for each ID.
     * @param start        Search window start.
     * @param end          Search window end.
     * @param max_mag      Stellar magnitude limit.
     * @param engine       AstDyn engine.
     * @return List of verified occultations.
     */
    static std::vector<OccultationCandidate> find_multi_asteroid_occultations(
        const std::vector<std::string>& asteroid_ids,
        class ChebyshevEphemerisManager& manager,
        time::EpochTDB start,
        time::EpochTDB end,
        const OccultationConfig& config,
        AstDynEngine& engine);

    /**
     * @brief Computes 1-sigma uncertainty parameters for an event.
     * 
     * @param params        The parameters to update.
     * @param covariance_t0 6x6 covariance at the initial epoch t0.
     * @param initial_state Nominal state at t0.
     * @param engine        AstDyn engine for STM propagation.
     */
    static void apply_uncertainty(
        OccultationParameters& params,
        const catalog::Star& star,
        const Eigen::Matrix<double, 6, 6>& covariance_t0,
        const physics::CartesianStateTyped<core::GCRF>& initial_state,
        AstDynEngine& engine);

    /**
     * @brief Convenience overload using SkyCoord.
     */
    template <typename Frame>
    static OccultationParameters compute_from_sky(
        const SkyCoord<Frame>& star,
        const SkyCoord<Frame>& asteroid,
        const physics::Distance& distance,
        const Angle& dra_dt,
        const Angle& ddec_dt,
        const time::EpochTDB& t_ca,
        const physics::Velocity& ddist_dt = physics::Velocity::zero()) 
    {
        return compute_parameters(
            star.ra(), star.dec(),
            asteroid.ra(), asteroid.dec(),
            distance,
            dra_dt, ddec_dt, ddist_dt, t_ca
        );
    }
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP
