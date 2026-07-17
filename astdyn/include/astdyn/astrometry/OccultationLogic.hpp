#ifndef ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP
#define ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP

#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/catalog/CatalogTypes.hpp"
#include "astdyn/astrometry/ClosestApproachFinder.hpp"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <memory>

namespace astdyn { 
    class AstDynEngine; 
    namespace ephemeris { class PlanetaryEphemeris; }
}

namespace astdyn::astrometry {

/**
 * @brief Apparent magnitude in the IAU HG system.
 *
 *     V = H + 5 log10(r * delta) - 2.5 log10[(1-G) Phi1 + G Phi2]
 *     Phi_i = exp(-A_i tan(alpha/2)^B_i)
 *
 * Reported as -5.00 by Occult4 when unavailable -- NOT as zero, which would
 * claim an object brighter than Vega and wreck any magnitude drop computed from
 * it. The same convention is used here.
 *
 * @param h_mag  Absolute magnitude H.
 * @param g      Slope parameter G (0.15 by convention when unmeasured).
 * @param r_au   Heliocentric distance.
 * @param d_au   Geocentric distance.
 * @param alpha  Sun-object-Earth phase angle.
 */
[[nodiscard]] double hg_magnitude(double h_mag, double g, double r_au, double d_au,
                                  Angle alpha);

struct OccultationParameters {
    physics::Distance xi_ca;              
    physics::Distance eta_ca;             
    physics::Distance impact_parameter;    
    physics::Velocity shadow_velocity;    
    physics::Velocity dxi_dt;             
    physics::Velocity deta_dt;            
    Angle position_angle;                 
    physics::Velocity relative_velocity_mag;
    time::EpochTDB t_ca;                  
    Angle sun_altitude;                  
    Angle moon_altitude;                 
    Angle moon_dist;                     
    Angle center_lat;                  
    Angle center_lon;                  
    /// Sub-star point: where the star is at the zenith. This is the reference
    /// point of the fundamental plane, and what Occult4 stores in <Earth>;
    /// the shadow centre above is displaced from it by the impact parameter.
    /// GEOCENTRIC latitude: this is exactly the star's apparent declination,
    /// which is the convention Occult4 uses in <Earth>.
    Angle substar_lat;                 
    Angle substar_lon;                 
    /// Sub-solar point (geocentric). Occult4 carries it in <Earth> fields 3-4
    /// so that the day/night terminator can be drawn.
    Angle subsolar_lat;                
    Angle subsolar_lon;                
    double star_mag;                      
    double mag_drop;                      
    bool is_daylight;                     
    double total_apparent_rate; // arcsec/hr
    /// Geocentric distance of the occulting object at t_ca.
    physics::Distance geocentric_distance;
    /// Heliocentric distance and Sun-object-Earth phase angle at t_ca. Both are
    /// needed for the HG apparent magnitude, and neither is recoverable from the
    /// geocentric distance alone.
    physics::Distance heliocentric_distance;
    Angle phase_angle;
    /// Apparent hourly rates of the object. Kept separately from
    /// total_apparent_rate because the occelmnt format reports them apart, and
    /// dRA there is in SECONDS OF TIME per hour, not arcsec.
    double d_ra_arcsec_hr = 0.0;
    double d_dec_arcsec_hr = 0.0;
    time::TimeDuration max_duration;
    std::string star_id;
    double moon_phase;                   // 0.0 (new) to 1.0 (full)
    
    physics::Distance cross_track_uncertainty;

    // ---- SCOPE: plane-of-sky uncertainty at the event ----
    /// A priori nonlinearity index N (Eq. 18): RMS of the second-order term over
    /// RMS of the first. Available WITHOUT a Monte Carlo. Around 1e-7 for a
    /// well-determined main-belt asteroid a few weeks out, i.e. the linear
    /// theory is exact and N says so; large for a short-arc NEO, and N says that
    /// too. Zero means it was never computed.
    double nonlinearity_index = 0.0;
    /// 1-sigma error ellipse in the plane of sky, semi-axes and position angle
    /// of the major axis. Includes the stellar contribution (Eq. 17).
    Angle err_major, err_minor, err_pa;
    /// Second-order bias of the shadow position (Eq. 14), cross-track. The
    /// linear theory predicts exactly zero, so this is a pure second-order term.
    physics::Distance bias_cross_track;
    Eigen::Vector3d shadow_velocity_vector;
    physics::Velocity shadow_velocity_fundamental_plane;
};

struct BodyOccultation {
    std::string name;
    OccultationParameters params;
    physics::Distance diameter;
};

struct OccultationCandidate {
    std::string asteroid_id;
    catalog::Star star;
    OccultationParameters params;
};

struct OccultationSystemCandidate {
    catalog::Star star;
    std::vector<BodyOccultation> bodies;
};

enum class OccultationRefinementMode {
    ChebyshevDaily, 
    Propagate,      
    FitObservations 
};

class OccultationLogic {
public:
    static OccultationParameters compute_parameters(
        const RightAscension& star_ra, const Declination& star_dec,
        const RightAscension& ast_ra, const Declination& ast_dec,
        const physics::Distance& ast_dist,
        const Angle& ast_dra_dt, const Angle& ast_ddec_dt,
        const physics::Velocity& ast_ddist_dt,
        const time::EpochTDB& t_ca,
        std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem);

    static std::vector<OccultationCandidate> find_occultations(
        const std::string& asteroid_id,
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
        time::EpochTDB start,
        time::EpochTDB end,
        const struct OccultationConfig& config,
        AstDynEngine& engine,
        OccultationRefinementMode mode = OccultationRefinementMode::ChebyshevDaily);

    static std::vector<OccultationSystemCandidate> find_system_occultations(
        const std::vector<std::string>& body_ids,
        const std::string& bsp_path,
        time::EpochTDB start,
        time::EpochTDB end,
        const struct OccultationConfig& config,
        AstDynEngine& engine);

    static std::vector<OccultationCandidate> find_multi_asteroid_occultations(
        const std::vector<std::string>& asteroid_ids,
        class ChebyshevEphemerisManager& manager,
        time::EpochTDB start,
        time::EpochTDB end,
        const struct OccultationConfig& config,
        AstDynEngine& engine);

    /**
     * @brief SCOPE: propagate the orbit covariance to the event and project it.
     *
     * @param covariance_t0  Orbit covariance at the epoch of @p initial_state,
     *        in AU and AU/day, ECLIPJ2000. From AstDyS via
     *        EquinoctialElements::jacobian_to_cartesian().
     * @param initial_state  State at t0. The frame tag is ECLIPJ2000 and it is
     *        meant literally: the previous signature said GCRF while being fed a
     *        heliocentric ecliptic state through cast_frame, which is exactly the
     *        kind of silent relabelling the frame tags exist to prevent.
     */
    static void apply_uncertainty(
        OccultationParameters& params,
        const catalog::Star& star,
        const Eigen::Matrix<double, 6, 6>& covariance_t0,
        const physics::CartesianStateTyped<core::ECLIPJ2000>& initial_state,
        AstDynEngine& engine);

private:
    static OccultationParameters compute_fundamental_plane_geometry(
        const RightAscension& star_ra, const Declination& star_dec,
        const RightAscension& ast_ra, const Declination& ast_dec,
        const physics::Distance& ast_dist);

    static void compute_shadow_velocity(
        OccultationParameters& params,
        const Declination& star_dec,
        const physics::Distance& ast_dist,
        const Angle& ast_dra_dt,
        const Angle& ast_ddec_dt);

    /// Geodetic lat/lon of the shadow centre at t_ca (uses the impact
    /// parameter xi/eta; NOT the sub-asteroid point).
    static void compute_shadow_centre(
        OccultationParameters& params,
        const time::EpochTDB& t_ca,
        const RightAscension& star_ra, const Declination& star_dec);

    static void compute_sky_conditions(
        OccultationParameters& params,
        const time::EpochTDB& t_ca,
        const RightAscension& ast_ra, const Declination& ast_dec,
        const physics::Distance& ast_distance,
        const RightAscension& star_ra, const Declination& star_dec,
        std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem);

    static void process_day_window(
        std::vector<OccultationCandidate>& results,
        double day_jd,
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
        const OccultationConfig& config,
        AstDynEngine& engine,
        double t_start_jd, double t_end_jd);

    /// @param diameter_km  Occulting object's diameter, needed for the maximum
    ///        duration. Zero leaves max_duration at zero, which the duration
    ///        filter reads: pass it whenever it is known.
    static void evaluate_candidate(
        std::vector<OccultationCandidate>& results,
        const catalog::ChebyshevSegment& segment,
        const catalog::Star& star,
        const OccultationConfig& config,
        AstDynEngine& engine,
        double t_start_jd, double t_end_jd,
        double diameter_km = 0.0);

    static void evaluate_candidate(
        const std::string& id,
        const catalog::Star& star,
        const catalog::ChebyshevSegment& segment,
        const std::vector<ClosestApproachResult>& candidates,
        std::vector<OccultationCandidate>& results,
        const OccultationConfig& config,
        AstDynEngine& engine);
};

struct OccultationConfig {
    double max_mag_star = 14.0;
    double min_mag_drop = 0.05;
    double min_sun_altitude = -6.0;
    double min_moon_dist = 5.0;
    double min_object_altitude = 10.0;
    bool use_proper_motion = true;
    bool use_parallax = true;
    bool filter_daylight = true;
    physics::Distance max_shadow_distance = physics::Distance::from_km(10000.0);
    bool compute_uncertainty = true;

    // Advanced Scientific Filters
    double min_duration_s = 0.0;
    double min_asteroid_diameter_km = 0.0;

    // Site-Specific Proximity Filters
    double obs_lat = 0.0;
    double obs_lon = 0.0;
    double max_obs_dist_km = 0.0; // 0.0 = disabled

    // Scientific Quality Filters
    double max_gaia_ruwe = 99.0;         // Rejects targets with RUWE > limit (e.g. 1.4)
    double max_moon_phase = 1.0;         // 0.0 to 1.0 (Full)
};

} // namespace astdyn::astrometry

#endif
