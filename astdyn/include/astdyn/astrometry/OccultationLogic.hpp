#ifndef ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP
#define ASTDYN_ASTROMETRY_OCCULTATION_LOGIC_HPP

#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/catalog/CatalogTypes.hpp"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <memory>

namespace astdyn { 
    class AstDynEngine; 
    namespace ephemeris { class PlanetaryEphemeris; }
}

namespace astdyn::astrometry {

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
    double star_mag;                      
    double mag_drop;                      
    bool is_daylight;                     
    double total_apparent_rate; // arcsec/hr
    time::TimeDuration max_duration;
    std::string star_id;
    
    physics::Distance cross_track_uncertainty;
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

    static void apply_uncertainty(
        OccultationParameters& params,
        const catalog::Star& star,
        const Eigen::Matrix<double, 6, 6>& covariance_t0,
        const physics::CartesianStateTyped<core::GCRF>& initial_state,
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

    static void compute_sky_conditions(
        OccultationParameters& params,
        const time::EpochTDB& t_ca,
        const RightAscension& ast_ra, const Declination& ast_dec,
        const physics::Distance& ast_distance,
        const RightAscension& star_ra, const Declination& star_dec,
        std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem);
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
};

} // namespace astdyn::astrometry

#endif
