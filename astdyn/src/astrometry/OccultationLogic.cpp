#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/io/MPCClient.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/io/SPKReader.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/astrometry/ChebyshevEphemerisManager.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/AstDynEngine.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace astdyn::astrometry {

OccultationParameters OccultationLogic::compute_parameters(
    const RightAscension& star_ra, const Declination& star_dec,
    const RightAscension& ast_ra, const Declination& ast_dec,
    const physics::Distance& ast_dist,
    const Angle& ast_dra_dt, const Angle& ast_ddec_dt,
    const physics::Velocity& ast_ddist_dt,
    const time::EpochTDB& t_ca,
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem)
{
    using namespace constants;

    double as = star_ra.to_rad();
    double ds = star_dec.to_rad();
    double a = ast_ra.to_rad();
    double d = ast_dec.to_rad();

    OccultationParameters params;
    params.t_ca = t_ca;

    Eigen::Vector3d k(std::cos(ds) * std::cos(as), std::cos(ds) * std::sin(as), std::sin(ds));
    Eigen::Vector3d r_ast = ast_dist.to_m() * Eigen::Vector3d(std::cos(d) * std::cos(a), std::cos(d) * std::sin(a), std::sin(d));
    
    // East/North basis in FP
    Eigen::Vector3d east(-std::sin(as), std::cos(as), 0.0);
    Eigen::Vector3d north = k.cross(east);

    // Shadow Vector in Fundamental Plane
    Eigen::Vector3d r_sub_gcrf = r_ast - (r_ast.dot(k)) * k;
    params.xi_ca = physics::Distance::from_m(r_sub_gcrf.dot(east));
    params.eta_ca = physics::Distance::from_m(r_sub_gcrf.dot(north));
    params.impact_parameter = physics::Distance::from_m(r_sub_gcrf.norm());

    // Shadow Velocity
    double dra = ast_dra_dt.to_rad();
    double ddec = ast_ddec_dt.to_rad();
    double r_val = ast_dist.to_m();
    
    // Dot dr/dt with basis
    params.dxi_dt = physics::Velocity::from_ms(r_val * dra * std::cos(ds));
    params.deta_dt = physics::Velocity::from_ms(r_val * ddec);
    params.shadow_velocity = physics::Velocity::from_ms(std::sqrt(params.dxi_dt.to_ms()*params.dxi_dt.to_ms() + params.deta_dt.to_ms()*params.deta_dt.to_ms()));
    params.position_angle = Angle::from_rad(std::atan2(params.dxi_dt.to_ms(), params.deta_dt.to_ms())).wrap_0_2pi();
    params.relative_velocity_mag = params.shadow_velocity;

    if (ephem) {
        auto sun_state = ephem->getState(ephemeris::CelestialBody::SUN, t_ca);
        Eigen::Vector3d n_sun = sun_state.position.to_eigen_si().normalized();
        params.sun_altitude = std::asin(std::clamp((r_ast.normalized()).dot(n_sun), -1.0, 1.0)) * 180.0 / M_PI;
        params.is_daylight = params.sun_altitude > -0.83;

        auto moon_ssb = ephem->getState(ephemeris::CelestialBody::MOON, t_ca);
        auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t_ca);
        Eigen::Vector3d n_moon = (moon_ssb.position.to_eigen_si() - earth_ssb.position.to_eigen_si()).normalized();
        params.moon_altitude = std::asin(std::clamp((r_ast.normalized()).dot(n_moon), -1.0, 1.0)) * 180.0 / M_PI;
        params.moon_dist = std::acos(std::clamp(k.dot(n_moon), -1.0, 1.0)) * 180.0 / M_PI;
    }

    return params;
}

std::vector<OccultationCandidate> OccultationLogic::find_occultations(
    const std::string& asteroid_id,
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
    time::EpochTDB start,
    time::EpochTDB end,
    const OccultationConfig& config,
    AstDynEngine& engine,
    OccultationRefinementMode mode)
{
    using namespace catalog;
    std::vector<OccultationCandidate> results;
    auto& catalog_inst = GaiaDR3Catalog::instance();

    // 1. Initial State Initialization
    physics::KeplerianStateTyped<core::ECLIPJ2000> working_elements = initial_elements;

    double t_start_jd = start.jd();
    double t_end_jd = end.jd();

    std::vector<double> days;
    for (double day_jd = std::floor(t_start_jd - 0.5) + 0.5; day_jd <= t_end_jd + 0.5; day_jd += 1.0) {
        if (day_jd >= t_start_jd - 0.5 && day_jd <= t_end_jd + 0.5) {
            days.push_back(day_jd);
        }
    }

    #pragma omp parallel
    {
        std::vector<OccultationCandidate> thread_results;
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < days.size(); ++i) {
            double day_jd = days[i];
            time::EpochTDB midnight = time::EpochTDB::from_jd(day_jd);

            // Per-thread engine state management if needed. 
            // In AstDyn 3.0, AstDynEngine computation methods are generally thread-safe 
            // as internal state (propagator) is accessed via shared_ptr.
            
            try {
                auto ephemeris = engine.compute_ephemeris(midnight - time::TimeDuration::from_hours(12), 
                                                        midnight + time::TimeDuration::from_hours(12), 1.0/24.0);
                
                // (Logic for star proximity search would go here)
            } catch (...) { }
        }
        
        #pragma omp critical
        results.insert(results.end(), thread_results.begin(), thread_results.end());
    }

    return results;
}

std::vector<OccultationSystemCandidate> OccultationLogic::find_system_occultations(
    const std::vector<std::string>& body_ids,
    const std::string& bsp_path,
    time::EpochTDB start,
    time::EpochTDB end,
    const OccultationConfig& config,
    AstDynEngine& engine)
{
    std::vector<OccultationSystemCandidate> results;
    // Implementation would use engine.getEphemeris()
    return results;
}

std::vector<OccultationCandidate> OccultationLogic::find_multi_asteroid_occultations(
    const std::vector<std::string>& asteroid_ids,
    ChebyshevEphemerisManager& manager,
    time::EpochTDB start,
    time::EpochTDB end,
    const OccultationConfig& config,
    AstDynEngine& engine)
{
    std::vector<OccultationCandidate> results;
    // Implementation would use engine.getEphemeris()
    return results;
}

void OccultationLogic::apply_uncertainty(
    OccultationParameters& params,
    const catalog::Star& star,
    const Eigen::Matrix<double, 6, 6>& covariance_t0,
    const physics::CartesianStateTyped<core::GCRF>& initial_state,
    AstDynEngine& engine)
{
    using namespace orbit_determination;
    StateTransitionMatrix<core::GCRF> stm_engine(engine.propagator());
    auto stm_res = stm_engine.compute(initial_state, params.t_ca);
    Eigen::Matrix<double, 6, 6> cov_tca = stm_res.map_covariance(covariance_t0);
    
    auto state_tca = stm_res.final_state;
    auto earth_tca = engine.getEphemeris()->getState(ephemeris::CelestialBody::EARTH, params.t_ca);
    Eigen::Vector3d rho_geo = state_tca.position.to_eigen_si() - earth_tca.position.to_eigen_si();

    params.cross_track_uncertainty = physics::Distance::from_km(
        AstrometryReducer::compute_cross_track_uncertainty(
            cov_tca, 
            rho_geo,
            state_tca.velocity.to_eigen_si(),
            star.ra.to_rad(), 
            star.dec.to_rad()
        )
    );
}

} // namespace astdyn::astrometry
