/**
 * @file export_kml_bundle.cpp
 * @brief Generates 4 KML files for Vesta and Nireus occultations.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/time/epoch.hpp"
#include <iostream>

using namespace astdyn;
using namespace astdyn::astrometry;

int main() {
    // 1. Vesta (March 22, 2026)
    OccultationParameters vesta_params;
    vesta_params.xi_ca = physics::Distance::from_km(-1242.0);
    vesta_params.eta_ca = physics::Distance::from_km(3699.0);
    vesta_params.dxi_dt = physics::Velocity::from_ms(17109.0);
    vesta_params.deta_dt = physics::Velocity::from_ms(5818.0);
    vesta_params.shadow_velocity = physics::Velocity::from_ms(std::sqrt(17109*17109 + 5818*5818));
    vesta_params.position_angle = Angle::from_rad(std::atan2(17110, 5819)).wrap_0_2pi();
    vesta_params.closest_approach_time_offset = time::TimeDuration::zero();
    vesta_params.cross_track_uncertainty = physics::Distance::from_km(15.0); 

    RightAscension vesta_star_ra = RightAscension(Angle::from_deg(22.56756 * 15.0));
    Declination vesta_star_dec = Declination(Angle::from_deg(-12.6708));
    time::EpochUTC vesta_t0 = time::EpochUTC::from_mjd(time::calendar_to_mjd(2026, 3, 22, 12.5826 / 24.0));

    auto vesta_path = OccultationMapper::compute_path(vesta_params, vesta_star_ra, vesta_star_dec, physics::Distance::from_km(529.0), vesta_t0, time::TimeDuration::from_seconds(3600.0));
    OccultationMapper::export_kml(vesta_path, "vesta_astdyn_refined.kml");

    OccultationParameters v_off;
    v_off.xi_ca = physics::Distance::from_km(-0.19653 * constants::R_EARTH);
    v_off.eta_ca = physics::Distance::from_km(0.58091 * constants::R_EARTH);
    v_off.dxi_dt = physics::Velocity::from_km_s(24.82 * constants::R_EARTH / 3600.0);
    v_off.deta_dt = physics::Velocity::from_km_s(8.44 * constants::R_EARTH / 3600.0);
    v_off.position_angle = Angle::from_rad(std::atan2(24.82, 8.44));
    v_off.closest_approach_time_offset = time::TimeDuration::zero();
    auto vesta_path_off = OccultationMapper::compute_path(v_off, vesta_star_ra, vesta_star_dec, physics::Distance::from_km(529.0), vesta_t0, time::TimeDuration::from_seconds(3600.0));
    OccultationMapper::export_kml(vesta_path_off, "vesta_official_refined.kml");

    // 2. Nireus (March 12, 2026)
    OccultationParameters nireus_params;
    nireus_params.xi_ca = physics::Distance::from_km(96.0);
    nireus_params.eta_ca = physics::Distance::from_km(4596.0);
    nireus_params.dxi_dt = physics::Velocity::from_ms(7310.0 * std::sin(271.0 * constants::DEG_TO_RAD));
    nireus_params.deta_dt = physics::Velocity::from_ms(7310.0 * std::cos(271.0 * constants::DEG_TO_RAD));
    nireus_params.shadow_velocity = physics::Velocity::from_km_s(7.31);
    nireus_params.position_angle = Angle::from_deg(271.0);
    nireus_params.closest_approach_time_offset = time::TimeDuration::zero();
    nireus_params.cross_track_uncertainty = physics::Distance::from_km(45.0);

    RightAscension nireus_star_ra = RightAscension(Angle::from_deg(204.7357));
    Declination nireus_star_dec = Declination(Angle::from_deg(-8.8142));
    time::EpochUTC nireus_t0 = time::EpochUTC::from_mjd(61112.115671);

    auto nireus_path = OccultationMapper::compute_path(nireus_params, nireus_star_ra, nireus_star_dec, physics::Distance::from_km(100.0), nireus_t0, time::TimeDuration::from_seconds(1200.0));
    OccultationMapper::export_kml(nireus_path, "nireus_astdyn_refined.kml");

    OccultationParameters n_off;
    n_off.xi_ca = physics::Distance::from_m(96000.0);
    n_off.eta_ca = physics::Distance::from_m(4596000.0);
    n_off.dxi_dt = physics::Velocity::from_ms(7310.0 * std::sin(271.0 * constants::DEG_TO_RAD));
    n_off.deta_dt = physics::Velocity::from_ms(7310.0 * std::cos(271.0 * constants::DEG_TO_RAD));
    n_off.position_angle = Angle::from_deg(271.0);
    n_off.closest_approach_time_offset = time::TimeDuration::zero();
    auto nireus_path_off = OccultationMapper::compute_path(n_off, nireus_star_ra, nireus_star_dec, physics::Distance::from_km(100.0), nireus_t0, time::TimeDuration::from_seconds(1200.0));
    OccultationMapper::export_kml(nireus_path_off, "nireus_official_refined.kml");

    std::cout << "Generated 4 KML files: vesta_astdyn_refined.kml, vesta_official_refined.kml, nireus_astdyn_refined.kml, nireus_official_refined.kml" << std::endl;
    return 0;
}
