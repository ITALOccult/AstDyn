/**
 * @file occultation_vesta_2026.cpp
 * @brief Benchmark for (4) Vesta stellar occultation on 2026-03-22.
 * Comparisons vs Occult4 Official Data using BSP ephemeris for Vesta.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/io/SPKReader.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::astrometry;

int main() {
    std::cout << "=== ASTDYN BENCHMARK: (4) VESTA OCCULTATION 2026-03-22 (BSP MODE) ===" << std::endl;

    // 1. Setup Ephemeris
    auto ephem_provider = std::make_shared<ephemeris::DE441Provider>("/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp");
    ephemeris::PlanetaryEphemeris::setProvider(ephem_provider);

    // 2. Setup Vesta BSP Reader (sb441-n16.bsp)
    io::SPKReader vesta_reader("/Users/michelebigi/.ioccultcalc/ephemerides/sb441-n16.bsp");
    const int VESTA_NAIF_ID = 2000004;

    // 3. Definizione Stella: TYC 5817-00509-1 (Mean J2000 from Occult4)
    double star_ra_h = 22.56756392;
    double star_dec_deg = -12.6708136;
    SkyCoord<core::GCRF> target_star = SkyCoord<core::GCRF>::from_deg(star_ra_h * 15.0, star_dec_deg);

    // 4. Time Reference (T0 = 12.58265 UTC from Occult4 <Elements> tag)
    double t0_hours = 12.5826525;
    auto t0_utc = time::EpochUTC::from_mjd(time::calendar_to_mjd(2026, 3, 22, t0_hours / 24.0));
    auto t0_tdb = time::to_tdb(t0_utc);

    // 5. AstDyn Nominal Computation using BSP for Vesta
    auto earth_ssb = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t0_tdb);

    double lt = 0.0;
    Eigen::Vector3d rho_vec;
    for(int i=0; i<3; ++i) {
        auto t_emit = t0_tdb - time::TimeDuration::from_seconds(lt);
        double et_emit = (t_emit.jd() - constants::JD2000) * 86400.0;
        
        Eigen::VectorXd v_ssb_emit_raw = vesta_reader.getState(VESTA_NAIF_ID, et_emit);
        Eigen::Vector3d v_ssb_emit = v_ssb_emit_raw.head<3>() * 1000.0; // to meters
        
        rho_vec = v_ssb_emit - earth_ssb.position.to_eigen_si();
        lt = rho_vec.norm() / (constants::C_LIGHT * 1000.0);
    }

    // Positions
    double ra_rad = std::atan2(rho_vec.y(), rho_vec.x());
    double dec_rad = std::asin(rho_vec.z() / rho_vec.norm());
    RightAscension vesta_ra(Angle::from_rad(ra_rad));
    Declination vesta_dec(Angle::from_rad(dec_rad));

    // Rate computation
    auto t1_tdb = t0_tdb + time::TimeDuration::from_seconds(1.0);
    auto t_emit1 = t1_tdb - time::TimeDuration::from_seconds(lt);
    double et_emit1 = (t_emit1.jd() - constants::JD2000) * 86400.0;
    Eigen::VectorXd v_ssb_emit1_raw = vesta_reader.getState(VESTA_NAIF_ID, et_emit1);
    Eigen::Vector3d v_ssb_emit1 = v_ssb_emit1_raw.head<3>() * 1000.0;
    
    auto earth_ssb1 = ephemeris::PlanetaryEphemeris::getState(ephemeris::CelestialBody::EARTH, t1_tdb);
    Eigen::Vector3d rho_vec1 = v_ssb_emit1 - earth_ssb1.position.to_eigen_si();

    double ra_rad1 = std::atan2(rho_vec1.y(), rho_vec1.x());
    double dec_rad1 = std::asin(rho_vec1.z() / rho_vec1.norm());

    // Fix RA jump by wrapping difference to [-pi, pi]
    double dra_val = ra_rad1 - ra_rad;
    while (dra_val > constants::PI) dra_val -= 2.0 * constants::PI;
    while (dra_val < -constants::PI) dra_val += 2.0 * constants::PI;

    Angle dra = Angle::from_rad(dra_val);
    Angle ddec = Angle::from_rad(dec_rad1 - dec_rad);
    physics::Velocity v_rad = physics::Velocity::from_ms(rho_vec1.norm() - rho_vec.norm());

    auto params_astdyn = OccultationLogic::compute_parameters(
        target_star.ra(), target_star.dec(),
        vesta_ra, vesta_dec,
        physics::Distance::from_m(rho_vec.norm()),
        dra, ddec, v_rad
    );

    // 6. Official Occult4 Elements
    double R_E = constants::R_EARTH;
    OccultationParameters params_official;
    params_official.xi_ca = physics::Distance::from_km(-0.1965300 * R_E);
    params_official.eta_ca = physics::Distance::from_km(0.5809135 * R_E);
    params_official.shadow_velocity = physics::Velocity::from_km_s(std::sqrt(std::pow(24.8204827, 2) + std::pow(8.4418234, 2)) * R_E / 3600.0);
    params_official.position_angle = Angle::from_rad(std::atan2(24.8204827, 8.4418234)).wrap_0_2pi();
    params_official.dxi_dt = physics::Velocity::from_ms(24.8204827 * R_E * 1000.0 / 3600.0);
    params_official.deta_dt = physics::Velocity::from_ms(8.4418234 * R_E * 1000.0 / 3600.0);
    params_official.closest_approach_time_offset = time::TimeDuration::zero();

    // 7. Output and Comparison
    std::cout << "\n--- RESULTS COMPARISON (Fundamental Plane) ---" << std::endl;
    std::cout << std::left << std::setw(20) << "Parameter" << std::setw(20) << "AstDyn" << std::setw(20) << "Occult4" << std::setw(20) << "Diff" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    auto print_row = [](const std::string& name, double val1, double val2, const std::string& unit) {
        std::cout << std::left << std::setw(20) << name 
                  << std::fixed << std::setprecision(3) << std::setw(20) << val1 
                  << std::setw(20) << val2 
                  << std::setw(20) << (val1 - val2) << " " << unit << std::endl;
    };

    print_row("xi_ca", params_astdyn.xi_ca.to_km(), params_official.xi_ca.to_km(), "km");
    print_row("eta_ca", params_astdyn.eta_ca.to_km(), params_official.eta_ca.to_km(), "km");
    print_row("Velocity", params_astdyn.shadow_velocity.to_km_s(), params_official.shadow_velocity.to_km_s(), "km/s");
    print_row("Pos Angle", params_astdyn.position_angle.to_deg(), params_official.position_angle.to_deg(), "deg");

    double dist_err = std::sqrt(std::pow((params_astdyn.xi_ca - params_official.xi_ca).to_km(), 2) + 
                               std::pow((params_astdyn.eta_ca - params_official.eta_ca).to_km(), 2));
    
    std::cout << "\nTotal Path Shift: " << dist_err << " km" << std::endl;

    // 8. Mapping
    physics::Distance diameter = physics::Distance::from_km(529.226);
    auto path_astdyn = OccultationMapper::compute_path(params_astdyn, target_star.ra(), target_star.dec(), diameter, t0_utc, time::TimeDuration::from_seconds(3600.0));
    auto path_official = OccultationMapper::compute_path(params_official, target_star.ra(), target_star.dec(), diameter, t0_utc, time::TimeDuration::from_seconds(3600.0));

    OccultationMapper::export_global_svg({path_astdyn, path_official}, {"AstDyn BSP", "Occult4 Official"}, {"#06b6d4", "#ef4444"}, "vesta_2026_map.svg");
    OccultationMapper::export_kml(path_astdyn, "vesta_astdyn.kml");
    OccultationMapper::export_kml(path_official, "vesta_official.kml");

    return 0;
}
