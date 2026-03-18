#include <astdyn/AstDynEngine.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/RelativeMultiBodyPropagator.hpp>
#include <astdyn/astrometry/OccultationLogic.hpp>
#include <astdyn/astrometry/OccultationMapper.hpp>
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::astrometry;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
        
        AstDynEngine engine(cfg);
        auto de441 = engine.getEphemeris();
        
        PropagatorSettings settings;
        settings.include_relativity = true;
        settings.include_asteroids = false;
        settings.include_earth_j2 = true;
        
        auto force_field = std::make_shared<ForceField>(settings, de441);
        auto integrator = std::make_shared<RKF78Integrator>(30.0, 1e-15);
        RelativeMultiBodyPropagator rel_prop(integrator, force_field);
        
        double s_ra = 220.2437083;
        double s_dec = 14.6737777;
        auto star_ra = RightAscension::from_deg(s_ra);
        auto star_dec = Declination::from_deg(s_dec);
        
        auto start_time = time::EpochTDB::from_mjd(61164.0);
        
        std::vector<MultiBodyState> initial;
        MultiBodyState h; h.name = "Haumea"; h.gm = physics::GravitationalParameter::from_km3_s2(267.4);
        h.position = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(-5.513013154703953e+09 * 1000.0, -3.563873466876123e+09 * 1000.0, +3.520538626605666e+09 * 1000.0);
        h.velocity = math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(2.411700979939033e+00 * 1000.0, -3.024644217319804e+00 * 1000.0, -2.449900086844126e-01 * 1000.0);
        initial.push_back(h);

        MultiBodyState na; na.name = "Namaka"; na.gm = physics::GravitationalParameter::from_km3_s2(0.12);
        Eigen::Vector3d na_rel(+2.117923575673045e+04, +5.998545766923913e+03, -1.132606650049148e+04);
        Eigen::Vector3d na_v_rel(-2.695122300908965e-02, -8.481368036922111e-02, -5.418276147615637e-02);
        na.position = h.position + math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(na_rel.x()*1000, na_rel.y()*1000, na_rel.z()*1000);
        na.velocity = h.velocity + math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(na_v_rel.x()*1000, na_v_rel.y()*1000, na_v_rel.z()*1000);
        initial.push_back(na);
        
        auto rot_ecl_to_icrf = coordinates::ReferenceFrame::get_rotation<core::ECLIPJ2000, core::GCRF>();
        
        auto get_path = [&](int body_idx, double tca_mjd, double diam_km) {
            auto ep_tca = time::EpochTDB::from_mjd(tca_mjd);
            auto states = rel_prop.propagate(initial, start_time, ep_tca);
            
            auto provider = de441->getProvider();
            auto p_earth_ssb = provider->getPosition(ephemeris::CelestialBody::EARTH, ep_tca) - provider->getPosition(ephemeris::CelestialBody::SUN, ep_tca);
            double dist_m = (rot_ecl_to_icrf * states[body_idx].position - p_earth_ssb).norm().to_m();
            double tau = dist_m / (constants::C_LIGHT * 1000.0);
            
            auto states_ret = rel_prop.propagate(initial, start_time, time::EpochTDB::from_mjd(ep_tca.mjd() - tau/86400.0));
            auto rho = (rot_ecl_to_icrf * states_ret[body_idx].position) - p_earth_ssb;
            auto sky0 = SkyCoord<core::GCRF>::from_vector(rho);
            
            auto ep_t1 = time::EpochTDB::from_mjd(tca_mjd - 0.1/24.0);
            auto ep_t2 = time::EpochTDB::from_mjd(tca_mjd + 0.1/24.0);
            auto s1 = rel_prop.propagate(initial, start_time, ep_t1);
            auto s2 = rel_prop.propagate(initial, start_time, ep_t2);
            
            auto sky1 = SkyCoord<core::GCRF>::from_vector( (rot_ecl_to_icrf * s1[body_idx].position) - (provider->getPosition(ephemeris::CelestialBody::EARTH, ep_t1) - provider->getPosition(ephemeris::CelestialBody::SUN, ep_t1)) );
            auto sky2 = SkyCoord<core::GCRF>::from_vector( (rot_ecl_to_icrf * s2[body_idx].position) - (provider->getPosition(ephemeris::CelestialBody::EARTH, ep_t2) - provider->getPosition(ephemeris::CelestialBody::SUN, ep_t2)) );
            
            double dra_dt = (sky2.ra().to_deg() - sky1.ra().to_deg()) / (0.2/24.0);
            double ddec_dt = (sky2.dec().to_deg() - sky1.dec().to_deg()) / (0.2/24.0);
            
            auto params = OccultationLogic::compute_parameters(
                star_ra, star_dec, sky0.ra(), sky0.dec(),
                physics::Distance::from_m(dist_m),
                Angle::from_deg(dra_dt), Angle::from_deg(ddec_dt),
                physics::Velocity::from_ms(0), ep_tca, de441
            );
            
            return OccultationMapper::compute_path(params, star_ra, star_dec, physics::Distance::from_km(diam_km), time::to_utc(ep_tca), de441);
        };
        
        auto path_h = get_path(0, 61164.84742, 1100.0);
        auto path_n = get_path(1, 61164.85141, 170.0);
        
        std::vector<OccultationPath> paths = {path_h, path_n};
        std::vector<std::string> labels = {"Haumea (136108)", "Namaka Shadow"};
        std::vector<std::string> colors = {"#f43f5e", "#0ea5e9"};
        
        OccultationMapper::export_global_svg(paths, labels, colors, "haumea_system_occultation.svg", de441, Angle::from_deg(15), Angle::from_deg(30), 2.5);
        OccultationMapper::export_kml(paths, labels, "haumea_system_occultation.kml");
        
        std::cout << "Relative high-precision map generated with corrected ForceField.\n";

    } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << "\n"; }
    return 0;
}
