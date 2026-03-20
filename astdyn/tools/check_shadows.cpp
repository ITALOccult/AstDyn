#include <astdyn/AstDynEngine.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/MultiBodyPropagator.hpp>
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
        auto provider = de441->getProvider();
        auto integrator = std::make_shared<RKF78Integrator>(60.0, 1e-14);
        MultiBodyPropagator mb_prop(integrator, de441);
        
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
        
        auto check_body = [&](int idx, double tca, double diam_km) {
            auto ep = time::EpochTDB::from_mjd(tca);
            auto st = mb_prop.propagate(initial, start_time, ep)[idx];
            auto e_ssb = provider->getPosition(ephemeris::CelestialBody::EARTH, ep) - provider->getPosition(ephemeris::CelestialBody::SUN, ep);
            auto rho = (rot_ecl_to_icrf * st.position) - e_ssb;
            auto sc = SkyCoord<core::GCRF>::from_vector(rho);
            
            double sep = std::sqrt(std::pow((sc.ra().to_deg() - s_ra)*std::cos(s_dec*constants::DEG_TO_RAD), 2) + std::pow(sc.dec().to_deg() - s_dec, 2)) * 3600.0;
            double dist_m = rho.norm().to_m();
            double earth_rad_arcsec = (6371000.0 / dist_m) * constants::RAD_TO_DEG * 3600.0;
            double body_rad_arcsec = (diam_km * 500.0 / dist_m) * constants::RAD_TO_DEG * 3600.0;
            
            std::cout << "Body " << idx << " (" << st.name << ") at MJD " << std::fixed << std::setprecision(5) << tca << ":\n";
            std::cout << "  Separation: " << std::setprecision(4) << sep << " arcsec\n";
            std::cout << "  Earth Radius: " << earth_rad_arcsec << " arcsec\n";
            std::cout << "  Body Radius:  " << body_rad_arcsec << " arcsec\n";
            
            if (sep < earth_rad_arcsec + body_rad_arcsec) {
                std::cout << "  STATUS: SHADOW ON EARTH\n";
            } else {
                std::cout << "  STATUS: SHADOW MISSES EARTH\n";
            }
        };

        check_body(0, 61164.84739, 1500.0);
        check_body(1, 61164.85138, 170.0);

    } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << "\n"; }
    return 0;
}
