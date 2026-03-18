#include <astdyn/AstDynEngine.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/MultiBodyPropagator.hpp>
#include <astdyn/astrometry/OccultationLogic.hpp>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace astdyn;
using namespace astdyn::propagation;

int main() {
    try {
        std::cout << "=== Haumea System High-Precision Search (MultiBody + Planets) ===\n";
        
        AstDynConfig cfg;
        cfg.integrator_type = IntegratorType::RKF78;
        cfg.tolerance = 1e-14;
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
        cfg.propagator_settings.include_asteroids = false;
        
        AstDynEngine engine(cfg);
        auto de441 = engine.getEphemeris();
        auto provider = de441->getProvider();
        auto integrator = std::make_shared<RKF78Integrator>(60.0, 1e-14);
        
        MultiBodyPropagator mb_prop(integrator, de441);
        auto start_time = time::EpochTDB::from_mjd(61164.0);
        
        std::vector<MultiBodyState> initial;
        MultiBodyState h;
        h.name = "Haumea";
        h.gm = physics::GravitationalParameter::from_km3_s2(267.4);
        h.position = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(
            -5.513013154703953e+09 * 1000.0, -3.563873466876123e+09 * 1000.0, +3.520538626605666e+09 * 1000.0);
        h.velocity = math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(
            +2.411700979939033e+00 * 1000.0, -3.024644217319804e+00 * 1000.0, -2.449900086844126e-01 * 1000.0);
        initial.push_back(h);

        MultiBodyState hi;
        hi.name = "Hi'iaka";
        hi.gm = physics::GravitationalParameter::from_km3_s2(1.2);
        Eigen::Vector3d hi_rel(+3.230691247696847e+04, +2.899365818987711e+04, +1.578694964659445e+04);
        Eigen::Vector3d hi_v_rel(+4.835737731380611e-02, -2.198936577535662e-02, -5.659653697591440e-02);
        hi.position = h.position + math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(hi_rel.x()*1000, hi_rel.y()*1000, hi_rel.z()*1000);
        hi.velocity = h.velocity + math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(hi_v_rel.x()*1000, hi_v_rel.y()*1000, hi_v_rel.z()*1000);
        initial.push_back(hi);

        MultiBodyState na;
        na.name = "Namaka";
        na.gm = physics::GravitationalParameter::from_km3_s2(0.12);
        Eigen::Vector3d na_rel(+2.117923575673045e+04, +5.998545766923913e+03, -1.132606650049148e+04);
        Eigen::Vector3d na_v_rel(-2.695122300908965e-02, -8.481368036922111e-02, -5.418276147615637e-02);
        na.position = h.position + math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(na_rel.x()*1000, na_rel.y()*1000, na_rel.z()*1000);
        na.velocity = h.velocity + math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(na_v_rel.x()*1000, na_v_rel.y()*1000, na_v_rel.z()*1000);
        initial.push_back(na);

        // Batch Integration (From Midnight to cover 6h light time)
        double t_scan_start = 61164.840;
        double t_scan_end = 61164.860;
        double dt_step = 10.0; // 10s steps for the broad cache
        
        std::cout << "Caching trajectory from Midnight to cover Light-Time (6h)...\n";
        struct Point { double mjd; std::vector<MultiBodyState> states; };
        std::vector<Point> cache;
        
        cache.push_back({61164.0, initial});
        
        double cur_mjd = 61164.0;
        while (cur_mjd < t_scan_end + 0.05) {
            double next_mjd = cur_mjd + dt_step / 86400.0;
            auto next_states = mb_prop.propagate(cache.back().states, time::EpochTDB::from_mjd(cur_mjd), time::EpochTDB::from_mjd(next_mjd));
            cache.push_back({next_mjd, next_states});
            cur_mjd = next_mjd;
        }

        auto get_states = [&](double mjd_val) {
            auto it = std::lower_bound(cache.begin(), cache.end(), mjd_val, [](const Point& a, double v){ return a.mjd < v; });
            if (it == cache.begin()) return it->states;
            auto it_prev = std::prev(it);
            double alpha = (mjd_val - it_prev->mjd) / (it->mjd - it_prev->mjd);
            auto res = it_prev->states;
            for (size_t i=0; i<3; ++i) {
                res[i].position = it_prev->states[i].position * (1.0-alpha) + it->states[i].position * alpha;
            }
            return res;
        };

        double s_ra = 220.2437083;
        double s_dec = 14.6737777;
        auto rot_ecl_to_icrf = coordinates::ReferenceFrame::get_rotation<core::ECLIPJ2000, core::GCRF>();
        
        std::cout << "\nScanning for Min Separation...\n";
        struct Best { double sep = 10.0; double t = 0; };
        Best winners[3];

        for (double t = t_scan_start; t <= t_scan_end; t += 0.1 / 86400.0) {
            auto ep_t = time::EpochTDB::from_mjd(t);
            auto p_earth_ssb = provider->getPosition(ephemeris::CelestialBody::EARTH, ep_t) - provider->getPosition(ephemeris::CelestialBody::SUN, ep_t);
            
            auto cur_states = get_states(t);
            for (int i=0; i<3; ++i) {
                double dist_m = (rot_ecl_to_icrf * cur_states[i].position - p_earth_ssb).norm().to_m();
                double tau = dist_m / (constants::C_LIGHT * 1000.0);
                
                auto ret_states = get_states(t - tau/86400.0);
                auto rho = (rot_ecl_to_icrf * ret_states[i].position) - p_earth_ssb;
                auto sky = astrometry::SkyCoord<core::GCRF>::from_vector(rho);
                
                double dra = (sky.ra().to_deg() - s_ra) * 3600.0 * std::cos(s_dec * constants::DEG_TO_RAD);
                double ddec = (sky.dec().to_deg() - s_dec) * 3600.0;
                double sep = std::sqrt(dra*dra + ddec*ddec);
                
                if (sep < winners[i].sep) {
                    winners[i].sep = sep;
                    winners[i].t = t;
                }
            }
        }

        std::cout << "\n=== FINAL RESULTS (Planets Included) ===\n";
        for (int i=0; i<3; ++i) {
            auto utc = time::to_utc(time::EpochTDB::from_mjd(winners[i].t));
            std::cout << std::left << std::setw(15) << initial[i].name 
                      << std::setw(15) << std::fixed << std::setprecision(5) << winners[i].sep 
                      << winners[i].t << " MJD (" << utc.mjd() << ")\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
    return 0;
}
