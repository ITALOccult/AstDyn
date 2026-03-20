#include <astdyn/AstDynEngine.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/time/TimeScale.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/coordinates/CartesianState.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <astdyn/math/frame_algebra.hpp>
#include <astdyn/core/frame_tags.hpp>
#include <astdyn/core/physics_state.hpp>
#include <astdyn/core/Enums.hpp>
#include <astdyn/ephemeris/CelestialBody.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

using namespace astdyn;

int main() {
    try {
        std::cout << "=== Haumea Occultation Prediction Verification (arXiv:2603.15049) ===\n";
        
        double km_to_m = 1000.0;
        
        // --- 1. CONFIG ---
        astdyn::AstDynConfig cfg;
        cfg.integrator_type = astdyn::IntegratorType::RKF78;
        cfg.initial_step_size = 0.1;
        cfg.tolerance = 1e-13;
        cfg.ephemeris_type = astdyn::EphemerisType::DE441;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
        cfg.propagator_settings.include_asteroids = false;

        astdyn::AstDynEngine engine(cfg);
        
        // 2. Haumea @ 00:00 TDB (Ecliptic J2000)
        double dX = -5.513011317885303E+09, dY = -3.563870643567566E+09, dZ = 3.520535555074772E+09;
        double dVX = 2.411699650465367E+00, dVY = -3.024638539222640E+00, dVZ = -2.449948126755586E-01;
        
        auto epoch0 = time::EpochTDB::from_jd(2461164.5);
        auto cart0_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>::from_si(
            epoch0, dX * km_to_m, dY * km_to_m, dZ * km_to_m,
            dVX * km_to_m, dVY * km_to_m, dVZ * km_to_m, constants::GM_SUN * 1e9
        );
        engine.set_initial_orbit(propagation::cartesian_to_keplerian(cart0_ecl));
        
        // 3. Star (ICRF) 14 40 58.49 +14 40 25.6
        double s_ra = (14.0 + 40.0/60.0 + 58.49/3600.0) * 15.0;
        double s_dec = 14.0 + 40.0/60.0 + 25.6/3600.0;
        Eigen::Vector3d star_icrf;
        star_icrf << std::cos(s_dec * constants::DEG_TO_RAD)*std::cos(s_ra * constants::DEG_TO_RAD),
                     std::cos(s_dec * constants::DEG_TO_RAD)*std::sin(s_ra * constants::DEG_TO_RAD),
                     std::sin(s_dec * constants::DEG_TO_RAD);

        auto rot_ecl_to_icrf = coordinates::ReferenceFrame::get_rotation<core::ECLIPJ2000, core::GCRF>();
        auto provider = engine.getEphemeris()->getProvider();
        
        std::cout << "Starting cross-frame scan (Haumea-Ecliptic, Earth-ICRF)...\n";
        
        double min_sep = 1.0; 
        double tca_mjd = 0;
        
        // MJD window around 20:17 UT (61164.845)
        for (double t = 61164.83; t <= 61164.86; t += 1.0 / 86400.0) { 
            auto epoch_t = time::EpochTDB::from_mjd(t);
            
            // 1. Get Earth & Sun geocentric/SSB ICRF from provider
            auto p_earth_icrf = provider->getPosition(ephemeris::CelestialBody::EARTH, epoch_t);
            auto p_sun_icrf = provider->getPosition(ephemeris::CelestialBody::SUN, epoch_t);
            // Earth Heliocentric (ICRF)
            Eigen::Vector3d earth_helio_icrf = p_earth_icrf.to_eigen_si() - p_sun_icrf.to_eigen_si();
            
            // 2. Get Haumea Heliocentric (Ecliptic)
            auto state = engine.propagate_to(epoch_t);
            auto cart_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(state);
            
            // 3. Initial Range estimation in Ecliptic
            // (Just to get a rough tau, frame doesn't matter for norm if both are same, but I'll use Ecliptic Earth for tau)
            // Wait, let's just use ICRF for everything.
            auto haumea_icrf = (rot_ecl_to_icrf * cart_ecl.position).to_eigen_si();
            double range = (haumea_icrf - earth_helio_icrf).norm();
            double tau = range / 299792458.0;
            
            // 4. Corrected Haumea (Retarded)
            auto state_ret = engine.propagate_to(time::EpochTDB::from_mjd(t - tau / 86400.0));
            auto cart_ret_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(state_ret);
            auto haumea_ret_icrf = (rot_ecl_to_icrf * cart_ret_ecl.position).to_eigen_si();
            
            // 5. Geocentric ICRF relative vector
            Eigen::Vector3d rho_icrf = haumea_ret_icrf - earth_helio_icrf;
            
            double sep = std::acos(rho_icrf.normalized().dot(star_icrf));
            if (sep < min_sep) {
                min_sep = sep;
                tca_mjd = t;
            }
        }
        
        auto [y, mon, d_cal, f] = time::mjd_to_calendar(tca_mjd);
        int hh = static_cast<int>(f * 24.0);
        int mm = static_cast<int>((f * 24.0 - hh) * 60.0);
        double ss = ((f * 24.0 - hh) * 60.0 - mm) * 60.0;
        
        std::cout << "\n=== BENCHMARK RESULTS ===\n";
        std::cout << "Found TCA (TDB):      " << y << "-" << mon << "-" << d_cal << " " << hh << ":" << mm << ":" << ss << "\n";
        std::cout << "Min Angular Separation: " << min_sep * constants::RAD_TO_DEG * 3600.0 << " arcsec\n";
        
        std::cout << "\narXiv:2603.15049 Reference UT: ~20:17:00\n";
        std::cout << "AstDyn Predicted UT:         " << hh << ":" << mm << ":" << (ss - 69.1) << " UT\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
    return 0;
}
