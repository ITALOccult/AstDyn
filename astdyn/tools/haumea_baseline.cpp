#include <astdyn/AstDynEngine.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/astrometry/OccultationLogic.hpp>
#include <astdyn/astrometry/OccultationMapper.hpp>
#include <astdyn/io/HorizonsClient.hpp>
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::astrometry;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.verbose = true;
        
        AstDynEngine engine(cfg);
        auto de441 = engine.getEphemeris();
        
        // 1. Initial State from Horizons (J2026.34 approx)
        astdyn::io::HorizonsClient horizons;
        auto start_time = time::EpochTDB::from_jd(2461165.347);
        auto state_h = horizons.query_vectors("136108", start_time);
        if (!state_h) throw std::runtime_error("Horizons query failed");

        engine.set_initial_orbit(propagation::cartesian_to_keplerian<core::ECLIPJ2000>(state_h->cast_frame<core::ECLIPJ2000>()));
        
        // 2. Star Coordinates + 8.5 mas Offset
        double s_ra_base = 220.2437083;
        double s_dec = 14.6737777;
        double s_ra = s_ra_base + (8.5 / 3600000.0) / std::cos(s_dec * constants::DEG_TO_RAD);
        auto star_ra = RightAscension::from_deg(s_ra);
        auto star_dec = Declination::from_deg(s_dec);
        
        // 3. Compute Path
        auto t_ca = start_time;
        auto app_ast = engine.compute_asteroid_apparent_place(t_ca, "");
        std::cout << "AstDyn RA: " << app_ast.ra.to_deg() << " deg, Dec: " << app_ast.dec.to_deg() << " deg\n";
        
        auto hor_obs = horizons.query_observation("136108", start_time);
        if (hor_obs) {
            std::cout << "Horizons RA: " << hor_obs->ra.to_deg() << " deg, Dec: " << hor_obs->dec.to_deg() << " deg\n";
        }
        
        std::cout << "Star RA: " << star_ra.to_deg() << " deg, Dec: " << star_dec.to_deg() << " deg\n";
        
        // Use logic to get precise TCA
        auto ep_t1 = time::EpochTDB::from_mjd(t_ca.mjd() - 0.1/24.0);
        auto ep_t2 = time::EpochTDB::from_mjd(t_ca.mjd() + 0.1/24.0);
        auto app1 = engine.compute_asteroid_apparent_place(ep_t1, "");
        auto app2 = engine.compute_asteroid_apparent_place(ep_t2, "");
        
        double dra_dt = (app2.ra.to_deg() - app1.ra.to_deg()) / (0.2/24.0);
        double ddec_dt = (app2.dec.to_deg() - app1.dec.to_deg()) / (0.2/24.0);
        
        auto params = OccultationLogic::compute_parameters(
            star_ra, star_dec, app_ast.ra, app_ast.dec,
            app_ast.dist,
            Angle::from_deg(dra_dt), Angle::from_deg(ddec_dt),
            physics::Velocity::from_ms(0), t_ca, de441
        );
        
        auto path = OccultationMapper::compute_path(params, star_ra, star_dec, physics::Distance::from_km(2224.0), time::to_utc(t_ca), de441);
        
        OccultationMapper::export_global_svg({path}, {"Haumea Baseline (8.5 mas)"}, {"#ef4444"}, "/Users/michelebigi/.gemini/antigravity/brain/df2d9270-f3cf-4f0d-b697-705997e6ac28/haumea_baseline_calib.svg", de441, Angle::from_deg(15), Angle::from_deg(30), 2.5);
        
        std::cout << "Baseline map generated at /Users/michelebigi/.gemini/antigravity/brain/df2d9270-f3cf-4f0d-b697-705997e6ac28/haumea_baseline_calib.svg\n";
        std::cout << "Center: " << params.center_lat.to_deg() << " N, " << params.center_lon.to_deg() << " E\n";
        std::cout << "Impact Parameter: " << params.impact_parameter.to_km() << " km\n";

    } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << "\n"; }
    return 0;
}
