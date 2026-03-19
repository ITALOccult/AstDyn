#include "astdyn/AstDyn.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.ephemeris_type = EphemerisType::DE441;
        
        static AstDynEngine engine(cfg);
        static io::HorizonsClient horizons;
        
        // Target: May 4, 2026 00:00 UTC (JD 2461164.5 UTC)
        // JD TDB = JD UTC + (deltaT + 32.184) / 86400.
        // For 2026, deltaT is ~69.1s. Total ~101.3s = 0.00117d.
        double jd_utc = 2461164.5;
        auto t_utc = time::EpochTDB::from_jd(jd_utc); 
        
        std::cout << "Querying Horizons for (234) Barbara VECTORS at 2026-05-04 00:00 UTC..." << std::endl;
        auto state_opt = horizons.query_vectors("234", t_utc);
        
        if (state_opt) {
            auto cart = state_opt.value(); // This is the state at 00:00 UTC according to HorizonsClient
            
            astrometry::AstrometricSettings a_settings;
            a_settings.light_time_correction = true;
            a_settings.aberrazione_differenziale = true;
            a_settings.deflessione_relativistica = true;
            a_settings.frame_conversion_to_equatorial = true;
            
            std::cout << "Running Truth Test: Astrometery on JPL Vector directly..." << std::endl;
            // Use compute_observation_from_cartesian to avoid integration errors
            auto obs_res = astrometry::AstrometryReducer::compute_observation_from_cartesian(cart, t_utc, t_utc, cfg, a_settings);
            
            if (obs_res) {
                auto obs = obs_res.value();
                std::cout << "\n=== TRUTH TEST RESULTS ===\n";
                // JPL Expected: 13 11 51.10 +13 47 03.2
                std::cout << "JPL EXPECTED: 13 11 51.100 +13 47 03.20\n";
                std::cout << "ASTDYN RESULT: " << obs.ra.to_hms() << " " << obs.dec.to_dms() << "\n";
                std::cout << "RA (deg):  " << std::fixed << std::setprecision(8) << obs.ra.to_deg() << "\n";
                std::cout << "Dec (deg): " << std::fixed << std::setprecision(8) << obs.dec.to_deg() << "\n";
                
                double dra_s = (obs.ra.to_deg() - 197.96291667) * 3600.0 / 15.0; // 13 11 51.10 = 197.96291667
                double ddec_as = (obs.dec.to_deg() - 13.78422222) * 3600.0;     // 13 47 03.20 = 13.78422222
                
                std::cout << "Residuals (AstDyn - JPL):\n";
                std::cout << "  Delta RA:  " << dra_s << " s\n";
                std::cout << "  Delta Dec: " << ddec_as << " arcsec\n";
            }
        }
    } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << "\n"; }
    return 0;
}
