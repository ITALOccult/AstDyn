#include "astdyn/AstDyn.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    try {
        AstDynConfig cfg;
        cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
        cfg.ephemeris_type = EphemerisType::DE441;
        
        // --- NEW HIGH PRECISION SETTINGS ---
        cfg.integrator_type = IntegratorType::AAS;
        cfg.aas_precision = 1e-4; // Standard for high precision
        cfg.propagator_settings.include_asteroids = true;
        cfg.propagator_settings.use_default_asteroid_set = true;
        cfg.propagator_settings.include_relativity = true;
        // ------------------------------------
        
        static AstDynEngine engine(cfg);
        static io::HorizonsClient horizons;
        
        // Target: May 4, 2026 00:00 UTC
        auto t_target = time::EpochTDB::from_jd(2461164.5);
        
        std::cout << "Querying Horizons for (234) Barbara state (ICRF/GCRF)..." << std::endl;
        auto state_opt = horizons.query_vectors("234", t_target);
        
        if (state_opt) {
            auto keplerian = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(state_opt.value().cast_frame<core::ECLIPJ2000>());
            
            astrometry::AstrometricSettings a_settings;
            a_settings.light_time_correction = true;
            a_settings.aberrazione_differenziale = true;
            a_settings.deflessione_relativistica = true;
            a_settings.frame_conversion_to_equatorial = true;
            
            std::cout << "Computing with AAS + Massive Asteroids for May 4, 2026..." << std::endl;
            auto obs_res = astrometry::AstrometryReducer::compute_observation(keplerian, t_target, t_target, cfg, a_settings);
            
            if (obs_res) {
                auto obs = obs_res.value();
                std::cout << "\n=== RESULTS (234) BARBARA (AAS + ASTEROIDS) ===\n";
                std::cout << "RA:  " << obs.ra.to_hms() << " (" << std::fixed << std::setprecision(8) << obs.ra.to_deg() << " deg)\n";
                std::cout << "Dec: " << obs.dec.to_dms() << " (" << std::fixed << std::setprecision(8) << obs.dec.to_deg() << " deg)\n";
                std::cout << "Dist: " << obs.distance.to_au() << " AU\n";
                std::cout << "=============================\n";
            } else {
                std::cerr << "Astrometry calculation failed!\n";
            }
        } else {
            std::cerr << "Could not get state from Horizons!\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
    return 0;
}
