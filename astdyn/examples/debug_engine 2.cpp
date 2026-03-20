#include "astdyn/AstDynEngine.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    // 1. Setup Engine with DE441
    AstDynConfig config;
    config.ephemeris_type = "DE441";
    config.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    config.propagator_settings.include_planets = true;
    config.propagator_settings.include_moon = true;
    config.propagator_settings.include_relativity = true;
    config.propagator_settings.integrate_in_ecliptic = false;
    
    AstDynEngine engine(config);
    
    // 2. Query Horizons for Baruffetti (70621) at specific epoch
    io::HorizonsClient horizons;
    time::EpochTDB t_ref = time::EpochTDB::from_jd(2460320.5); // 2024-01-11 00:00 TDB
    
    std::string asteroid = "70621";
    auto h_vec = horizons.query_vectors(asteroid, t_ref);
    auto h_obs = horizons.query_observation(asteroid, t_ref, "500@399"); // Earth Geocenter
    
    if (!h_vec || !h_obs) {
        std::cerr << "Horizons query failed!\n";
        return 1;
    }
    
    // 3. Set the state in the engine
    auto cart = *h_vec;
    auto kep = propagation::cartesian_to_keplerian<core::GCRF>(cart);
    engine.set_initial_orbit(kep);
    
    // 4. Compute residual (Engine way)
    observations::OpticalObservation obs;
    obs.time = time::to_utc(t_ref);
    obs.ra = h_obs->ra.value;
    obs.dec = h_obs->dec.value;
    obs.observatory_code = "500"; // Geocenter
    
    // We can't call compute_residual directly (it's private or in calculator)
    // But we can call engine.fit_orbit() with 1 iteration to see the first residuals
    engine.add_observation(obs);
    
    // Manually trigger residual calculation via public-ish if possible, 
    // or just let it fit.
    // Actually, fit_orbit prints it.
    std::cout << "\n=== Engine Residual Verification ===\n";
    auto fit_res = engine.fit_orbit();
    
    std::cout << "Horizons RA:  " << h_obs->ra.value * constants::RAD_TO_DEG << " deg\n";
    std::cout << "Horizons Dec: " << h_obs->dec.value * constants::RAD_TO_DEG << " deg\n";
    
    if (!fit_res.residuals_ra.empty()) {
        std::cout << "Engine RMS RA:  " << fit_res.rms_ra << " arcsec\n";
        std::cout << "Engine RMS Dec: " << fit_res.rms_dec << " arcsec\n";
    } else {
        std::cout << "No residuals (rejected?)\n";
        std::cout << "Num rejected: " << fit_res.num_rejected << "\n";
    }

    return 0;
}
