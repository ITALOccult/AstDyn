/**
 * @file astdyn_chebyshev_demo.cpp
 * @brief Demonstration of the AsteroidChebyshevEphemeris class.
 */

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/astrometry/AsteroidChebyshevEphemeris.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include <iostream>
#include <iomanip>
#include <filesystem>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: astdyn_chebyshev_demo <asteroid_id/name> [start_mjd] [duration_days]\n";
        return 1;
    }

    std::string target = argv[1];
    double start_mjd = (argc > 2) ? std::stod(argv[2]) : 61121.0; // Default March 2026
    double duration = (argc > 3) ? std::stod(argv[3]) : 7.0;

    // --- 1. System Setup ---
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de440s.bsp";
    if (!std::filesystem::exists(bsp_path)) {
        bsp_path = "de440s.bsp";
    }

    astdyn::AstDynEngine engine;
    astdyn::AstDynConfig cfg;
    cfg.ephemeris_type = astdyn::EphemerisType::DE441;
    cfg.ephemeris_file = bsp_path;
    cfg.verbose = false;
    engine.set_config(cfg);

    try {
        // Essential: Provider must be set for planetary ephemerides
        auto provider = std::make_shared<astdyn::ephemeris::DE441Provider>(bsp_path);
        astdyn::ephemeris::PlanetaryEphemeris::setGlobalProvider(provider);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not load DE441: " << e.what() << "\n";
    }

    // --- 2. Input Data Retrieval ---
    try {
        astdyn::time::EpochTDB start_epoch = astdyn::time::EpochTDB::from_mjd(start_mjd);
        
        std::cout << "Fetching Cartesian state for " << target << " from JPL Horizons...\n";
        astdyn::io::HorizonsClient horizons;
        auto state_vec = horizons.query_vectors(target, start_epoch);
        if (!state_vec) {
             throw std::runtime_error("Failed to fetch data from Horizons.");
        }
        
        // Convert to Keplerian for initial state bridge
        auto state_vec_ecl = state_vec->cast_frame<astdyn::core::ECLIPJ2000>();
        auto elements = astdyn::propagation::cartesian_to_keplerian<astdyn::core::ECLIPJ2000>(state_vec_ecl);
        
        astdyn::time::EpochTDB end_epoch = astdyn::time::EpochTDB::from_mjd(start_mjd + duration);

        std::cout << "Calculating daily Chebyshev polynomials for " << duration << " days...\n";
        // Polynomial degree 12 is usually enough for 1-day segments (sub-milliarcsec precision)
        astdyn::astrometry::AsteroidChebyshevEphemeris ephem(elements, start_epoch, end_epoch, engine.config(), 12);

        std::cout << "\nComparison with high-precision AstrometryReducer (inclusive of light-time & aberration):\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  Time [MJD]  |  Cheby RA [deg] |  Truth RA [deg] |  Error [arcsec]\n";
        std::cout << "--------------|-----------------|-----------------|----------------\n";

        astdyn::astrometry::AstrometricSettings a_settings;
        a_settings.light_time_correction = true;
        a_settings.aberrazione_differenziale = false; // Chebyshev typically fits astrometric positions
        a_settings.frame_conversion_to_equatorial = true;

        for (double t = start_mjd; t <= start_mjd + duration; t += 1.0) {
            auto t_epoch = astdyn::time::EpochTDB::from_mjd(t);
            
            // 1. Position from Chebyshev
            auto eval_res = ephem.evaluate(t_epoch);
            double c_ra = std::get<0>(eval_res);

            // 2. Position from high-precision reducer
            auto res = astdyn::astrometry::AstrometryReducer::compute_observation(
                elements, elements.epoch, t_epoch, engine.config(), a_settings
            );
            
            if (!res) continue;
            
            double truth_ra_deg = res->ra.to_deg();
            double truth_dec_deg = res->dec.to_deg();
            
            double cos_dec = std::cos(truth_dec_deg * M_PI / 180.0);
            double error_arcsec = std::abs(c_ra - truth_ra_deg) * 3600.0 * cos_dec;

            std::cout << "  " << t << " |    " << c_ra << "   |    " << truth_ra_deg << "   |    " << std::scientific << std::setprecision(3) << error_arcsec << "\n" << std::fixed << std::setprecision(6);
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
