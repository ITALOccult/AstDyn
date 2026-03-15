/**
 * @file astdyn_multi_chebyshev_demo.cpp
 * @brief Demonstration of the ChebyshevEphemerisManager class with multiple asteroids.
 */

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/astrometry/ChebyshevEphemerisManager.hpp"
#include "astdyn/astrometry/Astrometry.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include <iostream>
#include <iomanip>
#include <filesystem>

int main() {
    // --- 1. System Setup ---
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de440s.bsp";
    if (!std::filesystem::exists(bsp_path)) bsp_path = "de440s.bsp";

    astdyn::AstDynEngine engine;
    astdyn::AstDynConfig cfg;
    cfg.ephemeris_type = astdyn::EphemerisType::DE441;
    cfg.ephemeris_file = bsp_path;
    cfg.verbose = false;
    engine.set_config(cfg);

    try {
        auto provider = std::make_shared<astdyn::ephemeris::DE441Provider>(bsp_path);
        astdyn::ephemeris::PlanetaryEphemeris::setProvider(provider);
    } catch (...) {}

    // --- 2. Multiple Asteroids Configuration ---
    std::vector<std::string> targets = {"1", "2", "4"}; // Ceres, Pallas, Vesta
    double start_mjd = 61121.0;
    double duration = 5.0;
    astdyn::time::EpochTDB start_epoch = astdyn::time::EpochTDB::from_mjd(start_mjd);
    astdyn::time::EpochTDB end_epoch = astdyn::time::EpochTDB::from_mjd(start_mjd + duration);

    astdyn::astrometry::ChebyshevEphemerisManager manager(engine.config());
    astdyn::io::HorizonsClient horizons;

    std::cout << "Pre-calculating polynomials for asteroids: ";
    for (const auto& id : targets) {
        std::cout << id << " ";
        auto state_vec = horizons.query_vectors(id, start_epoch);
        if (state_vec) {
            auto elements = astdyn::propagation::cartesian_to_keplerian<astdyn::core::ECLIPJ2000>(
                state_vec->cast_frame<astdyn::core::ECLIPJ2000>()
            );
            manager.add_asteroid(id, elements, start_epoch, end_epoch);
        }
    }
    std::cout << "DONE.\n\n";

    // --- 3. Evaluate Position AND Velocity ---
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Comparison for Vesta (ID: 4) at " << start_mjd + 2.5 << " MJD:\n";
    
    auto t_query = astdyn::time::EpochTDB::from_mjd(start_mjd + 2.5);
    auto [pos, vel] = manager.evaluate_full("4", t_query);

    auto [ra, dec, dist] = pos;
    auto [vra, vdec, vdist] = vel;

    std::cout << "  Position: RA=" << ra << " deg, Dec=" << dec << " deg, Dist=" << dist << " AU\n";
    std::cout << "  Velocity: vRA=" << vra << " deg/day, vDec=" << vdec << " deg/day, vDist=" << vdist << " AU/day\n";

    // Numerical check for velocity
    double dt = 0.001; // 1.44 minutes
    auto [p1, v1] = manager.evaluate_full("4", astdyn::time::EpochTDB::from_mjd(start_mjd + 2.5 - dt/2.0));
    auto [p2, v2] = manager.evaluate_full("4", astdyn::time::EpochTDB::from_mjd(start_mjd + 2.5 + dt/2.0));
    
    double num_vra = (std::get<0>(p2) - std::get<0>(p1)) / dt;
    double num_vdec = (std::get<1>(p2) - std::get<1>(p1)) / dt;

    std::cout << "\nConsistency Check (Analytical vs Numerical Derivative):\n";
    std::cout << "  vRA err:  " << std::scientific << std::setprecision(2) << std::abs(vra - num_vra) << " deg/day\n";
    std::cout << "  vDec err: " << std::scientific << std::setprecision(2) << std::abs(vdec - num_vdec) << " deg/day\n";

    return 0;
}
