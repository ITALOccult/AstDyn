/**
 * @file occultation_catalog_benchmark.cpp
 * @brief Benchmark comparing local vs online stellar catalog performance.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/time/epoch.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>

using namespace astdyn;
using namespace astdyn::astrometry;

void run_test(const std::string& name, const std::string& config_json) {
    std::cout << "[INFO] Setup Engine..." << std::endl;
    AstDynConfig cfg;
    cfg.ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    cfg.verbose = true;
    AstDynEngine engine(cfg);

    std::cout << "[INFO] Initialize Catalog with: " << config_json << std::endl;
    try {
        catalog::GaiaDR3Catalog::initialize(config_json);
    } catch (const std::exception& e) {
        std::cerr << "FAILED_INIT " << name << ": " << e.what() << std::endl;
        return;
    }

    auto vesta_elements = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_jd(2461041.5), // 2026-01-01
        2.3619, 0.0888, 7.14, 103.85, 151.20, 318.52, 
        physics::GravitationalParameter::sun()
    );

    time::EpochTDB start = time::EpochTDB::from_jd(2461121.5); // 2026-03-21 12:00
    time::EpochTDB end   = time::EpochTDB::from_jd(2461122.5); // 2026-03-22 12:00

    std::cout << "[INFO] Running find_occultations..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    auto events = OccultationLogic::find_occultations(
        "4", vesta_elements, start, end, 13.0, engine, OccultationRefinementMode::ChebyshevDaily);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    std::cout << "NAME: " << name << std::endl;
    std::cout << "TIME: " << duration << " ms" << std::endl;
    std::cout << "COUNT: " << events.size() << std::endl;

    std::cout << "[INFO] Shutdown Catalog..." << std::endl;
    catalog::GaiaDR3Catalog::shutdown();
}

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    std::string mode = argv[1];
    
    if (mode == "disk") {
        std::string sqlite_cfg = R"({
            "catalog_type": "sqlite_dr3", 
            "database_path": "/Users/michelebigi/.catalog/crossreference/gaia_dr3_occult_pro.db"
        })";
        run_test("DISK_SQLITE", sqlite_cfg);
    } else if (mode == "online") {
        std::string online_cfg = R"({
            "catalog_type": "online_esa",
            "timeout_seconds": 60
        })";
        run_test("ONLINE_ESA", online_cfg);
    }
    return 0;
}
