/**
 * @file simple_catalog_benchmark.cpp
 * @brief Simple benchmark for Gaia catalog performance.
 */

#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>

using namespace astdyn;

void run_bench(const std::string& name, const std::string& config) {
    std::cout << "\nTesting " << name << "..." << std::endl;
    try {
        catalog::GaiaDR3Catalog::initialize(config);
        auto& cat = catalog::GaiaDR3Catalog::instance();
        
        catalog::ConeQuery q;
        q.ra = catalog::RightAscension(catalog::Angle::from_deg(100.0));
        q.dec = catalog::Declination(catalog::Angle::from_deg(20.0));
        q.radius = catalog::Angle::from_arcsec(300.0);
        q.max_magnitude = 14.0;

        // Warmup
        cat.query_cone(q);

        auto start = std::chrono::high_resolution_clock::now();
        int total_stars = 0;
        for(int i=0; i<5; ++i) {
            auto stars = cat.query_cone(q);
            total_stars = stars.size();
        }
        auto end = std::chrono::high_resolution_clock::now();
        double avg_ms = std::chrono::duration<double, std::milli>(end - start).count() / 5.0;

        std::cout << "  - Result: " << total_stars << " stars" << std::endl;
        std::cout << "  - Average Time: " << std::fixed << std::setprecision(2) << avg_ms << " ms" << std::endl;
        
        catalog::GaiaDR3Catalog::shutdown();
    } catch (const std::exception& e) {
        std::cerr << "  - ERROR: " << e.what() << std::endl;
    }
}

int main() {
    std::string online_cfg = R"({"catalog_type":"online_esa","timeout_seconds":30})";
    std::string disk_cfg = R"({"catalog_type":"sqlite_dr3", "database_path":"/Users/michelebigi/.catalog/crossreference/gaia_dr3_occult_pro.db"})";

    run_bench("DISK (SQLite)", disk_cfg);
    run_bench("ONLINE (ESA)", online_cfg);

    return 0;
}
