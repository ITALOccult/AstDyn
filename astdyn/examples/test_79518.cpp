#include <astdyn/AstDyn.hpp>
#include <iostream>
#include <iomanip>

int main() {
    using namespace astdyn;
    using namespace astdyn::io;

    HorizonsClient horizons;

    // Epoch: "Tomorrow morning" -> 2026-03-07 09:00:00 UTC
    // MJD for 2026-03-07 00:00:00 is approx 61106.0
    // 09:00:00 is 9/24 = 0.375 days
    double target_mjd = 61106.375;
    utils::Instant epoch = utils::Instant::from_tt(utils::ModifiedJulianDate(target_mjd));

    std::cout << "--- Testing JPL Horizons for Asteroid 79518 ---" << std::endl;
    std::cout << "Target Time: 2026-03-07 09:00:00 UTC (MJD " << target_mjd << ")\n" << std::endl;

    // 1. Query Orbital Elements
    auto elements = horizons.query_elements("79518", epoch);
    if (elements) {
        std::cout << "[Orbital Elements (ECLIPJ2000)]" << std::endl;
        std::cout << "  a   = " << elements->a.to_au() << " AU" << std::endl;
        std::cout << "  e   = " << elements->e << std::endl;
        std::cout << "  i   = " << elements->i.to_deg() << " deg" << std::endl;
        std::cout << "  OM  = " << elements->node.to_deg() << " deg" << std::endl;
        std::cout << "  W   = " << elements->omega.to_deg() << " deg" << std::endl;
        std::cout << "  MA  = " << elements->M.to_deg() << " deg" << std::endl;
    } else {
        std::cerr << "Error fetching elements." << std::endl;
    }

    // 2. Query Astrometric Position
    auto obs = horizons.query_observation("79518", epoch, "500@399");
    if (obs) {
        std::cout << "\n[Astrometric Position (Geocentric)]" << std::endl;
        // RA in hours
        double ra_h = obs->ra.value * constants::RAD_TO_DEG / 15.0;
        double dec_d = obs->dec.value * constants::RAD_TO_DEG;
        
        std::cout << "  RA  = " << std::fixed << std::setprecision(6) << ra_h << " h" << std::endl;
        std::cout << "  DEC = " << std::fixed << std::setprecision(6) << dec_d << " deg" << std::endl;
        std::cout << "  Dist= " << obs->distance.value / (constants::AU * 1000.0) << " AU" << std::endl;
    } else {
        std::cerr << "Error fetching observation." << std::endl;
    }

    return 0;
}
