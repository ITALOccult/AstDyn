#include <astdyn/AstDyn.hpp>
#include <iostream>
#include <iomanip>

int main() {
    using namespace astdyn;
    using namespace astdyn::io;

    // 1. Initialize Client
    HorizonsClient horizons;

    // Epoch: Today (approx)
    time::EpochTDB epoch = time::EpochTDB::from_mjd(60310.5); // 2024-01-01

    std::cout << "--- Querying JPL Horizons for Ceres (1) ---" << std::endl;

    // 2. Query Orbital Elements
    auto elements = horizons.query_elements("1", epoch);
    if (elements) {
        std::cout << "\n[Elements]" << std::endl;
        std::cout << "  a  = " << elements->a.to_au() << " AU" << std::endl;
        std::cout << "  e  = " << elements->e << std::endl;
        std::cout << "  i  = " << elements->i.to_deg() << " deg" << std::endl;
    } else {
        std::cerr << "Failed to query elements." << std::endl;
    }

    // 3. Query State Vectors (Barycentric)
    auto vectors = horizons.query_vectors("1", epoch, "500@0");
    if (vectors) {
        std::cout << "\n[Vectors (GCRF/SSB)]" << std::endl;
        std::cout << "  X  = " << vectors->position.x_si() / (constants::AU * 1000.0) << " AU" << std::endl;
        std::cout << "  Y  = " << vectors->position.y_si() / (constants::AU * 1000.0) << " AU" << std::endl;
        std::cout << "  Z  = " << vectors->position.z_si() / (constants::AU * 1000.0) << " AU" << std::endl;
    }

    // 4. Query Astrometric Observation (from Earth)
    auto obs = horizons.query_observation("1", epoch, "500@399");
    if (obs) {
        std::cout << "\n[Geocentric Observation]" << std::endl;
        std::cout << "  RA  = " << obs->ra.value * constants::RAD_TO_DEG / 15.0 << " h" << std::endl;
        std::cout << "  DEC = " << obs->dec.value * constants::RAD_TO_DEG << " deg" << std::endl;
        std::cout << "  Dist= " << obs->distance.value / constants::AU / 1000.0 << " AU" << std::endl;
    }

    return 0;
}
