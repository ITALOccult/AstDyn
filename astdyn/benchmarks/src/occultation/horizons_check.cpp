/**
 * @file horizons_check.cpp
 * @brief Check with JPL Horizons what the actual expected coordinates are.
 */

#include "astdyn/io/HorizonsClient.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    io::HorizonsClient client;
    time::EpochTDB t = time::EpochTDB::from_mjd(61071.61);
    
    std::cout << "Querying Horizons for (704) at MJD 61071.61 (TDB)..." << std::endl;
    auto obs = client.query_observation("704", t, "500");
    
    if (obs) {
        std::cout << "Horizons RA:  " << std::fixed << std::setprecision(6) << (*obs).ra.value * constants::RAD_TO_DEG << " deg" << std::endl;
        std::cout << "Horizons Dec: " << (*obs).dec.value * constants::RAD_TO_DEG << " deg" << std::endl;
    } else {
        std::cout << "Query failed!" << std::endl;
    }

    return 0;
}
