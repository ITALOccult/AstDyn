/**
 * @file fetch_horizons_vectors.cpp
 * @brief Scarica coordinate iniziali da JPL Horizons per il benchmark
 */

#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/core/physics_state.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace astdyn;

struct Target {
    std::string name;
    std::string id;
};

int main() {
    io::HorizonsClient client;
    auto epoch = time::EpochTDB::from_mjd(60310.0); // 2024-Jan-01 00:00:00 TDB
    
    std::vector<Target> targets = {
        {"Ceres", "1;"},
        {"Pallas", "2"},
        {"Baruffetti", "304719"},
        {"Icarus", "1566"}
    };

    std::cout << std::fixed << std::setprecision(16);
    std::cout << "--- DOWNLOAD COORDINATE JPL HORIZONS (MJD 60310.0) ---" << std::endl;
    std::cout << "Frame: ICRF/GCRF SSB (@ssb)" << std::endl << std::endl;

    for (const auto& t : targets) {
        std::cout << "Interrogazione " << t.name << " [" << t.id << "]..." << std::endl;
        
        auto result = client.query_vectors(t.id, epoch, "500@0");
        
        if (result) {
            auto state = *result;
            
            // Estraiamo i valori in AU e AU/d usando interop con Eigen e costanti IAU
            auto pos_eigen = state.position.to_eigen_si() / (constants::AU * 1000.0);
            auto vel_eigen = state.velocity.to_eigen_si() * 86400.0 / (constants::AU * 1000.0);

            std::cout << "    { " << std::endl;
            std::cout << "        \"" << t.name << "\", 60310.0," << std::endl;
            std::cout << "        " << pos_eigen.x() << ", " << pos_eigen.y() << ", " << pos_eigen.z() << "," << std::endl;
            std::cout << "        " << vel_eigen.x() << ", " << vel_eigen.y() << ", " << vel_eigen.z() << std::endl;
            std::cout << "    }," << std::endl;
        } else {
            std::cerr << "    ERRORE per " << t.name << ": Target non trovato o errore di rete." << std::endl;
        }
    }

    return 0;
}
