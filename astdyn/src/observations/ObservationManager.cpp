#include "astdyn/observations/ObservationManager.hpp"
#include "astdyn/observations/MPCReader.hpp"
#include "astdyn/core/Constants.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>

namespace astdyn::observations {

int ObservationManager::load_mpc(const std::string& filename) {
    auto new_obs = MPCReader::readFile(filename);
    for (const auto& obs : new_obs) {
        observations_.push_back(obs);
    }
    return static_cast<int>(new_obs.size());
}

void ObservationManager::add(const OpticalObservation& obs) {
    observations_.push_back(obs);
}

void ObservationManager::sort_by_time() {
    std::sort(observations_.begin(), observations_.end(),
        [](const OpticalObservation& a, const OpticalObservation& b) {
            return a.time.mjd() < b.time.mjd();
        });
}

void ObservationManager::load_biases(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open()) return;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string code;
        double ra_bias, dec_bias;
        if (std::getline(ss, code, ',') && (ss >> ra_bias) && ss.ignore() && (ss >> dec_bias)) {
            catalog_biases_[code] = {astrometry::Angle::from_arcsec(ra_bias), 
                                     astrometry::Angle::from_arcsec(dec_bias)};
        }
    }
}

void ObservationManager::apply_biases() {
    for (auto& obs : observations_) {
        auto it = catalog_biases_.find(obs.observatory_code);
        if (it != catalog_biases_.end()) {
            obs.ra = astrometry::RightAscension(obs.ra - it->second.first);
            obs.dec = astrometry::Declination(obs.dec - it->second.second);
        }
    }
}

} // namespace astdyn::observations
