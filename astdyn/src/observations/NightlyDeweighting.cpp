/**
 * @file NightlyDeweighting.cpp
 * @brief Veres et al. (2017) batch deweighting of same-night observations.
 */
#include "astdyn/observations/NightlyDeweighting.hpp"
#include <cmath>
#include <map>
#include <string>
#include <utility>

namespace astdyn::observations {

namespace {

/// One observatory, one UTC day: the batch whose errors are correlated.
using NightKey = std::pair<std::string, long>;

NightKey night_of(const OpticalObservation& obs) {
    return {obs.observatory_code, static_cast<long>(std::floor(obs.time.mjd()))};
}

} // namespace

int deweight_nightly_batches(std::vector<OpticalObservation>& observations) {
    std::map<NightKey, int> count_per_night;
    for (const auto& obs : observations) ++count_per_night[night_of(obs)];

    int inflated_count = 0;
    for (auto& obs : observations) {
        const int count = count_per_night[night_of(obs)];
        if (count <= OBSERVATIONS_PER_NIGHT_UNWEIGHTED) continue;
        const double factor = std::sqrt(static_cast<double>(count) / OBSERVATIONS_PER_NIGHT_UNWEIGHTED);
        obs.sigma_ra = astrometry::Angle::from_rad(obs.sigma_ra.to_rad() * factor);
        obs.sigma_dec = astrometry::Angle::from_rad(obs.sigma_dec.to_rad() * factor);
        ++inflated_count;
    }
    return inflated_count;
}

} // namespace astdyn::observations
