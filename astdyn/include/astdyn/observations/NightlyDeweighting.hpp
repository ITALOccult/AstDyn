/**
 * @file NightlyDeweighting.hpp
 * @brief Weights of observations taken in batches by one observatory in one night.
 */
#ifndef ASTDYN_OBSERVATIONS_NIGHTLY_DEWEIGHTING_HPP
#define ASTDYN_OBSERVATIONS_NIGHTLY_DEWEIGHTING_HPP

#include "astdyn/observations/Observation.hpp"
#include <vector>

namespace astdyn::observations {

/// Observations above which a night's batch from one observatory is deweighted.
inline constexpr int OBSERVATIONS_PER_NIGHT_UNWEIGHTED = 4;

/**
 * @brief Inflate sigma by sqrt(N/4) for the N > 4 observations an observatory
 *        made in the same night (Veres et al. 2017, as used by JPL and the MPC).
 *
 * Observations of one night from one instrument share their systematic errors,
 * so they do not carry N independent pieces of information. Without this,
 * (3830)'s 2425 TESS positions -- 20 nights, up to 144 each, sigma 0.6" --
 * outweighed every ground observation together and pulled the orbit 180-265 km
 * (on the shadow) away from JPL's; without TESS the gap was 45-68 km.
 *
 * @return number of observations whose sigma was inflated
 */
int deweight_nightly_batches(std::vector<OpticalObservation>& observations);

} // namespace astdyn::observations

#endif // ASTDYN_OBSERVATIONS_NIGHTLY_DEWEIGHTING_HPP
