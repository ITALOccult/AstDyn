/**
 * @file StationWeights.hpp
 * @brief Per-observatory weights estimated from the fit's own residuals.
 */
#ifndef ASTDYN_ORBIT_DETERMINATION_STATION_WEIGHTS_HPP
#define ASTDYN_ORBIT_DETERMINATION_STATION_WEIGHTS_HPP

#include "astdyn/observations/Observation.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include <map>
#include <string>
#include <vector>

namespace astdyn::orbit_determination {

/// Observatories with fewer residuals keep their weights.
inline constexpr int STATION_WEIGHT_MIN_OBSERVATIONS = 20;
inline constexpr double STATION_FACTOR_MIN = 0.5;
inline constexpr double STATION_FACTOR_MAX = 10.0;

/**
 * @brief Factor by which each observatory's sigmas should be scaled so that its
 *        normalised residuals have unit spread (variance-component estimation).
 *
 * This is what a station-by-station weighting table such as Veres et al. (2017)
 * encodes, measured here on the object's own data. The spread is robust:
 * median(z^2) / 0.4549, the median of a chi-square with one degree of freedom,
 * so a few outliers do not inflate it. Clipped to [0.5, 10].
 */
std::map<std::string, double> station_sigma_factors(const std::vector<ObservationResidual>& residuals);

/**
 * @brief Scale sigma_ra and sigma_dec by the observatory's factor. Observations
 *        carrying their own full covariance (Gaia) are left untouched.
 * @return number of observations rescaled
 */
int apply_station_factors(std::vector<observations::OpticalObservation>& observations,
                          const std::map<std::string, double>& factors);

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_STATION_WEIGHTS_HPP
