/**
 * @file StationWeights.cpp
 * @brief Per-observatory variance components from the fit residuals.
 */
#include "astdyn/orbit_determination/StationWeights.hpp"
#include <algorithm>
#include <cmath>

namespace astdyn::orbit_determination {

namespace {

/// Median of a chi-square distribution with one degree of freedom.
constexpr double CHI2_ONE_DOF_MEDIAN = 0.4549;

double median(std::vector<double> values) {
    const auto mid = values.begin() + values.size() / 2;
    std::nth_element(values.begin(), mid, values.end());
    return *mid;
}

} // namespace

std::map<std::string, double> station_sigma_factors(const std::vector<ObservationResidual>& residuals) {
    std::map<std::string, std::vector<double>> z_squared;
    for (const auto& r : residuals) {
        if (r.outlier) continue;
        z_squared[r.observatory_code].push_back(r.normalized_ra * r.normalized_ra);
        z_squared[r.observatory_code].push_back(r.normalized_dec * r.normalized_dec);
    }
    std::map<std::string, double> factors;
    for (const auto& [code, values] : z_squared) {
        if (values.size() < 2 * STATION_WEIGHT_MIN_OBSERVATIONS) continue;
        const double factor = std::sqrt(median(values) / CHI2_ONE_DOF_MEDIAN);
        factors[code] = std::clamp(factor, STATION_FACTOR_MIN, STATION_FACTOR_MAX);
    }
    return factors;
}

int apply_station_factors(std::vector<observations::OpticalObservation>& observations,
                          const std::map<std::string, double>& factors) {
    int rescaled = 0;
    for (auto& obs : observations) {
        if (obs.radec_covariance_arcsec2) continue;
        const auto it = factors.find(obs.observatory_code);
        if (it == factors.end()) continue;
        obs.sigma_ra = astrometry::Angle::from_rad(obs.sigma_ra.to_rad() * it->second);
        obs.sigma_dec = astrometry::Angle::from_rad(obs.sigma_dec.to_rad() * it->second);
        ++rescaled;
    }
    return rescaled;
}

} // namespace astdyn::orbit_determination
