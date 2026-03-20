/**
 * @file CommonTypes.hpp
 * @brief Common types for orbit determination
 * @date 2025-12-09
 */

#ifndef ASTDYN_OD_COMMON_TYPES_HPP
#define ASTDYN_OD_COMMON_TYPES_HPP

#include "astdyn/time/epoch.hpp"
#include <string>

namespace astdyn::orbit_determination {

struct Observation {
    time::EpochUTC time;
    double ra_deg;
    double dec_deg;
    double ra_sigma_arcsec;
    double dec_sigma_arcsec;
    std::string observatory_code;
    double weight;
    bool rejected;
};

struct Residual {
    time::EpochUTC time;
    double ra_computed_deg;
    double dec_computed_deg;
    double ra_residual_arcsec;
    double dec_residual_arcsec;
    bool rejected;
};

struct ObservationResidual : public Residual {
    double ra_obs_deg;
    double dec_obs_deg;
    double weight_ra;
    double weight_dec;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_OD_COMMON_TYPES_HPP
