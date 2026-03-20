/**
 * @file ResidualAnalysis.hpp
 * @brief Analysis of orbit estimation quality (Observation minus Calculation).
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_RESIDUAL_ANALYSIS_HPP
#define ASTDYN_ORBIT_DETERMINATION_RESIDUAL_ANALYSIS_HPP

#include "astdyn/core/physics_state.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include <vector>
#include <string>

namespace astdyn::orbit_determination {

/** @brief Detailed information for a single measurement residual. */
struct ResidualReport {
    time::EpochUTC time;
    double ra_obs, ra_calc, ra_resid;   // arcsec
    double dec_obs, dec_calc, dec_resid; // arcsec
    double total_resid;                // arcsec
    double pos_error_km = -1.0;         // If truth available
};

/**
 * @brief Summary of the orbit validation against observations.
 */
struct OrbitValidationSummary {
    double rms_ra, rms_dec, rms_total; // arcsec
    double mean_pos_error_km;          // If truth available
    int num_observations;
    std::string report_text;
};

/**
 * @brief Utility class for orbit quality analysis.
 */
class ResidualAnalysis {
public:
    /**
     * @brief Generate a full report of residuals.
     * 
     * @param initial_state Initial state guess.
     * @param obs Selection of observations.
     * @param propagator Propagator to use for predictions.
     * @return Orbit validation summary.
     */
    static OrbitValidationSummary analyze_orbit(
        const physics::CartesianStateTyped<core::GCRF>& initial_state,
        const std::vector<observations::OpticalObservation>& obs,
        std::shared_ptr<propagation::Propagator> propagator);

};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_RESIDUAL_ANALYSIS_HPP
