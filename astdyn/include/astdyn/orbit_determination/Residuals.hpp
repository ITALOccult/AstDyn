/**
 * @file Residuals.hpp
 * @brief Observation residual calculations (O-C)
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 * 
 * Computes residuals (Observed minus Computed) for various observation types.
 * This is the foundation for orbit determination via differential corrections.
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_RESIDUALS_HPP
#define ASTDYN_ORBIT_DETERMINATION_RESIDUALS_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <vector>
#include <optional>
#include <memory>

namespace astdyn {
    namespace propagation {
        class Propagator;
    }
}

namespace astdyn::orbit_determination {

/**
 * @brief Single observation residual
 */
struct ObservationResidual {
    double mjd_utc;                      ///< Observation epoch
    std::string observatory_code;        ///< Observatory
    
    // Angular residuals [rad]
    double residual_ra;                  ///< O-C in right ascension
    double residual_dec;                 ///< O-C in declination
    
    // Normalized residuals (dimensionless)
    double normalized_ra;                ///< (O-C) / sigma_ra
    double normalized_dec;               ///< (O-C) / sigma_dec
    
    // Computed values
    double computed_ra;                  ///< Computed RA [rad]
    double computed_dec;                 ///< Computed Dec [rad]
    
    // Geometry
    double range;                        ///< Topocentric distance [AU]
    double range_rate;                   ///< Topocentric range rate [AU/day]
    
    // Quality flags
    bool outlier;                        ///< Marked as outlier?
    double chi_squared;                  ///< χ² = (normalized_ra)² + (normalized_dec)²
    
    /**
     * @brief Check if residual is within tolerance
     */
    bool is_outlier(double sigma_threshold = 3.0) const {
        return (std::abs(normalized_ra) > sigma_threshold ||
                std::abs(normalized_dec) > sigma_threshold);
    }
};

/**
 * @brief Statistics for residual set
 */
struct ResidualStatistics {
    int num_observations;
    int num_outliers;
    int degrees_of_freedom;              ///< N_obs - N_params
    
    // RMS residuals [arcsec]
    double rms_ra;
    double rms_dec;
    double rms_total;
    
    // Weighted RMS (normalized)
    double weighted_rms;
    
    // Chi-squared test
    double chi_squared;
    double reduced_chi_squared;          ///< χ²/dof
    
    // Max residuals
    double max_abs_ra;
    double max_abs_dec;
};

/**
 * @brief Residual calculator for orbit determination
 * 
 * Computes O-C residuals for optical observations given an orbital state.
 * Handles light-time correction, aberration, and topocentric parallax.
 */
class ResidualCalculator {
public:
    /**
     * @brief Constructor
     * @param ephemeris Planetary ephemeris for Earth position
     * @param propagator Orbit propagator (optional, for automatic propagation)
     */
    explicit ResidualCalculator(
        std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris,
        std::shared_ptr<astdyn::propagation::Propagator> propagator = nullptr);
    
    /**
     * @brief Compute residuals for all observations
     * 
     * @param observations Optical observations
     * @param state Orbital state (Cartesian, heliocentric)
     * @return Vector of residuals
     */
    std::vector<ObservationResidual> compute_residuals(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const astdyn::propagation::CartesianElements& state) const;
    
    /**
     * @brief Compute residual for single observation
     * 
     * @param obs Single optical observation
     * @param state Orbital state at observation epoch
     * @return Residual, or nullopt if computation fails
     */
    std::optional<ObservationResidual> compute_residual(
        const astdyn::observations::OpticalObservation& obs,
        const astdyn::propagation::CartesianElements& state) const;
    
    /**
     * @brief Compute statistics from residuals
     * 
     * @param residuals Vector of observation residuals
     * @param num_parameters Number of fitted parameters (usually 6)
     * @return Statistical summary
     */
    static ResidualStatistics compute_statistics(
        const std::vector<ObservationResidual>& residuals,
        int num_parameters = 6);
    
    /**
     * @brief Identify and mark outliers
     * 
     * Uses iterative 3-sigma clipping.
     * 
     * @param residuals Residuals to check (modified in place)
     * @param sigma_threshold Threshold for outlier detection (default 3.0)
     * @return Number of outliers found
     */
    static int identify_outliers(
        std::vector<ObservationResidual>& residuals,
        double sigma_threshold = 3.0);
    
    /**
     * @brief Enable/disable light-time correction
     */
    void set_light_time_correction(bool enable) { 
        light_time_correction_ = enable; 
    }
    
    /**
     * @brief Enable/disable aberration correction
     */
    void set_aberration_correction(bool enable) { 
        aberration_correction_ = enable; 
    }
    
    /**
     * @brief Get observer position at observation time
     * 
     * @param obs Observation with time and observatory code
     * @return Observer position (heliocentric) [AU], or nullopt if failed
     */
    std::optional<astdyn::Vector3d> get_observer_position(
        const astdyn::observations::OpticalObservation& obs) const;
    
    /**
     * @brief Get observer velocity at observation time
     * 
     * @param obs Observation with time and observatory code
     * @return Observer velocity (heliocentric) [AU/day], or nullopt if failed
     */
    std::optional<astdyn::Vector3d> get_observer_velocity(
        const astdyn::observations::OpticalObservation& obs) const;
    
    /**
     * @brief Convert UTC to TDB time scale
     * 
     * @param mjd_utc Modified Julian Date in UTC
     * @return MJD in TDB time scale
     */
    static double utc_to_tdb(double mjd_utc);

private:
    /**
     * @brief Compute topocentric position of object
     * 
     * @param heliocentric_pos Object position (heliocentric) [AU]
     * @param observer_pos Observer position (heliocentric) [AU]
     * @param observer_vel Observer velocity (heliocentric) [AU/day]
     * @param[out] range Topocentric distance [AU]
     * @param[out] range_rate Topocentric range rate [AU/day]
     * @return Unit vector from observer to object
     */
    astdyn::Vector3d compute_topocentric_direction(
        const astdyn::Vector3d& heliocentric_pos,
        const astdyn::Vector3d& observer_pos,
        const astdyn::Vector3d& observer_vel,
        double& range,
        double& range_rate) const;
    
    /**
     * @brief Convert topocentric Cartesian to RA/Dec
     * 
     * @param direction Unit vector (topocentric)
     * @param[out] ra Right ascension [rad]
     * @param[out] dec Declination [rad]
     */
    void cartesian_to_radec(
        const astdyn::Vector3d& direction,
        const astdyn::Vector3d& rho_vec,
        const astdyn::Vector3d& observer_pos,
        double& ra,
        double& dec) const;

private:
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<astdyn::propagation::Propagator> propagator_;
    
    // Correction flags
    bool light_time_correction_ = true;
    bool aberration_correction_ = true;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_RESIDUALS_HPP
