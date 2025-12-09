/**
 * @file ResidualCalculator.hpp
 * @brief Calculate astrometric residuals (O-C) for orbit determination
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * Computes residuals between observed and computed positions:
 * - Propagates orbit to observation epoch
 * - Converts to topocentric coordinates
 * - Applies light-time, aberration corrections
 * - Computes RA/Dec
 * - Returns O-C residuals
 */

#ifndef ASTDYN_RESIDUAL_CALCULATOR_HPP
#define ASTDYN_RESIDUAL_CALCULATOR_HPP

#include <Eigen/Dense>
#include <vector>
#include <string>

namespace astdyn::orbit_determination {

/**
 * @brief Astrometric observation
 */
struct Observation {
    double epoch_mjd;           ///< Modified Julian Date (UTC)
    double ra_deg;              ///< Right Ascension [degrees]
    double dec_deg;             ///< Declination [degrees]
    double ra_sigma_arcsec;     ///< RA uncertainty [arcsec]
    double dec_sigma_arcsec;    ///< Dec uncertainty [arcsec]
    std::string observatory_code; ///< MPC observatory code
    double weight;              ///< Observation weight (1/σ²)
    bool rejected;              ///< Outlier flag
};

/**
 * @brief Observatory position
 */
struct Observatory {
    std::string code;           ///< MPC code
    Eigen::Vector3d position;   ///< Geocentric position [AU]
    double longitude_deg;       ///< Longitude [degrees]
    double latitude_deg;        ///< Latitude [degrees]
    double altitude_m;          ///< Altitude [meters]
};

/**
 * @brief Residual (O-C)
 */
struct Residual {
    double epoch_mjd;
    double ra_residual_arcsec;  ///< RA*(O-C) [arcsec]
    double dec_residual_arcsec; ///< Dec(O-C) [arcsec]
    double ra_computed_deg;     ///< Computed RA [degrees]
    double dec_computed_deg;    ///< Computed Dec [degrees]
    bool rejected;
};

/**
 * @brief Calculate astrometric residuals
 */
class ResidualCalculator {
public:
    /**
     * @brief Construct residual calculator
     */
    ResidualCalculator();
    
    /**
     * @brief Compute residuals for all observations
     * 
     * @param observations List of observations
     * @param state0 Initial state [r, v] at epoch0
     * @param epoch0_mjd Epoch of initial state (MJD)
     * @return Vector of residuals
     */
    std::vector<Residual> compute_residuals(
        const std::vector<Observation>& observations,
        const Eigen::Vector<double, 6>& state0,
        double epoch0_mjd
    );
    
    /**
     * @brief Compute single residual
     */
    Residual compute_residual(
        const Observation& obs,
        const Eigen::Vector<double, 6>& state0,
        double epoch0_mjd
    );
    
    /**
     * @brief Set planetary ephemeris function
     * 
     * Function signature: (body_id, mjd) -> position [AU]
     */
    using EphemerisFunction = std::function<Eigen::Vector3d(int, double)>;
    void set_ephemeris_function(EphemerisFunction func) {
        ephemeris_func_ = func;
    }
    
    /**
     * @brief Load observatory database
     */
    void load_observatories(const std::string& filename);
    
    /**
     * @brief Get observatory by MPC code
     */
    Observatory get_observatory(const std::string& code) const;
    
private:
    EphemerisFunction ephemeris_func_;
    std::map<std::string, Observatory> observatories_;
    
    /**
     * @brief Convert Cartesian to RA/Dec
     * 
     * @param r Position vector [AU]
     * @param ra_deg Output RA [degrees]
     * @param dec_deg Output Dec [degrees]
     */
    void cartesian_to_radec(
        const Eigen::Vector3d& r,
        double& ra_deg,
        double& dec_deg
    ) const;
    
    /**
     * @brief Apply light-time correction
     * 
     * Iteratively solve for light-time delay
     */
    Eigen::Vector3d apply_light_time(
        const Eigen::Vector3d& r_helio,
        const Eigen::Vector3d& r_obs,
        double& light_time_days
    ) const;
    
    /**
     * @brief Get Earth position at epoch
     */
    Eigen::Vector3d get_earth_position(double mjd) const;
    
    /**
     * @brief Get observatory position (topocentric)
     */
    Eigen::Vector3d get_observatory_position(
        const std::string& code,
        double mjd
    ) const;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_RESIDUAL_CALCULATOR_HPP
