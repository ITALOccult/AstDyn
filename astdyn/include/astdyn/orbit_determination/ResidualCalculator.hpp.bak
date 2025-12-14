/**
 * @file ResidualCalculator.hpp
 * @brief Header for residual calculator
 * @date 2025-12-09
 */

#ifndef ASTDYN_RESIDUAL_CALCULATOR_HPP
#define ASTDYN_RESIDUAL_CALCULATOR_HPP

#include "astdyn/ephemeris/EphemerisProvider.hpp"
#include "astdyn/orbit_determination/CommonTypes.hpp"
#include <Eigen/Core>
#include <vector>
#include <map>
#include <memory>
#include <string>

namespace astdyn::orbit_determination {

struct ObservatoryData {
    std::string code;
    double longitude_deg; // East positive
    double rho_cos_phi;   // Parallax constant
    double rho_sin_phi;   // Parallax constant
};

class ResidualCalculator {
public:
    ResidualCalculator();

    /**
     * @brief Set ephemeris provider (needed for Earth position)
     */
    void set_ephemeris_provider(std::shared_ptr<ephemeris::EphemerisProvider> provider);
    
    /**
     * @brief Load MPC observatory file (optional, built-in common codes are available)
     */
    void load_observatories(const std::string& filename);
    
    /**
     * @brief Compute residual with full precision (Topocentric + Light-time)
     * 
     * @param obs The observation
     * @param state_at_obs The state [r, v] propagated to the observation time
     * @param epoch_mjd Epoch of the state (must match obs.epoch_mjd)
     */
    Residual compute_residual(
        const Observation& obs,
        const Eigen::Vector<double, 6>& state_at_obs, 
        double epoch_mjd
    );

private:
    std::shared_ptr<ephemeris::EphemerisProvider> ephemeris_provider_;
    std::map<std::string, ObservatoryData> observatories_;
    
    void initialize_default_observatories();

    // Helper: Compute Earth position from Ephemeris
    Eigen::Vector3d get_earth_position(double mjd) const;
    
    // Helper: Compute Observatory position in Geocentric Equatorial Frame (J2000)
    // Requires computing Local Sidereal Time
    Eigen::Vector3d get_observatory_position(const std::string& code, double mjd_utc) const;
    
    // Helper: Iterative light-time correction
    // Returns topocentric vector (Observer -> Asteroid) at emission time
    Eigen::Vector3d apply_light_time_correction(
        const Eigen::Vector3d& r_ast_obs, // Asteroid position at obs time
        const Eigen::Vector3d& v_ast_obs, // Asteroid velocity at obs time
        const Eigen::Vector3d& r_observer, // Observer position at obs time
        double& light_time_days
    ) const;
    
    // Helper: Greenwich Mean Sidereal Time
    double calc_gmst(double mjd_utc) const;

    void cartesian_to_radec(
        const Eigen::Vector3d& r,
        double& ra_deg,
        double& dec_deg
    ) const;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_RESIDUAL_CALCULATOR_HPP
