/**
 * @file GaussIOD.hpp
 * @brief Gauss method for Initial Orbit Determination
 * @author ITALOccult AstDyn Team
 * @date 2025-11-25
 * 
 * Implements the classical Gauss method for determining a preliminary orbit
 * from three optical observations (RA/Dec).
 * 
 * Algorithm:
 * 1. Select three well-separated observations
 * 2. Compute slant ranges using 8th order polynomial (Gauss-Gibbs)
 * 3. Calculate geocentric state vectors
 * 4. Convert to heliocentric and derive orbital elements
 * 
 * Reference: Bate, Mueller & White - "Fundamentals of Astrodynamics" (1971)
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_GAUSS_IOD_HPP
#define ASTDYN_ORBIT_DETERMINATION_GAUSS_IOD_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <vector>
#include <optional>

namespace astdyn::orbit_determination {

/**
 * @brief Settings for Gauss IOD
 */
struct GaussIODSettings {
    int max_iterations = 50;             ///< Maximum iterations for slant range
    double tolerance = 1e-8;             ///< Convergence tolerance [AU]
    double min_separation_days = 1.0;    ///< Minimum separation between obs [days]
    double max_separation_days = 60.0;   ///< Maximum separation between obs [days]
    bool use_light_time = true;          ///< Apply light-time correction
    bool verbose = false;                ///< Print debug information
};

/**
 * @brief Result from Gauss IOD
 */
struct GaussIODResult {
    bool success;
    std::string error_message;
    
    // Resulting orbit
    astdyn::propagation::CartesianElements state;   ///< Heliocentric state at middle obs
    double epoch_mjd_tdb;                           ///< Epoch of solution [MJD TDB]
    
    // Quality indicators
    double slant_range_1;                ///< Distance to object [AU] at obs 1
    double slant_range_2;                ///< Distance to object [AU] at obs 2
    double slant_range_3;                ///< Distance to object [AU] at obs 3
    int iterations;                      ///< Iterations required
    
    // Indices of observations used
    int obs_index_1;
    int obs_index_2;
    int obs_index_3;
    
    /**
     * @brief Print summary
     */
    void print_summary() const;
};

/**
 * @brief Gauss method for Initial Orbit Determination
 * 
 * Determines a preliminary orbit from three optical observations using
 * the classical Gauss method with 8th order polynomial refinement.
 */
class GaussIOD {
public:
    /**
     * @brief Constructor
     */
    explicit GaussIOD(const GaussIODSettings& settings = GaussIODSettings());
    
    /**
     * @brief Determine orbit from observations
     * 
     * Automatically selects three well-separated observations and computes
     * a preliminary orbit.
     * 
     * @param observations Vector of optical observations (requires ≥3)
     * @return IOD result
     */
    GaussIODResult compute(const std::vector<observations::OpticalObservation>& observations);
    
    /**
     * @brief Determine orbit from three specific observations
     * 
     * @param obs1 First observation
     * @param obs2 Second observation (middle, epoch of solution)
     * @param obs3 Third observation
     * @return IOD result
     */
    GaussIODResult compute_from_three(
        const observations::OpticalObservation& obs1,
        const observations::OpticalObservation& obs2,
        const observations::OpticalObservation& obs3);
    
    /**
     * @brief Get/set settings
     */
    const GaussIODSettings& settings() const { return settings_; }
    void set_settings(const GaussIODSettings& s) { settings_ = s; }

private:
    GaussIODSettings settings_;
    
    /**
     * @brief Select three optimal observations
     * 
     * Chooses observations with good time separation and geometry.
     * 
     * @param observations All available observations
     * @return Indices of selected observations (or empty if selection fails)
     */
    std::optional<std::array<int, 3>> select_observations(
        const std::vector<observations::OpticalObservation>& observations);
    
    /**
     * @brief Compute line-of-sight unit vector
     * 
     * @param ra Right ascension [rad]
     * @param dec Declination [rad]
     * @return Unit vector in equatorial frame
     */
    Vector3d compute_line_of_sight(double ra, double dec) const;
    
    /**
     * @brief Solve Gauss polynomial for slant ranges
     * 
     * Uses 8th order polynomial to iteratively refine slant ranges.
     * 
     * @param tau1 Time interval t2-t1 [days]
     * @param tau3 Time interval t3-t2 [days]
     * @param los1 Line of sight vector at obs 1
     * @param los2 Line of sight vector at obs 2
     * @param los3 Line of sight vector at obs 3
     * @param R1 Earth position at obs 1 [AU]
     * @param R2 Earth position at obs 2 [AU]
     * @param R3 Earth position at obs 3 [AU]
     * @param[out] rho1 Slant range at obs 1 [AU]
     * @param[out] rho2 Slant range at obs 2 [AU]
     * @param[out] rho3 Slant range at obs 3 [AU]
     * @param[out] iterations Number of iterations
     * @return true if converged
     */
    bool solve_slant_ranges(
        double tau1, double tau3,
        const Vector3d& los1, const Vector3d& los2, const Vector3d& los3,
        const Vector3d& R1, const Vector3d& R2, const Vector3d& R3,
        double& rho1, double& rho2, double& rho3,
        int& iterations);
    
    /**
     * @brief Compute Lagrange coefficients
     * 
     * f and g functions for two-body propagation.
     * 
     * @param r Position vector [AU]
     * @param v Velocity vector [AU/day]
     * @param dt Time interval [days]
     * @param mu Gravitational parameter [AU³/day²]
     * @return (f, g) coefficients
     */
    std::pair<double, double> compute_f_g_coefficients(
        const Vector3d& r, const Vector3d& v, double dt, double mu) const;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_GAUSS_IOD_HPP
