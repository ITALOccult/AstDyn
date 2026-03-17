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
#include "astdyn/core/physics_state.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "astdyn/astrometry/sky_types.hpp"
#include "src/core/frame_tags.hpp"
#include <vector>
#include <optional>

namespace astdyn::orbit_determination {

/**
 * @brief Settings for Gauss IOD
 */
struct GaussIODSettings {
    int max_iterations = 50;
    physics::Distance tolerance = physics::Distance::from_au(1e-8);
    time::TimeDuration min_separation = time::TimeDuration::from_days(1.0);
    time::TimeDuration max_separation = time::TimeDuration::from_days(60.0);
    physics::GravitationalParameter mu = physics::GravitationalParameter::sun();
    bool use_light_time = true;
    bool verbose = false;
};

/**
 * @brief Result from Gauss IOD
 */
struct GaussIODResult {
    bool success;
    std::string error_message;
    
    // Resulting orbit
    physics::CartesianStateTyped<core::GCRF> state;   ///< Heliocentric state at middle obs
    time::EpochTDB epoch;                           ///< Epoch of solution
    
    // Quality indicators
    physics::Distance slant_range_1;     ///< Distance to object at obs 1
    physics::Distance slant_range_2;     ///< Distance to object at obs 2
    physics::Distance slant_range_3;     ///< Distance to object at obs 3
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
    explicit GaussIOD(std::shared_ptr<ephemeris::PlanetaryEphemeris> ephem = nullptr, 
                      const GaussIODSettings& settings = GaussIODSettings());
    
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
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    
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
    Eigen::Vector3d compute_line_of_sight(astrometry::RightAscension ra, astrometry::Declination dec) const;
    
    /**
     * @brief Solve Gauss polynomial for slant ranges
     * 
     * Uses 8th order polynomial to iteratively refine slant ranges.
     * 
     * @param tau1 Time interval t2-t1 [days]
     * @param tau3 Time interval t3-t2 [days]
     * @param l1 Line of sight vector at obs 1
     * @param l2 Line of sight vector at obs 2
     * @param l3 Line of sight vector at obs 3
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
        time::TimeDuration tau1, time::TimeDuration tau3,
        const Eigen::Vector3d& l1, 
        const Eigen::Vector3d& l2, 
        const Eigen::Vector3d& l3,
        const math::Vector3<core::GCRF, physics::Distance>& R1, 
        const math::Vector3<core::GCRF, physics::Distance>& R2, 
        const math::Vector3<core::GCRF, physics::Distance>& R3,
        physics::Distance& rho1, physics::Distance& rho2, physics::Distance& rho3,
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
        const math::Vector3<core::GCRF, physics::Distance>& r, 
        const math::Vector3<core::GCRF, physics::Velocity>& v, 
        double dt, double mu) const;
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_GAUSS_IOD_HPP
