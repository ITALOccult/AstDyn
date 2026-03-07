/**
 * @file GoodingIOD.hpp
 * @brief Gooding's method for Initial Orbit Determination from 3 RA/Dec observations.
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_GOODING_IOD_HPP
#define ASTDYN_ORBIT_DETERMINATION_GOODING_IOD_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/observations/Observation.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "src/core/frame_tags.hpp"
#include <vector>
#include <optional>

namespace astdyn::orbit_determination {

/** @brief Results from Gooding IOD. Gooding can yield multiple candidate orbits. */
struct GoodingIODResult {
    bool success = false;
    std::string error_message;
    
    struct Solution {
        physics::CartesianStateTyped<core::GCRF> state; ///< Heliocentric state at t1
        time::EpochTDB epoch;                           ///< Epoch t1
        double rho1, rho2, rho3;                         ///< Slant ranges [AU]
        double rms_error;                               ///< Post-fit RMS on obs 2 [arcsec]
    };
    
    std::vector<Solution> solutions;
};

/**
 * @brief Gooding's method for Initial Orbit Determination.
 * 
 * Iterative method for determining an orbit from 3 RA/Dec observations.
 * Generally more robust than Gauss for long-arc observations.
 * 
 * Reference: R.H. Gooding, "A New Procedure for Orbit Determination Based on Three Lines of Sight" (1990)
 */
class GoodingIOD {
public:
    struct Settings {
        int max_iterations = 80;
        double tolerance = 1e-12; // in AU
        double mu = constants::GM_SUN_AU; // AU^3/day^2
        bool verbose = false;
    };

    explicit GoodingIOD(const Settings& settings = Settings());

    /**
     * @brief Compute orbit using Gooding's algorithm.
     * 
     * @param obs1 First observation
     * @param obs2 Second observation
     * @param obs3 Third observation
     * @param rho1_guess Initial guess for slant range 1 [AU]
     * @param rho3_guess Initial guess for slant range 3 [AU]
     * @return Result with potential solutions.
     */
    GoodingIODResult compute(
        const observations::OpticalObservation& obs1,
        const observations::OpticalObservation& obs2,
        const observations::OpticalObservation& obs3,
        double rho1_guess,
        double rho3_guess);

private:
    Settings settings_;

    /** @brief Newton iteration core. */
    bool solve_iteration(
        const time::EpochTDB& t1, const time::EpochTDB& t2, const time::EpochTDB& t3,
        const math::Vector3<core::GCRF, physics::Distance>& L1,
        const math::Vector3<core::GCRF, physics::Distance>& L2,
        const math::Vector3<core::GCRF, physics::Distance>& L3,
        const math::Vector3<core::GCRF, physics::Distance>& R1,
        const math::Vector3<core::GCRF, physics::Distance>& R2,
        const math::Vector3<core::GCRF, physics::Distance>& R3,
        double& rho1, double& rho3,
        physics::CartesianStateTyped<core::GCRF>& final_state);
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_GOODING_IOD_HPP
