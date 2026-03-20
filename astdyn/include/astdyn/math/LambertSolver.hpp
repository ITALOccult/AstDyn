/**
 * @file LambertSolver.hpp
 * @brief Solving Lambert's problem (Boundary Value Problem) in astrodynamics.
 */

#ifndef ASTDYN_MATH_LAMBERT_SOLVER_HPP
#define ASTDYN_MATH_LAMBERT_SOLVER_HPP

#include <Eigen/Dense>
#include "astdyn/core/Constants.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/core/frame_tags.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/math/frame_algebra.hpp"

namespace astdyn::math {

/**
 * @brief Solver for Lambert's problem (2-point boundary value problem).
 * 
 * Given two positions r1, r2 and a time interval dt, finds the velocity v1.
 */
class LambertSolver {
public:
    /**
     * @brief Universal Variables Lambert solver.
     * 
     * @param r1_vec Initial position vector [AU]
     * @param r2_vec Final position vector [AU]
     * @param dt Time interval [days]
     * @param mu Gravitational parameter [AU^3/day^2]
     * @param retrograde Whether the motion is retrograde
     * @return Velocity vector v1 [AU/day]
     */
    static math::Vector3<core::GCRF, physics::Velocity> solve(
        const math::Vector3<core::GCRF, physics::Distance>& r1_vec,
        const math::Vector3<core::GCRF, physics::Distance>& r2_vec,
        time::TimeDuration dt,
        physics::GravitationalParameter mu = physics::GravitationalParameter::sun(),
        bool retrograde = false);

private:
    /** @brief Universal Variable C(z) function. */
    static double C(double z);
    /** @brief Universal Variable S(z) function. */
    static double S(double z);
};

} // namespace astdyn::math

#endif // ASTDYN_MATH_LAMBERT_SOLVER_HPP
