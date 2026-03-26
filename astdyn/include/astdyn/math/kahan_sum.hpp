/**
 * @file kahan_sum.hpp
 * @brief Kahan compensated summation for Eigen::Vector3d accumulators.
 *
 * Implements Zeebe 2023 (AJ 166:1) option (1): compensation applied after
 * each term additively, not accumulated in separate variables.
 *
 * Usage:
 *   Eigen::Vector3d acc  = Eigen::Vector3d::Zero();
 *   Eigen::Vector3d comp = Eigen::Vector3d::Zero();
 *   for (const auto& term : terms)
 *       astdyn::math::kahan_add(acc, comp, term);
 *
 * Reduces floating-point summation error from O(n * eps) to O(eps).
 */

#ifndef ASTDYN_MATH_KAHAN_SUM_HPP
#define ASTDYN_MATH_KAHAN_SUM_HPP

#include <Eigen/Core>

namespace astdyn::math {

/// Kahan compensated vector addition (Zeebe 2023 option 1).
/// Call once per term added to acc. comp must persist across calls.
inline void kahan_add(
    Eigen::Vector3d&       acc,
    Eigen::Vector3d&       compensation,
    const Eigen::Vector3d& term) noexcept
{
    const Eigen::Vector3d y = term - compensation;
    const Eigen::Vector3d t = acc + y;
    compensation = (t - acc) - y;
    acc = t;
}

} // namespace astdyn::math

#endif // ASTDYN_MATH_KAHAN_SUM_HPP
