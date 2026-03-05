#ifndef ASTDYN_MATH_TRIGONOMETRY_HPP
#define ASTDYN_MATH_TRIGONOMETRY_HPP

#include "astdyn/core/Constants.hpp"
#include "src/core/units.hpp"
#include <cmath>

namespace astdyn::math {

using core::Radian;
namespace constants = astdyn::constants;

/**
 * @brief Normalizes an angle to the range [0, 2pi).
 * @param angle Any input angle in Radians.
 * @return Radian Normalized angle.
 */
[[nodiscard]] constexpr Radian wrap_to_2pi(const Radian angle) noexcept {
    const double raw = angle.value;
    const double tau = constants::TAU;
    
    double wrapped = std::fmod(raw, tau);
    
    if (wrapped < 0.0) {
        wrapped += tau;
    }
    
    return Radian(wrapped);
}

/**
 * @brief Normalizes an angle to the range [-pi, pi).
 * @param angle Any input angle in Radians.
 * @return Radian Normalized angle.
 */
[[nodiscard]] constexpr Radian wrap_to_pi(const Radian angle) noexcept {
    const Radian wrapped_2pi = wrap_to_2pi(angle);
    double raw = wrapped_2pi.value;
    const double pi = constants::PI;
    const double tau = constants::TAU;
    
    if (raw >= pi) {
        raw -= tau;
    }
    
    return Radian(raw);
}

} // namespace astdyn::math

#endif // ASTDYN_MATH_TRIGONOMETRY_HPP
