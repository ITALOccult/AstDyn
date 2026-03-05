#ifndef ASTDYN_COORDINATES_ROTATION_MATRICES_HPP
#define ASTDYN_COORDINATES_ROTATION_MATRICES_HPP

#include "src/core/units.hpp"
#include "src/types/matrices.hpp"
#include <cmath>

namespace astdyn::coordinates {

using core::Radian;
using types::Matrix3x3;

/** @brief Rotation matrix around X axis (R1). */
[[nodiscard]] inline Matrix3x3 rotation_x(const Radian angle) noexcept {
    const double cos_a = std::cos(angle.value);
    const double sin_a = std::sin(angle.value);
    return Matrix3x3({1.0, 0.0, 0.0, 
                      0.0, cos_a, sin_a, 
                      0.0, -sin_a, cos_a});
}

/** @brief Rotation matrix around Y axis (R2). */
[[nodiscard]] inline Matrix3x3 rotation_y(const Radian angle) noexcept {
    const double cos_a = std::cos(angle.value);
    const double sin_a = std::sin(angle.value);
    return Matrix3x3({cos_a, 0.0, -sin_a, 
                      0.0, 1.0, 0.0, 
                      sin_a, 0.0, cos_a});
}

/** @brief Rotation matrix around Z axis (R3). */
[[nodiscard]] inline Matrix3x3 rotation_z(const Radian angle) noexcept {
    const double cos_a = std::cos(angle.value);
    const double sin_a = std::sin(angle.value);
    return Matrix3x3({cos_a, sin_a, 0.0, 
                      -sin_a, cos_a, 0.0, 
                      0.0, 0.0, 1.0});
}

} // namespace astdyn::coordinates

#endif // ASTDYN_COORDINATES_ROTATION_MATRICES_HPP
