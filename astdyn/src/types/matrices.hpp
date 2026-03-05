#ifndef ASTDYN_TYPES_MATRICES_HPP
#define ASTDYN_TYPES_MATRICES_HPP

#include <array>
#include "vectors.hpp"

namespace astdyn::types {

/**
 * @brief Strong-typed 3x3 Matrix for orbital transformations.
 * Implicitly assumes consistent units.
 * 
 * CTFYH adherence: Immutable, minimalist.
 */
struct Matrix3x3 {
    const std::array<double, 9> elements;

    explicit constexpr Matrix3x3(const std::array<double, 9>& el) noexcept : elements(el) {}

    /** @brief Matrix-Vector multiplication. */
    template <typename Frame, typename Unit>
    [[nodiscard]] constexpr Vector3<Frame, Unit> multiply(const Vector3<Frame, Unit>& vec) const noexcept {
        const double x = elements[0] * vec.x + elements[1] * vec.y + elements[2] * vec.z;
        const double y = elements[3] * vec.x + elements[4] * vec.y + elements[5] * vec.z;
        const double z = elements[6] * vec.x + elements[7] * vec.y + elements[8] * vec.z;
        return Vector3<Frame, Unit>(x, y, z);
    }

    /** @brief Matrix-Matrix multiplication. */
    [[nodiscard]] constexpr Matrix3x3 multiply(const Matrix3x3& other) const noexcept {
        return Matrix3x3({
            elements[0]*other.elements[0] + elements[1]*other.elements[3] + elements[2]*other.elements[6],
            elements[0]*other.elements[1] + elements[1]*other.elements[4] + elements[2]*other.elements[7],
            elements[0]*other.elements[2] + elements[1]*other.elements[5] + elements[2]*other.elements[8],

            elements[3]*other.elements[0] + elements[4]*other.elements[3] + elements[5]*other.elements[6],
            elements[3]*other.elements[1] + elements[4]*other.elements[4] + elements[5]*other.elements[7],
            elements[3]*other.elements[2] + elements[4]*other.elements[5] + elements[5]*other.elements[8],

            elements[6]*other.elements[0] + elements[7]*other.elements[3] + elements[8]*other.elements[6],
            elements[6]*other.elements[1] + elements[7]*other.elements[4] + elements[8]*other.elements[7],
            elements[6]*other.elements[2] + elements[7]*other.elements[5] + elements[8]*other.elements[8]
        });
    }
};

} // namespace astdyn::types

#endif // ASTDYN_TYPES_MATRICES_HPP
