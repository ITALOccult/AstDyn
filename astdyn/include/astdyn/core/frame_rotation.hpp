#ifndef ASTDYN_CORE_FRAME_ROTATION_HPP
#define ASTDYN_CORE_FRAME_ROTATION_HPP

#include "frame_tags.hpp"
#include "orbital_state.hpp"
#include <concepts>

namespace astdyn::core {

/**
 * @brief Trait system to define allowed frame rotations.
 * By default, no rotation is defined.
 */
template <typename From, typename To>
struct RotationMatrixProvider {
    static constexpr bool exists = false;
};

/**
 * @brief Concept to verify if a rotation between two frames is mathematically defined.
 */
template <typename From, typename To>
concept IsRotatable = RotationMatrixProvider<From, To>::exists;

/**
 * @brief Rotates an orbital state into a target reference frame.
 * Compiles ONLY if a valid RotationMatrixProvider specialization exists for the frames.
 * 
 * @tparam TargetFrame The desired destination frame tag.
 * @tparam Frame Initial frame tag.
 * @tparam Rep State representation tag.
 * @param state The immutable orbital state.
 * @return A new OrbitalState in the target frame.
 */
template <typename TargetFrame, typename Frame, typename Rep>
requires IsRotatable<Frame, TargetFrame>
[[nodiscard]] constexpr OrbitalState<TargetFrame, Rep> rotate(const OrbitalState<Frame, Rep>& state) {
    // Transformation logic would use RotationMatrixProvider<Frame, TargetFrame>::get_matrix()
    // For now, it enforces the type safety and template constraints required.
    const auto transformed_data = RotationMatrixProvider<Frame, TargetFrame>::apply(state.raw_data());
    return OrbitalState<TargetFrame, Rep>(transformed_data);
}

} // namespace astdyn::core

#endif // ASTDYN_CORE_FRAME_ROTATION_HPP
