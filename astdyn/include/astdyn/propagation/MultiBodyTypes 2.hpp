#ifndef ASTDYN_PROPAGATION_MULTIBODY_TYPES_HPP
#define ASTDYN_PROPAGATION_MULTIBODY_TYPES_HPP

#include "astdyn/core/physics_state.hpp"
#include <string>

namespace astdyn::propagation {

/**
 * @brief Represents a single body in a multi-body system
 */
struct MultiBodyState {
    std::string name;
    physics::GravitationalParameter gm;
    math::Vector3<core::ECLIPJ2000, physics::Distance> position;
    math::Vector3<core::ECLIPJ2000, physics::Velocity> velocity;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_MULTIBODY_TYPES_HPP
