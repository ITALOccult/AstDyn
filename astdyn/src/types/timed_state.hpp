#ifndef ASTDYN_TYPES_TIMED_STATE_HPP
#define ASTDYN_TYPES_TIMED_STATE_HPP

#include "src/types/orbital_state.hpp"
#include "src/utils/time_types.hpp"

namespace astdyn::types {

/**
 * @brief Couple of and OrbitalState with its Instant (MJD + Scale).
 * Semantic carrier for the propagation interface.
 */
template <typename FrameTag, typename RepresentationTag>
struct TimedState {
    const OrbitalState<FrameTag, RepresentationTag> state;
    const utils::Instant instant;

    explicit constexpr TimedState(const OrbitalState<FrameTag, RepresentationTag>& s, const utils::Instant& i) noexcept
        : state(s), instant(i) {}
};

} // namespace astdyn::types

#endif // ASTDYN_TYPES_TIMED_STATE_HPP
