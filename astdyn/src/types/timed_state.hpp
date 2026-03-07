#ifndef ASTDYN_TYPES_TIMED_STATE_HPP
#define ASTDYN_TYPES_TIMED_STATE_HPP

#include "src/types/orbital_state.hpp"
#include "astdyn/time/epoch.hpp"

namespace astdyn::types {

/**
 * @brief Couple of an OrbitalState with its Epoch (time::EpochTDB).
 * Semantic carrier for the propagation interface.
 */
template <typename FrameTag, typename RepresentationTag>
struct TimedState {
    const OrbitalState<FrameTag, RepresentationTag> state;
    const time::EpochTDB instant;

    explicit constexpr TimedState(const OrbitalState<FrameTag, RepresentationTag>& s, const time::EpochTDB& i) noexcept
        : state(s), instant(i) {}
};

} // namespace astdyn::types

#endif // ASTDYN_TYPES_TIMED_STATE_HPP
