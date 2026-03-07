#ifndef ASTDYN_PROPAGATION_PROPAGATOR_BASE_HPP
#define ASTDYN_PROPAGATION_PROPAGATOR_BASE_HPP

#include "src/types/timed_state.hpp"
#include "src/core/messages.hpp"
#include <optional>

/** 
 * Note: Switched to std::optional for C++17 compatibility.
 */
namespace astdyn::propagation {

using core::MessageKey;
using types::TimedState;

/**
 * @brief Abstract base for orbital propagators.
 * Defines the contract for time evolution of orbital states.
 */
template <typename Frame, typename Rep>
class AnalyticalPropagator {
public:
    virtual ~AnalyticalPropagator() = default;

    /**
     * @brief Evolves the initial TimedState to the target_time.
     * @param initial Coupled state and instant at epoch.
     * @param target_time Integration target instant.
     * @return Expected new TimedState or an error key.
     */
    [[nodiscard]] virtual std::optional<TimedState<Frame, Rep>> 
    propagate(const TimedState<Frame, Rep>& initial, const time::EpochTDB& target_time) const = 0;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_PROPAGATOR_BASE_HPP
