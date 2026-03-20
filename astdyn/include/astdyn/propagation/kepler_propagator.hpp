#ifndef ASTDYN_PROPAGATION_KEPLER_PROPAGATOR_HPP
#define ASTDYN_PROPAGATION_KEPLER_PROPAGATOR_HPP

#include "propagator_base.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>

namespace astdyn::propagation {

using types::KeplerianTag;
namespace constants = ::astdyn::constants;

/**
 * @brief Computes mean motion based on gravitational parameter and semi-major axis.
 */
[[nodiscard]] constexpr double compute_mean_motion(const double a, const double mu) noexcept {
    const double a3 = a * a * a;
    return std::sqrt(mu / a3);
}

/**
 * @brief Updates mean anomaly based on time delta in seconds.
 */
[[nodiscard]] constexpr double update_anomaly(const double m0, const double n, const double dt_seconds) noexcept {
    return m0 + (n * dt_seconds);
}

/**
 * @brief Analytical Keplerian Propagator.
 * Adheres to CTFYH Rule: Max 3 constructor parameters, max 15 lines per function.
 */
template <typename Frame>
class KeplerPropagator final : public AnalyticalPropagator<Frame, KeplerianTag> {
public:
    explicit constexpr KeplerPropagator(const double gm) noexcept : gravitational_parameter_(gm) {}

    [[nodiscard]] std::optional<TimedState<Frame, KeplerianTag>> 
    propagate(const TimedState<Frame, KeplerianTag>& initial, const time::EpochTDB& target_time) const override {
        const double dt_days = target_time.mjd() - initial.instant.mjd();
        if (dt_days < 0.0) return std::nullopt;
        
        const double dt_sec = dt_days * constants::SECONDS_PER_DAY;
        const double a = initial.state.a();
        const double n = compute_mean_motion(a, gravitational_parameter_);
        
        const double m_new = update_anomaly(initial.state.m_anomaly(), n, dt_sec);
        const auto new_data = rebuild_data(initial.state.raw_values(), m_new);

        return TimedState<Frame, KeplerianTag>(types::OrbitalState<Frame, KeplerianTag>(new_data), target_time);
    }

private:
    const double gravitational_parameter_;

    /** @brief Helper for updating only the mean anomaly index in raw data array. */
    [[nodiscard]] static constexpr std::array<double, 6> rebuild_data(std::array<double, 6> raw, const double m_new) noexcept {
        raw[5] = m_new; // Mean anomaly is index 5
        return raw;
    }
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_KEPLER_PROPAGATOR_HPP
