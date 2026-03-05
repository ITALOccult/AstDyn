#ifndef ASTDYN_PROPAGATION_J2_PROPAGATOR_HPP
#define ASTDYN_PROPAGATION_J2_PROPAGATOR_HPP

#include "propagator_base.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>

namespace astdyn::propagation {

using types::KeplerianTag;
namespace constants = ::astdyn::constants;

/** 
 * @brief Secular drift rates for J2 perturbation. 
 * CTFYH: Extracted logic for readability and complexity reduction.
 */
struct J2Rates {
    const double raan_dot;
    const double omega_dot;
    const double m_dot;

    /** @brief Computes secular rates based on Earth J2 model. */
    static constexpr J2Rates compute(const double a, const double e, const double i, const double mu) noexcept {
        const double J2 = constants::EARTH_J2;
        const double RE = constants::R_EARTH;
        
        const double p = a * (1.0 - e * e);
        const double n = std::sqrt(mu / (a * a * a));
        const double factor = 0.75 * J2 * std::pow(RE / p, 2.0) * n;

        const double cos_i = std::cos(i);
        const double rdot = -2.0 * factor * cos_i; // RAAN dot
        const double wdot = factor * (5.0 * cos_i * cos_i - 1.0); // Argument of Periapsis dot
        const double mdot = n + factor * std::sqrt(1.0 - e * e) * (3.0 * cos_i * cos_i - 1.0); // Mean Anomaly dot

        return J2Rates{rdot, wdot, mdot};
    }
};

/**
 * @brief Secular J2 Propagator for Keplerian elements.
 * CTFYH adherence: Max 3 constructor parameters, const methods, explicit.
 */
template <typename Frame>
class J2Propagator final : public AnalyticalPropagator<Frame, KeplerianTag> {
public:
    explicit constexpr J2Propagator(const double gm) noexcept : mu_(gm) {}

    [[nodiscard]] std::optional<TimedState<Frame, KeplerianTag>> 
    propagate(const TimedState<Frame, KeplerianTag>& initial, const Instant& target_time) const override {
        const double dt_days = target_time.mjd.value - initial.instant.mjd.value;
        if (dt_days < 0.0) return std::nullopt;
        
        const double dt_sec = dt_days * constants::SECONDS_PER_DAY;
        const auto rates = J2Rates::compute(initial.state.a(), initial.state.e(), initial.state.i(), mu_);

        auto data = initial.state.raw_values();
        data[3] += rates.raan_dot * dt_sec;
        data[4] += rates.omega_dot * dt_sec;
        data[5] += rates.m_dot * dt_sec;

        return TimedState<Frame, KeplerianTag>(types::OrbitalState<Frame, KeplerianTag>(data), target_time);
    }

private:
    const double mu_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_J2_PROPAGATOR_HPP
