#ifndef ASTDYN_MATH_ANOMALY_CONVERSIONS_HPP
#define ASTDYN_MATH_ANOMALY_CONVERSIONS_HPP

#include "src/core/units.hpp"
#include <cmath>

namespace astdyn::math {

using core::Radian;

/**
 * @brief Converts eccentric anomaly to true anomaly for elliptic orbits.
 * @param E Eccentric anomaly.
 * @param eccentricity Dimensionless eccentricity (0 <= e < 1).
 * @return True anomaly in range [0, 2pi).
 */
[[nodiscard]] constexpr Radian eccentric_to_true_anomaly(const Radian E, const double eccentricity) noexcept {
    const double e = eccentricity;
    const double half_E = E.value / 2.0;
    
    const double sqrt_factor = std::sqrt((1.0 + e) / (1.0 - e));
    const double tan_half_E = std::tan(half_E);
    
    const double nu = 2.0 * std::atan(sqrt_factor * tan_half_E);
    return Radian(nu); // Normalization done by consumer if needed
}

/**
 * @brief Converts true anomaly to eccentric anomaly for elliptic orbits.
 * @param nu True anomaly.
 * @param eccentricity Dimensionless eccentricity (0 <= e < 1).
 * @return Eccentric anomaly.
 */
[[nodiscard]] constexpr Radian true_to_eccentric_anomaly(const Radian nu, const double eccentricity) noexcept {
    const double e = eccentricity;
    const double half_nu = nu.value / 2.0;

    const double sqrt_factor = std::sqrt((1.0 - e) / (1.0 + e));
    const double tan_half_nu = std::tan(half_nu);
    
    const double E = 2.0 * std::atan(sqrt_factor * tan_half_nu);
    return Radian(E);
}

} // namespace astdyn::math

#endif // ASTDYN_MATH_ANOMALY_CONVERSIONS_HPP
