#ifndef ASTDYN_MATH_KEPLER_SOLVER_HPP
#define ASTDYN_MATH_KEPLER_SOLVER_HPP

#include "astdyn/core/units.hpp"
#include "trigonometry.hpp"
#include <cmath>
#include <optional>

namespace astdyn::math {

using core::Radian;

/**
 * @brief Parameters for Kepler equation solver to keep function arguments <= 3.
 */
struct KeplerSolverOptions {
    const double eccentricity;
    const double tolerance = 1e-12;
    const int max_iterations = 100;
};

/**
 * @brief Computes a Newton-Raphson step for Kepler's equation: E - e*sin(E) - M = 0.
 * @return delta The adjustment to be applied to E.
 */
[[nodiscard]] constexpr double kepler_newton_step(const double E_val, const double M_val, const double e_val) noexcept {
    const double f_val = E_val - e_val * std::sin(E_val) - M_val;
    const double f_prime = 1.0 - e_val * std::cos(E_val);
    return -f_val / f_prime;
}

/**
 * @brief Solves Kepler's equation for elliptic orbits (e < 1.0).
 * @param M Mean anomaly in Radians.
 * @param options Solver options (eccentricity, tolerance, max_iterations).
 * @return std::optional<Radian> Eccentric anomaly, or empty if failed.
 */
[[nodiscard]] constexpr std::optional<Radian> solve_kepler_elliptic(const Radian M, const KeplerSolverOptions& options) noexcept {
    const double M_val = wrap_to_2pi(M).value;
    const double e_val = options.eccentricity;
    double E_val = M_val; // Simple initial guess
    
    if (e_val >= 1.0) return std::nullopt; // Guard clause: non-elliptic not supported here

    for (int iter = 0; iter < options.max_iterations; ++iter) {
        const double delta = kepler_newton_step(E_val, M_val, e_val);
        E_val += delta;
        if (std::abs(delta) < options.tolerance) {
            return Radian(E_val);
        }
    }
    
    return std::nullopt;
}

/**
 * @brief Solves the hyperbolic Kepler equation for hyperbolic orbits (e > 1).
 *
 * Solves: e * sinh(H) - H = Mh,  where Mh is the hyperbolic mean anomaly.
 * Uses Newton-Raphson iteration: H_new = H - (e*sinh(H) - H - Mh) / (e*cosh(H) - 1)
 *
 * @param Mh Hyperbolic mean anomaly (dimensionless).
 * @param options Solver options (eccentricity must be > 1, tolerance, max_iterations).
 * @return std::optional<double> Hyperbolic eccentric anomaly H, or empty if failed / e <= 1.
 */
[[nodiscard]] inline std::optional<double> solve_kepler_hyperbolic(const double Mh, const KeplerSolverOptions& options) noexcept {
    const double e = options.eccentricity;
    if (e <= 1.0) return std::nullopt; // Guard: only for hyperbolic orbits

    // Initial guess: for small |Mh|, H ≈ Mh/(e-1); for large |Mh|, H ≈ sign(Mh)*log(2|Mh|/e).
    double H;
    if (std::abs(Mh) < 1.0) {
        H = Mh / (e - 1.0);
    } else {
        H = std::copysign(std::log(2.0 * std::abs(Mh) / e + 1.8), Mh);
    }

    for (int iter = 0; iter < options.max_iterations; ++iter) {
        const double f      = e * std::sinh(H) - H - Mh;
        const double f_prime = e * std::cosh(H) - 1.0;
        if (std::abs(f_prime) < 1e-15) return std::nullopt; // Degenerate
        const double delta  = -f / f_prime;
        H += delta;
        if (std::abs(delta) < options.tolerance) {
            return H;
        }
    }
    return std::nullopt; // Did not converge
}

} // namespace astdyn::math

#endif // ASTDYN_MATH_KEPLER_SOLVER_HPP
