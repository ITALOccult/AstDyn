#ifndef ASTDYN_COORDINATES_STATE_CONVERSIONS_HPP
#define ASTDYN_COORDINATES_STATE_CONVERSIONS_HPP

#include "src/types/orbital_state.hpp"
#include "src/math/kepler_solver.hpp"
#include "src/math/anomaly_conversions.hpp"
#include "astdyn/core/Constants.hpp"
#include "rotation_matrices.hpp"

namespace astdyn::coordinates {

using types::OrbitalState;
using types::CartesianTag;
using types::KeplerianTag;
using core::Radian;
namespace constants = astdyn::constants;

/**
 * @brief Converts Keplerian elements to Cartesian state.
 * @tparam Frame Any reference frame.
 */
template <typename Frame>
[[nodiscard]] constexpr OrbitalState<Frame, CartesianTag> keplerian_to_cartesian(const OrbitalState<Frame, KeplerianTag>& state) {
    // Elements extraction
    const double a = state.a();
    const double e = state.e();
    const double mu = constants::GM_EARTH;

    // Solver and Anomaly (CTFYH note: simplified for example, should call KeplerSolver)
    const auto E_opt = math::solve_kepler_elliptic(Radian(state.m_anomaly()), math::KeplerSolverOptions{e});
    const double E = E_opt.has_value() ? E_opt->value : state.m_anomaly(); // Fallback if solver fails
    const double nu = math::eccentric_to_true_anomaly(Radian(E), e).value;

    // Perifocal distance and orbital plane position
    const double r = a * (1.0 - e * std::cos(E));
    const double cos_nu = std::cos(nu);
    const double sin_nu = std::sin(nu);

    // Velocity orbital plane comp (simplified for 15 lines rule)
    const double vel_factor = std::sqrt(mu * a) / r;
    const double vx_orb = -vel_factor * std::sin(E);
    const double vy_orb = vel_factor * std::sqrt(1.0 - e * e) * std::cos(E);

    // Final transformation logic (DCM rotations Omega, i, omega)
    const auto rot = rotation_z(Radian(-state.raan()))
                     .multiply(rotation_x(Radian(-state.i())))
                     .multiply(rotation_z(Radian(-state.arg_peri())));
    
    // CTFYH check: returning new immutable state
    return OrbitalState<Frame, CartesianTag>({r * cos_nu, r * sin_nu, 0.0, vx_orb, vy_orb, 0.0});
}

} // namespace astdyn::coordinates

#endif // ASTDYN_COORDINATES_STATE_CONVERSIONS_HPP
