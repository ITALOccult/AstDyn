#ifndef ASTDYN_COORDINATES_STATE_CONVERSIONS_HPP
#define ASTDYN_COORDINATES_STATE_CONVERSIONS_HPP

#include "astdyn/types/orbital_state.hpp"
#include "astdyn/math/kepler_solver.hpp"
#include "astdyn/math/anomaly_conversions.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "rotation_matrices.hpp"

namespace astdyn::coordinates {

using types::OrbitalState;
using types::CartesianTag;
using types::KeplerianTag;
using types::AuDayTag;
using core::Radian;
namespace constants = astdyn::constants;

/**
 * @brief Converts Keplerian elements to Cartesian state.
 * @tparam Frame Any reference frame.
 */
template <typename Frame>
[[nodiscard]] constexpr OrbitalState<Frame, CartesianTag, AuDayTag> keplerian_to_cartesian(
    const OrbitalState<Frame, KeplerianTag>& state, 
    const double mu = constants::GM_EARTH) {
    const double a = state.a(), e = state.e();
    const auto E_opt = math::solve_kepler_elliptic(Radian(state.m_anomaly()), math::KeplerSolverOptions{e});
    const double E = E_opt.has_value() ? E_opt->value : state.m_anomaly();
    const double nu = math::eccentric_to_true_anomaly(Radian(E), e).value;

    const double r = a * (1.0 - e * std::cos(E)), v_f = std::sqrt(mu * a) / r;
    const double cos_n = std::cos(nu), sin_n = std::sin(nu);
    
    // Perifocal state (Planar)
    math::Vector3<Frame, core::Meter> pos_p = math::Vector3<Frame, core::Meter>::from_si(r * cos_n, r * sin_n, 0.0);
    math::Vector3<Frame, core::Meter> vel_p = math::Vector3<Frame, core::Meter>::from_si(-v_f * std::sin(E), v_f * std::sqrt(1.0 - e * e) * std::cos(E), 0.0);

    // Rotation: Perifocal -> Local Frame
    const auto rot = rotation_z(Radian(-state.raan()))
                     .multiply(rotation_x(Radian(-state.i())))
                     .multiply(rotation_z(Radian(-state.arg_peri())));
    
    const auto final_pos = rot.multiply(pos_p);
    const auto final_vel = rot.multiply(vel_p);
    
    return OrbitalState<Frame, CartesianTag, AuDayTag>({
        final_pos.x_si(), final_pos.y_si(), final_pos.z_si(), 
        final_vel.x_si(), final_vel.y_si(), final_vel.z_si()
    });
}

} // namespace astdyn::coordinates

#endif // ASTDYN_COORDINATES_STATE_CONVERSIONS_HPP
