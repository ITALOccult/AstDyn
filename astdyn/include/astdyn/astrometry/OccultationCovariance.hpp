/**
 * @file OccultationCovariance.hpp
 * @brief The SCOPE composite event map: dynamics + projection, to second order.
 *
 * The shadow position on the fundamental plane is a composition of two maps: the
 * dynamical transport of the state from t0 to the event, and the projection of
 * that state onto the plane of the sky. Linearising the composition is what
 * every occultation prediction does; carrying it to second order is what lets
 * the bias, the skewness and -- above all -- the nonlinearity index be computed
 * WITHOUT a Monte Carlo.
 *
 * Notation follows the paper:
 *
 *   A_aj    = A_h Phi                                   composite Jacobian   (Eq. 5)
 *   B_ajk   = A_h Psi + H_h (Phi x Phi)                 composite Hessian    (Eq. 13)
 *   bias_a  = 1/2 B_ajk C0_jk                                                (Eq. 14)
 *   C_ab    = A_aj A_bk C0_jk + 1/2 B_ajk B_blm C0_jl C0_km                  (Eq. 15)
 *   kappa   = A_aj A_bk B_clm C0_jl C0_km + (cyclic)                         (Eq. 16)
 *   N       = 1/2 sqrt(2 B_n:C0:B_n) / sqrt(A_n C0 A_n)                      (Eq. 18)
 *
 * The index N is the point of the whole exercise: it is the ratio of the RMS
 * second-order term to the RMS first-order term along the cross-track direction,
 * and it is available a priori. For a well-determined main-belt asteroid 48 days
 * out it comes out around 1e-7, i.e. the linear theory is exact for all practical
 * purposes and N says so. For a NEO with a short arc it does not, and N says that
 * too -- before anyone has to run a Monte Carlo to find out.
 *
 * Validated against the NumPy reference in projection.py, itself validated
 * against Monte Carlo.
 */
#ifndef ASTDYN_ASTROMETRY_OCCULTATION_COVARIANCE_HPP
#define ASTDYN_ASTROMETRY_OCCULTATION_COVARIANCE_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/math/Tensor3.hpp"
#include <Eigen/Dense>
#include <array>

namespace astdyn::astrometry {

/// Orthonormal basis of the plane of the sky: e[0] and e[1] span it, both
/// perpendicular to the line of sight.
using TangentBasis = std::array<Vector3d, 2>;

/// Projection maps at the event (Eqs. 4 and 7).
struct ProjectionMaps {
    Eigen::Matrix<double, 2, 6> A_h;          ///< d xi_a / d x
    std::array<Matrix6d, 2>     H_h;          ///< d2 xi_a / dx dx, one 6x6 per component
};

/// Composite event map (Eqs. 5 and 13).
struct CompositeMap {
    Eigen::Matrix<double, 2, 6> A;            ///< A_h Phi
    std::array<Matrix6d, 2>     B;            ///< A_h Psi + H_h (Phi x Phi)
};

/// Leading moments of the projected observable (Eqs. 14-16).
struct Moments {
    Eigen::Vector2d bias;                     ///< 1/2 B_ajk C0_jk
    Eigen::Matrix2d C_xi;                     ///< linear + second-order
    Eigen::Matrix2d C_lin;                    ///< linear part alone
    /// kappa_abc, third cumulant; [a][b][c] with a,b,c in {0,1}.
    std::array<std::array<Eigen::Vector2d, 2>, 2> kappa;
};

/**
 * @brief Orthonormal basis of the plane of the sky about @p s_hat.
 *
 * @param s_hat  Unit line of sight.
 * @param pole   Reference direction used to fix the roll; must not be parallel
 *               to @p s_hat.
 */
[[nodiscard]] TangentBasis tangent_basis(const Vector3d& s_hat,
                                         const Vector3d& pole = Vector3d::UnitZ());

/**
 * @brief Projection Jacobian and Hessian at the event (Eqs. 4 and 7).
 *
 * Exact, not only at closest approach: the position block of A_h is
 * (1/rho)(e_a - (e_a . rho_hat) rho_hat), which reduces to e_a/rho when e_a is
 * perpendicular to the line of sight. The velocity blocks vanish because the
 * observable depends on the position alone.
 *
 * @param r      Object position relative to the observer's frame origin.
 * @param basis  Plane-of-sky basis.
 * @param r_obs  Observer position, if not at the origin.
 */
[[nodiscard]] ProjectionMaps projection_maps(const Vector3d& r,
                                             const TangentBasis& basis,
                                             const Vector3d& r_obs = Vector3d::Zero());

/**
 * @brief Compose the projection with the dynamics (Eqs. 5 and 13).
 *
 * @param pm   Projection maps at the event.
 * @param Phi  State transition matrix from t0 to the event.
 * @param Psi  State transition tensor from t0 to the event.
 */
[[nodiscard]] CompositeMap composite_map(const ProjectionMaps& pm,
                                         const Matrix6d& Phi,
                                         const math::Tensor6& Psi);

/**
 * @brief Leading moments for delta x0 ~ N(0, C0), by Isserlis (Eqs. 14-16).
 *
 * The linear theory predicts bias = 0 and kappa = 0 identically; both are pure
 * second-order effects, which is why they are the sharpest test of the map.
 */
[[nodiscard]] Moments moments(const CompositeMap& cm, const Matrix6d& C0);

/**
 * @brief A priori nonlinearity index along @p n_hat (Eq. 18).
 *
 * Ratio of the RMS second-order term to the RMS first-order term of the
 * projected observable. No Monte Carlo required.
 *
 * @param n_hat  Unit direction in the plane of sky, normally cross-track.
 */
[[nodiscard]] double nonlinearity_index(const CompositeMap& cm, const Matrix6d& C0,
                                        const Eigen::Vector2d& n_hat);

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_OCCULTATION_COVARIANCE_HPP
