/**
 * @file PotentialDerivatives.hpp
 * @brief Closed-form second and third spatial derivatives of the potential.
 *
 * Provides U_ij (Hessian, 3x3) and U_ijk (third derivative, 3x3x3) for the
 * full force model used by the STT: Keplerian central term, solar J2, and a
 * sum of direct N-body perturbers.  These are the analytic quantities of
 * Appendix B; the AAS integrator's compute_hessian() implements only the
 * Keplerian primary term (J2/N-body are stubbed there), so these functions
 * both complete and generalise it.
 *
 * Units: AU, day; gravitational parameters in AU^3/day^2.  Sign convention:
 * U = -mu/r (attractive), so a = -grad U and the acceleration gradient is
 * G = da/dr = -U_ij, the acceleration Hessian T3_ijk = -U_ijk.
 */

#ifndef ASTDYN_PROPAGATION_POTENTIAL_DERIVATIVES_HPP
#define ASTDYN_PROPAGATION_POTENTIAL_DERIVATIVES_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/math/Tensor3.hpp"
#include <Eigen/Dense>
#include <vector>

namespace astdyn::propagation {

/**
 * @brief Static (frozen-ephemeris) force model for the STT derivatives.
 *
 * central_gm    : central-body GM [AU^3/day^2].
 * j2, r_eq      : solar oblateness (dimensionless J2, radius in AU); the pole
 *                 is the +z axis of the working frame.  j2 == 0 disables it.
 * perturber_gm  : perturber GMs [AU^3/day^2].
 * perturber_pos : perturber positions [AU], index-aligned with perturber_gm;
 *                 refreshed by the caller each macro-step from the ephemeris.
 */
struct ForceModel {
    double central_gm = 0.0002959122082855911;   // GM_sun, AU^3/day^2
    double j2 = 0.0;
    double r_eq = 0.004650467261;                  // R_sun in AU
    std::vector<double> perturber_gm;
    std::vector<astdyn::Vector3d> perturber_pos;
};

/// Acceleration a = -grad U (Keplerian + J2 + direct N-body) [AU/day^2].
astdyn::Vector3d acceleration(const astdyn::Vector3d& r, const ForceModel& m);

/// Potential Hessian U_ij (3x3), Keplerian + J2 + N-body.
astdyn::Matrix3d potential_hessian(const astdyn::Vector3d& r, const ForceModel& m);

/// Potential third derivative U_ijk (3x3x3), Keplerian + J2 + N-body (App. B).
astdyn::math::Tensor3
potential_third_derivative(const astdyn::Vector3d& r, const ForceModel& m);

// --- individual analytic pieces (exposed for testing/reuse) --------------
astdyn::Matrix3d kepler_hessian(const astdyn::Vector3d& u, double mu);
astdyn::math::Tensor3 kepler_third(const astdyn::Vector3d& u, double mu);
astdyn::Matrix3d j2_hessian(const astdyn::Vector3d& r, double mu, double j2, double r_eq);
astdyn::math::Tensor3 j2_third(const astdyn::Vector3d& r, double mu, double j2, double r_eq);

}  // namespace astdyn::propagation

#endif  // ASTDYN_PROPAGATION_POTENTIAL_DERIVATIVES_HPP
