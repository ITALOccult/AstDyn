/**
 * @file StateTransitionTensor.hpp
 * @brief Joint propagation of state, STM and STT through the AAS kernel.
 *
 * Extends the analytical State Transition Matrix (STM) of AstDyn to the
 * second-order State Transition Tensor (STT), integrating (x, Phi, Psi)
 * together through the same symplectic Drift-Kick-Drift (Yoshida-4) kernel
 * and, crucially, with the *same adaptive step metric* as the AAS
 * integrator (step from the local potential-curvature / force-gradient).
 *
 * The STT is driven by the closed-form third derivative U_ijk
 * (PotentialDerivatives); each Kick carries the source -h U_ijk Phi Phi
 * (Appendix B, Eq. B8).  Public boundary uses strong frame/unit types;
 * internals are Eigen in au, day, mu in au^3/day^2 (ECLIPJ2000).
 */

#ifndef ASTDYN_PROPAGATION_STATE_TRANSITION_TENSOR_HPP
#define ASTDYN_PROPAGATION_STATE_TRANSITION_TENSOR_HPP

#include "astdyn/propagation/ForceField.hpp"
#include "astdyn/math/Tensor3.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/core/frame_tags.hpp"
#include "astdyn/time/epoch.hpp"
#include <Eigen/Dense>
#include <functional>
#include <memory>

namespace astdyn::propagation {

/**
 * @brief Second-order variational propagator (STM + STT) with AAS step control.
 */
class StateTransitionTensor {
public:
    using Frame = astdyn::core::ECLIPJ2000;
    using State = astdyn::physics::CartesianStateTyped<Frame>;

    /// Optional hook to refresh perturber positions at a given epoch

    struct Result {
        State final_state;              ///< propagated mean state at t1
        astdyn::Matrix6d phi;           ///< STM  Phi(t1, t0)   [6x6]
        astdyn::math::Tensor6 psi;      ///< STT  Psi(t1, t0)   [6x6x6]
    };

    /**
     * @param model     force model (central GM, J2, perturbers) in au/day units
     * @param precision AAS step-metric precision (as in AASIntegrator)
     */
    explicit StateTransitionTensor(std::shared_ptr<const ForceField> force, double precision = 1e-4);

    /**
     * @brief Propagate (state, Phi, Psi) from x0.epoch to t1 (forward in time).
     * @param refresh optional per-step ephemeris update of the perturbers.
     */
    [[nodiscard]] Result propagate(const State& x0, time::EpochTDB t1) const;

    /**
     * @brief Fixed-step variant (uniform n_steps, no adaptive metric).
     * Reproducible discrete map — used for the finite-difference validation
     * of Phi and Psi (the nominal and perturbed flows share the same steps).
     */
    [[nodiscard]] Result propagate_fixed(const State& x0, time::EpochTDB t1,
                                         int n_steps) const;

private:
    [[nodiscard]] State to_state(const astdyn::Vector3d& q, const astdyn::Vector3d& p,
                                 time::EpochTDB t,
                                 const astdyn::physics::GravitationalParameter& gm) const;

    // one adaptive macro-step (7-stage Yoshida DKD) advancing dt days
    void step(astdyn::Vector3d& q, astdyn::Vector3d& p,
              astdyn::Matrix6d& phi, astdyn::math::Tensor6& psi,
              double dt, time::EpochTDB t) const;
    void kick(const astdyn::Vector3d& q, astdyn::Vector3d& p,
              astdyn::Matrix6d& phi, astdyn::math::Tensor6& psi,
              double h, time::EpochTDB t) const;
    void drift(astdyn::Vector3d& q, const astdyn::Vector3d& p,
               astdyn::Matrix6d& phi, astdyn::math::Tensor6& psi, double h) const;
    // AAS adaptive metric (spectral radius via ForceField::gradient_spectral_radius)
    [[nodiscard]] double estimate_step(const astdyn::Vector3d& q,
                                       const astdyn::Vector3d& p,
                                       double target_dt, time::EpochTDB t) const;

    std::shared_ptr<const ForceField> force_;
    double precision_;
    double w1_, w0_, c1_, c2_, d1_, d2_;   ///< Yoshida-4 coefficients
};

}  // namespace astdyn::propagation

#endif  // ASTDYN_PROPAGATION_STATE_TRANSITION_TENSOR_HPP
