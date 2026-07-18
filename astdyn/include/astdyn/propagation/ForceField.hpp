#ifndef ASTDYN_PROPAGATION_FORCE_FIELD_HPP
#define ASTDYN_PROPAGATION_FORCE_FIELD_HPP

#include "astdyn/core/Constants.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/AsteroidPerturbations.hpp"
#include "astdyn/propagation/PropagatorSettings.hpp"
#include "astdyn/math/Tensor3.hpp"
#include <Eigen/Core>
#include <memory>
#include <stdexcept>
#include <vector>

namespace astdyn::propagation {

/**
 * @brief Base interface for force model components (CTFYH Standard)
 */
class ForceModel {
public:
    virtual ~ForceModel() = default;
    
    /**
     * @brief Compute acceleration [AU/day^2]
     * @param t Time (TDB)
     * @param pos_au Particle position [AU] (Heliocentric or Barycentric depending on field)
     * @param vel_au_d Particle velocity [AU/day]
     * @return Acceleration [AU/day^2]
     */
    virtual Eigen::Vector3d compute_acceleration(
        time::EpochTDB t, 
        const Eigen::Vector3d& pos_au, 
        const Eigen::Vector3d& vel_au_d) const = 0;

    /**
     * @brief Gradient of the acceleration wrt position: dA_i/dr_j [1/day^2].
     *
     * This is the Hessian of the potential (with a sign), and it is what the
     * state transition tensor needs to propagate the covariance. It is NOT the
     * same information as the acceleration: an integrator can advance the
     * trajectory with compute_acceleration() alone, but the variational
     * equations of a covariance cannot.
     *
     * Default: throws. A force term supplies this only when its analytic
     * derivative is worth carrying. Terms whose contribution to the tensor is
     * negligible (relativity, Yarkovsky: ~1e-8 of the Keplerian gradient) may
     * legitimately leave it unimplemented and contribute to the acceleration
     * only -- that omission is physics, and belongs in the paper, not a silent
     * shortcut.
     */
    [[nodiscard]] virtual Eigen::Matrix3d acceleration_gradient(
        time::EpochTDB /*t*/,
        const Eigen::Vector3d& /*pos_au*/,
        const Eigen::Vector3d& /*vel_au_d*/) const {
        throw std::logic_error(
            "acceleration_gradient not implemented for this force term");
    }

    /**
     * @brief Whether acceleration_gradient() is available.
     *
     * Lets the tensor assemble the gradient from only the terms that provide
     * one, and record which terms were skipped, instead of discovering the gap
     * through a thrown exception mid-propagation.
     */
    [[nodiscard]] virtual bool has_acceleration_gradient() const { return false; }

    /**
     * @brief Second gradient of acceleration: d^2 A_i/dr_j dr_k [1/(AU day^2)].
     *
     * Third spatial derivative of the potential, U_ijk. The state transition
     * tensor needs it for the second-order variational equations (Psi): bias,
     * skewness, and the non-Gaussianity index N. Phi (first-order STM) needs
     * only the Hessian above; Psi cannot be propagated without this.
     *
     * Default: throws, mirroring acceleration_gradient().
     */
    [[nodiscard]] virtual math::Tensor3 acceleration_second_gradient(
        time::EpochTDB /*t*/,
        const Eigen::Vector3d& /*pos_au*/,
        const Eigen::Vector3d& /*vel_au_d*/) const {
        throw std::logic_error(
            "acceleration_second_gradient not implemented for this force term");
    }

    [[nodiscard]] virtual bool has_acceleration_second_gradient() const { return false; }

    /**
     * @brief Spectral radius of the acceleration gradient at pos_au [1/day^2].
     *
     * AAS adaptive step control needs max local |2 mu / r^3| (plus J2) across
     * all attractors, not the full Hessian. Kept separate so the step metric
     * does not force a full gradient assembly at every substep.
     */
    [[nodiscard]] virtual double gradient_spectral_radius(
        time::EpochTDB /*t*/,
        const Eigen::Vector3d& /*pos_au*/) const {
        throw std::logic_error(
            "gradient_spectral_radius not implemented for this force term");
    }
};

/**
 * @brief Standard Force Field aggregator
 */
class ForceField {
public:
    explicit ForceField(const PropagatorSettings& settings, 
                        std::shared_ptr<ephemeris::PlanetaryEphemeris> ephem = nullptr);

    /**
     * @brief Compute total acceleration including all active models
     * Output is in AU/day^2.
     */
    Eigen::Vector3d total_acceleration(
        time::EpochTDB t, 
        const Eigen::Vector3d& pos_au, 
        const Eigen::Vector3d& vel_au_d) const;

    const PropagatorSettings& settings() const { return settings_; }
    void update_settings(const PropagatorSettings& s) { settings_ = s; }

private:
    PropagatorSettings settings_;
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<ephemeris::AsteroidPerturbations> asteroids_;

    // Individual force calculations (internal refactoring for reuse)
    Eigen::Vector3d n_body_perturbation(time::EpochTDB t, const Eigen::Vector3d& pos_au) const;
    Eigen::Vector3d relativistic_correction(const Eigen::Vector3d& r, const Eigen::Vector3d& v) const;
    Eigen::Vector3d j2_correction(time::EpochTDB t, const Eigen::Vector3d& pos_au) const;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_FORCE_FIELD_HPP
