#ifndef ASTDYN_PROPAGATION_FORCE_FIELD_HPP
#define ASTDYN_PROPAGATION_FORCE_FIELD_HPP

#include "astdyn/core/Constants.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/AsteroidPerturbations.hpp"
#include "astdyn/propagation/PropagatorSettings.hpp"
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
