#ifndef ASTDYN_PROPAGATION_RELATIVE_MULTI_BODY_PROPAGATOR_HPP
#define ASTDYN_PROPAGATION_RELATIVE_MULTI_BODY_PROPAGATOR_HPP

#include "astdyn/propagation/MultiBodyTypes.hpp"
#include "astdyn/propagation/ForceField.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include <memory>
#include <vector>
#include <Eigen/Dense>

namespace astdyn::propagation {

/**
 * @brief High-precision multi-body dynamic model using relative vectors.
 * y format: [r0, v0, rho1, drho1, ... rhoN, drhoN]
 */
struct RelativeDynamics {
    std::shared_ptr<ForceField> force_field;
    std::vector<double> gms; // [gm0, gm1, ... gmN]
    size_t n_bodies;
    time::EpochTDB t0;

    Eigen::VectorXd operator()(double t_sec, const Eigen::VectorXd& y) const;
};

/**
 * @brief High-precision multi-body propagator using relative vectors (CTFYH Standard)
 * 
 * Instead of absolute heliocentric coordinates, this class integrates:
 * 1. The primary body's absolute position.
 * 2. The relative vectors of satellites relative to the primary.
 */
class RelativeMultiBodyPropagator {
public:
    explicit RelativeMultiBodyPropagator(std::shared_ptr<Integrator> integrator,
                                         std::shared_ptr<ForceField> force_field);

    /**
     * @brief Propagate a system of bodies
     * @param initial_states States in heliocentric ECLIPJ2000
     * @param start_time Epoch
     * @param target_time Target epoch
     */
    std::vector<MultiBodyState> propagate(
        const std::vector<MultiBodyState>& initial_states,
        time::EpochTDB start_time,
        time::EpochTDB target_time);

private:
    std::shared_ptr<Integrator> integrator_;
    std::shared_ptr<ForceField> force_field_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_RELATIVE_MULTI_BODY_PROPAGATOR_HPP
