#ifndef ASTDYN_PROPAGATION_MULTI_BODY_PROPAGATOR_HPP
#define ASTDYN_PROPAGATION_MULTI_BODY_PROPAGATOR_HPP

#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <vector>
#include <string>
#include <memory>

namespace astdyn::propagation {

/**
 * @brief State of a single body in a multi-body system
 */
struct MultiBodyState {
    std::string name;
    physics::GravitationalParameter gm;
    math::Vector3<core::ECLIPJ2000, physics::Distance> position;
    math::Vector3<core::ECLIPJ2000, physics::Velocity> velocity;
};

/**
 * @brief Propagator for systems of mutually interacting bodies
 * 
 * Unlike the standard Propagator (which handles a test particle in a 
 * fixed planetary background), this class evolves N bodies 
 * simultaneously, accounting for their mutual gravitational attractions
 * and the Sun's influence.
 */
class MultiBodyPropagator {
public:
    explicit MultiBodyPropagator(std::shared_ptr<Integrator> integrator, 
                                std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris = nullptr);

    /**
     * @brief Propagate the system to a target epoch
     * 
     * @param initial_states Initial states of all bodies (must be at the same epoch)
     * @param start_time Epoch of initial states
     * @param target_time Target epoch
     * @return Vector of evolved states at the target epoch
     */
    std::vector<MultiBodyState> propagate(
        const std::vector<MultiBodyState>& initial_states,
        time::EpochTDB start_time,
        time::EpochTDB target_time);

private:
    struct Dynamics {
        double gm_sun;
        std::vector<double> gms;
        size_t n_bodies;
        time::EpochTDB t0_ephemeris;
        std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris;

        Eigen::VectorXd operator()(double t_sec, const Eigen::VectorXd& y) const;
    };

    std::shared_ptr<Integrator> integrator_;
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris_;
};

} // namespace astdyn::propagation

#endif // ASTDYN_PROPAGATION_MULTI_BODY_PROPAGATOR_HPP
