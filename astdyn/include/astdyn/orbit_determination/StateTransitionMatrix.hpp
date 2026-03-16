/**
 * @file StateTransitionMatrix.hpp
 * @brief State transition matrix (STM) computation
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 * 
 * Computes the state transition matrix Φ(t,t₀) = ∂x(t)/∂x₀
 * by integrating variational equations alongside the equations of motion.
 * 
 * Essential for differential corrections in orbit determination.
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_STATE_TRANSITION_MATRIX_HPP
#define ASTDYN_ORBIT_DETERMINATION_STATE_TRANSITION_MATRIX_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "src/core/frame_tags.hpp"
#include <memory>

namespace astdyn::orbit_determination {

/**
 * @brief STM computation settings
 */
struct STMSettings {
    bool use_numerical_jacobian = true; ///< Default to OrbFit style
    double differentiation_step = 1e-7;  ///< Default step [AU]
};

/**
 * @brief State transition matrix result
 */
template <typename Frame>
struct STMResult {
    astdyn::Matrix6d phi;                        ///< Φ(t,t₀) - 6x6 state transition matrix
    physics::CartesianStateTyped<Frame> final_state;       ///< Propagated state at time t
    
    /**
     * @brief Map covariance from t₀ to t
     * 
     * Cov(t) = Φ(t,t₀) * Cov(t₀) * Φ(t,t₀)ᵀ
     * 
     * @param cov_t0 Covariance at initial epoch
     * @return Mapped covariance at final epoch
     */
    astdyn::Matrix6d map_covariance(const astdyn::Matrix6d& cov_t0) const {
        return phi * cov_t0 * phi.transpose();
    }
};

/**
 * @brief Computes state transition matrix via variational equations
 * 
 * @tparam Frame Reference frame for integration
 */
template <typename Frame>
class StateTransitionMatrix {
public:
    /**
     * @brief Constructor
     * 
     * @param propagator Orbit propagator (provides force model)
     */
    explicit StateTransitionMatrix(std::shared_ptr<astdyn::propagation::Propagator> propagator,
                                   std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem = nullptr);
    
    /**
     * @brief Compute STM from t₀ to t
     * 
     * @param initial Initial Cartesian state at t₀
     * @param target_mjd_tdb Target time t
     * @return STM result with Φ(t,t₀) and final state
     */
    STMResult<Frame> compute(
        const physics::CartesianStateTyped<Frame>& initial,
        time::EpochTDB target_time);
    
    /**
     * @brief Compute STM and partials w.r.t observations
     * 
     * For orbit determination, we need:
     * - Φ(t,t₀): state transition matrix
     * - ∂(RA,Dec)/∂x: observation partials (2x6 matrix)
     * 
     * @param initial Initial state
     * @param target_mjd_tdb Observation time
     * @param observer_pos Observer position [AU]
     * @return STM and observation partials
     */
    struct ObservationPartials {
        astdyn::Matrix6d phi;                    ///< State transition matrix
        Eigen::Matrix<double, 2, 6> partial_radec; ///< ∂(RA,Dec)/∂x
        physics::CartesianStateTyped<Frame> final_state;
    };
    
    ObservationPartials compute_with_partials(
        const physics::CartesianStateTyped<Frame>& initial,
        time::EpochTDB target_time,
        const math::Vector3<core::GCRF, physics::Distance>& observer_pos);
    
    /**
     * @brief Compute STM at multiple observation epochs
     * 
     * Uses optimized batch integration hit targets without resetting the integrator.
     * 
     * @param initial Initial state at t0
     * @param target_times Vector of target times
     * @param observer_positions Optional observer positions for calculating observation partials
     * @return Vector of results
     */
    std::vector<ObservationPartials> compute_batch(
        const physics::CartesianStateTyped<Frame>& initial,
        const std::vector<time::EpochTDB>& target_times,
        const std::vector<math::Vector3<core::GCRF, physics::Distance>>& observer_positions);
    
    /**
     * @brief Set integrator for variational equations
     */
    void set_integrator(std::shared_ptr<astdyn::propagation::Integrator> integrator) {
        integrator_ = integrator;
    }

    /**
     * @brief Apply STM settings
     */
    void apply_settings(const STMSettings& settings) {
        use_numerical_jacobian_ = settings.use_numerical_jacobian;
        diff_step_ = settings.differentiation_step;
    }
    
    /**
     * @brief Set numerical differentiation step size
     * 
     * For computing Jacobian ∂f/∂x numerically (if not analytical).
     * 
     * @param step Step size [AU for position, AU/day for velocity]
     */
    void set_differentiation_step(double step) {
        diff_step_ = step;
    }

    /**
     * @brief Enable or disable numerical Jacobian computation
     * 
     * If true, uses finite differences on the full propagator derivatives (matching OrbFit).
     * If false, uses analytical gravity partials (traditional AstDyn behavior).
     */
    void set_use_numerical_jacobian(bool use_num) {
        use_numerical_jacobian_ = use_num;
    }

private:
    /**
     * @brief Compute Jacobian matrix A(t) = ∂f/∂x
     * 
     * For two-body problem:
     * A = [  0₃ₓ₃    I₃ₓ₃  ]
     *     [ ∂a/∂r  ∂a/∂v  ]
     * 
     * where a = acceleration = -μr/r³ (+ perturbations)
     * 
     * @param t Time
     * @param state State vector [x,y,z,vx,vy,vz]
     * @return 6x6 Jacobian matrix
     */
    astdyn::Matrix6d compute_jacobian(time::EpochTDB t, const Eigen::VectorXd& state);
    
    /**
     * @brief Compute ∂a/∂r for gravitational acceleration
     * 
     * For a = -μr/r³:
     * ∂a/∂r = -μ/r³ [I - 3(r⊗r)/r²]
     * 
     * @param r Position vector [AU]
     * @param mu Gravitational parameter [AU³/day²]
     * @return 3x3 partial derivative matrix
     */
    Eigen::Matrix3d compute_acceleration_position_partial(
        const astdyn::Vector3d& r,
        double mu) const;
    
    /**
     * @brief Propagate state and STM together
     * 
     * Integrates augmented state [x(6), Φ(36)] as a single 42-element vector.
     * 
     * @param initial Initial state
     * @param target_mjd Target time
     * @return Final state and STM
     */
    STMResult<Frame> propagate_with_stm(
        const physics::CartesianStateTyped<Frame>& initial,
        time::EpochTDB target_time);

    /**
     * @brief Internal batch propagator for augmented state
     */
    std::vector<STMResult<Frame>> propagate_with_stm_batch(
        const physics::CartesianStateTyped<Frame>& initial,
        const std::vector<time::EpochTDB>& target_times);
    
    /**
     * @brief Compute observation partials ∂(RA,Dec)/∂x
     * 
     * Chain rule: ∂obs/∂x₀ = ∂obs/∂x * ∂x/∂x₀ = ∂obs/∂x * Φ
     * 
     * @param state Final Cartesian state
     * @param observer_pos Observer position
     * @return 2x6 matrix of partials
     */
    Eigen::Matrix<double, 2, 6> compute_observation_partials(
        const physics::CartesianStateTyped<Frame>& state,
        const math::Vector3<core::GCRF, physics::Distance>& observer_pos) const;

private:
    std::shared_ptr<astdyn::propagation::Propagator> propagator_;
    std::shared_ptr<astdyn::propagation::Integrator> integrator_;
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephemeris_;
    
    double diff_step_ = 1e-7;            ///< Numerical differentiation step (OrbFit default)
    bool use_numerical_jacobian_ = true; ///< Default to true to match OrbFit logic
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_STATE_TRANSITION_MATRIX_HPP
