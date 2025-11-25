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
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include <memory>

namespace astdyn::orbit_determination {

/**
 * @brief State transition matrix result
 */
struct STMResult {
    astdyn::Matrix6d phi;                        ///< Φ(t,t₀) - 6x6 state transition matrix
    astdyn::propagation::CartesianElements final_state;       ///< Propagated state at time t
    
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
 * Integrates both the state equations and variational equations:
 * 
 * dx/dt = f(x, t)
 * dΦ/dt = A(t) * Φ
 * 
 * where A(t) = ∂f/∂x is the Jacobian of the dynamics.
 * 
 * Initial condition: Φ(t₀,t₀) = I (identity matrix)
 */
class StateTransitionMatrix {
public:
    /**
     * @brief Constructor
     * 
     * @param propagator Orbit propagator (provides force model)
     */
    explicit StateTransitionMatrix(std::shared_ptr<astdyn::propagation::Propagator> propagator);
    
    /**
     * @brief Compute STM from t₀ to t
     * 
     * @param initial Initial Cartesian state at t₀
     * @param target_mjd_tdb Target time t
     * @return STM result with Φ(t,t₀) and final state
     */
    STMResult compute(
        const astdyn::propagation::CartesianElements& initial,
        double target_mjd_tdb);
    
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
        astdyn::propagation::CartesianElements final_state;
    };
    
    ObservationPartials compute_with_partials(
        const astdyn::propagation::CartesianElements& initial,
        double target_mjd_tdb,
        const astdyn::Vector3d& observer_pos);
    
    /**
     * @brief Set integrator for variational equations
     */
    void set_integrator(std::shared_ptr<astdyn::propagation::Integrator> integrator) {
        integrator_ = integrator;
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
    astdyn::Matrix6d compute_jacobian(double t, const Eigen::VectorXd& state);
    
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
    STMResult propagate_with_stm(
        const astdyn::propagation::CartesianElements& initial,
        double target_mjd_tdb);
    
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
        const astdyn::propagation::CartesianElements& state,
        const astdyn::Vector3d& observer_pos) const;

private:
    std::shared_ptr<astdyn::propagation::Propagator> propagator_;
    std::shared_ptr<astdyn::propagation::Integrator> integrator_;
    
    double diff_step_ = 1e-8;            ///< Numerical differentiation step
};

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_STATE_TRANSITION_MATRIX_HPP
