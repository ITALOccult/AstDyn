/**
 * @file CovariancePropagator.hpp
 * @brief Propagates covariance matrix P using STM (Φ P Φᵀ).
 */

#ifndef ASTDYN_ORBIT_DETERMINATION_COVARIANCE_PROPAGATOR_HPP
#define ASTDYN_ORBIT_DETERMINATION_COVARIANCE_PROPAGATOR_HPP

#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include <Eigen/Dense>

namespace astdyn::orbit_determination {

/**
 * @brief Manages propagation of orbital uncertainty (Covariance).
 * @tparam Frame Reference frame of the state and covariance.
 */
template <typename Frame = core::GCRF>
class CovariancePropagator {
public:
    using Matrix6d = Eigen::Matrix<double, 6, 6>;

    /**
     * @brief Constructor
     * @param propagator The shared propagator model.
     */
    explicit CovariancePropagator(std::shared_ptr<propagation::Propagator> propagator)
        : stm_engine_(propagator) {}

    /**
     * @brief Set the initial state and covariance at t0.
     */
    void set_initial(const physics::CartesianStateTyped<Frame>& state, 
                    const Matrix6d& covariance) {
        initial_state_ = state;
        current_state_ = state;
        initial_covariance_ = covariance;
        current_covariance_ = covariance;
        current_stm_ = Matrix6d::Identity();
    }

    /**
     * @brief Propagate state and covariance to target time.
     * @param target_time Target epoch (TDB).
     */
    void propagate(time::EpochTDB target_time) {
        // 1. Compute state and STM from t0 (initial) to t (target)
        auto stm_res = stm_engine_.compute(initial_state_, target_time);
        
        // 2. Update internal state
        current_state_ = stm_res.final_state;
        current_stm_ = stm_res.phi;
        
        // 3. Map covariance: P(t) = Φ(t, t0) * P(t0) * Φ(t, t0)ᵀ
        current_covariance_ = stm_res.map_covariance(initial_covariance_);
    }

    /**
     * @brief Get the current propagated covariance matrix.
     */
    Matrix6d get_covariance() const { return current_covariance_; }

    /**
     * @brief Get the current propagated state.
     */
    physics::CartesianStateTyped<Frame> get_state() const { return current_state_; }

    /**
     * @brief Get the current State Transition Matrix Φ(t, t0).
     */
    Matrix6d get_stm() const { return current_stm_; }

private:
    StateTransitionMatrix<Frame> stm_engine_;
    physics::CartesianStateTyped<Frame> initial_state_;
    physics::CartesianStateTyped<Frame> current_state_;
    Matrix6d initial_covariance_;
    Matrix6d current_covariance_;
    Matrix6d current_stm_;
};

// Alias for common GCRF case
using CovariancePropagatorGCRF = CovariancePropagator<core::GCRF>;

} // namespace astdyn::orbit_determination

#endif // ASTDYN_ORBIT_DETERMINATION_COVARIANCE_PROPAGATOR_HPP
