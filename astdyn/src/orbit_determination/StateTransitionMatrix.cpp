/**
 * @file StateTransitionMatrix.cpp
 * @brief Implementation of state transition matrix computation
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include <stdexcept>

namespace astdyn::orbit_determination {

using namespace astdyn::propagation;

// ============================================================================
// StateTransitionMatrix Implementation
// ============================================================================

template <typename Frame>
StateTransitionMatrix<Frame>::StateTransitionMatrix(
    std::shared_ptr<propagation::Propagator> propagator)
    : propagator_(propagator) {
    
    if (propagator_) {
        integrator_ = propagator_->get_integrator();
    }
    
    if (!integrator_) {
        integrator_ = std::make_shared<RKF78Integrator>(
            0.1,      
            1e-12     
        );
    }
}

template <typename Frame>
STMResult<Frame> StateTransitionMatrix<Frame>::compute(
    const physics::CartesianStateTyped<Frame>& initial,
    time::EpochTDB target_time) {
    
    return propagate_with_stm(initial, target_time);
}

template <typename Frame>
typename StateTransitionMatrix<Frame>::ObservationPartials 
StateTransitionMatrix<Frame>::compute_with_partials(
    const physics::CartesianStateTyped<Frame>& initial,
    time::EpochTDB target_time,
    const math::Vector3<core::GCRF, physics::Distance>& observer_pos) {
    
    auto stm_result = propagate_with_stm(initial, target_time);
    auto obs_partials = compute_observation_partials(
        stm_result.final_state, observer_pos);
    
    return {stm_result.phi, obs_partials, stm_result.final_state};
}

template <typename Frame>
STMResult<Frame> StateTransitionMatrix<Frame>::propagate_with_stm(
    const physics::CartesianStateTyped<Frame>& initial,
    time::EpochTDB target_time) {
    
    Eigen::VectorXd y0(42);
    
    y0.segment<3>(0) = initial.position.to_eigen_si() / (constants::AU * 1000.0);
    y0.segment<3>(3) = initial.velocity.to_eigen_si() * 86400.0 / (constants::AU * 1000.0);
    
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            y0[6 + i * 6 + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    auto f_augmented = [this](double t_val, const Eigen::VectorXd& y) -> Eigen::VectorXd {
        const auto t_tdb = time::EpochTDB::from_mjd(t_val);
        Eigen::VectorXd state = y.head<6>();
        
        if (y.size() < 42) {
            // Internal integrator step (like AAS symplectic kick) only needs the state derivatives
            return propagator_->compute_derivatives(t_tdb, state);
        }

        astdyn::Matrix6d phi;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                phi(i, j) = y[6 + i * 6 + j];
            }
        }
        
        Eigen::VectorXd f_state = propagator_->compute_derivatives(t_tdb, state);
        astdyn::Matrix6d A = compute_jacobian(t_tdb, state);
        astdyn::Matrix6d phi_dot = A * phi;
        
        Eigen::VectorXd dy(42);
        dy.segment<6>(0) = f_state;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                dy[6 + i * 6 + j] = phi_dot(i, j);
            }
        }
        return dy;
    };
    
    Eigen::VectorXd yf = integrator_->integrate(
        f_augmented, y0, initial.epoch.mjd(), target_time.mjd());
    
    auto final_state = physics::CartesianStateTyped<Frame>::from_si(
        target_time,
        yf[0] * constants::AU * 1000.0, yf[1] * constants::AU * 1000.0, yf[2] * constants::AU * 1000.0,
        yf[3] * constants::AU * 1000.0 / 86400.0, yf[4] * constants::AU * 1000.0 / 86400.0, yf[5] * constants::AU * 1000.0 / 86400.0,
        initial.gm.to_m3_s2()
    );
    
    astdyn::Matrix6d phi;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            phi(i, j) = yf[6 + i * 6 + j];
        }
    }
    
    return {phi, final_state};
}

template <typename Frame>
astdyn::Matrix6d StateTransitionMatrix<Frame>::compute_jacobian(time::EpochTDB t, const Eigen::VectorXd& state) {
    astdyn::Matrix6d A = astdyn::Matrix6d::Zero();
    A.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    
    const double h = diff_step_;
    Eigen::VectorXd state_plus = state;
    Eigen::VectorXd state_minus = state;
    
    for (int i = 0; i < 6; ++i) {
        state_plus[i] += h;
        state_minus[i] -= h;
        
        Eigen::VectorXd f_plus = propagator_->compute_derivatives(t, state_plus);
        Eigen::VectorXd f_minus = propagator_->compute_derivatives(t, state_minus);
        
        A.block<3, 1>(3, i) = (f_plus.segment<3>(3) - f_minus.segment<3>(3)) / (2.0 * h);
        
        state_plus[i] = state[i];
        state_minus[i] = state[i];
    }
    
    return A;
}

template <typename Frame>
Eigen::Matrix3d StateTransitionMatrix<Frame>::compute_acceleration_position_partial(
    const astdyn::Vector3d& r,
    double mu) const {
    
    double r_mag = r.norm();
    double r3 = r_mag * r_mag * r_mag;
    double r5 = r3 * r_mag * r_mag;
    
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d r_outer_r = r * r.transpose();
    
    return (-mu / r3) * I + (3.0 * mu / r5) * r_outer_r;
}

template <typename Frame>
Eigen::Matrix<double, 2, 6> StateTransitionMatrix<Frame>::compute_observation_partials(
    const physics::CartesianStateTyped<Frame>& state,
    const math::Vector3<core::GCRF, physics::Distance>& observer_pos) const {
    
    // 1. Must use GCRF for RA/Dec partials
    auto pos_gcrf = coordinates::ReferenceFrame::transform_pos<Frame, core::GCRF>(state.position, state.epoch);
    
    Eigen::Vector3d rho(pos_gcrf.x_si() - observer_pos.x_si(),
                        pos_gcrf.y_si() - observer_pos.y_si(),
                        pos_gcrf.z_si() - observer_pos.z_si());
    double range = rho.norm();
    
    double x = rho[0] / range;
    double y = rho[1] / range;
    double z = rho[2] / range;
    
    double rho_xy2 = x * x + y * y;
    double sqrt_1mz2 = std::sqrt(1.0 - z * z);
    
    Eigen::Matrix<double, 2, 3> partial_radec_rho_gcrf;
    partial_radec_rho_gcrf(0, 0) = -y / (range * rho_xy2);
    partial_radec_rho_gcrf(0, 1) = x / (range * rho_xy2);
    partial_radec_rho_gcrf(0, 2) = 0.0;
    
    partial_radec_rho_gcrf(1, 0) = -x * z / (range * sqrt_1mz2);
    partial_radec_rho_gcrf(1, 1) = -y * z / (range * sqrt_1mz2);
    partial_radec_rho_gcrf(1, 2) = sqrt_1mz2 / range;
    
    // 2. Chain rule: ∂obs/∂x_frame = ∂obs/∂x_gcrf * ∂x_gcrf/∂x_frame
    // ∂x_gcrf/∂x_frame = R (rotation matrix)
    Eigen::Matrix3d rotation;
    if constexpr (std::is_same_v<Frame, core::GCRF>) {
        rotation = Eigen::Matrix3d::Identity();
    } else {
        rotation = coordinates::ReferenceFrame::get_rotation<Frame, core::GCRF>(state.epoch).to_eigen();
    }
    
    Eigen::Matrix<double, 2, 3> partial_radec_rho_frame = partial_radec_rho_gcrf * rotation;

    Eigen::Matrix<double, 2, 6> partial_radec_state_au;
    // AU scaling
    partial_radec_state_au.block<2, 3>(0, 0) = partial_radec_rho_frame * (constants::AU * 1000.0);
    partial_radec_state_au.block<2, 3>(0, 3).setZero();
    
    return partial_radec_state_au;
}

// Explicit instantiations
template class StateTransitionMatrix<core::GCRF>;
template class StateTransitionMatrix<core::ECLIPJ2000>;

} // namespace astdyn::orbit_determination
