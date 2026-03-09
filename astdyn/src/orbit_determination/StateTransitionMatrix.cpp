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
    
    // Convert to AU/Day for numeric stability during integration
    y0.segment<3>(0) = initial.position.to_eigen_si() / physics::Distance::from_au(1.0).to_m();
    y0.segment<3>(3) = initial.velocity.to_eigen_si() * 86400.0 / physics::Distance::from_au(1.0).to_m();
    
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            y0[6 + i * 6 + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    auto f_augmented = [this](double t_val, const Eigen::VectorXd& y) -> Eigen::VectorXd {
        const auto t_tdb = time::EpochTDB::from_mjd(t_val);
        Eigen::VectorXd state = y.head<6>();
        
        if (y.size() < 42) {
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
    
    const double au_m = physics::Distance::from_au(1.0).to_m();

    auto final_state = physics::CartesianStateTyped<Frame>::from_si(
        target_time,
        yf[0] * au_m, yf[1] * au_m, yf[2] * au_m,
        yf[3] * au_m / 86400.0, yf[4] * au_m / 86400.0, yf[5] * au_m / 86400.0,
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
    
    // Numerical differentiation step (diff_step_ is relative or absolute AU)
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
    
    // rho = r_obj - r_obs (Heliocentric vectors in GCRF)
    Eigen::Vector3d rho_vec(pos_gcrf.x_si() - observer_pos.x_si(),
                            pos_gcrf.y_si() - observer_pos.y_si(),
                            pos_gcrf.z_si() - observer_pos.z_si());
    
    double rho = rho_vec.norm();
    double rho2 = rho * rho;
    double rho_xy2 = rho_vec.x() * rho_vec.x() + rho_vec.y() * rho_vec.y();
    double cos_dec = std::sqrt(rho_xy2) / rho;

    // Handle singularity at poles
    if (rho_xy2 < 1e-12) rho_xy2 = 1e-12;

    Eigen::Matrix<double, 2, 3> d_radec_d_rho_gcrf;
    
    // Partially projected RA: d(alpha * cos(dec)) / d(rho)
    // alpha = atan2(y, x) -> d(alpha)/dx = -y/(x^2+y^2), d(alpha)/dy = x/(x^2+y^2)
    // d(alpha * cos(dec))/dx = (-y / rho_xy2) * cos(dec)
    d_radec_d_rho_gcrf(0, 0) = -rho_vec.y() / rho_xy2 * cos_dec;
    d_radec_d_rho_gcrf(0, 1) =  rho_vec.x() / rho_xy2 * cos_dec;
    d_radec_d_rho_gcrf(0, 2) =  0.0;
    
    // Dec: sin(dec) = z/rho -> cos(dec) d(dec) = d(z/rho)
    // d(dec)/dx = -x*z / (rho^2 * sqrt(rho^2 - z^2))
    double sqrt_rho2_mz2 = std::max(1e-6, std::sqrt(rho2 - rho_vec.z() * rho_vec.z()));
    d_radec_d_rho_gcrf(1, 0) = -rho_vec.x() * rho_vec.z() / (rho2 * sqrt_rho2_mz2);
    d_radec_d_rho_gcrf(1, 1) = -rho_vec.y() * rho_vec.z() / (rho2 * sqrt_rho2_mz2);
    d_radec_d_rho_gcrf(1, 2) = sqrt_rho2_mz2 / rho2;
    
    // 2. Chain rule: ∂obs/∂x_frame = ∂obs/∂x_gcrf * ∂x_gcrf/∂x_frame
    // ∂x_gcrf/∂x_frame = R (rotation matrix)
    Eigen::Matrix3d rotation;
    if constexpr (std::is_same_v<Frame, core::GCRF>) {
        rotation = Eigen::Matrix3d::Identity();
    } else {
        rotation = coordinates::ReferenceFrame::get_rotation<Frame, core::GCRF>(state.epoch).to_eigen();
    }
    
    Eigen::Matrix<double, 2, 3> d_radec_d_pos_frame = d_radec_d_rho_gcrf * rotation;

    Eigen::Matrix<double, 2, 6> partial_radec_state_si;
    // Position partials: [rad/m]
    partial_radec_state_si.block<2, 3>(0, 0) = d_radec_d_pos_frame;
    
    // Velocity partials: 
    // The STM 'phi' is [dx_AU/dx0_AU, dx_AU/dv0_AUd, ...]
    // We need d_radec / dv0_si.
    // d_radec / dv0_si = (d_radec / dx_si) * (dx_si / dx_AU) * (dx_AU / dv0_AUd) * (dv0_AUd / dv0_si)
    // dv0_AUd / dv0_si = 86400.0 / (constants::AU * 1000.0)
    double au_to_m = constants::AU * 1000.0;
    double v_scale = 86400.0 / au_to_m;
    
    // We don't have direct access to velocity partials of RA/Dec (they are zero at epsilon=0),
    // but the chain rule for the STM velocity part is:
    // A_vel = (d_radec / dx_si) * (au_to_m * phi_pos_vel) * v_scale 
    //       = (d_radec / dx_si) * phi_pos_vel * 86400.0
    // Actually, when we multiply A * Phi later in DifferentialCorrector,
    // we just need to ensure A is consistent.
    
    // Let's simplify: DifferentialCorrector does: A_obs = partials.partial_radec * cumulative_phi;
    // If we want correction in SI:
    // partial_radec_state_si should be [ d_radec/dx_si, d_radec/dv_si ] at time t.
    // Since RA/Dec only depends on position: d_radec/dv_si = 0.
    partial_radec_state_si.block<2, 3>(0, 3).setZero();
    
    return partial_radec_state_si;
}

// Explicit instantiations
template class StateTransitionMatrix<core::GCRF>;
template class StateTransitionMatrix<core::ECLIPJ2000>;

} // namespace astdyn::orbit_determination
