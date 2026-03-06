/**
 * @file StateTransitionMatrix.cpp
 * @brief Implementation of state transition matrix computation
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/core/Constants.hpp"
#include <stdexcept>

namespace astdyn::orbit_determination {

using namespace astdyn::propagation;

// ============================================================================
// StateTransitionMatrix Implementation
// ============================================================================

StateTransitionMatrix::StateTransitionMatrix(
    std::shared_ptr<propagation::Propagator> propagator)
    : propagator_(propagator) {
    
    // Create default integrator if not set
    if (!integrator_) {
        integrator_ = std::make_shared<RKF78Integrator>(
            0.1,      // initial step
            1e-12,    // tolerance
            1e-6,     // min step
            10.0      // max step
        );
    }
}

STMResult StateTransitionMatrix::compute(
    const physics::CartesianStateTyped<core::GCRF>& initial,
    time::EpochTDB target_time) {
    
    return propagate_with_stm(initial, target_time);
}

StateTransitionMatrix::ObservationPartials 
StateTransitionMatrix::compute_with_partials(
    const physics::CartesianStateTyped<core::GCRF>& initial,
    time::EpochTDB target_time,
    const math::Vector3<core::GCRF, physics::Distance>& observer_pos) {
    
    // Compute STM
    auto stm_result = propagate_with_stm(initial, target_time);
    
    // Compute observation partials ∂(RA,Dec)/∂x
    auto obs_partials = compute_observation_partials(
        stm_result.final_state, observer_pos);
    
    ObservationPartials result{stm_result.phi, obs_partials, stm_result.final_state};
    return result;
}

STMResult StateTransitionMatrix::propagate_with_stm(
    const physics::CartesianStateTyped<core::GCRF>& initial,
    time::EpochTDB target_time) {
    
    // Augmented state vector: [x(6), Φ(36)]
    // Total dimension: 42
    Eigen::VectorXd y0(42);
    
    // Initial state in AU and AU/day
    y0.segment<3>(0) = initial.position.to_eigen_si() / (constants::AU * 1000.0);
    y0.segment<3>(3) = initial.velocity.to_eigen_si() * 86400.0 / (constants::AU * 1000.0);
    
    // Initial STM: Φ(t₀,t₀) = I
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            int idx = 6 + i * 6 + j;
            y0[idx] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // Define augmented derivative function
    auto f_augmented = [this](double t_val, const Eigen::VectorXd& y) -> Eigen::VectorXd {
        const auto t_tdb = time::EpochTDB::from_mjd(t_val);
        // Extract state
        Eigen::VectorXd state = y.segment<6>(0);
        
        // Extract STM (as 6x6 matrix)
        Matrix6d phi;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                phi(i, j) = y[6 + i * 6 + j];
            }
        }
        
        // Compute state derivative dx/dt = f(x,t)
        Eigen::VectorXd f_state = propagator_->compute_derivatives(t_tdb, state);
        
        // Compute Jacobian A = ∂f/∂x
        Matrix6d A = compute_jacobian(t_tdb, state);
        
        // STM derivative: dΦ/dt = A * Φ
        Matrix6d phi_dot = A * phi;
        
        // Pack into augmented derivative
        Eigen::VectorXd dy(42);
        dy.segment<6>(0) = f_state;
        
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                dy[6 + i * 6 + j] = phi_dot(i, j);
            }
        }
        
        return dy;
    };
    
    // Integrate augmented system
    // Integrate augmented system
    Eigen::VectorXd yf = integrator_->integrate(
        f_augmented, y0, initial.epoch.mjd(), target_time.mjd());
    
    // Extract final state
    auto target_time_tdb = target_time;
    auto final_state = physics::CartesianStateTyped<core::GCRF>::from_si(
        target_time_tdb,
        yf[0] * constants::AU * 1000.0, yf[1] * constants::AU * 1000.0, yf[2] * constants::AU * 1000.0,
        yf[3] * constants::AU * 1000.0 / 86400.0, yf[4] * constants::AU * 1000.0 / 86400.0, yf[5] * constants::AU * 1000.0 / 86400.0,
        initial.gm.to_m3_s2()
    );
    
    // Extract STM
    astdyn::Matrix6d phi;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            phi(i, j) = yf[6 + i * 6 + j];
        }
    }
    
    STMResult result{phi, final_state};
    return result;
}

Matrix6d StateTransitionMatrix::compute_jacobian(time::EpochTDB t, const Eigen::VectorXd& state) {
    
    Vector3d r = state.segment<3>(0);
    Vector3d v = state.segment<3>(3);
    
    double mu = propagator_->settings().central_body_gm;
    
    Matrix6d A = Matrix6d::Zero();
    
    // Upper-right block: ∂(dx/dt)/∂v = I
    A.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    
    // Lower-left block: ∂(dv/dt)/∂r = ∂a/∂r
    A.block<3, 3>(3, 0) = compute_acceleration_position_partial(r, mu);
    
    // Lower-right block: ∂(dv/dt)/∂v
    // For two-body problem, acceleration doesn't depend on velocity
    // For drag/radiation pressure, would have ∂a/∂v terms
    A.block<3, 3>(3, 3) = Eigen::Matrix3d::Zero();
    
    return A;
}

Eigen::Matrix3d StateTransitionMatrix::compute_acceleration_position_partial(
    const Vector3d& r,
    double mu) const {
    
    double r_mag = r.norm();
    double r3 = r_mag * r_mag * r_mag;
    double r5 = r3 * r_mag * r_mag;
    
    // For a = -μr/r³:
    // ∂a/∂r = -μ/r³ [I - 3(r⊗r)/r²]
    //       = -μ/r³ I + 3μ(r⊗r)/r⁵
    
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d r_outer_r = r * r.transpose();
    
    Eigen::Matrix3d partial = (-mu / r3) * I + (3.0 * mu / r5) * r_outer_r;
    
    return partial;
}

Eigen::Matrix<double, 2, 6> StateTransitionMatrix::compute_observation_partials(
    const physics::CartesianStateTyped<core::GCRF>& state,
    const math::Vector3<core::GCRF, physics::Distance>& observer_pos) const {
    
    // Topocentric position vector in meters
    Eigen::Vector3d rho(state.position.x_si() - observer_pos.x_si(),
                        state.position.y_si() - observer_pos.y_si(),
                        state.position.z_si() - observer_pos.z_si());
    double range = rho.norm();
    
    // Unit direction vector
    double x = rho[0] / range;
    double y = rho[1] / range;
    double z = rho[2] / range;
    
    // Partials of RA and Dec w.r.t. topocentric Cartesian coordinates
    // 
    // RA = atan2(y, x)
    // Dec = asin(z/r) = asin(z)  [since r=1 for unit vector]
    //
    // ∂RA/∂x = -y/(x²+y²)
    // ∂RA/∂y = x/(x²+y²)
    // ∂RA/∂z = 0
    //
    // ∂Dec/∂x = -xz/sqrt(1-z²)
    // ∂Dec/∂y = -yz/sqrt(1-z²)
    // ∂Dec/∂z = sqrt(1-z²)
    
    double rho_xy2 = x * x + y * y;
    double sqrt_1mz2 = std::sqrt(1.0 - z * z);
    
    Eigen::Matrix<double, 2, 3> partial_radec_rho;
    
    // ∂RA/∂ρ (first row)
    partial_radec_rho(0, 0) = -y / (range * rho_xy2);  // ∂RA/∂x
    partial_radec_rho(0, 1) = x / (range * rho_xy2);   // ∂RA/∂y
    partial_radec_rho(0, 2) = 0.0;                      // ∂RA/∂z
    
    // ∂Dec/∂ρ (second row)
    partial_radec_rho(1, 0) = -x * z / (range * sqrt_1mz2); // ∂Dec/∂x
    partial_radec_rho(1, 1) = -y * z / (range * sqrt_1mz2); // ∂Dec/∂y
    partial_radec_rho(1, 2) = sqrt_1mz2 / range;             // ∂Dec/∂z
    
    // ∂ρ/∂r_helio = I (topocentric = heliocentric - observer)
    // Observer position is fixed at observation time
    
    // State is in SI internally (for external APIs), but phi relates variations in AU!
    // We want partials of RA, Dec w.r.t initial SI positions/velocities?
    // Wait... if STM `phi` relates [x, y, z, vx, vy, vz]_AU to [x0, y0, z0, vx0, vy0, vz0]_AU
    // But `compute_with_partials` needs to output `partial_radec` that multiplies `phi`.
    // It's much simpler if the differential corrector keeps taking `state` in SI...
    // Actually, Differential Corrector handles its own units via STM, wait.
    // If we multiply A_obs = partial_radec_state * phi, 
    // partial_radec_state has units rad / AU. So A_obs is rad / AU(0).
    // And `correction` will naturally be in AU / AU/d. 
    // We must apply a chain rule if we want SI. Let's provide partials w.r.t AU, 
    // so the DifferentialCorrector works entirely with AU for `correction`.

    Eigen::Matrix<double, 2, 6> partial_radec_state_au;
    // ρ is in meters, so partial_radec_rho is rad / meter.
    // We multiply by meters per AU to get rad / AU.
    partial_radec_state_au.block<2, 3>(0, 0) = partial_radec_rho * (constants::AU * 1000.0);
    partial_radec_state_au.block<2, 3>(0, 3).setZero();
    
    return partial_radec_state_au;
}

} // namespace astdyn::orbit_determination
