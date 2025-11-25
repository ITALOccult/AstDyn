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
    const CartesianElements& initial,
    double target_mjd_tdb) {
    
    return propagate_with_stm(initial, target_mjd_tdb);
}

StateTransitionMatrix::ObservationPartials 
StateTransitionMatrix::compute_with_partials(
    const CartesianElements& initial,
    double target_mjd_tdb,
    const Vector3d& observer_pos) {
    
    // Compute STM
    auto stm_result = propagate_with_stm(initial, target_mjd_tdb);
    
    // Compute observation partials ∂(RA,Dec)/∂x
    auto obs_partials = compute_observation_partials(
        stm_result.final_state, observer_pos);
    
    ObservationPartials result;
    result.phi = stm_result.phi;
    result.partial_radec = obs_partials;
    result.final_state = stm_result.final_state;
    
    return result;
}

STMResult StateTransitionMatrix::propagate_with_stm(
    const CartesianElements& initial,
    double target_mjd_tdb) {
    
    // Augmented state vector: [x(6), Φ(36)]
    // Total dimension: 42
    Eigen::VectorXd y0(42);
    
    // Initial state
    y0.segment<3>(0) = initial.position;
    y0.segment<3>(3) = initial.velocity;
    
    // Initial STM: Φ(t₀,t₀) = I
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            int idx = 6 + i * 6 + j;
            y0[idx] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // Define augmented derivative function
    auto f_augmented = [this](double t, const Eigen::VectorXd& y) -> Eigen::VectorXd {
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
        Eigen::VectorXd f_state = propagator_->compute_derivatives(t, state);
        
        // Compute Jacobian A = ∂f/∂x
        Matrix6d A = compute_jacobian(t, state);
        
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
    Eigen::VectorXd yf = integrator_->integrate(
        f_augmented, y0, initial.epoch_mjd_tdb, target_mjd_tdb);
    
    // Extract final state
    STMResult result;
    result.final_state.epoch_mjd_tdb = target_mjd_tdb;
    result.final_state.gravitational_parameter = initial.gravitational_parameter;
    result.final_state.position = yf.segment<3>(0);
    result.final_state.velocity = yf.segment<3>(3);
    
    // Extract STM
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            result.phi(i, j) = yf[6 + i * 6 + j];
        }
    }
    
    return result;
}

Matrix6d StateTransitionMatrix::compute_jacobian(
    double t, 
    const Eigen::VectorXd& state) {
    
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
    const CartesianElements& state,
    const Vector3d& observer_pos) const {
    
    // Topocentric position vector
    Vector3d rho = state.position - observer_pos;
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
    
    // Full partials: ∂(RA,Dec)/∂(r,v)
    // Only position affects direction, not velocity
    Eigen::Matrix<double, 2, 6> partial_radec_state;
    partial_radec_state.block<2, 3>(0, 0) = partial_radec_rho;
    partial_radec_state.block<2, 3>(0, 3).setZero();
    
    return partial_radec_state;
}

} // namespace astdyn::orbit_determination
