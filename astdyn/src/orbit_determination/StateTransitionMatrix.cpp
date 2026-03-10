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

template <typename Frame>
StateTransitionMatrix<Frame>::StateTransitionMatrix(
    std::shared_ptr<propagation::Propagator> propagator)
    : propagator_(propagator) {
    
    // Always use an explicit RKF78 integrator for the augmented ODE (6 state + 36 STM = 42
    // components).
    integrator_ = std::make_shared<RKF78Integrator<Eigen::VectorXd>>(
        0.05,    // initial step: 0.05 day
        1e-12,   // relative tolerance
        1e-12    // absolute tolerance (min_step)
    );
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
    y0.segment<3>(3) = initial.velocity.to_eigen_si() / physics::Velocity::from_au_d(1.0).to_ms();
    
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            y0[6 + i * 6 + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    auto f_augmented = [this](time::EpochTDB t_tdb, const Eigen::VectorXd& y) -> Eigen::VectorXd {
        StateAU state_au;
        state_au.epoch = t_tdb;
        state_au.x = y[0]; state_au.y = y[1]; state_au.z = y[2];
        state_au.vx = y[3]; state_au.vy = y[4]; state_au.vz = y[5];
        state_au.gm = physics::GravitationalParameter::from_au3_d2(propagator_->settings().central_body_gm);

        if (y.size() < 42) {
            auto dot = propagator_->compute_derivatives_au(t_tdb, state_au);
            Eigen::VectorXd res(6);
            res << dot.dx.to_au_d(), dot.dy.to_au_d(), dot.dz.to_au_d(),
                   dot.dvx.to_au_d2(), dot.dvy.to_au_d2(), dot.dvz.to_au_d2();
            return res;
        }

        astdyn::Matrix6d phi;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                phi(i, j) = y[6 + i * 6 + j];
            }
        }
        
        auto dot = propagator_->compute_derivatives_au(t_tdb, state_au);
        Eigen::VectorXd f_state(6);
        f_state << dot.dx.to_au_d(), dot.dy.to_au_d(), dot.dz.to_au_d(),
                   dot.dvx.to_au_d2(), dot.dvy.to_au_d2(), dot.dvz.to_au_d2();
                   
        astdyn::Matrix6d A = compute_jacobian(t_tdb, y.head<6>());
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
    
    const double au_m = constants::AU * 1000.0;
    const double aud_ms = au_m / 86400.0;
    
    auto final_state = physics::CartesianStateTyped<Frame>::from_si(
        target_time,
        yf[0] * au_m, yf[1] * au_m, yf[2] * au_m,
        yf[3] * aud_ms, yf[4] * aud_ms, yf[5] * aud_ms,
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
    
    // Top-right 3x3: Velocity identity
    A.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    
    // Bottom-left 3x3: Acceleration-Position partials (Analytical for Central Body)
    Eigen::Vector3d r = state.head<3>();
    double r_mag = r.norm();
    double mu = propagator_->settings().central_body_gm;
    
    // Grad(a) = -mu/r^3 * [I - 3(r*r^T)/r^2]
    A.block<3, 3>(3, 0) = compute_acceleration_position_partial(r, mu);

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
    
    // rho = r_obj - r_obs (Heliocentric vectors in GCRF, in AU)
    const double au_m = constants::AU * 1000.0;
    Eigen::Vector3d rho_vec(pos_gcrf.x_si() / au_m - observer_pos.x_si() / au_m,
                            pos_gcrf.y_si() / au_m - observer_pos.y_si() / au_m,
                            pos_gcrf.z_si() / au_m - observer_pos.z_si() / au_m);
    
    double rho2 = rho_vec.squaredNorm();
    double rho = std::sqrt(rho2);
    double rho_xy2 = rho_vec.x() * rho_vec.x() + rho_vec.y() * rho_vec.y();
    double cos_dec = std::sqrt(rho_xy2) / rho;

    // Handle singularity at poles
    if (rho_xy2 < 1e-12) rho_xy2 = 1e-12;

    double rho_xy = std::sqrt(rho_xy2);
    
    Eigen::Matrix<double, 2, 3> d_radec_d_rho_gcrf;
    
    // RA: ra = atan2(y, x).  d(ra) = (x dy - y dx) / (x^2 + y^2) [rad/AU]
    // Properly project RA partial with cos(dec) to match residuals: d(RA*cos(dec))/dx
    d_radec_d_rho_gcrf(0, 0) = (-rho_vec.y() / rho_xy2) * cos_dec;
    d_radec_d_rho_gcrf(0, 1) = ( rho_vec.x() / rho_xy2) * cos_dec;
    d_radec_d_rho_gcrf(0, 2) =  0.0;
    
    // Dec: dec = asin(z/rho). d(dec) = [ (x^2+y^2).dz - z(x.dx + y.dy) ] / [ rho^2 * sqrt(x^2+y^2) ]
    d_radec_d_rho_gcrf(1, 0) = -rho_vec.x() * rho_vec.z() / (rho2 * rho_xy);
    d_radec_d_rho_gcrf(1, 1) = -rho_vec.y() * rho_vec.z() / (rho2 * rho_xy);
    d_radec_d_rho_gcrf(1, 2) = rho_xy / rho2;
    
    // 2. Chain rule: ∂obs/∂x_frame = ∂obs/∂x_gcrf * ∂x_gcrf/∂x_frame
    // ∂x_gcrf/∂x_frame = R (rotation matrix)
    Eigen::Matrix3d rotation;
    if constexpr (std::is_same_v<Frame, core::GCRF>) {
        rotation = Eigen::Matrix3d::Identity();
    } else {
        rotation = coordinates::ReferenceFrame::get_rotation<Frame, core::GCRF>(state.epoch).to_eigen();
    }
    
    Eigen::Matrix<double, 2, 3> d_radec_d_pos_frame = d_radec_d_rho_gcrf * rotation;

    Eigen::Matrix<double, 2, 6> partial_radec_state_au;
    // Position partials: [rad/AU]
    // rotation is [AU_frame/AU_gcrf], d_radec_d_rho_gcrf is [rad/AU] if we use AU distances
    partial_radec_state_au.block<2, 3>(0, 0) = d_radec_d_pos_frame;
    
    // Velocity partials: zero for pure geometric RA/Dec
    partial_radec_state_au.block<2, 3>(0, 3).setZero();
    
    return partial_radec_state_au;
}

template <typename Frame>
std::vector<typename StateTransitionMatrix<Frame>::ObservationPartials> StateTransitionMatrix<Frame>::compute_batch(
    const physics::CartesianStateTyped<Frame>& initial,
    const std::vector<std::pair<time::EpochTDB, math::Vector3<core::GCRF, physics::Distance>>>& observations) 
{
    std::vector<ObservationPartials> results;
    results.reserve(observations.size());

    if (observations.empty()) return results;

    // 1. Initial augmented state [x(6), phi(36)]
    // Integration is done in AU and AU/day for numerical stability
    Eigen::VectorXd y = Eigen::VectorXd::Zero(42);
    
    // Explicit conversion from SI to AU/day
    const double au_m = constants::AU * 1000.0;
    const double aud_ms = au_m / 86400.0;
    
    auto state_si = initial.to_eigen_si();
    y.template head<3>() = state_si.template head<3>() / au_m;
    y.template segment<3>(3) = state_si.template tail<3>() / aud_ms;
    
    // Identity mapping for I(t0)
    for (int i = 0; i < 6; ++i) y[6 + i * 6 + i] = 1.0;

    double t_current = initial.epoch.mjd();

    // Augmented ODE system: dy/dt = [v, a, dPhi/dt]
    auto f_augmented = [&](time::EpochTDB t, const Eigen::VectorXd& y_aug) -> Eigen::VectorXd {
        Eigen::VectorXd dy = Eigen::VectorXd::Zero(42);
        
        StateAU state_au;
        state_au.epoch = t;
        state_au.x = y_aug[0]; state_au.y = y_aug[1]; state_au.z = y_aug[2];
        state_au.vx = y_aug[3]; state_au.vy = y_aug[4]; state_au.vz = y_aug[5];
        state_au.gm = physics::GravitationalParameter::from_au3_d2(propagator_->settings().central_body_gm);

        // State derivatives [v, a]
        auto dot = propagator_->compute_derivatives_au(t, state_au);
        dy[0] = dot.dx.to_au_d(); dy[1] = dot.dy.to_au_d(); dy[2] = dot.dz.to_au_d();
        dy[3] = dot.dvx.to_au_d2(); dy[4] = dot.dvy.to_au_d2(); dy[5] = dot.dvz.to_au_d2();
        
        // Phi derivatives: dPhi/dt = A(t) * Phi
        astdyn::Matrix6d Phi;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) Phi(i, j) = y_aug[6 + i * 6 + j];
        }
        
        astdyn::Matrix6d A = compute_jacobian(t, y_aug.template head<6>());
        astdyn::Matrix6d Phi_dot = A * Phi;
        
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) dy[6 + i * 6 + j] = Phi_dot(i, j);
        }
        
        return dy;
    };

    // 2. Continuous Propagation through all sorted observations
    std::vector<double> target_times;
    target_times.reserve(observations.size());
    for (const auto& obs : observations) {
        target_times.push_back(obs.first.mjd());
    }

    // Single batch integration call - allows the integrator (e.g. Adams) 
    // to maintain internal state/history across all points.
    auto states_batch = integrator_->integrate_batch(f_augmented, y, t_current, target_times);

    // 3. Extract results
    for (size_t i = 0; i < observations.size(); ++i) {
        const auto& obs = observations[i];
        const auto& y_res = states_batch[i];
        
        ObservationPartials res;
        
        // 5. Final state in SI: Position [m], Velocity [m/s]
        const double au_m = constants::AU * 1000.0;
        const double aud_ms = au_m / 86400.0;
        
        res.final_state = physics::CartesianStateTyped<Frame>::from_si(
            obs.first,
            y_res[0] * au_m, y_res[1] * au_m, y_res[2] * au_m,
            y_res[3] * aud_ms, y_res[4] * aud_ms, y_res[5] * aud_ms,
            initial.gm.to_m3_s2());
        
        // STM Units: The integrator works in AU and AU/day.
        // res.phi is used directly in AU units for normalized design matrix
        for (int r = 0; r < 6; ++r) {
            for (int r2 = 0; r2 < 6; ++r2) res.phi(r, r2) = y_res[6 + r * 6 + r2];
        }
        
        // RA/Dec partials w.r.t integration frame (now specifically in rad/AU)
        res.partial_radec = compute_observation_partials(res.final_state, obs.second);
        results.push_back(res);
    }

    return results;
}

// Explicit instantiations
template class StateTransitionMatrix<core::GCRF>;
template class StateTransitionMatrix<core::ECLIPJ2000>;

} // namespace astdyn::orbit_determination
