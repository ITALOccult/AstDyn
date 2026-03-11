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
    
    // STM always uses a dedicated high-precision RKF78 integrator
    // to avoid interference with the propagator's integrator and ensure stability
    integrator_ = std::make_shared<RKF78Integrator>(
        0.1,      // initial step [days]
        1e-13     // high precision for variational equations
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
std::vector<typename StateTransitionMatrix<Frame>::ObservationPartials> 
StateTransitionMatrix<Frame>::compute_batch(
    const physics::CartesianStateTyped<Frame>& initial,
    const std::vector<time::EpochTDB>& target_times,
    const std::vector<math::Vector3<core::GCRF, physics::Distance>>& observer_positions) {
    
    if (target_times.empty()) return {};
    
    // 1. Bulk propagate STM (Φ(t_i, t_0))
    auto stm_results = propagate_with_stm_batch(initial, target_times);
    
    // 2. Compute observation partials at each point
    std::vector<ObservationPartials> results;
    results.reserve(target_times.size());
    
    for (size_t i = 0; i < target_times.size(); ++i) {
        auto obs_partials = compute_observation_partials(
            stm_results[i].final_state, observer_positions[i]);
        
        results.push_back({stm_results[i].phi, obs_partials, stm_results[i].final_state});
    }
    
    return results;
}

template <typename Frame>
STMResult<Frame> StateTransitionMatrix<Frame>::propagate_with_stm(
    const physics::CartesianStateTyped<Frame>& initial,
    time::EpochTDB target_time) {
    
    // Unified Strategy: Always integrate in ECLIPJ2000 AU/day internally
    // if propagator is configured for Ecliptic.
    const bool use_ecl = (propagator_ && propagator_->settings().integrate_in_ecliptic);
    
    Eigen::VectorXd y0(42);
    
    if (use_ecl) {
        auto initial_ecl = initial.template cast_frame<core::ECLIPJ2000>();
        y0.head<6>() = initial_ecl.to_eigen_au_aud();
    } else {
        y0.head<6>() = initial.to_eigen_au_aud();
    }
    
    // Initial identity for STM
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            y0[6 + i * 6 + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    auto f_augmented = [this](double t_val, const Eigen::VectorXd& y) -> Eigen::VectorXd {
        const auto t_tdb = time::EpochTDB::from_mjd(t_val);
        const Eigen::VectorXd state = y.head<6>();
        
        astdyn::Matrix6d phi;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) phi(i, j) = y[6 + i * 6 + j];
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
    
    physics::CartesianStateTyped<Frame> final_state;
    if (use_ecl) {
        auto final_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>::from_au_aud(
            target_time, yf.head<6>(), initial.gm);
        final_state = final_ecl.template cast_frame<Frame>();
    } else {
        final_state = physics::CartesianStateTyped<Frame>::from_au_aud(
            target_time, yf.head<6>(), initial.gm);
    }
    
    astdyn::Matrix6d phi_integrated;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            phi_integrated(i, j) = yf[6 + i * 6 + j];
        }
    }
    
    astdyn::Matrix6d phi = phi_integrated;
    if (use_ecl && !std::is_same_v<Frame, core::ECLIPJ2000>) {
        // STM was integrated in ECLIPJ2000. Frame is GCRF (typically).
        // Phi_GCRF = R_E2G * Phi_ECL * R_G2E
        auto rot_g2e = coordinates::ReferenceFrame::get_rotation<core::GCRF, core::ECLIPJ2000>().to_eigen();
        auto rot_e2g = rot_g2e.transpose();
        
        astdyn::Matrix6d R6 = astdyn::Matrix6d::Zero();
        R6.block<3,3>(0,0) = rot_g2e;
        R6.block<3,3>(3,3) = rot_g2e;
        
        astdyn::Matrix6d R6_inv = astdyn::Matrix6d::Zero();
        R6_inv.block<3,3>(0,0) = rot_e2g;
        R6_inv.block<3,3>(3,3) = rot_e2g;
        
        phi = R6_inv * phi_integrated * R6;
    }
    
    return {phi, final_state};
}

template <typename Frame>
std::vector<STMResult<Frame>> StateTransitionMatrix<Frame>::propagate_with_stm_batch(
    const physics::CartesianStateTyped<Frame>& initial,
    const std::vector<time::EpochTDB>& target_times) {
    
    const bool use_ecl = (propagator_ && propagator_->settings().integrate_in_ecliptic);
    
    Eigen::VectorXd y0(42);
    if (use_ecl) {
        auto initial_ecl = initial.template cast_frame<core::ECLIPJ2000>();
        y0.head<6>() = initial_ecl.to_eigen_au_aud();
    } else {
        y0.head<6>() = initial.to_eigen_au_aud();
    }
    
    // Initial identity for STM
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            y0[6 + i * 6 + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // Shared augmented derivative function
    auto f_augmented = [this](double t_val, const Eigen::VectorXd& y) -> Eigen::VectorXd {
        const auto t_tdb = time::EpochTDB::from_mjd(t_val);
        const Eigen::VectorXd state = y.head<6>();
        
        astdyn::Matrix6d phi;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) phi(i, j) = y[6 + i * 6 + j];
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
    
    std::vector<double> mjd_targets;
    mjd_targets.reserve(target_times.size());
    for (const auto& t : target_times) mjd_targets.push_back(t.mjd());
    
    std::vector<Eigen::VectorXd> y_results = integrator_->integrate_at(
        f_augmented, y0, initial.epoch.mjd(), mjd_targets);
    
    std::vector<STMResult<Frame>> results;
    results.reserve(y_results.size());
    
    const bool needs_rotation = (use_ecl && !std::is_same_v<Frame, core::ECLIPJ2000>);
    Eigen::Matrix3d rot_g2e, rot_e2g;
    astdyn::Matrix6d R6, R6_inv;
    
    if (needs_rotation) {
        rot_g2e = coordinates::ReferenceFrame::get_rotation<core::GCRF, core::ECLIPJ2000>().to_eigen();
        rot_e2g = rot_g2e.transpose();
        R6.setZero();
        R6.block<3,3>(0,0) = rot_g2e;
        R6.block<3,3>(3,3) = rot_g2e;
        R6_inv.setZero();
        R6_inv.block<3,3>(0,0) = rot_e2g;
        R6_inv.block<3,3>(3,3) = rot_e2g;
    }
    
    for (size_t i = 0; i < y_results.size(); ++i) {
        const auto& y_f = y_results[i];
        const auto& t_f = target_times[i];
        
        physics::CartesianStateTyped<Frame> final_state;
        if (use_ecl) {
            auto final_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>::from_au_aud(
                t_f, y_f.head<6>(), initial.gm);
            final_state = final_ecl.template cast_frame<Frame>();
        } else {
            final_state = physics::CartesianStateTyped<Frame>::from_au_aud(
                t_f, y_f.head<6>(), initial.gm);
        }
        
        astdyn::Matrix6d phi_integrated;
        for (int row = 0; row < 6; ++row) {
            for (int col = 0; col < 6; ++col) {
                phi_integrated(row, col) = y_f[6 + row * 6 + col];
            }
        }
        
        astdyn::Matrix6d phi = phi_integrated;
        if (needs_rotation) {
            phi = R6_inv * phi_integrated * R6;
        }

        results.push_back({phi, final_state});
    }
    
    return results;
}

template <typename Frame>
astdyn::Matrix6d StateTransitionMatrix<Frame>::compute_jacobian(time::EpochTDB t, const Eigen::VectorXd& state) {
    astdyn::Matrix6d A = astdyn::Matrix6d::Zero();
    A.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    
    const Eigen::Vector3d r_ast = state.head<3>();
    const double mu_sun = propagator_->settings().central_body_gm;
    
    // 1. Central body (Sun) partials
    A.block<3, 3>(3, 0) = compute_acceleration_position_partial(r_ast, mu_sun);
    
    // 2. Planetary perturbations partials (Analytical)
    if (propagator_->settings().include_planets) {
        auto provider = ephemeris::PlanetaryEphemeris::getProvider();
        auto sun_pos_bary = ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(t);
        Eigen::Vector3d sun_pos_bary_au = sun_pos_bary.to_eigen_si() / (constants::AU * 1000.0);
        
        static const math::RotationMatrix<core::GCRF, core::ECLIPJ2000> rot_ecl = 
            coordinates::ReferenceFrame::get_rotation<core::GCRF, core::ECLIPJ2000>();

        if (propagator_->settings().integrate_in_ecliptic) {
            sun_pos_bary_au = rot_ecl.to_eigen() * sun_pos_bary_au;
        }

        auto add_planetary_jacobian = [&](ephemeris::CelestialBody body, double gm_au) {
            auto p_pos_bary = provider->getPosition(body, t);
            Eigen::Vector3d p_pos_bary_au = p_pos_bary.to_eigen_si() / (constants::AU * 1000.0);
            
            if (propagator_->settings().integrate_in_ecliptic) {
                p_pos_bary_au = rot_ecl.to_eigen() * p_pos_bary_au;
            }
            
            Eigen::Vector3d r_planet_helio = p_pos_bary_au - sun_pos_bary_au;
            // Analytical gradient of GM * (r_p - r) / |r_p - r|^3 w.r.t r
            // This is exactly compute_acceleration_position_partial with relative vector (r - r_p)
            A.block<3, 3>(3, 0) += compute_acceleration_position_partial(r_ast - r_planet_helio, gm_au);
        };

        const auto& s = propagator_->settings();
        if (s.perturb_mercury) add_planetary_jacobian(ephemeris::CelestialBody::MERCURY, constants::GM_MERCURY_AU);
        if (s.perturb_venus)   add_planetary_jacobian(ephemeris::CelestialBody::VENUS,   constants::GM_VENUS_AU);
        if (s.perturb_earth)   add_planetary_jacobian(ephemeris::CelestialBody::EARTH,   constants::GM_EARTH_AU);
        if (s.perturb_mars)    add_planetary_jacobian(ephemeris::CelestialBody::MARS,    constants::GM_MARS_AU);
        if (s.perturb_jupiter) add_planetary_jacobian(ephemeris::CelestialBody::JUPITER, constants::GM_JUPITER_AU);
        if (s.perturb_saturn)  add_planetary_jacobian(ephemeris::CelestialBody::SATURN,  constants::GM_SATURN_AU);
        if (s.perturb_uranus)  add_planetary_jacobian(ephemeris::CelestialBody::URANUS,  constants::GM_URANUS_AU);
        if (s.perturb_neptune) add_planetary_jacobian(ephemeris::CelestialBody::NEPTUNE, constants::GM_NEPTUNE_AU);
    }
    
    // Note: Relativity and Non-gravitational partials (Yarkovsky) are tiny and neglected in Jacobian
    // for performance. They will be captured by the residuals but not the STM geometry.

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

    Eigen::Matrix<double, 2, 6> partial_radec_state_au;
    
    double au_to_m = constants::AU * 1000.0;

    // Position partials: [rad/m] * [m/AU] = [rad/AU]
    partial_radec_state_au.block<2, 3>(0, 0) = d_radec_d_pos_frame * au_to_m;
    
    // Velocity partials are identically zero at epsilon=0 for geometric observation
    partial_radec_state_au.block<2, 3>(0, 3).setZero();
    
    return partial_radec_state_au;
}

// Explicit instantiations
template class StateTransitionMatrix<core::GCRF>;
template class StateTransitionMatrix<core::ECLIPJ2000>;

} // namespace astdyn::orbit_determination
