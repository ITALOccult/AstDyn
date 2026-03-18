#include "astdyn/propagation/RelativeMultiBodyPropagator.hpp"
#include "astdyn/core/Constants.hpp"

namespace astdyn::propagation {

RelativeMultiBodyPropagator::RelativeMultiBodyPropagator(std::shared_ptr<Integrator> integrator,
                                                         std::shared_ptr<ForceField> force_field)
    : integrator_(std::move(integrator)), force_field_(std::move(force_field)) {}

Eigen::VectorXd RelativeMultiBodyPropagator::RelativeDynamics::operator()(double t_sec, const Eigen::VectorXd& y) const {
    // y structure: [r0, v0, rho1, d_rho1, ... rhoN, d_rhoN]
    // r0, v0 (6): Primary Absolute Helio/SSB
    // rho_i, d_rho_i (6): Satellite relative to Primary
    
    Eigen::VectorXd dy(n_bodies * 6);
    time::EpochTDB t = time::EpochTDB::from_mjd(t0.mjd() + t_sec / 86400.0);
    
    // 1. Primary absolute acceleration
    Eigen::Vector3d r0 = y.segment<3>(0);
    Eigen::Vector3d v0 = y.segment<3>(3);
    Eigen::Vector3d a0 = force_field->total_acceleration(t, r0, v0);
    
    // Apply mutual attraction from satellites on primary (Reflex motion)
    for (size_t i = 1; i < n_bodies; ++i) {
        Eigen::Vector3d rho_i = y.segment<3>(i * 6);
        double dist = rho_i.norm();
        if (dist > 1e-12) {
            a0 += gms[i] * rho_i / (dist * dist * dist);
        }
    }
    
    dy.segment<3>(0) = v0 / 86400.0;
    dy.segment<3>(3) = a0 / (86400.0 * 86400.0);
    
    //  satellite relative accelerations
    for (size_t i = 1; i < n_bodies; ++i) {
        Eigen::Vector3d rho_i = y.segment<3>(i * 6);
        Eigen::Vector3d v_rho_i = y.segment<3>(i * 6 + 3);
        
        Eigen::Vector3d ri = r0 + rho_i; // Absolute satellite
        Eigen::Vector3d vi = v0 + v_rho_i;
        
        // Total absolute acceleration of satellite
        Eigen::Vector3d ai = force_field->total_acceleration(t, ri, vi);
        
        // Add mutual attraction from other satellites
        for (size_t j = 1; j < n_bodies; ++j) {
            if (i == j) continue;
            Eigen::Vector3d rho_j = y.segment<3>(j * 6);
            Eigen::Vector3d dist_vec = rho_i - rho_j; // ri - rj
            double dist = dist_vec.norm();
            if (dist > 1e-12) {
                ai += -gms[j] * dist_vec / (dist * dist * dist);
            }
        }
        
        // Equation: d^2(ri - r0)/dt^2 = ai - a0
        dy.segment<3>(i * 6) = v_rho_i / 86400.0;
        dy.segment<3>(i * 6 + 3) = (ai - a0) / (86400.0 * 86400.0);
    }
    
    return dy;
}

std::vector<MultiBodyState> RelativeMultiBodyPropagator::propagate(
    const std::vector<MultiBodyState>& initial_states,
    time::EpochTDB start_time,
    time::EpochTDB target_time) 
{
    if (initial_states.empty()) return {};
    size_t N = initial_states.size();
    
    Eigen::VectorXd y0(N * 6);
    RelativeDynamics dynamics;
    dynamics.force_field = force_field_;
    dynamics.n_bodies = N;
    dynamics.t0 = start_time;
    
    // Primary - Always absolute
    y0(0) = physics::Distance::from_si(initial_states[0].position.x_si()).to_au();
    y0(1) = physics::Distance::from_si(initial_states[0].position.y_si()).to_au();
    y0(2) = physics::Distance::from_si(initial_states[0].position.z_si()).to_au();
    y0(3) = physics::Velocity::from_si(initial_states[0].velocity.x_si()).to_au_d();
    y0(4) = physics::Velocity::from_si(initial_states[0].velocity.y_si()).to_au_d();
    y0(5) = physics::Velocity::from_si(initial_states[0].velocity.z_si()).to_au_d();
    dynamics.gms.push_back(initial_states[0].gm.to_au3_d2());
    
    for (size_t i = 1; i < N; ++i) {
        auto pos_rel = initial_states[i].position - initial_states[0].position;
        auto vel_rel = initial_states[i].velocity - initial_states[0].velocity;
        
        y0(i*6 + 0) = physics::Distance::from_si(pos_rel.x_si()).to_au();
        y0(i*6 + 1) = physics::Distance::from_si(pos_rel.y_si()).to_au();
        y0(i*6 + 2) = physics::Distance::from_si(pos_rel.z_si()).to_au();
        y0(i*6 + 3) = physics::Velocity::from_si(vel_rel.x_si()).to_au_d();
        y0(i*6 + 4) = physics::Velocity::from_si(vel_rel.y_si()).to_au_d();
        y0(i*6 + 5) = physics::Velocity::from_si(vel_rel.z_si()).to_au_d();
        dynamics.gms.push_back(initial_states[i].gm.to_au3_d2());
    }
    
    double dt_sec = (target_time.mjd() - start_time.mjd()) * 86400.0;
    Eigen::VectorXd yf = integrator_->integrate(dynamics, y0, 0.0, dt_sec);
    
    std::vector<MultiBodyState> results;
    results.reserve(N);
    
    // Primary Absolute
    MultiBodyState p;
    p.name = initial_states[0].name;
    p.gm = initial_states[0].gm;
    p.position = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(
        physics::Distance::from_au(yf(0)).to_m(),
        physics::Distance::from_au(yf(1)).to_m(),
        physics::Distance::from_au(yf(2)).to_m()
    );
    p.velocity = math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(
        physics::Velocity::from_au_d(yf(3)).to_ms(),
        physics::Velocity::from_au_d(yf(4)).to_ms(),
        physics::Velocity::from_au_d(yf(5)).to_ms()
    );
    results.push_back(p);
    
    for (size_t i = 1; i < N; ++i) {
        MultiBodyState s;
        s.name = initial_states[i].name;
        s.gm = initial_states[i].gm;
        // Absolute = Primary + Relative
        auto pos_rel = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(
            physics::Distance::from_au(yf(i*6+0)).to_m(),
            physics::Distance::from_au(yf(i*6+1)).to_m(),
            physics::Distance::from_au(yf(i*6+2)).to_m()
        );
        auto vel_rel = math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(
            physics::Velocity::from_au_d(yf(i*6+3)).to_ms(),
            physics::Velocity::from_au_d(yf(i*6+4)).to_ms(),
            physics::Velocity::from_au_d(yf(i*6+5)).to_ms()
        );
        s.position = p.position + pos_rel;
        s.velocity = p.velocity + vel_rel;
        results.push_back(s);
    }
    
    return results;
}

} // namespace astdyn::propagation
