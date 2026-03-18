#include "astdyn/propagation/MultiBodyPropagator.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>

namespace astdyn::propagation {

MultiBodyPropagator::MultiBodyPropagator(std::shared_ptr<Integrator> integrator,
                                         std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris)
    : integrator_(std::move(integrator)), ephemeris_(std::move(ephemeris)) {}

Eigen::VectorXd MultiBodyPropagator::Dynamics::operator()(double t_sec, const Eigen::VectorXd& y) const {
    Eigen::VectorXd dy(n_bodies * 6);
    double t_mjd = t0_ephemeris.mjd() + t_sec / 86400.0;
    time::EpochTDB epoch_t = time::EpochTDB::from_mjd(t_mjd);
    
    // Units of y: AU and AU/day
    // We need derivatives in AU/sec and AU/sec^2
    
    // Pre-calculate planetary positions relative to Sun in AU (ECLIPTIC J2000)
    // To match the Heliocentric integration
    std::vector<Eigen::Vector3d> planet_pos;
    if (ephemeris) {
        auto rot_icrf_to_ecl = coordinates::ReferenceFrame::get_rotation<core::GCRF, core::ECLIPJ2000>();
        auto sun_ssb = ephemeris->getSunBarycentricPosition(epoch_t).to_eigen_si() / (constants::AU * 1000.0);
        
        for (int i = 1; i <= 9; ++i) { // Mercury to Pluto
            auto p_ssb = ephemeris->getPosition(static_cast<ephemeris::CelestialBody>(i), epoch_t).to_eigen_si() / (constants::AU * 1000.0);
            Eigen::Vector3d diff = p_ssb - sun_ssb;
            auto v_icrf = math::Vector3<core::GCRF, physics::Distance>::from_si(diff.x(), diff.y(), diff.z());
            planet_pos.push_back((rot_icrf_to_ecl * v_icrf).to_eigen_si());
        }
    }

    for (size_t i = 0; i < n_bodies; ++i) {
        Eigen::Vector3d r1 = y.segment<3>(i * 6);
        Eigen::Vector3d v1 = y.segment<3>(i * 6 + 3);
        
        // 1. Solar acceleration
        double r1_norm = r1.norm();
        Eigen::Vector3d acc = -gm_sun * r1 / (r1_norm * r1_norm * r1_norm);
        
        // 2. Mutual accelerations (among N bodies)
        for (size_t j = 0; j < n_bodies; ++j) {
            if (i == j) continue;
            Eigen::Vector3d r2 = y.segment<3>(j * 6);
            Eigen::Vector3d dist = r1 - r2;
            double d_norm = dist.norm();
            if (d_norm > 1e-12) {
                acc += -gms[j] * dist / (d_norm * d_norm * d_norm);
            }
        }
        
        // 3. Planetary Perturbations (in AU/day^2)
        if (ephemeris) {
            for (int p = 1; p <= 9; ++p) {
                double gm_p = ephemeris::PlanetaryEphemeris::planet_gm(static_cast<ephemeris::CelestialBody>(p));
                Eigen::Vector3d rp = planet_pos[p-1];
                Eigen::Vector3d r1p = r1 - rp;
                double r1p_norm = r1p.norm();
                double rp_norm = rp.norm();
                // N-body acceleration: gm * ( (rp-r1)/|r1-p|^3 - rp/|rp|^3 )
                acc += gm_p * ( -r1p / (r1p_norm * r1p_norm * r1p_norm) - rp / (rp_norm * rp_norm * rp_norm) );
            }
        }
        
        dy.segment<3>(i * 6) = v1 / 86400.0;
        dy.segment<3>(i * 6 + 3) = acc / (86400.0 * 86400.0);
    }
    return dy;
}

std::vector<MultiBodyState> MultiBodyPropagator::propagate(
    const std::vector<MultiBodyState>& initial_states,
    time::EpochTDB start_time,
    time::EpochTDB target_time) 
{
    if (initial_states.empty()) return {};
    
    size_t N = initial_states.size();
    Eigen::VectorXd y0(N * 6);
    Dynamics dynamics;
    dynamics.n_bodies = N;
    dynamics.gm_sun = constants::GMS;
    dynamics.t0_ephemeris = start_time;
    dynamics.ephemeris = ephemeris_;
    dynamics.gms.reserve(N);
    
    for (size_t i = 0; i < N; ++i) {
        const auto& s = initial_states[i];
        y0(i*6 + 0) = physics::Distance::from_si(s.position.x_si()).to_au();
        y0(i*6 + 1) = physics::Distance::from_si(s.position.y_si()).to_au();
        y0(i*6 + 2) = physics::Distance::from_si(s.position.z_si()).to_au();
        y0(i*6 + 3) = physics::Velocity::from_si(s.velocity.x_si()).to_au_d();
        y0(i*6 + 4) = physics::Velocity::from_si(s.velocity.y_si()).to_au_d();
        y0(i*6 + 5) = physics::Velocity::from_si(s.velocity.z_si()).to_au_d();
        dynamics.gms.push_back(s.gm.to_au3_d2());
    }
    
    double t0_mjd = start_time.mjd();
    double tf_mjd = target_time.mjd();
    
    Eigen::VectorXd yf = integrator_->integrate(dynamics, y0, 0.0, (tf_mjd - t0_mjd) * 86400.0);
    
    std::vector<MultiBodyState> results;
    results.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        MultiBodyState res;
        res.name = initial_states[i].name;
        res.gm = initial_states[i].gm;
        res.position = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(
            physics::Distance::from_au(yf(i*6+0)).to_m(),
            physics::Distance::from_au(yf(i*6+1)).to_m(),
            physics::Distance::from_au(yf(i*6+2)).to_m()
        );
        res.velocity = math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(
            physics::Velocity::from_au_d(yf(i*6+3)).to_ms(),
            physics::Velocity::from_au_d(yf(i*6+4)).to_ms(),
            physics::Velocity::from_au_d(yf(i*6+5)).to_ms()
        );
        results.push_back(res);
    }
    return results;
}

} // namespace astdyn::propagation
