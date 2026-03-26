#include "astdyn/propagation/MultiBodyPropagator.hpp"
#include "astdyn/math/kahan_sum.hpp"
#include "astdyn/core/Constants.hpp"

namespace astdyn::propagation {

// PPN 1PN correction with beta=gamma=1 (standard GR), inputs and output in AU/day units
static Eigen::Vector3d relativistic_correction_1pn(
    const Eigen::Vector3d& r_au,
    const Eigen::Vector3d& v_au_d,
    double gm_au3_d2)
{
    const double r_mag = r_au.norm();
    const double c2    = constants::SPEED_OF_LIGHT_AU_PER_DAY * constants::SPEED_OF_LIGHT_AU_PER_DAY;
    // PPN terms with beta=gamma=1: term1 = 4*mu/r - v^2, term2 = 4*(r·v)
    const double term1     = 4.0 * gm_au3_d2 / r_mag - v_au_d.squaredNorm();
    const double term2     = 4.0 * r_au.dot(v_au_d);
    return (gm_au3_d2 / (c2 * r_mag * r_mag * r_mag)) * (term1 * r_au + term2 * v_au_d);
}

MultiBodyPropagator::MultiBodyPropagator(std::shared_ptr<Integrator> integrator,
                                         std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris)
    : integrator_(std::move(integrator)), ephemeris_(std::move(ephemeris)) {}

Eigen::VectorXd MultiBodyPropagator::Dynamics::operator()(double t_sec, const Eigen::VectorXd& y) const {
    Eigen::VectorXd dy(n_bodies * 6);
    const double t_mjd = t0_ephemeris.mjd() + t_sec / constants::SECONDS_PER_DAY;
    const time::EpochTDB epoch_t = time::EpochTDB::from_mjd(t_mjd);

    // Units of y: AU and AU/day
    // Pre-calculate planetary positions relative to Sun in AU (ECLIPJ2000)
    // Bug 3 fix: use same body set as ForceField::n_body_perturbation()
    // {Mercury..Neptune, Moon} — no Pluto, Moon conditional on include_moon
    using ephemeris::CelestialBody;
    static constexpr CelestialBody kPerturbingBodies[] = {
        CelestialBody::MERCURY, CelestialBody::VENUS,   CelestialBody::EARTH,
        CelestialBody::MARS,    CelestialBody::JUPITER, CelestialBody::SATURN,
        CelestialBody::URANUS,  CelestialBody::NEPTUNE, CelestialBody::MOON
    };
    const int n_perturbers = include_moon ? 9 : 8;

    std::vector<Eigen::Vector3d> planet_pos;
    if (ephemeris) {
        const auto sun_ssb_au = ephemeris->getSunBarycentricPosition(epoch_t).to_eigen_si()
                                / (constants::AU * 1000.0);

        for (int i = 0; i < n_perturbers; ++i) {
            const auto p_ssb_au = ephemeris->getPosition(kPerturbingBodies[i], epoch_t).to_eigen_si()
                                  / (constants::AU * 1000.0);
            planet_pos.push_back(p_ssb_au - sun_ssb_au);
        }
    }

    for (size_t i = 0; i < n_bodies; ++i) {
        const Eigen::Vector3d r1 = y.segment<3>(i * 6);
        const Eigen::Vector3d v1 = y.segment<3>(i * 6 + 3);

        // 1. Solar acceleration
        const double r1_norm = r1.norm();
        Eigen::Vector3d acc = -gm_sun * r1 / (r1_norm * r1_norm * r1_norm);

        // 2 & 3. Mutual + planetary accelerations with Kahan compensation
        Eigen::Vector3d kahan_comp = Eigen::Vector3d::Zero();

        for (size_t j = 0; j < n_bodies; ++j) {
            if (i == j) continue;
            const Eigen::Vector3d r2     = y.segment<3>(j * 6);
            const Eigen::Vector3d dist   = r1 - r2;
            const double          d_norm = dist.norm();
            if (d_norm > 1e-12) {
                const Eigen::Vector3d term = -gms[j] * dist / (d_norm * d_norm * d_norm);
                math::kahan_add(acc, kahan_comp, term);
            }
        }

        for (int p = 0; ephemeris && p < n_perturbers; ++p) {
            const double          gm_p     = ephemeris::PlanetaryEphemeris::planet_gm(kPerturbingBodies[p]);
            const Eigen::Vector3d rp       = planet_pos[p];
            const Eigen::Vector3d r1p      = r1 - rp;
            const double          r1p_norm = r1p.norm();
            const double          rp_norm  = rp.norm();
            const Eigen::Vector3d term     = gm_p * (-r1p / (r1p_norm * r1p_norm * r1p_norm)
                                                      - rp  / (rp_norm  * rp_norm  * rp_norm));
            math::kahan_add(acc, kahan_comp, term);
        }

        // 4. Relativistic correction (PPN 1PN, beta=gamma=1)
        if (include_relativity) {
            acc += relativistic_correction_1pn(r1, v1, gm_sun);
        }

        dy.segment<3>(i * 6)     = v1 / constants::SECONDS_PER_DAY;
        dy.segment<3>(i * 6 + 3) = acc / (constants::SECONDS_PER_DAY * constants::SECONDS_PER_DAY);
    }
    return dy;
}

std::vector<MultiBodyState> MultiBodyPropagator::propagate(
    const std::vector<MultiBodyState>& initial_states,
    time::EpochTDB start_time,
    time::EpochTDB target_time)
{
    if (initial_states.empty()) return {};

    const size_t n_bodies = initial_states.size();
    Eigen::VectorXd y0(n_bodies * 6);
    Dynamics dynamics;
    dynamics.n_bodies       = n_bodies;
    dynamics.gm_sun         = constants::GMS;
    dynamics.t0_ephemeris   = start_time;
    dynamics.ephemeris      = ephemeris_;
    dynamics.gms.reserve(n_bodies);

    for (size_t i = 0; i < n_bodies; ++i) {
        const auto& s = initial_states[i];
        y0(i*6 + 0) = physics::Distance::from_si(s.position.x_si()).to_au();
        y0(i*6 + 1) = physics::Distance::from_si(s.position.y_si()).to_au();
        y0(i*6 + 2) = physics::Distance::from_si(s.position.z_si()).to_au();
        y0(i*6 + 3) = physics::Velocity::from_si(s.velocity.x_si()).to_au_d();
        y0(i*6 + 4) = physics::Velocity::from_si(s.velocity.y_si()).to_au_d();
        y0(i*6 + 5) = physics::Velocity::from_si(s.velocity.z_si()).to_au_d();
        dynamics.gms.push_back(s.gm.to_au3_d2());
    }

    const double t0_sec = 0.0;
    const double tf_sec = (target_time.mjd() - start_time.mjd()) * constants::SECONDS_PER_DAY;
    const Eigen::VectorXd yf = integrator_->integrate(dynamics, y0, t0_sec, tf_sec);

    std::vector<MultiBodyState> results;
    results.reserve(n_bodies);
    for (size_t i = 0; i < n_bodies; ++i) {
        MultiBodyState res;
        res.name = initial_states[i].name;
        res.gm   = initial_states[i].gm;
        res.position = math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(
            physics::Distance::from_au(yf(i*6+0)).to_m(),
            physics::Distance::from_au(yf(i*6+1)).to_m(),
            physics::Distance::from_au(yf(i*6+2)).to_m());
        res.velocity = math::Vector3<core::ECLIPJ2000, physics::Velocity>::from_si(
            physics::Velocity::from_au_d(yf(i*6+3)).to_ms(),
            physics::Velocity::from_au_d(yf(i*6+4)).to_ms(),
            physics::Velocity::from_au_d(yf(i*6+5)).to_ms());
        results.push_back(res);
    }
    return results;
}

} // namespace astdyn::propagation
