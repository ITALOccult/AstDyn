/**
 * @file Propagator.cpp
 * @brief Implementation of orbital propagation
 */

#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/core/frame_tags.hpp"
#include "astdyn/propagation/kepler_propagator.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include <cmath>
#include <iostream>

namespace astdyn::propagation {
using namespace astdyn::constants;

using ephemeris::CelestialBody;

// ============================================================================
// Propagator Implementation
// ============================================================================

Propagator::Propagator(std::shared_ptr<Integrator> integrator,
                      std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris,
                      const PropagatorSettings& settings)
    : integrator_(std::move(integrator)), ephemeris_(std::move(ephemeris)), settings_(settings) {
    mat_ecl_ = coordinates::ReferenceFrame::j2000_to_ecliptic();
    if (settings_.include_asteroids) {
        asteroids_ = std::make_shared<ephemeris::AsteroidPerturbations>();
        if (settings_.use_default_asteroid_set) asteroids_->loadAstDynDefaultSet();
        else if (settings_.use_default_30_set) asteroids_->loadDefault30Asteroids();
        if (!settings_.asteroid_ephemeris_file.empty()) asteroids_->loadSPK(settings_.asteroid_ephemeris_file);
    }
    auto aas = std::dynamic_pointer_cast<AASIntegrator>(integrator_);
    if (aas) {
        bool sun = std::abs(settings_.central_body_gm - constants::GMS) < 1e-5;
        double j2 = sun ? (settings_.include_sun_j2 ? constants::SUN_J2 : 0.0) : (settings_.include_earth_j2 ? constants::EARTH_J2 : 0.0);
        double r_eq = sun ? constants::R_SUN_AU : (constants::R_EARTH/constants::AU);
        aas->set_central_body(settings_.central_body_gm, j2, r_eq);
    }
}

void Propagator::update_force_cache(time::EpochTDB t) {
    if (cache_valid_ && std::abs(t.mjd() - last_t_cache_.mjd()) < 1e-13) return;
    last_t_cache_ = t;
    auto sun_pos_bary = ephemeris_->getSunBarycentricPosition(t);
    last_sun_pos_bary_cache_ = sun_pos_bary.to_eigen_si() / (constants::AU * 1000.0);
    if (settings_.integrate_in_ecliptic) last_sun_pos_bary_cache_ = mat_ecl_ * last_sun_pos_bary_cache_;
    num_planets_cached_ = 0;
    auto provider = ephemeris_->getProvider();
    if (settings_.include_planets && provider) {
        static const CelestialBody b[] = { CelestialBody::MERCURY, CelestialBody::VENUS, CelestialBody::EARTH, CelestialBody::MARS, CelestialBody::JUPITER, CelestialBody::SATURN, CelestialBody::URANUS, CelestialBody::NEPTUNE, CelestialBody::MOON };
        int count = settings_.include_moon ? 9 : 8;
        for (int i = 0; i < count; ++i) {
            Eigen::Vector3d p_bary = provider->getPosition(b[i], t).to_eigen_si() / (constants::AU * 1000.0);
            if (settings_.integrate_in_ecliptic) p_bary = mat_ecl_ * p_bary;
            planet_cache_[num_planets_cached_++] = {b[i], p_bary};
        }
    }
    cache_valid_ = true;
}

Eigen::Vector3d Propagator::compute_n_body_acceleration(const Eigen::Vector3d& position) {
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    if (settings_.baricentric_integration) {
        Eigen::Vector3d ds = last_sun_pos_bary_cache_ - position;
        acc += (constants::GMS / std::pow(ds.norm(), 3)) * ds;
    } else {
        double r = position.norm();
        acc = -settings_.central_body_gm * position / (r * r * r);
    }
    for (int i = 0; i < num_planets_cached_; ++i) {
        const auto& p = planet_cache_[i];
        double gm = ephemeris::PlanetaryEphemeris::planet_gm(p.body);
        Eigen::Vector3d target = settings_.baricentric_integration ? p.pos_bary_au : (p.pos_bary_au - last_sun_pos_bary_cache_);
        Eigen::Vector3d delta = target - position;
        acc += gm * (delta / std::pow(delta.norm(), 3));
        if (!settings_.baricentric_integration) acc -= gm * target / std::pow(target.norm(), 3);
    }
    return acc;
}

Eigen::Vector3d Propagator::compute_harmonic_acceleration(const Eigen::Vector3d& position, time::EpochTDB t) {
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    if (settings_.include_earth_j2) {
        auto provider = ephemeris_->getProvider();
        if (provider) {
            Eigen::Vector3d e_bary = provider->getPosition(CelestialBody::EARTH, t).to_eigen_si() / (constants::AU * 1000.0);
            if (settings_.integrate_in_ecliptic) e_bary = mat_ecl_ * e_bary;
            Eigen::Vector3d r_rel = position - e_bary; double r = r_rel.norm();
            if (r > 1e-6) {
                Eigen::Vector3d r_eq = settings_.integrate_in_ecliptic ? (coordinates::ReferenceFrame::ecliptic_to_j2000() * r_rel) : r_rel;
                double j2_c = 1.5 * constants::EARTH_J2 * constants::GM_EARTH_AU * std::pow(constants::R_EARTH/constants::AU, 2);
                double z_r = r_eq.z()/r, r5 = std::pow(r, 5);
                Eigen::Vector3d a_j2; a_j2 << r_eq.x()*(5*z_r*z_r-1), r_eq.y()*(5*z_r*z_r-1), r_eq.z()*(5*z_r*z_r-3);
                acc += (settings_.integrate_in_ecliptic ? coordinates::ReferenceFrame::j2000_to_ecliptic() : Eigen::Matrix3d::Identity()) * (j2_c/r5 * a_j2);
            }
        }
    }
    if (settings_.include_sun_j2) {
        Eigen::Vector3d r_h = position - last_sun_pos_bary_cache_; double r = r_h.norm();
        if (r > 1e-6) {
            static const double a0 = 286.13 * constants::DEG_TO_RAD, d0 = 63.87 * constants::DEG_TO_RAD;
            Eigen::Vector3d pole_f = Eigen::Vector3d(std::cos(d0)*std::cos(a0), std::cos(d0)*std::sin(a0), std::sin(d0));
            if (settings_.integrate_in_ecliptic) pole_f = mat_ecl_ * pole_f;
            double j2_c = 1.5 * constants::SUN_J2 * constants::GMS * std::pow(constants::R_SUN_AU, 2);
            double z_n = r_h.dot(pole_f)/r;
            acc += (j2_c / std::pow(r, 5)) * ((5.0*z_n*z_n - 1.0)*pole_f - (2.0*z_n)*(r_h/r));
        }
    }
    return acc;
}

Eigen::Vector3d Propagator::compute_non_gravitational_acceleration(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, time::EpochTDB t) {
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    if (settings_.include_asteroids && asteroids_) acc += asteroids_->computePerturbationRaw(position, t.mjd(), last_sun_pos_bary_cache_, settings_.integrate_in_ecliptic);
    if (settings_.include_relativity) acc += relativistic_correction(position, velocity);
    if (settings_.include_yarkovsky && std::abs(settings_.yarkovsky_a2) > 1e-15) {
        double r = position.norm(), v = velocity.norm();
        if (r > 1e-4 && v > 1e-6) acc += (settings_.yarkovsky_a2 / (r*r)) * (velocity/v);
    }
    return acc;
}

Eigen::VectorXd Propagator::compute_derivatives(time::EpochTDB t, const Eigen::VectorXd& state) {
    Eigen::Vector3d position = state.head<3>(), velocity = state.tail<3>();
    update_force_cache(t);
    Eigen::Vector3d acc = compute_n_body_acceleration(position);
    acc += compute_harmonic_acceleration(position, t);
    acc += compute_non_gravitational_acceleration(position, velocity, t);
    Eigen::VectorXd xdot(6); xdot << velocity, acc;
    return xdot;
}


Eigen::Vector3d Propagator::relativistic_correction(const Eigen::Vector3d& position,
                                                   const Eigen::Vector3d& velocity) const {
    double r = position.norm();
    double v2 = velocity.squaredNorm();
    double c_au_d = physics::SpeedOfLight::to_au_d();
    double c2 = c_au_d * c_au_d;
    double mu = settings_.central_body_gm;
    
    double beta = settings_.ppn_beta;
    double gamma = settings_.ppn_gamma;
    
    double term1_coeff = 2.0 * (beta + gamma) * mu / r - gamma * v2;
    double rv_dot = position.dot(velocity);
    
    Eigen::Vector3d term1 = term1_coeff * position;
    Eigen::Vector3d term2 = 2.0 * (1.0 + gamma) * rv_dot * velocity;
    
    return (mu / (r * r * r * c2)) * (term1 + term2);
}

Eigen::VectorXd Propagator::integrate_raw_au(const Eigen::VectorXd& y0_au, double t0_mjd, double tf_mjd) {
    DerivativeFunction f = [this](double t_val, const Eigen::VectorXd& y) {
        return compute_derivatives(time::EpochTDB::from_mjd(t_val), y);
    };
    
    return integrator_->integrate(f, y0_au, t0_mjd, tf_mjd);
}

std::vector<Eigen::VectorXd> Propagator::integrate_raw_au_batch(const Eigen::VectorXd& y0_au, double t0_mjd, const std::vector<double>& tf_mjds) {
    DerivativeFunction f = [this](double t_val, const Eigen::VectorXd& y) {
        return compute_derivatives(time::EpochTDB::from_mjd(t_val), y);
    };
    
    return integrator_->integrate_at(f, y0_au, t0_mjd, tf_mjds);
}

} // namespace astdyn::propagation
