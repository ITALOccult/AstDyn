/**
 * @file Propagator.cpp
 * @brief Implementation of orbital propagation
 */

#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "src/core/frame_tags.hpp"
#include "src/propagation/kepler_propagator.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
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
    : integrator_(std::move(integrator)),
      ephemeris_(std::move(ephemeris)),
      settings_(settings) {
    
    mat_ecl_ = coordinates::ReferenceFrame::j2000_to_ecliptic();

    if (settings_.include_asteroids) {
        asteroids_ = std::make_shared<ephemeris::AsteroidPerturbations>();
        if (!settings_.asteroid_ephemeris_file.empty()) {
            asteroids_->loadSPK(settings_.asteroid_ephemeris_file);
        }
    }
}

Eigen::VectorXd Propagator::compute_derivatives(time::EpochTDB t, const Eigen::VectorXd& state) {
    Eigen::Vector3d position = state.head<3>();
    Eigen::Vector3d velocity = state.tail<3>();
    
    // 1. Update/Check Cache
    if (!cache_valid_ || std::abs(t.mjd() - last_t_cache_.mjd()) > 1e-13) {
        last_t_cache_ = t;
        
        auto sun_pos_bary = ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(t);
        last_sun_pos_bary_cache_ = sun_pos_bary.to_eigen_si() / (constants::AU * 1000.0);
        if (settings_.integrate_in_ecliptic) {
            last_sun_pos_bary_cache_ = mat_ecl_ * last_sun_pos_bary_cache_;
        }
        
        // Refresh planets
        num_planets_cached_ = 0;
        auto provider = ephemeris::PlanetaryEphemeris::getProvider();
        if (settings_.include_planets && provider) {
            static const CelestialBody bodies[] = {
                CelestialBody::MERCURY, CelestialBody::VENUS, CelestialBody::EARTH,
                CelestialBody::MARS, CelestialBody::JUPITER, CelestialBody::SATURN,
                CelestialBody::URANUS, CelestialBody::NEPTUNE, CelestialBody::MOON
            };
            int count = settings_.include_moon ? 9 : 8;
            
            for (int i = 0; i < count; ++i) {
                auto p_pos_state = provider->getPosition(bodies[i], t);
                Eigen::Vector3d p_pos_bary_au = p_pos_state.to_eigen_si() / (constants::AU * 1000.0);
                if (settings_.integrate_in_ecliptic) p_pos_bary_au = mat_ecl_ * p_pos_bary_au;
                planet_cache_[num_planets_cached_++] = {bodies[i], p_pos_bary_au};
            }
        }
        cache_valid_ = true;
    }
    
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();

    if (settings_.baricentric_integration) {
        // SSB: Acc = sum [ GM_k * (r_k - r) / |r_k - r|^3 ] including Sun
        // Sun contribution
        Eigen::Vector3d delta_sun = last_sun_pos_bary_cache_ - position;
        double d_sun = delta_sun.norm();
        acc += (constants::GMS / (d_sun * d_sun * d_sun)) * delta_sun;
        
        // Planets
        for (int i = 0; i < num_planets_cached_; ++i) {
            const auto& p = planet_cache_[i];
            double gm = 0;
            switch(p.body) {
                case CelestialBody::MERCURY: gm = constants::GM_MERCURY_AU; break;
                case CelestialBody::VENUS:   gm = constants::GM_VENUS_AU; break;
                case CelestialBody::EARTH:   gm = constants::GM_EARTH_AU; break;
                case CelestialBody::MARS:    gm = constants::GM_MARS_AU; break;
                case CelestialBody::JUPITER: gm = constants::GM_JUPITER_AU; break;
                case CelestialBody::SATURN:  gm = constants::GM_SATURN_AU; break;
                case CelestialBody::URANUS:  gm = constants::GM_URANUS_AU; break;
                case CelestialBody::NEPTUNE: gm = constants::GM_NEPTUNE_AU; break;
                case CelestialBody::MOON:    gm = constants::GM_MOON_AU; break;
                default: break;
            }
            Eigen::Vector3d delta = p.pos_bary_au - position;
            double r2 = delta.squaredNorm();
            acc += (gm / (r2 * std::sqrt(r2))) * delta;
        }
    } else {
        // Heliocentric: Main term is -mu*r/r^3
        acc = two_body_acceleration(position);
        
        if (settings_.include_planets) {
            for (int i = 0; i < num_planets_cached_; ++i) {
                const auto& p = planet_cache_[i];
                double gm = 0;
                switch(p.body) {
                    case CelestialBody::MERCURY: gm = constants::GM_MERCURY_AU; break;
                    case CelestialBody::VENUS:   gm = constants::GM_VENUS_AU; break;
                    case CelestialBody::EARTH:   gm = constants::GM_EARTH_AU; break;
                    case CelestialBody::MARS:    gm = constants::GM_MARS_AU; break;
                    case CelestialBody::JUPITER: gm = constants::GM_JUPITER_AU; break;
                    case CelestialBody::SATURN:  gm = constants::GM_SATURN_AU; break;
                    case CelestialBody::URANUS:  gm = constants::GM_URANUS_AU; break;
                    case CelestialBody::NEPTUNE: gm = constants::GM_NEPTUNE_AU; break;
                    case CelestialBody::MOON:    gm = constants::GM_MOON_AU; break;
                    default: break;
                }
                Eigen::Vector3d planet_pos_helio_au = p.pos_bary_au - last_sun_pos_bary_cache_;
                Eigen::Vector3d delta = planet_pos_helio_au - position;
                
                double r2_delta = delta.squaredNorm();
                double r2_p = planet_pos_helio_au.squaredNorm();
                acc += gm * (delta / (r2_delta * std::sqrt(r2_delta)) - planet_pos_helio_au / (r2_p * std::sqrt(r2_p)));
            }
        }
    }
    
    if (settings_.include_asteroids && asteroids_) {
        acc += asteroids_->computePerturbationRaw(position, t.mjd(), last_sun_pos_bary_cache_, settings_.integrate_in_ecliptic);
    }
    
    if (settings_.include_relativity) {
        acc += relativistic_correction(position, velocity);
    }
    
    // Yarkovsky Effect
    if (settings_.include_yarkovsky && std::abs(settings_.yarkovsky_a2) > 0.0) {
        double r_au = position.norm();
        double v_au_d = velocity.norm();
        if (r_au > 1e-4 && v_au_d > 1e-6) {
            double acc_val = settings_.yarkovsky_a2 / (r_au * r_au);
            acc += acc_val * (velocity / v_au_d);
        }
    }
    
    Eigen::VectorXd xdot(6);
    xdot.head<3>() = velocity;
    xdot.tail<3>() = acc;
    return xdot;
}

Eigen::Vector3d Propagator::two_body_acceleration(const Eigen::Vector3d& position) const {
    double r = position.norm();
    double r3 = r * r * r;
    // central_body_gm is configured in AU^3/d^2
    double mu = settings_.central_body_gm;
    return -mu * position / r3;
}

Eigen::Vector3d Propagator::planetary_perturbations(const Eigen::Vector3d& position,
                                                  time::EpochTDB t,
                                                  const Eigen::Vector3d& sun_pos_bary) {
    Eigen::Vector3d perturbation = Eigen::Vector3d::Zero();
    auto provider = ephemeris::PlanetaryEphemeris::getProvider();
    const double au_m = constants::AU * 1000.0;

    auto add_planet_perturbation = [&]( CelestialBody planet, double planet_gm_au) {
        // Direct AU barycentric position
        auto p_pos_state = provider->getPosition(planet, t);
        Eigen::Vector3d p_pos_bary_au = p_pos_state.to_eigen_si() / (constants::AU * 1000.0);

        if (settings_.integrate_in_ecliptic) {
            p_pos_bary_au = mat_ecl_ * p_pos_bary_au;
        }
        
        Eigen::Vector3d planet_pos_helio_au = p_pos_bary_au - sun_pos_bary;
        Eigen::Vector3d delta = planet_pos_helio_au - position;
        
        const double d_norm = delta.norm();
        const double d3 = d_norm * d_norm * d_norm;
        const double p_dist = planet_pos_helio_au.norm();
        const double p_dist3 = p_dist * p_dist * p_dist;
        
        perturbation += planet_gm_au * (delta / d3 - planet_pos_helio_au / p_dist3);
    };
    
    if (settings_.perturb_mercury) add_planet_perturbation(CelestialBody::MERCURY, constants::GM_MERCURY_AU);
    if (settings_.perturb_venus)   add_planet_perturbation(CelestialBody::VENUS,   constants::GM_VENUS_AU);
    if (settings_.perturb_earth)   add_planet_perturbation(CelestialBody::EARTH,   constants::GM_EARTH_AU);
    if (settings_.perturb_mars)    add_planet_perturbation(CelestialBody::MARS,    constants::GM_MARS_AU);
    if (settings_.perturb_jupiter) add_planet_perturbation(CelestialBody::JUPITER, constants::GM_JUPITER_AU);
    if (settings_.perturb_saturn)  add_planet_perturbation(CelestialBody::SATURN,  constants::GM_SATURN_AU);
    if (settings_.perturb_uranus)  add_planet_perturbation(CelestialBody::URANUS,  constants::GM_URANUS_AU);
    if (settings_.perturb_neptune) add_planet_perturbation(CelestialBody::NEPTUNE, constants::GM_NEPTUNE_AU);
    
    if (settings_.include_moon) {
        double moon_gm_au = physics::GravitationalParameter::from_km3_s2(constants::GM_MOON).to_au3_d2();
        add_planet_perturbation(CelestialBody::MOON, moon_gm_au);
    }
    
    return perturbation;
}

Eigen::Vector3d Propagator::asteroid_perturbations(const Eigen::Vector3d& position,
                                                 time::EpochTDB t,
                                                 const Eigen::Vector3d& sun_pos_bary) {
    // Asteroids library assumes AU and returns AU/day^2.
    // Use the raw interface to avoid state creation overhead during integration.
    return asteroids_->computePerturbationRaw(position, t.mjd(), sun_pos_bary, settings_.integrate_in_ecliptic);
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
