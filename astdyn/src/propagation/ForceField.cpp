#include "astdyn/propagation/ForceField.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>

namespace astdyn::propagation {

ForceField::ForceField(const PropagatorSettings& settings, 
                      std::shared_ptr<ephemeris::PlanetaryEphemeris> ephem)
    : settings_(settings), ephemeris_(std::move(ephem)) 
{
    if (settings_.include_asteroids) {
        asteroids_ = std::make_shared<ephemeris::AsteroidPerturbations>();
        if (settings_.use_default_asteroid_set) asteroids_->loadAstDynDefaultSet();
        else if (settings_.use_default_30_set) asteroids_->loadDefault30Asteroids();
        if (!settings_.asteroid_ephemeris_file.empty()) {
            asteroids_->loadSPK(settings_.asteroid_ephemeris_file);
        }
    }
}

Eigen::Vector3d ForceField::total_acceleration(
    time::EpochTDB t, 
    const Eigen::Vector3d& pos_au, 
    const Eigen::Vector3d& vel_au_d) const 
{
    // 1. Central Body Gravity
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    
    if (settings_.baricentric_integration && ephemeris_) {
        // Absolute SSB integration: attraction from Sun
        auto sun_ssb = ephemeris_->getSunBarycentricPosition(t).to_eigen_si() / (constants::AU * 1000.0);
        Eigen::Vector3d ds = sun_ssb - pos_au;
        acc = (constants::GMS / std::pow(ds.norm(), 3)) * ds;
    } else {
        // Heliocentric (default)
        double r_mag = pos_au.norm();
        acc = -settings_.central_body_gm * pos_au / (r_mag * r_mag * r_mag);
    }

    // 2. N-Body Perturbations (Planets + Moon)
    if (settings_.include_planets && ephemeris_) {
        acc += n_body_perturbation(t, pos_au);
    }

    // 3. Asteroid Perturbations
    if (settings_.include_asteroids && asteroids_ && ephemeris_) {
        auto sun_ssb = ephemeris_->getSunBarycentricPosition(t).to_eigen_si() / (constants::AU * 1000.0);
        acc += asteroids_->computePerturbationRaw(pos_au, t.mjd(), sun_ssb, settings_.integrate_in_ecliptic);
    }

    // 4. Relativistic Correction
    if (settings_.include_relativity) {
        acc += relativistic_correction(pos_au, vel_au_d);
    }

    // 5. Harmonic Corrections (J2)
    acc += j2_correction(t, pos_au);

    return acc;
}

Eigen::Vector3d ForceField::n_body_perturbation(time::EpochTDB t, const Eigen::Vector3d& pos_au) const {
    Eigen::Vector3d acc_p = Eigen::Vector3d::Zero();
    auto provider = ephemeris_->getProvider();
    if (!provider) return acc_p;

    auto sun_ssb = ephemeris_->getSunBarycentricPosition(t).to_eigen_si() / (constants::AU * 1000.0);
    auto rot_icrf_to_ecl = coordinates::ReferenceFrame::get_rotation<core::GCRF, core::ECLIPJ2000>();
    
    using ephemeris::CelestialBody;
    static const CelestialBody bodies[] = {
        CelestialBody::MERCURY, CelestialBody::VENUS, CelestialBody::EARTH, CelestialBody::MARS,
        CelestialBody::JUPITER, CelestialBody::SATURN, CelestialBody::URANUS, CelestialBody::NEPTUNE,
        CelestialBody::MOON
    };
    
    int count = settings_.include_moon ? 9 : 8;
    
    for (int i = 0; i < count; ++i) {
        double gm = ephemeris::PlanetaryEphemeris::planet_gm(bodies[i]);
        auto p_ssb = provider->getPosition(bodies[i], t).to_eigen_si() / (constants::AU * 1000.0);
        
        if (settings_.baricentric_integration) {
            // Absolute acceleration in SSB frame
            Eigen::Vector3d delta = p_ssb - pos_au;
            acc_p += gm * delta / std::pow(delta.norm(), 3);
        } else {
            // Relative acceleration in Heliocentric frame
            Eigen::Vector3d rp = p_ssb - sun_ssb;
            Eigen::Vector3d r1p = pos_au - rp;
            acc_p += gm * ( -r1p / std::pow(r1p.norm(), 3) - rp / std::pow(rp.norm(), 3) );
        }
    }
    return acc_p;
}

Eigen::Vector3d ForceField::relativistic_correction(const Eigen::Vector3d& r, const Eigen::Vector3d& v) const {
    double mu = settings_.central_body_gm;
    double r_mag = r.norm();
    double v2 = v.squaredNorm();
    double r_dot_v = r.dot(v);
    double c2 = constants::SPEED_OF_LIGHT_AU_PER_DAY * constants::SPEED_OF_LIGHT_AU_PER_DAY;
    
    double term1 = (2.0 * (settings_.ppn_beta + settings_.ppn_gamma) * mu / r_mag - settings_.ppn_gamma * v2);
    double term2 = 2.0 * (1.0 + settings_.ppn_gamma) * r_dot_v;
    
    return (mu / (c2 * r_mag * r_mag * r_mag)) * (term1 * r + term2 * v);
}

Eigen::Vector3d ForceField::j2_correction(time::EpochTDB t, const Eigen::Vector3d& pos_au) const {
    Eigen::Vector3d acc_j2 = Eigen::Vector3d::Zero();
    
    // Simplifed Sun J2 in Ecliptic
    if (settings_.include_sun_j2) {
        double r_mag = pos_au.norm();
        double factor = -1.5 * settings_.central_body_gm * constants::SUN_J2 * std::pow(constants::R_SUN_AU / r_mag, 2) / (r_mag * r_mag * r_mag);
        acc_j2.x() = factor * pos_au.x() * (1.0 - 5.0 * std::pow(pos_au.z() / r_mag, 2));
        acc_j2.y() = factor * pos_au.y() * (1.0 - 5.0 * std::pow(pos_au.z() / r_mag, 2));
        acc_j2.z() = factor * pos_au.z() * (3.0 - 5.0 * std::pow(pos_au.z() / r_mag, 2));
    }
    
    // Earth J2 (Simplified - needs proper frame for ultimate precision)
    if (settings_.include_earth_j2 && ephemeris_) {
        Eigen::Vector3d e_ssb = ephemeris_->getPosition(ephemeris::CelestialBody::EARTH, t).to_eigen_si() / (constants::AU * 1000.0);
        Eigen::Vector3d sun_ssb = ephemeris_->getSunBarycentricPosition(t).to_eigen_si() / (constants::AU * 1000.0);
        Eigen::Vector3d r_earth = settings_.baricentric_integration ? e_ssb : Eigen::Vector3d(e_ssb - sun_ssb);
        Eigen::Vector3d rel = pos_au - r_earth;
        double r = rel.norm();
        if (r > 1e-6) {
            double factor = -1.5 * constants::GM_EARTH_AU * constants::EARTH_J2 * std::pow(constants::R_EARTH_EQUATORIAL / constants::AU / r, 2) / (r * r * r);
            // Assuming Earth spin axis is J2000 Z (GCRF)
            acc_j2 += factor * rel; // Placeholder for full rotation-based J2
        }
    }
    
    return acc_j2;
}

} // namespace astdyn::propagation
