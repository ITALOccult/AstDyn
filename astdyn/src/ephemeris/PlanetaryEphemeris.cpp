/**
 * @file PlanetaryEphemeris.cpp
 * @brief Implementation of planetary ephemeris calculations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Low-precision planetary positions using Simon et al. (1994) formulae.
 * Accuracy: 1-20 arcsec over 1800-2050.
 */

#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/math/MathUtils.hpp"
#include <cmath>
#include "astdyn/coordinates/ReferenceFrame.hpp"

namespace astdyn {
namespace ephemeris {

using Eigen::Vector3d;
using coordinates::CartesianState;

// Convert GM from km³/s² to AU³/day²
static constexpr double GM_SUN_AU_DAY2 = 1.32712440018e11 / (1.49597870700e8 * 1.49597870700e8 * 1.49597870700e8) * (86400.0 * 86400.0);

// ============================================================================
// Public Interface
// ============================================================================

std::shared_ptr<EphemerisProvider> PlanetaryEphemeris::global_provider_ = nullptr;

void PlanetaryEphemeris::setProvider(std::shared_ptr<EphemerisProvider> provider) {
    global_provider_ = provider;
}

Vector3d PlanetaryEphemeris::getPosition(CelestialBody body, double jd_tdb) {
    if (global_provider_ && global_provider_->isAvailable()) {
        return global_provider_->getPosition(body, jd_tdb);
    }

    if (body == CelestialBody::SUN) {
        return Vector3d::Zero();
    }
    
    double T = julianCenturies(jd_tdb);
    double elements[6];
    computeOrbitalElements(body, T, elements);
    
    Vector3d pos = elementsToPosition(elements);
    Vector3d pert = computePerturbations(body, T);
    
    return pos + pert;
}

Vector3d PlanetaryEphemeris::getVelocity(CelestialBody body, double jd_tdb) {
    if (global_provider_ && global_provider_->isAvailable()) {
        return global_provider_->getVelocity(body, jd_tdb);
    }

    if (body == CelestialBody::SUN) {
        return Vector3d::Zero();
    }
    
    double T = julianCenturies(jd_tdb);
    double elements[6];
    computeOrbitalElements(body, T, elements);
    
    return elementsToVelocity(elements, GM_SUN_AU_DAY2);
}

CartesianState PlanetaryEphemeris::getState(CelestialBody body, double jd_tdb) {
    // If provider is available, use it (optimized if getState() exists in interface, but getPos/getVel work too)
    // EphemerisProvider interface typically has getPosition/getVelocity separately.
    if (global_provider_ && global_provider_->isAvailable()) {
        Vector3d pos = global_provider_->getPosition(body, jd_tdb);
        Vector3d vel = global_provider_->getVelocity(body, jd_tdb);
        return CartesianState(pos, vel);
    }

    Vector3d pos = getPosition(body, jd_tdb);
    Vector3d vel = getVelocity(body, jd_tdb);
    return CartesianState(pos, vel);
}

Vector3d PlanetaryEphemeris::getSunBarycentricPosition(double jd_tdb) {
    // Simplified: Sun offset from barycenter due to planets
    // Dominated by Jupiter (~5e-4 AU) and Saturn (~3e-4 AU)
    
    Vector3d r_sun = Vector3d::Zero();
    
    // Major contributors
    double total_mass = PlanetaryData::MASS_SUN;
    
    // Jupiter
    Vector3d r_jup = getPosition(CelestialBody::JUPITER, jd_tdb);
    r_sun -= r_jup * (PlanetaryData::MASS_JUPITER / total_mass);
    
    // Saturn
    Vector3d r_sat = getPosition(CelestialBody::SATURN, jd_tdb);
    r_sun -= r_sat * (PlanetaryData::MASS_SATURN / total_mass);
    
    // Uranus
    Vector3d r_ura = getPosition(CelestialBody::URANUS, jd_tdb);
    r_sun -= r_ura * (PlanetaryData::MASS_URANUS / total_mass);
    
    // Neptune
    Vector3d r_nep = getPosition(CelestialBody::NEPTUNE, jd_tdb);
    r_sun -= r_nep * (PlanetaryData::MASS_NEPTUNE / total_mass);
    
    return r_sun;
}

CartesianState PlanetaryEphemeris::heliocentricToBarycentric(
    const CartesianState& heliocentric_state, 
    double jd_tdb) 
{
    Vector3d r_sun = getSunBarycentricPosition(jd_tdb);
    
    // Barycentric position = heliocentric position - Sun position
    Vector3d r_bary = heliocentric_state.position() - r_sun;
    
    // Approximate velocity (ignore Sun's velocity for now)
    Vector3d v_bary = heliocentric_state.velocity();
    
    return CartesianState(r_bary, v_bary);
}

// ============================================================================
// Private Implementation
// ============================================================================

void PlanetaryEphemeris::computeOrbitalElements(
    CelestialBody body, 
    double T,
    double elements[6]) 
{
    // Simon et al. (1994) formulae for mean orbital elements
    // Elements: [a, e, i, L, omega_bar, Omega]
    // a [AU], e [-], i [rad], L [rad], omega_bar [rad], Omega [rad]
    
    double a, e, i, L, omega_bar, Omega;
    
    switch (body) {
        case CelestialBody::MERCURY:
            a = 0.38709927;
            e = 0.20563593 + 0.00001906 * T;
            i = (7.00497902 - 0.00594749 * T) * constants::DEG_TO_RAD;
            L = (252.25032350 + 149472.67411175 * T) * constants::DEG_TO_RAD;
            omega_bar = (77.45779628 + 0.16047689 * T) * constants::DEG_TO_RAD;
            Omega = (48.33076593 - 0.12534081 * T) * constants::DEG_TO_RAD;
            break;
            
        case CelestialBody::VENUS:
            a = 0.72333566;
            e = 0.00677672 - 0.00004107 * T;
            i = (3.39467605 - 0.00078890 * T) * constants::DEG_TO_RAD;
            L = (181.97909950 + 58517.81538729 * T) * constants::DEG_TO_RAD;
            omega_bar = (131.60246718 + 0.00268329 * T) * constants::DEG_TO_RAD;
            Omega = (76.67984255 - 0.27769418 * T) * constants::DEG_TO_RAD;
            break;
            
        case CelestialBody::EARTH:
            a = 1.00000261 + 0.00000562 * T;
            e = 0.01671123 - 0.00004392 * T;
            i = (-0.00001531 - 0.01294668 * T) * constants::DEG_TO_RAD;
            L = (100.46457166 + 35999.37244981 * T) * constants::DEG_TO_RAD;
            omega_bar = (102.93768193 + 0.32327364 * T) * constants::DEG_TO_RAD;
            Omega = 0.0 * constants::DEG_TO_RAD;
            break;
            
        case CelestialBody::MARS:
            a = 1.52371034;
            e = 0.09339410 + 0.00007882 * T;
            i = (1.84969142 - 0.00813131 * T) * constants::DEG_TO_RAD;
            L = (-4.55343205 + 19140.30268499 * T) * constants::DEG_TO_RAD;
            omega_bar = (-23.94362959 + 0.44441088 * T) * constants::DEG_TO_RAD;
            Omega = (49.55953891 - 0.29257343 * T) * constants::DEG_TO_RAD;
            break;
            
        case CelestialBody::JUPITER:
            a = 5.20288700;
            e = 0.04838624 - 0.00013537 * T;
            i = (1.30439695 - 0.00183714 * T) * constants::DEG_TO_RAD;
            L = (34.39644051 + 3034.74612775 * T) * constants::DEG_TO_RAD;
            omega_bar = (14.72847983 + 0.21252668 * T) * constants::DEG_TO_RAD;
            Omega = (100.47390909 + 0.20469106 * T) * constants::DEG_TO_RAD;
            break;
            
        case CelestialBody::SATURN:
            a = 9.53667594;
            e = 0.05386179 - 0.00050991 * T;
            i = (2.48599187 + 0.00193609 * T) * constants::DEG_TO_RAD;
            L = (49.95424423 + 1222.49362201 * T) * constants::DEG_TO_RAD;
            omega_bar = (92.59887831 - 0.41897216 * T) * constants::DEG_TO_RAD;
            Omega = (113.66242448 - 0.28867794 * T) * constants::DEG_TO_RAD;
            break;
            
        case CelestialBody::URANUS:
            a = 19.18916464;
            e = 0.04725744 - 0.00004397 * T;
            i = (0.77263783 - 0.00242939 * T) * constants::DEG_TO_RAD;
            L = (313.23810451 + 428.48202785 * T) * constants::DEG_TO_RAD;
            omega_bar = (170.95427630 + 0.40805281 * T) * constants::DEG_TO_RAD;
            Omega = (74.01692503 + 0.04240589 * T) * constants::DEG_TO_RAD;
            break;
            
        case CelestialBody::NEPTUNE:
            a = 30.06992276;
            e = 0.00859048 + 0.00005105 * T;
            i = (1.77004347 + 0.00035372 * T) * constants::DEG_TO_RAD;
            L = (-55.12002969 + 218.45945325 * T) * constants::DEG_TO_RAD;
            omega_bar = (44.96476227 - 0.32241464 * T) * constants::DEG_TO_RAD;
            Omega = (131.78422574 - 0.00508664 * T) * constants::DEG_TO_RAD;
            break;
            
        case CelestialBody::PLUTO:
            // Simplified Pluto (low accuracy)
            a = 39.48211675;
            e = 0.24882730;
            i = 17.14001206 * constants::DEG_TO_RAD;
            L = (238.92903833 + 145.20780515 * T) * constants::DEG_TO_RAD;
            omega_bar = 224.06891629 * constants::DEG_TO_RAD;
            Omega = 110.30393684 * constants::DEG_TO_RAD;
            break;
            
        case CelestialBody::MOON:
            // Simplified lunar orbit (geocentric)
            // For heliocentric Moon, add Earth position
            a = 0.00256955529; // ~384,400 km
            e = 0.0549;
            i = 5.145 * constants::DEG_TO_RAD;
            L = (218.316 + 481267.881 * T) * constants::DEG_TO_RAD;
            omega_bar = 318.15 * constants::DEG_TO_RAD;
            Omega = (125.08 - 1934.136 * T) * constants::DEG_TO_RAD;
            break;
            
        default:
            // Unknown body
            a = 1.0;
            e = 0.0;
            i = 0.0;
            L = 0.0;
            omega_bar = 0.0;
            Omega = 0.0;
            break;
    }
    
    elements[0] = a;
    elements[1] = e;
    elements[2] = i;
    elements[3] = L;
    elements[4] = omega_bar;
    elements[5] = Omega;
}

Vector3d PlanetaryEphemeris::elementsToPosition(const double elements[6]) {
    double a = elements[0];
    double e = elements[1];
    double i = elements[2];
    double L = elements[3];
    double omega_bar = elements[4];
    double Omega = elements[5];
    
    // Mean anomaly
    double omega = omega_bar - Omega; // argument of perihelion
    double M = L - omega_bar; // mean anomaly
    M = math::normalize_angle_positive(M);
    
    // Solve Kepler's equation for eccentric anomaly E
    double E = M;
    for (int iter = 0; iter < 10; ++iter) {
        double dE = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-12) break;
    }
    
    // True anomaly
    double nu = 2.0 * std::atan2(
        std::sqrt(1.0 + e) * std::sin(E / 2.0),
        std::sqrt(1.0 - e) * std::cos(E / 2.0)
    );
    
    // Distance
    double r = a * (1.0 - e * std::cos(E));
    
    // Position in orbital plane
    double cos_nu = std::cos(nu);
    double sin_nu = std::sin(nu);
    double x_orb = r * cos_nu;
    double y_orb = r * sin_nu;
    
    // Rotation to ecliptic frame
    double cos_omega = std::cos(omega);
    double sin_omega = std::sin(omega);
    double cos_i = std::cos(i);
    double sin_i = std::sin(i);
    double cos_Omega = std::cos(Omega);
    double sin_Omega = std::sin(Omega);
    
    // Rotation matrices: R3(-Omega) * R1(-i) * R3(-omega)
    double x = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * x_orb +
               (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * y_orb;
    
    double y = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * x_orb +
               (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * y_orb;
    
    double z = (sin_omega * sin_i) * x_orb +
               (cos_omega * sin_i) * y_orb;
    
    Vector3d pos_ecliptic(x, y, z);
    
    // Convert Ecliptic J2000 -> Equatorial J2000
    // This rotation is critical for consistency with JPL ephemerides and the propagator frame
    return coordinates::ReferenceFrame::transform_position(
        pos_ecliptic, 
        coordinates::FrameType::ECLIPTIC, 
        coordinates::FrameType::J2000
    );
}

Vector3d PlanetaryEphemeris::elementsToVelocity(const double elements[6], double gm) {
    double a = elements[0];
    double e = elements[1];
    double i = elements[2];
    double L = elements[3];
    double omega_bar = elements[4];
    double Omega = elements[5];
    
    // Mean motion [rad/day]
    double n = std::sqrt(gm / (a * a * a));
    
    // Mean anomaly
    double omega = omega_bar - Omega;
    double M = L - omega_bar;
    M = math::normalize_angle_positive(M);
    
    // Solve Kepler's equation
    double E = M;
    for (int iter = 0; iter < 10; ++iter) {
        double dE = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-12) break;
    }
    
    // True anomaly
    double nu = 2.0 * std::atan2(
        std::sqrt(1.0 + e) * std::sin(E / 2.0),
        std::sqrt(1.0 - e) * std::cos(E / 2.0)
    );
    
    // Distance and velocity magnitude in orbital plane
    double r = a * (1.0 - e * std::cos(E));
    double h = std::sqrt(gm * a * (1.0 - e * e)); // specific angular momentum
    
    double cos_nu = std::cos(nu);
    double sin_nu = std::sin(nu);
    
    // Velocity in orbital plane [AU/day]
    double vx_orb = -(gm / h) * sin_nu;
    double vy_orb = (gm / h) * (e + cos_nu);
    
    // Rotation to ecliptic frame
    double cos_omega = std::cos(omega);
    double sin_omega = std::sin(omega);
    double cos_i = std::cos(i);
    double sin_i = std::sin(i);
    double cos_Omega = std::cos(Omega);
    double sin_Omega = std::sin(Omega);
    
    double vx = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * vx_orb +
                (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * vy_orb;
    
    double vy = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * vx_orb +
                (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * vy_orb;
    
    double vz = (sin_omega * sin_i) * vx_orb +
                (cos_omega * sin_i) * vy_orb;
    
    Vector3d vel_ecliptic(vx, vy, vz);
    
    // Convert Ecliptic J2000 -> Equatorial J2000
    // (Velocity rotation is same as position for inertial frames)
    return coordinates::ReferenceFrame::transform_velocity(
        Vector3d::Zero(), // Position doesn't matter for inertial rotation
        vel_ecliptic,
        coordinates::FrameType::ECLIPTIC,
        coordinates::FrameType::J2000
    );
}

Vector3d PlanetaryEphemeris::computePerturbations(CelestialBody body, double T) {
    // Placeholder for higher-order perturbations
    // Can be extended with VSOP87 series terms for improved accuracy
    return Vector3d::Zero();
}

} // namespace ephemeris
} // namespace astdyn
