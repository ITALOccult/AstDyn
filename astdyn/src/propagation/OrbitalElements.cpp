/**
 * @file OrbitalElements.cpp
 * @brief Implementation of orbital element conversions
 */

#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/propagation/PlanetaryPeriodicPerturbations.hpp"

// New AstDyn Architecture (CTFYH)
#include "astdyn/core/units.hpp"
#include "astdyn/types/orbital_state.hpp"
#include "astdyn/math/kepler_solver.hpp"
#include "astdyn/math/anomaly_conversions.hpp"
#include "astdyn/coordinates/state_conversions.hpp"
#include "astdyn/core/Constants.hpp"

#include <cmath>
#include <stdexcept>

namespace astdyn::propagation {
using namespace astdyn::constants;

// ============================================================================
// KeplerianElements methods
// ============================================================================

double KeplerianElements::period() const {
    if (!is_elliptic()) {
        throw std::domain_error("Period undefined for non-elliptic orbit");
    }
    return constants::TWO_PI * std::sqrt(semi_major_axis * semi_major_axis * semi_major_axis 
                                        / gravitational_parameter);
}

double KeplerianElements::mean_motion() const {
    if (!is_elliptic()) {
        throw std::domain_error("Mean motion undefined for non-elliptic orbit");
    }
    return std::sqrt(gravitational_parameter / 
                    (semi_major_axis * semi_major_axis * semi_major_axis));
}

double KeplerianElements::perihelion_distance() const {
    return semi_major_axis * (1.0 - eccentricity);
}

double KeplerianElements::aphelion_distance() const {
    if (!is_elliptic()) {
        throw std::domain_error("Aphelion undefined for non-elliptic orbit");
    }
    return semi_major_axis * (1.0 + eccentricity);
}

// ============================================================================
// CartesianElements methods
// ============================================================================

double CartesianElements::energy() const {
    double v = velocity.norm().to_au_d();
    double r = position.norm().to_au();
    double AU_M = constants::AU * 1000.0;
    double m3s2_to_au3d2 = (86400.0 * 86400.0) / (AU_M * AU_M * AU_M);
    double mu = gravitational_parameter * m3s2_to_au3d2;
    return 0.5 * v * v - mu / r;
}

Eigen::Vector3d CartesianElements::angular_momentum() const {
    auto r_au = position.to_eigen_si() / (constants::AU * 1000.0);
    auto v_aud = velocity.to_eigen_si() * (86400.0 / (constants::AU * 1000.0));
    return r_au.cross(v_aud);
}

// ============================================================================
// CometaryElements methods
// ============================================================================

double CometaryElements::semi_major_axis() const {
    if (is_parabolic()) {
        throw std::domain_error("Semi-major axis infinite for parabolic orbit");
    }
    return perihelion_distance / (1.0 - eccentricity);
}

// ============================================================================
// Kepler's Equation Solver
// ============================================================================

double solve_kepler_equation(double M, double e, double tolerance, int max_iter) {
    using astdyn::core::Radian;
    using astdyn::math::KeplerSolverOptions;

    const auto options = KeplerSolverOptions{e, tolerance, max_iter};
    const auto result = astdyn::math::solve_kepler_elliptic(Radian(M), options);

    if (!result.has_value()) {
        throw std::runtime_error("Kepler equation failed to converge (AstDyn Kernel)");
    }

    return result->value;
}

double eccentric_to_true_anomaly(double E, double e) {
    using astdyn::core::Radian;
    const auto nu = astdyn::math::eccentric_to_true_anomaly(Radian(E), e);
    return nu.value;
}

double true_to_eccentric_anomaly(double nu, double e) {
    using astdyn::core::Radian;
    const auto E = astdyn::math::true_to_eccentric_anomaly(Radian(nu), e);
    return E.value;
}

// ============================================================================
// Keplerian <-> Cartesian Conversions
// ============================================================================

static CartesianElements finalize_cartesian(const KeplerianElements& kep, const std::array<double, 6>& raw) {
    CartesianElements cart; cart.epoch = kep.epoch;
    double AU_M = constants::AU * 1000.0, AUd_MS = AU_M / 86400.0;
    double au3d2_to_m3s2 = (AU_M * AU_M * AU_M) / (86400.0 * 86400.0);
    cart.gravitational_parameter = kep.gravitational_parameter * au3d2_to_m3s2;
    cart.position = math::Vector3<core::GCRF, physics::Distance>::from_si(raw[0]*AU_M, raw[1]*AU_M, raw[2]*AU_M);
    cart.velocity = math::Vector3<core::GCRF, physics::Velocity>::from_si(raw[3]*AUd_MS, raw[4]*AUd_MS, raw[5]*AUd_MS);
    return cart;
}

CartesianElements keplerian_to_cartesian(const KeplerianElements& kep) {
    using namespace astdyn::types;
    const auto state_kep = OrbitalState<core::GCRF, KeplerianTag>({
        kep.semi_major_axis, kep.eccentricity, kep.inclination,
        kep.longitude_ascending_node, kep.argument_perihelion, kep.mean_anomaly
    });
    const auto state_cart = astdyn::coordinates::keplerian_to_cartesian(state_kep, kep.gravitational_parameter);
    return finalize_cartesian(kep, state_cart.raw_values());
}

KeplerianElements cartesian_to_keplerian(const CartesianElements& cart) {
    KeplerianElements kep;
    kep.epoch = cart.epoch;
    // CartesianElements.gravitational_parameter is in m³/s² (SI).
    // KeplerianElements.gravitational_parameter is in AU³/day².
    constexpr double AU_M = AU * 1000.0;
    constexpr double M3s2_TO_AU3d2 = (86400.0 * 86400.0) / (AU_M * AU_M * AU_M);
    kep.gravitational_parameter = cart.gravitational_parameter * M3s2_TO_AU3d2;

    double mu = cart.gravitational_parameter; // m³/s²

    // CartesianElements is always in SI (m, m/s, m³/s²) — no heuristic needed.
    Eigen::Vector3d r(cart.position.x_si(), cart.position.y_si(), cart.position.z_si()); // [m]
    Eigen::Vector3d v(cart.velocity.x_si(), cart.velocity.y_si(), cart.velocity.z_si()); // [m/s]
    
    double r_mag = r.norm();
    double v_mag = v.norm();
    
    // Angular momentum
    Eigen::Vector3d h = r.cross(v);
    double h_mag = h.norm();
    
    // Eccentricity vector
    double rdot = r.dot(v);
    Eigen::Vector3d e_vec = ((v_mag * v_mag - mu / r_mag) * r - rdot * v) / mu;
    double e = e_vec.norm();
    
    // Semi-major axis from specific energy
    // ε = v²/2 - μ/r = -μ/(2a)  =>  a = -μ/(2ε)
    double specific_energy = 0.5 * v_mag * v_mag - mu / r_mag;
    double a = -mu / (2.0 * specific_energy);
    // CartesianElements is always SI (m, m/s, m³/s²), so 'a' is in meters.
    // KeplerianElements.semi_major_axis is always in AU. AU_M defined above.
    double a_au = a / AU_M;
    
    kep.semi_major_axis = a_au;
    kep.eccentricity = e;
    
    // Inclination
    double i = std::acos(h(2) / h_mag);
    kep.inclination = i;
    
    // Node vector (k × h)
    Eigen::Vector3d k(0, 0, 1);
    Eigen::Vector3d n = k.cross(h);
    double n_mag = n.norm();
    
    // Longitude of ascending node
    double Omega = 0.0;
    if (n_mag > 1e-10) {
        Omega = std::acos(n(0) / n_mag);
        if (n(1) < 0.0) {
            Omega = constants::TWO_PI - Omega;
        }
    }
    kep.longitude_ascending_node = Omega;
    
    // Argument of perihelion
    double omega = 0.0;
    if (n_mag > 1e-10 && e > 1e-10) {
        omega = std::acos(n.dot(e_vec) / (n_mag * e));
        if (e_vec(2) < 0.0) {
            omega = constants::TWO_PI - omega;
        }
    }
    kep.argument_perihelion = omega;
    
    // True anomaly
    double nu = 0.0;
    if (e > 1e-10) {
        nu = std::acos(e_vec.dot(r) / (e * r_mag));
        if (r.dot(v) < 0.0) {
            nu = constants::TWO_PI - nu;
        }
    }
    
    // Eccentric anomaly and mean anomaly
    double E = true_to_eccentric_anomaly(nu, e);
    double M = E - e * std::sin(E);
    kep.mean_anomaly = M;
    
    return kep;
}

// ============================================================================
// Keplerian <-> Equinoctial Conversions
// ============================================================================

EquinoctialElements keplerian_to_equinoctial(const KeplerianElements& kep) {
    EquinoctialElements eq;
    eq.epoch = kep.epoch;
    eq.gravitational_parameter = kep.gravitational_parameter;
    
    double e = kep.eccentricity;
    double i = kep.inclination;
    double Omega = kep.longitude_ascending_node;
    double omega = kep.argument_perihelion;
    double M = kep.mean_anomaly;
    
    eq.a = kep.semi_major_axis;
    eq.h = e * std::sin(omega + Omega);
    eq.k = e * std::cos(omega + Omega);
    eq.p = std::tan(i / 2.0) * std::sin(Omega);
    eq.q = std::tan(i / 2.0) * std::cos(Omega);
    eq.lambda = M + omega + Omega;
    
    return eq;
}

KeplerianElements equinoctial_to_keplerian(const EquinoctialElements& eq) {
    KeplerianElements kep;
    kep.epoch = eq.epoch;
    kep.gravitational_parameter = eq.gravitational_parameter;
    
    kep.semi_major_axis = eq.a;
    kep.eccentricity = std::sqrt(eq.h * eq.h + eq.k * eq.k);
    
    double i = 2.0 * std::atan(std::sqrt(eq.p * eq.p + eq.q * eq.q));
    kep.inclination = i;
    
    double Omega = std::atan2(eq.p, eq.q);
    kep.longitude_ascending_node = Omega;
    
    double omega_plus_Omega = std::atan2(eq.h, eq.k);
    kep.argument_perihelion = std::fmod(omega_plus_Omega - Omega, constants::TWO_PI);
    if (kep.argument_perihelion < 0) kep.argument_perihelion += constants::TWO_PI;
    
    kep.mean_anomaly = std::fmod(eq.lambda - omega_plus_Omega, constants::TWO_PI);
    if (kep.mean_anomaly < 0) kep.mean_anomaly += constants::TWO_PI;
    
    kep.longitude_ascending_node = std::fmod(kep.longitude_ascending_node, constants::TWO_PI);
    if (kep.longitude_ascending_node < 0) kep.longitude_ascending_node += constants::TWO_PI;
    
    return kep;
}

// ============================================================================
// Keplerian <-> Cometary Conversions
// ============================================================================

CometaryElements keplerian_to_cometary(const KeplerianElements& kep) {
    CometaryElements com;
    com.epoch = kep.epoch;
    com.gravitational_parameter = kep.gravitational_parameter;
    
    com.perihelion_distance = kep.perihelion_distance();
    com.eccentricity = kep.eccentricity;
    com.inclination = kep.inclination;
    com.longitude_ascending_node = kep.longitude_ascending_node;
    com.argument_perihelion = kep.argument_perihelion;
    
    // Compute time of perihelion from M and epoch
    double n = kep.mean_motion();
    double time_since_perihelion = -kep.mean_anomaly / n;
    com.time_perihelion = time::EpochTDB::from_mjd(kep.epoch.mjd() + time_since_perihelion);
    
    return com;
}

KeplerianElements cometary_to_keplerian(const CometaryElements& com) {
    KeplerianElements kep;
    kep.epoch = com.epoch;
    kep.gravitational_parameter = com.gravitational_parameter;
    
    kep.semi_major_axis = com.semi_major_axis();
    kep.eccentricity = com.eccentricity;
    kep.inclination = com.inclination;
    kep.longitude_ascending_node = com.longitude_ascending_node;
    kep.argument_perihelion = com.argument_perihelion;
    
    // Compute mean anomaly from time of perihelion
    double n = std::sqrt(com.gravitational_parameter / 
                        (kep.semi_major_axis * kep.semi_major_axis * kep.semi_major_axis));
    double time_since_perihelion = com.epoch.mjd() - com.time_perihelion.mjd();
    kep.mean_anomaly = n * time_since_perihelion;
    
    return kep;
}

static void apply_planetary_corrections(KeplerianElements& osc, const KeplerianElements& mean) {
    PlanetaryPeriodicPerturbations mk;
    if (mean.gravitational_parameter < 1e-10) const_cast<KeplerianElements&>(mean).gravitational_parameter = GMS;
    auto corr = mk.calculateCorrections(mean, mean.epoch.mjd(), false);
    osc.semi_major_axis += corr[0];
}

KeplerianElements mean_to_osculating(const KeplerianElements& mean, double j2, double R) {
    KeplerianElements osc = mean;
    double a = mean.semi_major_axis, e = mean.eccentricity, i = mean.inclination, M = mean.mean_anomaly;
    double k = j2 * (R/a) * (R/a), eta = std::sqrt(1.0 - e*e);
    double nu = eccentric_to_true_anomaly(solve_kepler_equation(M, e), e);
    double sin_i = std::sin(i), cos_i = std::cos(i), sin2nu = std::sin(2.0*nu), cos2nu = std::cos(2.0*nu);
    osc.eccentricity += (k/8.0) * e * eta * (1.0 - 11.0*cos_i*cos_i) * sin2nu;
    osc.inclination += -(k/2.0) * e * sin_i * cos_i * cos2nu;
    osc.argument_perihelion += (k/8.0) * (4.0 - 5.0*sin_i*sin_i) * (2.0 + e*std::cos(nu)) * sin2nu / eta;
    if (a > 1.8 && a < 4.0) apply_planetary_corrections(osc, mean);
    return osc;
}

KeplerianElements osculating_to_mean(
    const KeplerianElements& osc_elements,
    double j2,
    double central_body_radius)
{
    KeplerianElements mean = osc_elements;
    
    // Iterate to find mean elements that produce the given osculating elements
    // This handles any complex perturbations implemented in mean_to_osculating
    for (int iter = 0; iter < 10; ++iter) {
        KeplerianElements osc_test = mean_to_osculating(mean, j2, central_body_radius);
        
        mean.semi_major_axis += (osc_elements.semi_major_axis - osc_test.semi_major_axis);
        mean.eccentricity     += (osc_elements.eccentricity     - osc_test.eccentricity);
        mean.inclination      += (osc_elements.inclination      - osc_test.inclination);
        
        double d_Omega = osc_elements.longitude_ascending_node - osc_test.longitude_ascending_node;
        while (d_Omega >  constants::PI) d_Omega -= constants::TWO_PI;
        while (d_Omega < -constants::PI) d_Omega += constants::TWO_PI;
        mean.longitude_ascending_node += d_Omega;
        
        double d_omega = osc_elements.argument_perihelion - osc_test.argument_perihelion;
        while (d_omega >  constants::PI) d_omega -= constants::TWO_PI;
        while (d_omega < -constants::PI) d_omega += constants::TWO_PI;
        mean.argument_perihelion += d_omega;
        
        double d_M = osc_elements.mean_anomaly - osc_test.mean_anomaly;
        while (d_M >  constants::PI) d_M -= constants::TWO_PI;
        while (d_M < -constants::PI) d_M += constants::TWO_PI;
        mean.mean_anomaly += d_M;
    }
    
    // Normalize angles
    mean.longitude_ascending_node = std::fmod(mean.longitude_ascending_node, constants::TWO_PI);
    if (mean.longitude_ascending_node < 0.0) mean.longitude_ascending_node += constants::TWO_PI;
    
    mean.argument_perihelion = std::fmod(mean.argument_perihelion, constants::TWO_PI);
    if (mean.argument_perihelion < 0.0) mean.argument_perihelion += constants::TWO_PI;
    
    mean.mean_anomaly = std::fmod(mean.mean_anomaly, constants::TWO_PI);
    if (mean.mean_anomaly < 0.0) mean.mean_anomaly += constants::TWO_PI;

    return mean;
}

} // namespace astdyn::propagation
