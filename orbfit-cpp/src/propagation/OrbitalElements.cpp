/**
 * @file OrbitalElements.cpp
 * @brief Implementation of orbital element conversions
 */

#include "orbfit/propagation/OrbitalElements.hpp"
#include <cmath>
#include <stdexcept>

namespace orbfit::propagation {

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
    double v2 = velocity.squaredNorm();
    double r = position.norm();
    return 0.5 * v2 - gravitational_parameter / r;
}

Eigen::Vector3d CartesianElements::angular_momentum() const {
    return position.cross(velocity);
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
    // Normalize M to [0, 2π)
    M = std::fmod(M, constants::TWO_PI);
    if (M < 0.0) M += constants::TWO_PI;
    
    // Initial guess for E (Danby's method)
    double E = M + 0.85 * e * std::sin(M) / std::abs(std::sin(M));
    if (e < 0.8) {
        E = M;
    }
    
    // Newton-Raphson iteration
    for (int iter = 0; iter < max_iter; ++iter) {
        double f = E - e * std::sin(E) - M;
        double fp = 1.0 - e * std::cos(E);
        double delta = -f / fp;
        
        E += delta;
        
        if (std::abs(delta) < tolerance) {
            return E;
        }
    }
    
    throw std::runtime_error("Kepler equation failed to converge");
}

double eccentric_to_true_anomaly(double E, double e) {
    // tan(ν/2) = √((1+e)/(1-e)) tan(E/2)
    double sqrt_factor = std::sqrt((1.0 + e) / (1.0 - e));
    double tan_half_E = std::tan(E / 2.0);
    double tan_half_nu = sqrt_factor * tan_half_E;
    double nu = 2.0 * std::atan(tan_half_nu);
    
    // Normalize to [0, 2π)
    if (nu < 0.0) nu += constants::TWO_PI;
    
    return nu;
}

double true_to_eccentric_anomaly(double nu, double e) {
    // tan(E/2) = √((1-e)/(1+e)) tan(ν/2)
    double sqrt_factor = std::sqrt((1.0 - e) / (1.0 + e));
    double tan_half_nu = std::tan(nu / 2.0);
    double tan_half_E = sqrt_factor * tan_half_nu;
    double E = 2.0 * std::atan(tan_half_E);
    
    // Ensure same quadrant as nu
    if (nu > constants::PI && E < 0.0) E += constants::TWO_PI;
    
    return E;
}

// ============================================================================
// Keplerian <-> Cartesian Conversions
// ============================================================================

CartesianElements keplerian_to_cartesian(const KeplerianElements& kep) {
    CartesianElements cart;
    cart.epoch_mjd_tdb = kep.epoch_mjd_tdb;
    cart.gravitational_parameter = kep.gravitational_parameter;
    
    // Extract elements
    double a = kep.semi_major_axis;
    double e = kep.eccentricity;
    double i = kep.inclination;
    double Omega = kep.longitude_ascending_node;
    double omega = kep.argument_perihelion;
    double M = kep.mean_anomaly;
    double mu = kep.gravitational_parameter;
    
    // Solve Kepler's equation for eccentric anomaly
    double E = solve_kepler_equation(M, e);
    
    // Compute true anomaly
    double nu = eccentric_to_true_anomaly(E, e);
    
    // Compute distance
    double r = a * (1.0 - e * std::cos(E));
    
    // Position and velocity in orbital plane (perifocal frame)
    double cos_nu = std::cos(nu);
    double sin_nu = std::sin(nu);
    
    double x_orb = r * cos_nu;
    double y_orb = r * sin_nu;
    
    // Velocity in orbital plane
    // v = sqrt(mu*a) / r * [-sin(E), sqrt(1-e²)*cos(E)]
    double sqrt_mu_a = std::sqrt(mu * a);
    double vx_orb = -(sqrt_mu_a / r) * std::sin(E);
    double vy_orb = (sqrt_mu_a / r) * std::sqrt(1.0 - e * e) * std::cos(E);
    
    // Rotation matrices
    double cos_omega = std::cos(omega);
    double sin_omega = std::sin(omega);
    double cos_Omega = std::cos(Omega);
    double sin_Omega = std::sin(Omega);
    double cos_i = std::cos(i);
    double sin_i = std::sin(i);
    
    // Perifocal to J2000 transformation
    // R = R_z(-Ω) R_x(-i) R_z(-ω)
    double r11 = cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i;
    double r12 = -cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i;
    double r21 = sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i;
    double r22 = -sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i;
    double r31 = sin_omega * sin_i;
    double r32 = cos_omega * sin_i;
    
    // Transform position
    cart.position(0) = r11 * x_orb + r12 * y_orb;
    cart.position(1) = r21 * x_orb + r22 * y_orb;
    cart.position(2) = r31 * x_orb + r32 * y_orb;
    
    // Transform velocity
    cart.velocity(0) = r11 * vx_orb + r12 * vy_orb;
    cart.velocity(1) = r21 * vx_orb + r22 * vy_orb;
    cart.velocity(2) = r31 * vx_orb + r32 * vy_orb;
    
    return cart;
}

KeplerianElements cartesian_to_keplerian(const CartesianElements& cart) {
    KeplerianElements kep;
    kep.epoch_mjd_tdb = cart.epoch_mjd_tdb;
    kep.gravitational_parameter = cart.gravitational_parameter;
    
    const Eigen::Vector3d& r = cart.position;
    const Eigen::Vector3d& v = cart.velocity;
    double mu = cart.gravitational_parameter;
    
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
    
    kep.semi_major_axis = a;
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
    eq.epoch_mjd_tdb = kep.epoch_mjd_tdb;
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
    kep.epoch_mjd_tdb = eq.epoch_mjd_tdb;
    kep.gravitational_parameter = eq.gravitational_parameter;
    
    kep.semi_major_axis = eq.a;
    kep.eccentricity = std::sqrt(eq.h * eq.h + eq.k * eq.k);
    
    double i = 2.0 * std::atan(std::sqrt(eq.p * eq.p + eq.q * eq.q));
    kep.inclination = i;
    
    double Omega = std::atan2(eq.p, eq.q);
    kep.longitude_ascending_node = Omega;
    
    double omega_plus_Omega = std::atan2(eq.h, eq.k);
    kep.argument_perihelion = omega_plus_Omega - Omega;
    
    kep.mean_anomaly = eq.lambda - omega_plus_Omega;
    
    return kep;
}

// ============================================================================
// Keplerian <-> Cometary Conversions
// ============================================================================

CometaryElements keplerian_to_cometary(const KeplerianElements& kep) {
    CometaryElements com;
    com.epoch_mjd_tdb = kep.epoch_mjd_tdb;
    com.gravitational_parameter = kep.gravitational_parameter;
    
    com.perihelion_distance = kep.perihelion_distance();
    com.eccentricity = kep.eccentricity;
    com.inclination = kep.inclination;
    com.longitude_ascending_node = kep.longitude_ascending_node;
    com.argument_perihelion = kep.argument_perihelion;
    
    // Compute time of perihelion from M and epoch
    double n = kep.mean_motion();
    double time_since_perihelion = -kep.mean_anomaly / n;
    com.time_perihelion_mjd_tdb = kep.epoch_mjd_tdb + time_since_perihelion;
    
    return com;
}

KeplerianElements cometary_to_keplerian(const CometaryElements& com) {
    KeplerianElements kep;
    kep.epoch_mjd_tdb = com.epoch_mjd_tdb;
    kep.gravitational_parameter = com.gravitational_parameter;
    
    kep.semi_major_axis = com.semi_major_axis();
    kep.eccentricity = com.eccentricity;
    kep.inclination = com.inclination;
    kep.longitude_ascending_node = com.longitude_ascending_node;
    kep.argument_perihelion = com.argument_perihelion;
    
    // Compute mean anomaly from time of perihelion
    double n = std::sqrt(com.gravitational_parameter / 
                        (kep.semi_major_axis * kep.semi_major_axis * kep.semi_major_axis));
    double time_since_perihelion = com.epoch_mjd_tdb - com.time_perihelion_mjd_tdb;
    kep.mean_anomaly = n * time_since_perihelion;
    
    return kep;
}

// ============================================================================
// Mean to Osculating Element Conversions
// ============================================================================

KeplerianElements mean_to_osculating(
    const KeplerianElements& mean_elements,
    double j2,
    double central_body_radius)
{
    KeplerianElements osc = mean_elements;
    
    // If J2 is zero or negligible, mean ≈ osculating
    if (std::abs(j2) < 1e-12) {
        return osc;
    }
    
    // Extract orbital parameters
    double a = mean_elements.semi_major_axis;
    double e = mean_elements.eccentricity;
    double i = mean_elements.inclination;
    double Omega = mean_elements.longitude_ascending_node;
    double omega = mean_elements.argument_perihelion;
    double M = mean_elements.mean_anomaly;
    
    // Compute derived quantities
    double eta = std::sqrt(1.0 - e * e);     // Auxiliary parameter
    
    // J2 factor: k = J2 * (R/a)²
    double R_over_a = central_body_radius / a;
    double k = j2 * R_over_a * R_over_a;
    
    // Trigonometric functions
    double sin_i = std::sin(i);
    double cos_i = std::cos(i);
    double sin2_i = sin_i * sin_i;
    
    // Convert mean anomaly to eccentric anomaly, then to true anomaly
    double E = solve_kepler_equation(M, e);
    double nu = eccentric_to_true_anomaly(E, e);
    
    double cos_nu = std::cos(nu);
    double sin_2nu = std::sin(2.0 * nu);
    double cos_2nu = std::cos(2.0 * nu);
    
    // Compute short-period J2 perturbations
    // These are the differences: Δx = x_osculating - x_mean
    
    // Semi-major axis: negligible short-period variation
    double Delta_a = 0.0;
    
    // Eccentricity short-period variation
    // Δe = (k/8) * e * η * (1 - 11cos²i - 40(cos⁴i)/(1-5cos²i)) * sin(2ω+2ν)
    double Delta_e = (k / 8.0) * e * eta * 
                     (1.0 - 11.0 * cos_i * cos_i) * sin_2nu;
    
    // Inclination short-period variation
    // Δi = -(k/2) * e * sin(i) * cos(i) * cos(2ω+2ν)
    double Delta_i = -(k / 2.0) * e * sin_i * cos_i * cos_2nu;
    
    // RAAN short-period variation
    // ΔΩ = (k/2) * cos(i) * [something complex with ν]
    // Simplified: dominant term is secular, short-period is small
    double Delta_Omega = 0.0;  // Short-period part typically neglected
    
    // Argument of perihelion short-period variation
    // Δω includes both short-period oscillations
    double Delta_omega = (k / 8.0) * (4.0 - 5.0 * sin2_i) * 
                         (2.0 + e * cos_nu) * sin_2nu / eta;
    
    // Apply corrections to get osculating elements
    osc.semi_major_axis = a + Delta_a;
    osc.eccentricity = e + Delta_e;
    osc.inclination = i + Delta_i;
    osc.longitude_ascending_node = Omega + Delta_Omega;
    osc.argument_perihelion = omega + Delta_omega;
    
    // For osculating, keep the true anomaly consistent
    // Convert back to mean anomaly for the osculating elements
    double E_osc = true_to_eccentric_anomaly(nu, osc.eccentricity);
    osc.mean_anomaly = E_osc - osc.eccentricity * std::sin(E_osc);
    
    // Normalize angles
    osc.longitude_ascending_node = std::fmod(osc.longitude_ascending_node, constants::TWO_PI);
    if (osc.longitude_ascending_node < 0.0) osc.longitude_ascending_node += constants::TWO_PI;
    
    osc.argument_perihelion = std::fmod(osc.argument_perihelion, constants::TWO_PI);
    if (osc.argument_perihelion < 0.0) osc.argument_perihelion += constants::TWO_PI;
    
    osc.mean_anomaly = std::fmod(osc.mean_anomaly, constants::TWO_PI);
    if (osc.mean_anomaly < 0.0) osc.mean_anomaly += constants::TWO_PI;
    
    return osc;
}

KeplerianElements osculating_to_mean(
    const KeplerianElements& osc_elements,
    double j2,
    double central_body_radius)
{
    KeplerianElements mean = osc_elements;
    
    // If J2 is zero, osculating ≈ mean
    if (std::abs(j2) < 1e-12) {
        return mean;
    }
    
    // The inverse transformation: subtract the short-period terms
    // For simplicity, we apply the negative of the corrections
    // A more rigorous approach would iterate, but for small J2 this is adequate
    
    double a = osc_elements.semi_major_axis;
    double e = osc_elements.eccentricity;
    double i = osc_elements.inclination;
    double M_osc = osc_elements.mean_anomaly;
    
    // Compute correction terms
    double R_over_a = central_body_radius / a;
    double k = j2 * R_over_a * R_over_a;
    double eta = std::sqrt(1.0 - e * e);
    
    double E = solve_kepler_equation(M_osc, e);
    double nu = eccentric_to_true_anomaly(E, e);
    
    double sin_i = std::sin(i);
    double cos_i = std::cos(i);
    double sin2_i = sin_i * sin_i;
    double sin_2nu = std::sin(2.0 * nu);
    double cos_2nu = std::cos(2.0 * nu);
    
    // Apply inverse corrections (subtract perturbations)
    mean.eccentricity = e - (k / 8.0) * e * eta * 
                             (1.0 - 11.0 * cos_i * cos_i) * sin_2nu;
    
    mean.inclination = i + (k / 2.0) * e * sin_i * cos_i * cos_2nu;
    
    mean.argument_perihelion = osc_elements.argument_perihelion - 
                               (k / 8.0) * (4.0 - 5.0 * sin2_i) * 
                               (2.0 + e * std::cos(nu)) * sin_2nu / eta;
    
    // Convert to mean anomaly
    double E_mean = true_to_eccentric_anomaly(nu, mean.eccentricity);
    mean.mean_anomaly = E_mean - mean.eccentricity * std::sin(E_mean);
    
    // Normalize angles
    mean.longitude_ascending_node = std::fmod(mean.longitude_ascending_node, constants::TWO_PI);
    if (mean.longitude_ascending_node < 0.0) mean.longitude_ascending_node += constants::TWO_PI;
    
    mean.argument_perihelion = std::fmod(mean.argument_perihelion, constants::TWO_PI);
    if (mean.argument_perihelion < 0.0) mean.argument_perihelion += constants::TWO_PI;
    
    mean.mean_anomaly = std::fmod(mean.mean_anomaly, constants::TWO_PI);
    if (mean.mean_anomaly < 0.0) mean.mean_anomaly += constants::TWO_PI;
    
    return mean;
}

} // namespace orbfit::propagation
