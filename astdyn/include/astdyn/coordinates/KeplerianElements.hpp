/**
 * @file KeplerianElements.hpp
 * @brief Keplerian orbital elements representation
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Classical Keplerian orbital elements: a, e, i, Ω, ω, M (or ν, E).
 * Includes conversions to/from Cartesian coordinates.
 */

#ifndef ASTDYN_COORDINATES_KEPLERIANELEMENTS_HPP
#define ASTDYN_COORDINATES_KEPLERIANELEMENTS_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include <cmath>
#include <optional>

namespace astdyn {
namespace coordinates {

/**
 * @brief Classical Keplerian orbital elements
 * 
 * Elements:
 * - a: Semi-major axis [km]
 * - e: Eccentricity [dimensionless]
 * - i: Inclination [rad]
 * - Ω (RAAN): Right ascension of ascending node [rad]
 * - ω (omega): Argument of periapsis [rad]
 * - M: Mean anomaly [rad] (or ν for true anomaly, E for eccentric)
 */
class KeplerianElements {
public:
    // ========================================================================
    // Constructors
    // ========================================================================
    
    /**
     * @brief Default constructor (circular equatorial orbit at 1 AU)
     */
    KeplerianElements()
        : a_(constants::AU),
          e_(0.0),
          i_(0.0),
          Omega_(0.0),
          omega_(0.0),
          M_(0.0),
          mu_(constants::GM_SUN) {}
    
    /**
     * @brief Construct from elements
     * @param a Semi-major axis [km]
     * @param e Eccentricity
     * @param i Inclination [rad]
     * @param Omega RAAN [rad]
     * @param omega Argument of periapsis [rad]
     * @param M Mean anomaly [rad]
     * @param mu Gravitational parameter [km³/s²]
     */
    KeplerianElements(double a, double e, double i,
                     double Omega, double omega, double M,
                     double mu = constants::GM_SUN)
        : a_(a), e_(e), i_(i),
          Omega_(Omega), omega_(omega), M_(M),
          mu_(mu) {}
    
    // ========================================================================
    // Accessors
    // ========================================================================
    
    double semi_major_axis() const { return a_; }
    double eccentricity() const { return e_; }
    double inclination() const { return i_; }
    double RAAN() const { return Omega_; }
    double argument_of_periapsis() const { return omega_; }
    double mean_anomaly() const { return M_; }
    double mu() const { return mu_; }
    
    void set_semi_major_axis(double a) { a_ = a; }
    void set_eccentricity(double e) { e_ = e; }
    void set_inclination(double i) { i_ = i; }
    void set_RAAN(double Omega) { Omega_ = Omega; }
    void set_argument_of_periapsis(double omega) { omega_ = omega; }
    void set_mean_anomaly(double M) { M_ = M; }
    void set_mu(double mu) { mu_ = mu; }
    
    // ========================================================================
    // Derived Quantities
    // ========================================================================
    
    /**
     * @brief Periapsis distance q = a(1 - e)
     * @return Periapsis [km]
     */
    double periapsis_distance() const {
        return a_ * (1.0 - e_);
    }
    
    /**
     * @brief Apoapsis distance Q = a(1 + e)
     * @return Apoapsis [km] (infinity if e >= 1)
     */
    double apoapsis_distance() const {
        if (e_ >= 1.0) {
            return std::numeric_limits<double>::infinity();
        }
        return a_ * (1.0 + e_);
    }
    
    /**
     * @brief Orbital period T = 2π√(a³/μ)
     * @return Period [s] (infinity if e >= 1)
     */
    double period() const {
        if (e_ >= 1.0 || a_ <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        return 2.0 * constants::PI * std::sqrt(a_ * a_ * a_ / mu_);
    }
    
    /**
     * @brief Mean motion n = √(μ/a³)
     * @return Mean motion [rad/s]
     */
    double mean_motion() const {
        return std::sqrt(mu_ / (a_ * a_ * a_));
    }
    
    /**
     * @brief Specific orbital energy ε = -μ/(2a)
     * @return Energy [km²/s²]
     */
    double specific_energy() const {
        return -mu_ / (2.0 * a_);
    }
    
    /**
     * @brief Semi-latus rectum p = a(1 - e²)
     * @return Semi-latus rectum [km]
     */
    double semi_latus_rectum() const {
        return a_ * (1.0 - e_ * e_);
    }
    
    // ========================================================================
    // Anomaly Conversions
    // ========================================================================
    
    /**
     * @brief Convert mean anomaly to eccentric anomaly (Kepler's equation)
     * @param M Mean anomaly [rad]
     * @param e Eccentricity
     * @param tolerance Convergence tolerance
     * @param max_iter Maximum iterations
     * @return Eccentric anomaly [rad]
     */
    static double mean_to_eccentric_anomaly(double M, double e,
                                           double tolerance = 1e-12,
                                           int max_iter = 50) {
        // Newton-Raphson solver for E - e*sin(E) = M
        double E = M; // Initial guess
        
        for (int i = 0; i < max_iter; ++i) {
            double f = E - e * std::sin(E) - M;
            double fp = 1.0 - e * std::cos(E);
            double delta = f / fp;
            E -= delta;
            
            if (std::abs(delta) < tolerance) {
                break;
            }
        }
        
        return E;
    }
    
    /**
     * @brief Convert eccentric anomaly to mean anomaly
     * @param E Eccentric anomaly [rad]
     * @param e Eccentricity
     * @return Mean anomaly [rad]
     */
    static double eccentric_to_mean_anomaly(double E, double e) {
        return E - e * std::sin(E);
    }
    
    /**
     * @brief Convert eccentric anomaly to true anomaly
     * @param E Eccentric anomaly [rad]
     * @param e Eccentricity
     * @return True anomaly [rad]
     */
    static double eccentric_to_true_anomaly(double E, double e) {
        double sin_nu = std::sqrt(1.0 - e * e) * std::sin(E) / (1.0 - e * std::cos(E));
        double cos_nu = (std::cos(E) - e) / (1.0 - e * std::cos(E));
        return std::atan2(sin_nu, cos_nu);
    }
    
    /**
     * @brief Convert true anomaly to eccentric anomaly
     * @param nu True anomaly [rad]
     * @param e Eccentricity
     * @return Eccentric anomaly [rad]
     */
    static double true_to_eccentric_anomaly(double nu, double e) {
        double cos_E = (e + std::cos(nu)) / (1.0 + e * std::cos(nu));
        double sin_E = std::sqrt(1.0 - e * e) * std::sin(nu) / (1.0 + e * std::cos(nu));
        return std::atan2(sin_E, cos_E);
    }
    
    /**
     * @brief Convert mean anomaly to true anomaly
     * @param M Mean anomaly [rad]
     * @param e Eccentricity
     * @return True anomaly [rad]
     */
    static double mean_to_true_anomaly(double M, double e) {
        double E = mean_to_eccentric_anomaly(M, e);
        return eccentric_to_true_anomaly(E, e);
    }
    
    /**
     * @brief Convert true anomaly to mean anomaly
     * @param nu True anomaly [rad]
     * @param e Eccentricity
     * @return Mean anomaly [rad]
     */
    static double true_to_mean_anomaly(double nu, double e) {
        double E = true_to_eccentric_anomaly(nu, e);
        return eccentric_to_mean_anomaly(E, e);
    }
    
    /**
     * @brief Get current eccentric anomaly
     * @return E [rad]
     */
    double eccentric_anomaly() const {
        return mean_to_eccentric_anomaly(M_, e_);
    }
    
    /**
     * @brief Get current true anomaly
     * @return ν [rad]
     */
    double true_anomaly() const {
        return mean_to_true_anomaly(M_, e_);
    }
    
    // ========================================================================
    // Conversions to/from Cartesian
    // ========================================================================
    
    /**
     * @brief Convert to Cartesian state
     * @return Cartesian state (position, velocity)
     */
    CartesianState to_cartesian() const;
    
    /**
     * @brief Construct from Cartesian state
     * @param state Cartesian state
     * @return Keplerian elements
     */
    static KeplerianElements from_cartesian(const CartesianState& state);
    
    // ========================================================================
    // Jacobian Matrices (for covariance propagation)
    // ========================================================================
    
    /**
     * @brief Compute Jacobian matrix ∂(Cartesian)/∂(Keplerian)
     * 
     * Returns the 6x6 partial derivative matrix of Cartesian state
     * with respect to Keplerian elements at current element values.
     * Used for covariance propagation: P_cart = J * P_kep * J^T
     * 
     * @return 6x6 Jacobian matrix [∂(x,y,z,vx,vy,vz)/∂(a,e,i,Ω,ω,M)]
     */
    Matrix6d jacobian_to_cartesian() const;
    
    /**
     * @brief Compute Jacobian matrix ∂(Keplerian)/∂(Cartesian)
     * 
     * Returns the 6x6 partial derivative matrix of Keplerian elements
     * with respect to Cartesian state at current element values.
     * 
     * @param state Cartesian state (must correspond to current elements)
     * @return 6x6 Jacobian matrix [∂(a,e,i,Ω,ω,M)/∂(x,y,z,vx,vy,vz)]
     */
    static Matrix6d jacobian_from_cartesian(const CartesianState& state);
    
    // ========================================================================
    // Orbit Classification
    // ========================================================================
    
    bool is_elliptic() const { return e_ < 1.0; }
    bool is_parabolic() const { return std::abs(e_ - 1.0) < 1e-6; }
    bool is_hyperbolic() const { return e_ > 1.0; }
    bool is_circular(double tol = 1e-6) const { return e_ < tol; }
    bool is_equatorial(double tol = 1e-6) const { 
        return i_ < tol || std::abs(i_ - constants::PI) < tol; 
    }
    
    // ========================================================================
    // String Representation
    // ========================================================================
    
    std::string to_string() const {
        std::ostringstream oss;
        oss << "KeplerianElements:\n"
            << "  a [km]: " << a_ << "\n"
            << "  e: " << e_ << "\n"
            << "  i [deg]: " << i_ * constants::RAD_TO_DEG << "\n"
            << "  Ω [deg]: " << Omega_ * constants::RAD_TO_DEG << "\n"
            << "  ω [deg]: " << omega_ * constants::RAD_TO_DEG << "\n"
            << "  M [deg]: " << M_ * constants::RAD_TO_DEG << "\n"
            << "  Period [days]: " << period() / constants::DAY << "\n"
            << "  q [km]: " << periapsis_distance();
        return oss.str();
    }

private:
    double a_;      ///< Semi-major axis [km]
    double e_;      ///< Eccentricity
    double i_;      ///< Inclination [rad]
    double Omega_;  ///< RAAN [rad]
    double omega_;  ///< Argument of periapsis [rad]
    double M_;      ///< Mean anomaly [rad]
    double mu_;     ///< Gravitational parameter [km³/s²]
};

} // namespace coordinates
} // namespace astdyn

#endif // ASTDYN_COORDINATES_KEPLERIANELEMENTS_HPP
