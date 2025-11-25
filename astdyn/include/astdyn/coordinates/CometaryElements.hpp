/**
 * @file CometaryElements.hpp
 * @brief Cometary orbital elements (perihelion-based representation)
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Cometary elements are particularly suited for:
 * - Parabolic orbits (e ≈ 1)
 * - Hyperbolic orbits (e > 1)
 * - Long-period comets
 * - Interstellar objects
 * 
 * Elements:
 * - q: Perihelion distance [km]
 * - e: Eccentricity
 * - i: Inclination [rad]
 * - Ω: Right ascension of ascending node [rad]
 * - ω: Argument of periapsis [rad]
 * - T: Time of perihelion passage [JD or MJD]
 * 
 * Reference: Murray & Dermott, "Solar System Dynamics"
 */

#ifndef ASTDYN_COORDINATES_COMETARYELEMENTS_HPP
#define ASTDYN_COORDINATES_COMETARYELEMENTS_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/coordinates/KeplerianElements.hpp"
#include <cmath>
#include <limits>

namespace astdyn {
namespace coordinates {

/**
 * @brief Cometary orbital elements (perihelion-based)
 * 
 * Advantages:
 * - Natural for parabolic/hyperbolic orbits
 * - Well-defined for all eccentricities
 * - Direct observational connection (perihelion passage time)
 * - Used in Minor Planet Center orbit catalogs
 */
class CometaryElements {
public:
    // ========================================================================
    // Constructors
    // ========================================================================
    
    /**
     * @brief Default constructor (parabolic comet at 1 AU perihelion)
     */
    CometaryElements()
        : q_(constants::AU),
          e_(1.0),
          i_(0.0),
          Omega_(0.0),
          omega_(0.0),
          T_(0.0),
          mu_(constants::GM_SUN) {}
    
    /**
     * @brief Construct from cometary elements
     * @param q Perihelion distance [km]
     * @param e Eccentricity
     * @param i Inclination [rad]
     * @param Omega Right ascension of ascending node [rad]
     * @param omega Argument of periapsis [rad]
     * @param T Time of perihelion passage [JD or MJD]
     * @param mu Gravitational parameter [km³/s²]
     */
    CometaryElements(double q, double e, double i,
                    double Omega, double omega, double T,
                    double mu = constants::GM_SUN)
        : q_(q), e_(e), i_(i), Omega_(Omega), omega_(omega), 
          T_(T), mu_(mu) {}
    
    // ========================================================================
    // Accessors
    // ========================================================================
    
    double perihelion_distance() const { return q_; }
    double eccentricity() const { return e_; }
    double inclination() const { return i_; }
    double RAAN() const { return Omega_; }
    double argument_of_periapsis() const { return omega_; }
    double time_of_perihelion() const { return T_; }
    double mu() const { return mu_; }
    
    void set_perihelion_distance(double q) { q_ = q; }
    void set_eccentricity(double e) { e_ = e; }
    void set_inclination(double i) { i_ = i; }
    void set_RAAN(double Omega) { Omega_ = Omega; }
    void set_argument_of_periapsis(double omega) { omega_ = omega; }
    void set_time_of_perihelion(double T) { T_ = T; }
    void set_mu(double mu) { mu_ = mu; }
    
    // ========================================================================
    // Derived Quantities
    // ========================================================================
    
    /**
     * @brief Semi-major axis a = q/(1-e)
     * @note Returns infinity for parabolic orbits (e=1)
     *       Returns negative for hyperbolic orbits (e>1)
     */
    double semi_major_axis() const {
        if (std::abs(e_ - 1.0) < 1e-10) {
            return std::numeric_limits<double>::infinity();
        }
        return q_ / (1.0 - e_);
    }
    
    /**
     * @brief Aphelion distance Q = q(1+e)/(1-e) = a(1+e)
     * @note Only valid for elliptic orbits (e < 1)
     */
    double aphelion_distance() const {
        if (e_ >= 1.0) {
            return std::numeric_limits<double>::infinity();
        }
        return q_ * (1.0 + e_) / (1.0 - e_);
    }
    
    /**
     * @brief Period T = 2π√(a³/μ)
     * @note Only valid for elliptic orbits
     */
    double period() const {
        if (e_ >= 1.0) {
            return std::numeric_limits<double>::infinity();
        }
        double a = semi_major_axis();
        return 2.0 * constants::PI * std::sqrt(a * a * a / mu_);
    }
    
    /**
     * @brief Semi-latus rectum p = q(1+e)
     */
    double semi_latus_rectum() const {
        return q_ * (1.0 + e_);
    }
    
    /**
     * @brief Mean motion n = √(μ/a³)
     * @note Only valid for elliptic orbits
     */
    double mean_motion() const {
        if (e_ >= 1.0) {
            return 0.0;
        }
        double a = semi_major_axis();
        return std::sqrt(mu_ / (a * a * a));
    }
    
    // ========================================================================
    // Anomaly Computations
    // ========================================================================
    
    /**
     * @brief Compute mean anomaly at given time
     * @param t Time [same units as T]
     * @return Mean anomaly [rad] (for elliptic) or hyperbolic mean anomaly (for hyperbolic)
     */
    double mean_anomaly_at_time(double t) const {
        double dt = t - T_;
        
        if (e_ < 1.0) {
            // Elliptic: M = n·Δt
            return mean_motion() * dt;
        } else if (e_ > 1.0) {
            // Hyperbolic: M_h = √(μ/(-a)³)·Δt
            double a = semi_major_axis(); // negative for hyperbolic
            return std::sqrt(mu_ / (-a * a * a)) * dt;
        } else {
            // Parabolic: use Barker's equation (not implemented here)
            // For now, return a placeholder
            return 0.0;
        }
    }
    
    // ========================================================================
    // Conversions
    // ========================================================================
    
    /**
     * @brief Convert to Keplerian elements at given time
     * @param t Time [same units as T]
     * @return Keplerian elements
     */
    KeplerianElements to_keplerian(double t) const {
        double M = mean_anomaly_at_time(t);
        double a = semi_major_axis();
        
        return KeplerianElements(a, e_, i_, Omega_, omega_, M, mu_);
    }
    
    /**
     * @brief Construct from Keplerian elements
     * @param kep Keplerian elements
     * @param T Time of perihelion passage
     * @return Cometary elements
     */
    static CometaryElements from_keplerian(const KeplerianElements& kep, double T) {
        double a = kep.semi_major_axis();
        double e = kep.eccentricity();
        double q = a * (1.0 - e);
        
        return CometaryElements(q, e, kep.inclination(),
                               kep.RAAN(), kep.argument_of_periapsis(),
                               T, kep.mu());
    }
    
    /**
     * @brief Convert to Cartesian state at given time
     * @param t Time [same units as T]
     * @return Cartesian state
     */
    CartesianState to_cartesian(double t) const {
        return to_keplerian(t).to_cartesian();
    }
    
    /**
     * @brief Construct from Cartesian state (requires time of perihelion)
     * @param state Cartesian state
     * @param T Time of perihelion passage
     * @return Cometary elements
     */
    static CometaryElements from_cartesian(const CartesianState& state, double T) {
        KeplerianElements kep = KeplerianElements::from_cartesian(state);
        return from_keplerian(kep, T);
    }
    
    // ========================================================================
    // Orbit Classification
    // ========================================================================
    
    bool is_elliptic() const { return e_ < 1.0; }
    bool is_parabolic(double tol = 1e-6) const { 
        return std::abs(e_ - 1.0) < tol; 
    }
    bool is_hyperbolic() const { return e_ > 1.0; }
    bool is_circular(double tol = 1e-6) const { return e_ < tol; }
    
    /**
     * @brief Check if orbit is bound (elliptic)
     */
    bool is_bound() const { return is_elliptic(); }
    
    /**
     * @brief Check if orbit is unbound (parabolic or hyperbolic)
     */
    bool is_unbound() const { return e_ >= 1.0; }
    
    // ========================================================================
    // String Representation
    // ========================================================================
    
    std::string to_string() const {
        std::ostringstream oss;
        oss << "CometaryElements:\n"
            << "  q [AU]: " << q_ / constants::AU << "\n"
            << "  e: " << e_ << "\n"
            << "  i [deg]: " << i_ * constants::RAD_TO_DEG << "\n"
            << "  Ω [deg]: " << Omega_ * constants::RAD_TO_DEG << "\n"
            << "  ω [deg]: " << omega_ * constants::RAD_TO_DEG << "\n"
            << "  T [JD]: " << T_ << "\n"
            << "  Type: ";
        
        if (is_elliptic()) {
            oss << "Elliptic\n"
                << "  a [AU]: " << semi_major_axis() / constants::AU << "\n"
                << "  T [years]: " << period() / constants::YEAR;
        } else if (is_parabolic()) {
            oss << "Parabolic";
        } else {
            oss << "Hyperbolic\n"
                << "  a [AU]: " << semi_major_axis() / constants::AU;
        }
        
        return oss.str();
    }

private:
    double q_;      ///< Perihelion distance [km]
    double e_;      ///< Eccentricity
    double i_;      ///< Inclination [rad]
    double Omega_;  ///< Right ascension of ascending node [rad]
    double omega_;  ///< Argument of periapsis [rad]
    double T_;      ///< Time of perihelion passage [JD or MJD]
    double mu_;     ///< Gravitational parameter [km³/s²]
};

} // namespace coordinates
} // namespace astdyn

#endif // ASTDYN_COORDINATES_COMETARYELEMENTS_HPP
