/**
 * @file EquinoctialElements.hpp
 * @brief Equinoctial orbital elements (non-singular representation)
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Equinoctial elements avoid singularities for circular and equatorial orbits.
 * Particularly useful for low eccentricity/inclination orbits and perturbation theory.
 * 
 * Elements:
 * - a: Semi-major axis [km]
 * - h = e·sin(ω + Ω): eccentricity x-component
 * - k = e·cos(ω + Ω): eccentricity y-component
 * - p = tan(i/2)·sin(Ω): inclination x-component
 * - q = tan(i/2)·cos(Ω): inclination y-component
 * - λ = M + ω + Ω: mean longitude [rad]
 * 
 * Reference: Battin, "An Introduction to the Mathematics and Methods of Astrodynamics"
 */

#ifndef ASTDYN_COORDINATES_EQUINOCTIALELEMENTS_HPP
#define ASTDYN_COORDINATES_EQUINOCTIALELEMENTS_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/coordinates/KeplerianElements.hpp"
#include <cmath>

namespace astdyn {
namespace coordinates {

/**
 * @brief Equinoctial orbital elements (non-singular)
 * 
 * Advantages:
 * - No singularities for e=0 (circular orbits)
 * - No singularities for i=0 (equatorial orbits)
 * - Better numerical stability for near-circular/near-equatorial orbits
 * - Smooth behavior for orbit propagation with perturbations
 */
class EquinoctialElements {
public:
    // ========================================================================
    // Constructors
    // ========================================================================
    
    /**
     * @brief Default constructor (circular equatorial orbit at 1 AU)
     */
    EquinoctialElements()
        : a_(constants::AU),
          h_(0.0),
          k_(0.0),
          p_(0.0),
          q_(0.0),
          lambda_(0.0),
          mu_(constants::GM_SUN) {}
    
    /**
     * @brief Construct from equinoctial elements
     * @param a Semi-major axis [km]
     * @param h e·sin(ω + Ω)
     * @param k e·cos(ω + Ω)
     * @param p tan(i/2)·sin(Ω)
     * @param q tan(i/2)·cos(Ω)
     * @param lambda Mean longitude [rad]
     * @param mu Gravitational parameter [km³/s²]
     */
    EquinoctialElements(double a, double h, double k,
                       double p, double q, double lambda,
                       double mu = constants::GM_SUN)
        : a_(a), h_(h), k_(k), p_(p), q_(q), 
          lambda_(lambda), mu_(mu) {}
    
    // ========================================================================
    // Accessors
    // ========================================================================
    
    double semi_major_axis() const { return a_; }
    double h() const { return h_; }
    double k() const { return k_; }
    double p() const { return p_; }
    double q() const { return q_; }
    double mean_longitude() const { return lambda_; }
    double mu() const { return mu_; }
    
    void set_semi_major_axis(double a) { a_ = a; }
    void set_h(double h) { h_ = h; }
    void set_k(double k) { k_ = k; }
    void set_p(double p) { p_ = p; }
    void set_q(double q) { q_ = q; }
    void set_mean_longitude(double lambda) { lambda_ = lambda; }
    void set_mu(double mu) { mu_ = mu; }
    
    // ========================================================================
    // Derived Classical Elements
    // ========================================================================
    
    /**
     * @brief Eccentricity e = √(h² + k²)
     */
    double eccentricity() const {
        return std::sqrt(h_ * h_ + k_ * k_);
    }
    
    /**
     * @brief Inclination i = 2·atan(√(p² + q²))
     * @return Inclination [rad]
     */
    double inclination() const {
        double s = std::sqrt(p_ * p_ + q_ * q_);
        return 2.0 * std::atan(s);
    }
    
    /**
     * @brief RAAN Ω = atan2(p, q)
     * @return Right ascension of ascending node [rad]
     */
    double RAAN() const {
        return std::atan2(p_, q_);
    }
    
    /**
     * @brief Argument of periapsis + RAAN: ω + Ω = atan2(h, k)
     * @return ω + Ω [rad]
     */
    double longitude_of_periapsis() const {
        return std::atan2(h_, k_);
    }
    
    /**
     * @brief Argument of periapsis ω = (ω + Ω) - Ω
     * @return Argument of periapsis [rad]
     */
    double argument_of_periapsis() const {
        return longitude_of_periapsis() - RAAN();
    }
    
    /**
     * @brief Mean anomaly M = λ - (ω + Ω)
     * @return Mean anomaly [rad]
     */
    double mean_anomaly() const {
        return lambda_ - longitude_of_periapsis();
    }
    
    /**
     * @brief Period T = 2π√(a³/μ)
     */
    double period() const {
        if (a_ <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        return 2.0 * constants::PI * std::sqrt(a_ * a_ * a_ / mu_);
    }
    
    /**
     * @brief Periapsis distance q = a(1 - e)
     */
    double periapsis_distance() const {
        return a_ * (1.0 - eccentricity());
    }
    
    // ========================================================================
    // Conversions
    // ========================================================================
    
    /**
     * @brief Convert to Keplerian elements
     * @return Classical Keplerian elements
     */
    KeplerianElements to_keplerian() const {
        double e = eccentricity();
        double i = inclination();
        double Omega = RAAN();
        double omega = argument_of_periapsis();
        double M = mean_anomaly();
        
        // Normalize angles
        auto normalize = [](double angle) {
            angle = std::fmod(angle, 2.0 * constants::PI);
            if (angle < 0.0) angle += 2.0 * constants::PI;
            return angle;
        };
        
        return KeplerianElements(a_, e, i, 
                                normalize(Omega), 
                                normalize(omega), 
                                normalize(M), mu_);
    }
    
    /**
     * @brief Construct from Keplerian elements
     * @param kep Classical Keplerian elements
     * @return Equinoctial elements
     */
    static EquinoctialElements from_keplerian(const KeplerianElements& kep) {
        double a = kep.semi_major_axis();
        double e = kep.eccentricity();
        double i = kep.inclination();
        double Omega = kep.RAAN();
        double omega = kep.argument_of_periapsis();
        double M = kep.mean_anomaly();
        
        // Compute equinoctial elements
        double lon_peri = omega + Omega; // ω + Ω
        double h = e * std::sin(lon_peri);
        double k = e * std::cos(lon_peri);
        double tan_half_i = std::tan(i / 2.0);
        double p = tan_half_i * std::sin(Omega);
        double q = tan_half_i * std::cos(Omega);
        double lambda = M + lon_peri;
        
        return EquinoctialElements(a, h, k, p, q, lambda, kep.mu());
    }
    
    /**
     * @brief Convert to Cartesian state
     * @return Cartesian position and velocity
     */
    CartesianState to_cartesian() const {
        // Convert via Keplerian (most straightforward)
        return to_keplerian().to_cartesian();
    }
    
    /**
     * @brief Construct from Cartesian state
     * @param state Cartesian state
     * @return Equinoctial elements
     */
    static EquinoctialElements from_cartesian(const CartesianState& state) {
        KeplerianElements kep = KeplerianElements::from_cartesian(state);
        return from_keplerian(kep);
    }
    
    // ========================================================================
    // Orbit Classification
    // ========================================================================
    
    bool is_circular(double tol = 1e-6) const {
        return eccentricity() < tol;
    }
    
    bool is_equatorial(double tol = 1e-6) const {
        return inclination() < tol;
    }
    
    bool is_elliptic() const {
        return eccentricity() < 1.0;
    }
    
    // ========================================================================
    // String Representation
    // ========================================================================
    
    std::string to_string() const {
        std::ostringstream oss;
        oss << "EquinoctialElements:\n"
            << "  a [km]: " << a_ << "\n"
            << "  h: " << h_ << "\n"
            << "  k: " << k_ << "\n"
            << "  p: " << p_ << "\n"
            << "  q: " << q_ << "\n"
            << "  λ [deg]: " << lambda_ * constants::RAD_TO_DEG << "\n"
            << "  Derived:\n"
            << "    e: " << eccentricity() << "\n"
            << "    i [deg]: " << inclination() * constants::RAD_TO_DEG << "\n"
            << "    Ω [deg]: " << RAAN() * constants::RAD_TO_DEG;
        return oss.str();
    }

private:
    double a_;       ///< Semi-major axis [km]
    double h_;       ///< e·sin(ω + Ω)
    double k_;       ///< e·cos(ω + Ω)
    double p_;       ///< tan(i/2)·sin(Ω)
    double q_;       ///< tan(i/2)·cos(Ω)
    double lambda_;  ///< Mean longitude [rad]
    double mu_;      ///< Gravitational parameter [km³/s²]
};

} // namespace coordinates
} // namespace astdyn

#endif // ASTDYN_COORDINATES_EQUINOCTIALELEMENTS_HPP
