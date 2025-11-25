/**
 * @file OrbitalElements.hpp
 * @brief Orbital element representations and conversions
 * 
 * This module defines various representations of orbital elements:
 * - Keplerian elements (a, e, i, Ω, ω, M)
 * - Cartesian state vectors (position, velocity)
 * - Equinoctial elements (singularity-free for near-circular/equatorial orbits)
 * - Cometary elements (q, e, i, Ω, ω, Tp)
 * 
 * Conversions between representations are provided.
 */

#ifndef ASTDYN_ORBITAL_ELEMENTS_HPP
#define ASTDYN_ORBITAL_ELEMENTS_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/Constants.hpp"
#include <Eigen/Dense>
#include <string>
#include <optional>

namespace astdyn::propagation {

/**
 * @brief Keplerian orbital elements
 * 
 * Standard osculating elements at a given epoch.
 * Units: km, radians, days (MJD)
 */
struct KeplerianElements {
    double epoch_mjd_tdb;           ///< Epoch (MJD TDB)
    double semi_major_axis;          ///< Semi-major axis [AU]
    double eccentricity;             ///< Eccentricity [dimensionless]
    double inclination;              ///< Inclination [rad]
    double longitude_ascending_node; ///< Longitude of ascending node (Ω) [rad]
    double argument_perihelion;      ///< Argument of perihelion (ω) [rad]
    double mean_anomaly;             ///< Mean anomaly (M) [rad]
    
    double gravitational_parameter;  ///< GM of central body [AU³/day²]
    
    // Optional uncertainties (covariance matrix diagonal)
    std::optional<Eigen::Matrix<double, 6, 6>> covariance;
    
    KeplerianElements() = default;
    
    /**
     * @brief Compute orbital period
     * @return Period in days
     */
    double period() const;
    
    /**
     * @brief Compute mean motion
     * @return Mean motion in rad/day
     */
    double mean_motion() const;
    
    /**
     * @brief Compute perihelion distance
     * @return q = a(1-e) [AU]
     */
    double perihelion_distance() const;
    
    /**
     * @brief Compute aphelion distance
     * @return Q = a(1+e) [AU]
     */
    double aphelion_distance() const;
    
    /**
     * @brief Check if orbit is hyperbolic
     */
    bool is_hyperbolic() const { return eccentricity >= 1.0; }
    
    /**
     * @brief Check if orbit is parabolic
     */
    bool is_parabolic() const { return std::abs(eccentricity - 1.0) < 1e-8; }
    
    /**
     * @brief Check if orbit is elliptic
     */
    bool is_elliptic() const { return eccentricity < 1.0; }
};

/**
 * @brief Cartesian state vector
 * 
 * Position and velocity in J2000 equatorial frame.
 * Units: AU, AU/day
 */
struct CartesianElements {
    double epoch_mjd_tdb;              ///< Epoch (MJD TDB)
    Eigen::Vector3d position;          ///< Position [AU]
    Eigen::Vector3d velocity;          ///< Velocity [AU/day]
    
    double gravitational_parameter;    ///< GM of central body [AU³/day²]
    
    // Optional covariance (6x6: position, velocity)
    std::optional<Eigen::Matrix<double, 6, 6>> covariance;
    
    CartesianElements() = default;
    
    /**
     * @brief Compute orbital energy
     * @return Specific energy [AU²/day²]
     */
    double energy() const;
    
    /**
     * @brief Compute angular momentum vector
     * @return h = r × v [AU²/day]
     */
    Eigen::Vector3d angular_momentum() const;
    
    /**
     * @brief Compute distance from central body
     * @return |r| [AU]
     */
    double distance() const { return position.norm(); }
    
    /**
     * @brief Compute speed
     * @return |v| [AU/day]
     */
    double speed() const { return velocity.norm(); }
};

/**
 * @brief Equinoctial orbital elements (singularity-free)
 * 
 * Useful for near-circular and near-equatorial orbits.
 * See Walker et al. (1985) Celestial Mechanics.
 */
struct EquinoctialElements {
    double epoch_mjd_tdb;     ///< Epoch (MJD TDB)
    double a;                 ///< Semi-major axis [AU]
    double h;                 ///< h = e sin(ω + Ω)
    double k;                 ///< k = e cos(ω + Ω)
    double p;                 ///< p = tan(i/2) sin(Ω)
    double q;                 ///< q = tan(i/2) cos(Ω)
    double lambda;            ///< Mean longitude λ = M + ω + Ω [rad]
    
    double gravitational_parameter;
    
    std::optional<Eigen::Matrix<double, 6, 6>> covariance;
    
    EquinoctialElements() = default;
};

/**
 * @brief Cometary orbital elements
 * 
 * Defined by perihelion distance instead of semi-major axis.
 * Useful for parabolic/hyperbolic orbits.
 */
struct CometaryElements {
    double epoch_mjd_tdb;              ///< Epoch (MJD TDB)
    double perihelion_distance;        ///< q [AU]
    double eccentricity;               ///< e [dimensionless]
    double inclination;                ///< i [rad]
    double longitude_ascending_node;   ///< Ω [rad]
    double argument_perihelion;        ///< ω [rad]
    double time_perihelion_mjd_tdb;    ///< Tp [MJD TDB]
    
    double gravitational_parameter;    ///< GM [AU³/day²]
    
    std::optional<Eigen::Matrix<double, 6, 6>> covariance;
    
    CometaryElements() = default;
    
    /**
     * @brief Compute semi-major axis (if elliptic)
     * @return a = q/(1-e) [AU]
     */
    double semi_major_axis() const;
    
    /**
     * @brief Check orbit type
     */
    bool is_hyperbolic() const { return eccentricity > 1.0; }
    bool is_parabolic() const { return std::abs(eccentricity - 1.0) < 1e-8; }
    bool is_elliptic() const { return eccentricity < 1.0; }
};

// ============================================================================
// Conversion Functions
// ============================================================================

/**
 * @brief Convert Keplerian to Cartesian elements
 * 
 * Uses standard two-body orbital mechanics.
 * Solves Kepler's equation M = E - e sin(E) for eccentric anomaly.
 * 
 * @param kep Keplerian elements
 * @return Cartesian state vector
 */
CartesianElements keplerian_to_cartesian(const KeplerianElements& kep);

/**
 * @brief Convert Cartesian to Keplerian elements
 * 
 * @param cart Cartesian state vector
 * @return Keplerian elements
 */
KeplerianElements cartesian_to_keplerian(const CartesianElements& cart);

/**
 * @brief Convert Keplerian to Equinoctial
 * 
 * @param kep Keplerian elements
 * @return Equinoctial elements
 */
EquinoctialElements keplerian_to_equinoctial(const KeplerianElements& kep);

/**
 * @brief Convert Equinoctial to Keplerian
 * 
 * @param eq Equinoctial elements
 * @return Keplerian elements
 */
KeplerianElements equinoctial_to_keplerian(const EquinoctialElements& eq);

/**
 * @brief Convert Keplerian to Cometary
 * 
 * Computes time of perihelion passage from mean anomaly and epoch.
 * 
 * @param kep Keplerian elements
 * @return Cometary elements
 */
CometaryElements keplerian_to_cometary(const KeplerianElements& kep);

/**
 * @brief Convert Cometary to Keplerian
 * 
 * Computes mean anomaly at epoch from time of perihelion.
 * 
 * @param com Cometary elements
 * @return Keplerian elements
 */
KeplerianElements cometary_to_keplerian(const CometaryElements& com);

/**
 * @brief Solve Kepler's equation M = E - e sin(E) for eccentric anomaly
 * 
 * Uses Newton-Raphson iteration with excellent convergence.
 * 
 * @param M Mean anomaly [rad]
 * @param e Eccentricity
 * @param tolerance Convergence tolerance (default 1e-12)
 * @param max_iter Maximum iterations (default 50)
 * @return Eccentric anomaly E [rad]
 */
double solve_kepler_equation(double M, double e, 
                             double tolerance = 1e-12, 
                             int max_iter = 50);

/**
 * @brief Compute true anomaly from eccentric anomaly
 * 
 * Uses exact formula: tan(ν/2) = √((1+e)/(1-e)) tan(E/2)
 * 
 * @param E Eccentric anomaly [rad]
 * @param e Eccentricity
 * @return True anomaly ν [rad]
 */
double eccentric_to_true_anomaly(double E, double e);

/**
 * @brief Compute eccentric anomaly from true anomaly
 * 
 * @param nu True anomaly [rad]
 * @param e Eccentricity
 * @return Eccentric anomaly E [rad]
 */
double true_to_eccentric_anomaly(double nu, double e);

// ============================================================================
// Mean to Osculating Element Conversions
// ============================================================================

/**
 * @brief Convert mean Keplerian elements to osculating elements
 * 
 * Applies first-order J2 perturbation corrections to convert from mean
 * elements (averaged over short-period perturbations) to osculating elements
 * (instantaneous orbital elements).
 * 
 * This is primarily useful for planetary orbits where mean elements are
 * provided (e.g., from JPL ephemeris approximations) but osculating elements
 * are needed for accurate MOID calculations or propagation.
 * 
 * **Theory:**
 * Mean elements are averaged over short-period perturbations (one orbital period).
 * Osculating elements are instantaneous and include these short-period variations.
 * 
 * The main perturbation for Solar System bodies is the J2 oblateness term of the Sun,
 * though it's very small. For Earth satellites, J2 corrections are much larger.
 * 
 * **First-order corrections:**
 * - Δa = 0 (semi-major axis secular term only)
 * - Δe ≈ J2 terms (small for planets)
 * - Δi ≈ J2 terms (small for planets)
 * - ΔΩ = -3/2 * n * J2 * (R/a)² * cos(i) * t (secular + short-period)
 * - Δω = 3/4 * n * J2 * (R/a)² * (4 - 5sin²i) * t (secular + short-period)
 * - ΔM = 3/4 * n * J2 * (R/a)² * √(1-e²) * (2 - 3sin²i) * t
 * 
 * where:
 * - n = mean motion
 * - J2 = second zonal harmonic (~1e-7 for Sun, ~1e-3 for Earth)
 * - R = equatorial radius of central body
 * - a = semi-major axis
 * - t = time since epoch
 * 
 * **Note:** For Solar System planets around the Sun, J2 effects are negligible
 * (~1 arcsecond over decades). This function includes the formalism for completeness
 * and can be extended for Earth satellites where J2 is significant.
 * 
 * @param mean_elements Mean Keplerian elements
 * @param j2 Second zonal harmonic coefficient (default: 0 for heliocentric, use 0.00108263 for geocentric)
 * @param central_body_radius Equatorial radius of central body [AU] (default: Sun radius)
 * @return Osculating Keplerian elements
 */
KeplerianElements mean_to_osculating(
    const KeplerianElements& mean_elements,
    double j2 = 0.0,
    double central_body_radius = 0.00465047); // Sun radius in AU

/**
 * @brief Convert osculating Keplerian elements to mean elements
 * 
 * Inverse of mean_to_osculating(). Removes short-period perturbations to
 * produce averaged mean elements.
 * 
 * This is less commonly needed than the forward conversion, but useful when:
 * - Comparing orbits with published mean elements
 * - Long-term orbital evolution studies
 * - Generating mean element ephemerides
 * 
 * The conversion applies inverse J2 corrections. For heliocentric orbits,
 * the corrections are tiny and osculating ≈ mean to high precision.
 * 
 * @param osc_elements Osculating Keplerian elements
 * @param j2 Second zonal harmonic coefficient (default: 0)
 * @param central_body_radius Equatorial radius of central body [AU]
 * @return Mean Keplerian elements
 */
KeplerianElements osculating_to_mean(
    const KeplerianElements& osc_elements,
    double j2 = 0.0,
    double central_body_radius = 0.00465047);

} // namespace astdyn::propagation

#endif // ASTDYN_ORBITAL_ELEMENTS_HPP
