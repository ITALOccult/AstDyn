/**
 * @file Types.hpp
 * @brief Common type definitions for OrbFit
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#ifndef ASTDYN_CORE_TYPES_HPP
#define ASTDYN_CORE_TYPES_HPP

#include <Eigen/Dense>
#include <string>
#include <cstdint>
#include <limits>
#include <cmath>

namespace astdyn {

// ============================================================================
// Eigen Type Aliases
// ============================================================================

/// 3D vector (position or velocity)
using Vector3d = Eigen::Vector3d;

/// 6D state vector (position + velocity)
using Vector6d = Eigen::Matrix<double, 6, 1>;

/// Generic dynamic vector
using VectorXd = Eigen::VectorXd;

/// 3x3 matrix
using Matrix3d = Eigen::Matrix3d;

/// 6x6 matrix (state transition matrix, covariance)
using Matrix6d = Eigen::Matrix<double, 6, 6>;

/// Generic dynamic matrix
using MatrixXd = Eigen::MatrixXd;

// ============================================================================
// Time Types
// ============================================================================

/// Modified Julian Date (days)
using MJD = double;

/// Julian Date (days)
using JD = double;

/// Time interval (days)
using TimeSpan = double;

// ============================================================================
// Scalar Types
// ============================================================================

/// Angle in radians
using Radians = double;

/// Angle in degrees
using Degrees = double;

/// Distance in AU
using AU_Distance = double;

/// Distance in km
using KM_Distance = double;

/// Velocity in AU/day
using AU_Per_Day = double;

/// Velocity in km/s
using KM_Per_Second = double;

// ============================================================================
// Enumeration Types
// ============================================================================

/**
 * @brief Coordinate system types
 */
enum class CoordinateSystem {
    ECLIPTIC_J2000,     ///< Ecliptic coordinates, J2000.0
    EQUATORIAL_J2000,   ///< Equatorial coordinates, J2000.0
    BARYCENTRIC,        ///< Solar system barycentric
    HELIOCENTRIC        ///< Heliocentric
};

/**
 * @brief Orbital element representation types
 */
enum class ElementType {
    KEPLERIAN,          ///< Classical Keplerian elements (a,e,i,Ω,ω,M)
    CARTESIAN,          ///< Cartesian state vector (x,y,z,vx,vy,vz)
    COMETARY,           ///< Cometary elements (q,e,i,Ω,ω,T)
    EQUINOCTIAL,        ///< Equinoctial elements (for small e, i)
    DELAUNAY            ///< Delaunay canonical elements
};

/**
 * @brief Time scale types
 */
enum class TimeScale {
    UTC,                ///< Coordinated Universal Time
    UT1,                ///< Universal Time 1
    TAI,                ///< International Atomic Time
    TT,                 ///< Terrestrial Time
    TDB,                ///< Barycentric Dynamical Time
    GPS                 ///< GPS Time
};

/**
 * @brief Observation types
 */
enum class ObservationType {
    OPTICAL_RA_DEC,     ///< Optical astrometry (RA, Dec)
    RADAR_RANGE,        ///< Radar ranging
    RADAR_DOPPLER,      ///< Radar Doppler
    SATELLITE,          ///< Satellite-based observation
    OCCULTATION,        ///< Stellar occultation
    ROVING_OBSERVER     ///< Roving observer
};

/**
 * @brief Integration method types
 */
enum class IntegratorType {
    RADAU15,            ///< Radau 15th order implicit
    RUNGE_KUTTA_GAUSS,  ///< Runge-Kutta-Gauss
    ADAMS_BASHFORTH,    ///< Adams-Bashforth-Moulton
    BULIRSCH_STOER      ///< Bulirsch-Stoer
};

/**
 * @brief Force model components
 */
enum class ForceComponent {
    POINT_MASS,         ///< N-body point mass gravity
    RELATIVISTIC,       ///< General relativistic corrections
    YARKOVSKY,          ///< Yarkovsky thermal effect
    YORP,               ///< YORP (spin effect)
    SOLAR_PRESSURE,     ///< Solar radiation pressure
    EARTH_HARMONICS,    ///< Earth spherical harmonics
    ASTEROID_PERT       ///< Asteroid perturbations
};

// ============================================================================
// Special Values
// ============================================================================

/// Not-a-Number for indicating missing/invalid data
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

/// Positive infinity
constexpr double Infinity = std::numeric_limits<double>::infinity();

/// Check if value is NaN
inline bool isNaN(double value) {
    return std::isnan(value);
}

/// Check if value is finite
inline bool isFinite(double value) {
    return std::isfinite(value);
}

// ============================================================================
// Result Types (for error handling)
// ============================================================================

/**
 * @brief Generic result type with success/failure status
 * @tparam T Result value type
 */
template<typename T>
struct Result {
    T value;                    ///< Result value (if successful)
    bool success;               ///< Success flag
    std::string error_message;  ///< Error message (if failed)
    
    /// Constructor for success
    static Result<T> Success(const T& val) {
        return Result<T>{val, true, ""};
    }
    
    /// Constructor for failure
    static Result<T> Failure(const std::string& msg) {
        return Result<T>{T{}, false, msg};
    }
    
    /// Check if result is successful
    explicit operator bool() const {
        return success;
    }
};

// ============================================================================
// ID Types
// ============================================================================

/// Observatory code (MPC format)
using ObservatoryCode = std::string;

/// Object designation/name
using ObjectName = std::string;

/// Planet ID (1-9 for planets, 10 for Moon, 0 for Sun)
using PlanetID = int;

} // namespace astdyn

#endif // ASTDYN_CORE_TYPES_HPP
