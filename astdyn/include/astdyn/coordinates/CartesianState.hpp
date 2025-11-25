/**
 * @file CartesianState.hpp
 * @brief Cartesian state representation (position and velocity)
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Represents orbital state in Cartesian coordinates (r, v).
 * Provides utilities for computing orbital properties and conversions.
 */

#ifndef ASTDYN_COORDINATES_CARTESIANSTATE_HPP
#define ASTDYN_COORDINATES_CARTESIANSTATE_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>
#include <optional>

namespace astdyn {
namespace coordinates {

/**
 * @brief Cartesian state vector (position and velocity)
 * 
 * Represents the state of an object in 3D Cartesian coordinates.
 * Position in km, velocity in km/s (default units).
 */
class CartesianState {
public:
    // ========================================================================
    // Constructors
    // ========================================================================
    
    /**
     * @brief Default constructor (zero state)
     */
    CartesianState() 
        : position_(Vector3d::Zero()), 
          velocity_(Vector3d::Zero()),
          mu_(constants::GM_SUN) {}
    
    /**
     * @brief Construct from position and velocity
     * @param position Position vector [km]
     * @param velocity Velocity vector [km/s]
     * @param mu Gravitational parameter [km^3/s^2] (default: Sun)
     */
    CartesianState(const Vector3d& position, 
                   const Vector3d& velocity,
                   double mu = constants::GM_SUN)
        : position_(position), 
          velocity_(velocity),
          mu_(mu) {}
    
    /**
     * @brief Construct from 6D state vector [rx, ry, rz, vx, vy, vz]
     * @param state 6D state vector
     * @param mu Gravitational parameter [km^3/s^2]
     */
    explicit CartesianState(const Vector6d& state, 
                           double mu = constants::GM_SUN)
        : position_(state.head<3>()), 
          velocity_(state.tail<3>()),
          mu_(mu) {}
    
    // ========================================================================
    // Accessors
    // ========================================================================
    
    const Vector3d& position() const { return position_; }
    const Vector3d& velocity() const { return velocity_; }
    double mu() const { return mu_; }
    
    void set_position(const Vector3d& pos) { position_ = pos; }
    void set_velocity(const Vector3d& vel) { velocity_ = vel; }
    void set_mu(double mu) { mu_ = mu; }
    
    /**
     * @brief Get full 6D state vector
     * @return [rx, ry, rz, vx, vy, vz]
     */
    Vector6d state_vector() const {
        Vector6d state;
        state << position_, velocity_;
        return state;
    }
    
    // ========================================================================
    // Orbital Properties
    // ========================================================================
    
    /**
     * @brief Position magnitude (distance from origin)
     * @return |r| [km]
     */
    double radius() const {
        return position_.norm();
    }
    
    /**
     * @brief Velocity magnitude (speed)
     * @return |v| [km/s]
     */
    double speed() const {
        return velocity_.norm();
    }
    
    /**
     * @brief Specific angular momentum vector h = r × v
     * @return Angular momentum [km^2/s]
     */
    Vector3d angular_momentum() const {
        return position_.cross(velocity_);
    }
    
    /**
     * @brief Specific angular momentum magnitude
     * @return |h| [km^2/s]
     */
    double angular_momentum_magnitude() const {
        return angular_momentum().norm();
    }
    
    /**
     * @brief Specific orbital energy ε = v²/2 - μ/r
     * @return Energy per unit mass [km^2/s^2]
     */
    double specific_energy() const {
        double v2 = velocity_.squaredNorm();
        double r = radius();
        return 0.5 * v2 - mu_ / r;
    }
    
    /**
     * @brief Semi-major axis a = -μ/(2ε)
     * @return Semi-major axis [km], or infinity for parabolic/hyperbolic
     */
    double semi_major_axis() const {
        double energy = specific_energy();
        if (std::abs(energy) < 1e-10) {
            return std::numeric_limits<double>::infinity(); // Parabolic
        }
        return -mu_ / (2.0 * energy);
    }
    
    /**
     * @brief Eccentricity vector e = (v×h)/μ - r/|r|
     * @return Eccentricity vector (dimensionless)
     */
    Vector3d eccentricity_vector() const {
        Vector3d h = angular_momentum();
        Vector3d r_hat = position_.normalized();
        return velocity_.cross(h) / mu_ - r_hat;
    }
    
    /**
     * @brief Orbital eccentricity magnitude
     * @return e (dimensionless)
     */
    double eccentricity() const {
        return eccentricity_vector().norm();
    }
    
    /**
     * @brief Inclination angle (angle between h and z-axis)
     * @return Inclination [rad]
     */
    double inclination() const {
        Vector3d h = angular_momentum();
        double h_mag = h.norm();
        if (h_mag < 1e-10) {
            return 0.0; // Rectilinear orbit
        }
        return std::acos(h.z() / h_mag);
    }
    
    /**
     * @brief Node vector N = z × h
     * @return Node vector
     */
    Vector3d node_vector() const {
        Vector3d z_axis(0, 0, 1);
        return z_axis.cross(angular_momentum());
    }
    
    /**
     * @brief Check if orbit is elliptic (e < 1, energy < 0)
     */
    bool is_elliptic() const {
        return eccentricity() < 1.0 && specific_energy() < 0.0;
    }
    
    /**
     * @brief Check if orbit is parabolic (e ≈ 1, energy ≈ 0)
     */
    bool is_parabolic() const {
        double e = eccentricity();
        double energy = specific_energy();
        return std::abs(e - 1.0) < 1e-6 && std::abs(energy) < 1e-6;
    }
    
    /**
     * @brief Check if orbit is hyperbolic (e > 1, energy > 0)
     */
    bool is_hyperbolic() const {
        return eccentricity() > 1.0 && specific_energy() > 0.0;
    }
    
    /**
     * @brief Check if orbit is circular (e ≈ 0)
     */
    bool is_circular(double tolerance = 1e-6) const {
        return eccentricity() < tolerance;
    }
    
    /**
     * @brief Check if orbit is equatorial (i ≈ 0 or i ≈ π)
     */
    bool is_equatorial(double tolerance = 1e-6) const {
        double i = inclination();
        return i < tolerance || std::abs(i - constants::PI) < tolerance;
    }
    
    // ========================================================================
    // Unit Conversions
    // ========================================================================
    
    /**
     * @brief Convert position units (in-place)
     * @param factor Conversion factor (e.g., constants::AU for km → AU)
     */
    void convert_position_units(double factor) {
        position_ *= factor;
    }
    
    /**
     * @brief Convert velocity units (in-place)
     * @param factor Conversion factor
     */
    void convert_velocity_units(double factor) {
        velocity_ *= factor;
    }
    
    /**
     * @brief Convert to AU and AU/day units
     * @return New state in AU, AU/day
     */
    CartesianState to_AU_per_day() const {
        CartesianState result = *this;
        result.convert_position_units(1.0 / constants::AU); // km → AU
        result.convert_velocity_units(constants::DAY / constants::AU); // km/s → AU/day
        result.set_mu(mu_ * std::pow(constants::DAY / constants::AU, 2) / constants::AU); // μ adjustment
        return result;
    }
    
    /**
     * @brief Convert from AU and AU/day to km and km/s
     * @return New state in km, km/s
     */
    CartesianState to_km_per_s() const {
        CartesianState result = *this;
        result.convert_position_units(constants::AU); // AU → km
        result.convert_velocity_units(constants::AU / constants::DAY); // AU/day → km/s
        result.set_mu(mu_ * constants::AU / std::pow(constants::DAY / constants::AU, 2));
        return result;
    }
    
    // ========================================================================
    // Operators
    // ========================================================================
    
    /**
     * @brief Equality operator
     */
    bool operator==(const CartesianState& other) const {
        return position_.isApprox(other.position_, 1e-10) &&
               velocity_.isApprox(other.velocity_, 1e-10) &&
               std::abs(mu_ - other.mu_) < 1e-10;
    }
    
    /**
     * @brief Inequality operator
     */
    bool operator!=(const CartesianState& other) const {
        return !(*this == other);
    }
    
    // ========================================================================
    // String Representation
    // ========================================================================
    
    /**
     * @brief String representation for debugging
     */
    std::string to_string() const {
        std::ostringstream oss;
        oss << "CartesianState:\n"
            << "  Position [km]: [" << position_.transpose() << "]\n"
            << "  Velocity [km/s]: [" << velocity_.transpose() << "]\n"
            << "  μ [km³/s²]: " << mu_ << "\n"
            << "  |r| [km]: " << radius() << "\n"
            << "  |v| [km/s]: " << speed() << "\n"
            << "  a [km]: " << semi_major_axis() << "\n"
            << "  e: " << eccentricity() << "\n"
            << "  i [deg]: " << inclination() * constants::RAD_TO_DEG;
        return oss.str();
    }

private:
    Vector3d position_;  ///< Position vector [km]
    Vector3d velocity_;  ///< Velocity vector [km/s]
    double mu_;          ///< Gravitational parameter [km³/s²]
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Compute relative state between two objects
 * @param state1 First state
 * @param state2 Second state
 * @return Relative state (state1 - state2)
 */
inline CartesianState relative_state(const CartesianState& state1,
                                     const CartesianState& state2) {
    return CartesianState(
        state1.position() - state2.position(),
        state1.velocity() - state2.velocity(),
        state1.mu()
    );
}

/**
 * @brief Compute distance between two states
 * @param state1 First state
 * @param state2 Second state
 * @return Distance [km]
 */
inline double distance(const CartesianState& state1,
                      const CartesianState& state2) {
    return (state1.position() - state2.position()).norm();
}

/**
 * @brief Compute relative velocity between two states
 * @param state1 First state
 * @param state2 Second state
 * @return Relative velocity magnitude [km/s]
 */
inline double relative_velocity(const CartesianState& state1,
                                const CartesianState& state2) {
    return (state1.velocity() - state2.velocity()).norm();
}

} // namespace coordinates
} // namespace astdyn

#endif // ASTDYN_COORDINATES_CARTESIANSTATE_HPP
