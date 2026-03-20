#ifndef ASTDYN_CORE_PHYSICS_TYPES_HPP
#define ASTDYN_CORE_PHYSICS_TYPES_HPP

#include "astdyn/core/Constants.hpp"
#include <cmath>

namespace astdyn::physics {

// ============================================================================
// Distance
// ============================================================================
class Distance {
public:
    // Constructors
    constexpr Distance() noexcept : value_m_(0.0) {}
    
    // Factory methods
    [[nodiscard]] static constexpr Distance from_si(double m)  noexcept { return Distance(m); }
    [[nodiscard]] static constexpr Distance from_m(double m)   noexcept { return Distance(m); }
    [[nodiscard]] static constexpr Distance from_km(double km) noexcept { return Distance(km * 1000.0); }
    [[nodiscard]] static constexpr Distance from_au(double au) noexcept { return Distance(au * constants::AU * 1000.0); }
    [[nodiscard]] static constexpr Distance zero()             noexcept { return Distance(0.0); }

    // Extraction methods (explicit casts)
    [[nodiscard]] constexpr double to_m()  const noexcept { return value_m_; }
    [[nodiscard]] constexpr double to_km() const noexcept { return value_m_ / 1000.0; }
    [[nodiscard]] constexpr double to_au() const noexcept { return value_m_ / (constants::AU * 1000.0); }

    // Arithmetic operations
    [[nodiscard]] constexpr Distance operator+(const Distance& other) const noexcept { return Distance(value_m_ + other.value_m_); }
    [[nodiscard]] constexpr Distance operator-(const Distance& other) const noexcept { return Distance(value_m_ - other.value_m_); }
    constexpr Distance& operator+=(const Distance& other) noexcept { value_m_ += other.value_m_; return *this; }
    constexpr Distance& operator-=(const Distance& other) noexcept { value_m_ -= other.value_m_; return *this; }
    
    // Scalar multiplication/division
    [[nodiscard]] constexpr Distance operator*(double scalar) const noexcept { return Distance(value_m_ * scalar); }
    [[nodiscard]] constexpr Distance operator/(double scalar) const noexcept { return Distance(value_m_ / scalar); }
    
    // Comparisons
    [[nodiscard]] constexpr bool operator==(const Distance& other) const noexcept { return value_m_ == other.value_m_; }
    [[nodiscard]] constexpr bool operator!=(const Distance& other) const noexcept { return value_m_ != other.value_m_; }
    [[nodiscard]] constexpr bool operator<(const Distance& other) const noexcept { return value_m_ < other.value_m_; }
    [[nodiscard]] constexpr bool operator>(const Distance& other) const noexcept { return value_m_ > other.value_m_; }
    [[nodiscard]] constexpr bool operator<=(const Distance& other) const noexcept { return value_m_ <= other.value_m_; }
    [[nodiscard]] constexpr bool operator>=(const Distance& other) const noexcept { return value_m_ >= other.value_m_; }

private:
    explicit constexpr Distance(double m) noexcept : value_m_(m) {}
    double value_m_; // Internal representation is always meters (SI)
};

/**
 * @brief Distance in Astronomical Units (AU)
 * Specialized for numerically stable orbital integration.
 */
class DistanceAU {
public:
    constexpr DistanceAU() noexcept : value_au_(0.0) {}
    [[nodiscard]] static constexpr DistanceAU from_au(double au) noexcept { return DistanceAU(au); }
    [[nodiscard]] static constexpr DistanceAU from_si(Distance dist) noexcept { return DistanceAU(dist.to_au()); }
    
    [[nodiscard]] constexpr double to_au() const noexcept { return value_au_; }
    [[nodiscard]] constexpr Distance to_si() const noexcept { return Distance::from_au(value_au_); }
    
    [[nodiscard]] constexpr DistanceAU operator+(const DistanceAU& other) const noexcept { return DistanceAU(value_au_ + other.value_au_); }
    [[nodiscard]] constexpr DistanceAU operator-(const DistanceAU& other) const noexcept { return DistanceAU(value_au_ - other.value_au_); }
    [[nodiscard]] constexpr DistanceAU operator*(double scalar) const noexcept { return DistanceAU(value_au_ * scalar); }
    [[nodiscard]] constexpr DistanceAU operator/(double scalar) const noexcept { return DistanceAU(value_au_ / scalar); }
    
    constexpr DistanceAU& operator+=(const DistanceAU& other) noexcept { value_au_ += other.value_au_; return *this; }
    
    [[nodiscard]] constexpr bool operator<(const DistanceAU& other) const noexcept { return value_au_ < other.value_au_; }
    [[nodiscard]] constexpr bool operator>(const DistanceAU& other) const noexcept { return value_au_ > other.value_au_; }

private:
    explicit constexpr DistanceAU(double au) noexcept : value_au_(au) {}
    double value_au_;
};

// ============================================================================
// Velocity
// ============================================================================
class Velocity {
public:
    // Constructors
    constexpr Velocity() noexcept : value_ms_(0.0) {}
    
    // Factory methods
    [[nodiscard]] static constexpr Velocity from_si(double ms)     noexcept { return Velocity(ms); }
    [[nodiscard]] static constexpr Velocity from_ms(double ms)     noexcept { return Velocity(ms); }
    [[nodiscard]] static constexpr Velocity from_km_s(double kms)  noexcept { return Velocity(kms * 1000.0); }
    [[nodiscard]] static constexpr Velocity from_au_d(double aud)  noexcept { return Velocity(aud * (constants::AU * 1000.0) / 86400.0); }
    [[nodiscard]] static constexpr Velocity zero()                 noexcept { return Velocity(0.0); }

    // Extraction methods
    [[nodiscard]] constexpr double to_ms()  const noexcept { return value_ms_; }
    [[nodiscard]] constexpr double to_km_s() const noexcept { return value_ms_ / 1000.0; }
    [[nodiscard]] constexpr double to_au_d() const noexcept { return value_ms_ * 86400.0 / (constants::AU * 1000.0); }

    // Arithmetic
    [[nodiscard]] constexpr Velocity operator+(const Velocity& other) const noexcept { return Velocity(value_ms_ + other.value_ms_); }
    [[nodiscard]] constexpr Velocity operator-(const Velocity& other) const noexcept { return Velocity(value_ms_ - other.value_ms_); }
    [[nodiscard]] constexpr Velocity operator*(double scalar) const noexcept { return Velocity(value_ms_ * scalar); }
    [[nodiscard]] constexpr Velocity operator/(double scalar) const noexcept { return Velocity(value_ms_ / scalar); }
    constexpr Velocity& operator+=(const Velocity& other) noexcept { value_ms_ += other.value_ms_; return *this; }
    constexpr Velocity& operator-=(const Velocity& other) noexcept { value_ms_ -= other.value_ms_; return *this; }

private:
    explicit constexpr Velocity(double ms) noexcept : value_ms_(ms) {}
    double value_ms_; // Internal representation is always m/s (SI)
};

/**
 * @brief Velocity in AU per day (AUD)
 */
class VelocityAUD {
public:
    constexpr VelocityAUD() noexcept : value_aud_(0.0) {}
    [[nodiscard]] static constexpr VelocityAUD from_au_d(double aud) noexcept { return VelocityAUD(aud); }
    [[nodiscard]] static constexpr VelocityAUD from_si(Velocity v) noexcept { return VelocityAUD(v.to_au_d()); }
    
    [[nodiscard]] constexpr double to_au_d() const noexcept { return value_aud_; }
    [[nodiscard]] constexpr Velocity to_si() const noexcept { return Velocity::from_au_d(value_aud_); }
    
    [[nodiscard]] constexpr VelocityAUD operator+(const VelocityAUD& other) const noexcept { return VelocityAUD(value_aud_ + other.value_aud_); }
    [[nodiscard]] constexpr VelocityAUD operator-(const VelocityAUD& other) const noexcept { return VelocityAUD(value_aud_ - other.value_aud_); }
    [[nodiscard]] constexpr VelocityAUD operator*(double scalar) const noexcept { return VelocityAUD(value_aud_ * scalar); }

private:
    explicit constexpr VelocityAUD(double aud) noexcept : value_aud_(aud) {}
    double value_aud_;
};

// ============================================================================
// Acceleration
// ============================================================================
class Acceleration {
public:
    // Constructors
    constexpr Acceleration() noexcept : value_ms2_(0.0) {}
    
    // Factory methods
    [[nodiscard]] static constexpr Acceleration from_si(double ms2)    noexcept { return Acceleration(ms2); }
    [[nodiscard]] static constexpr Acceleration from_ms2(double ms2)   noexcept { return Acceleration(ms2); }
    [[nodiscard]] static constexpr Acceleration from_au_d2(double aud2) noexcept { 
        return Acceleration(aud2 * (constants::AU * 1000.0) / (86400.0 * 86400.0)); 
    }
    [[nodiscard]] static constexpr Acceleration zero()                 noexcept { return Acceleration(0.0); }

    // Extraction
    [[nodiscard]] constexpr double to_ms2()  const noexcept { return value_ms2_; }
    [[nodiscard]] constexpr double to_au_d2() const noexcept { 
        return value_ms2_ / ((constants::AU * 1000.0) / (86400.0 * 86400.0)); 
    }

    // Arithmetic
    [[nodiscard]] constexpr Acceleration operator+(const Acceleration& other) const noexcept { return Acceleration(value_ms2_ + other.value_ms2_); }
    [[nodiscard]] constexpr Acceleration operator-(const Acceleration& other) const noexcept { return Acceleration(value_ms2_ - other.value_ms2_); }
    [[nodiscard]] constexpr Acceleration operator*(double scalar) const noexcept { return Acceleration(value_ms2_ * scalar); }
    constexpr Acceleration& operator+=(const Acceleration& other) noexcept { value_ms2_ += other.value_ms2_; return *this; }

private:
    explicit constexpr Acceleration(double ms2) noexcept : value_ms2_(ms2) {}
    double value_ms2_; // Internal representation is always m/s² (SI)
};

/**
 * @brief Acceleration in AU per day^2
 */
class AccelerationAUD2 {
public:
    constexpr AccelerationAUD2() noexcept : value_aud2_(0.0) {}
    [[nodiscard]] static constexpr AccelerationAUD2 from_au_d2(double aud2) noexcept { return AccelerationAUD2(aud2); }
    
    [[nodiscard]] constexpr double to_au_d2() const noexcept { return value_aud2_; }
    [[nodiscard]] constexpr Acceleration to_si() const noexcept { return Acceleration::from_au_d2(value_aud2_); }
    
    [[nodiscard]] constexpr AccelerationAUD2 operator+(const AccelerationAUD2& other) const noexcept { return AccelerationAUD2(value_aud2_ + other.value_aud2_); }
    [[nodiscard]] constexpr AccelerationAUD2 operator-(const AccelerationAUD2& other) const noexcept { return AccelerationAUD2(value_aud2_ - other.value_aud2_); }
    [[nodiscard]] constexpr AccelerationAUD2 operator*(double scalar) const noexcept { return AccelerationAUD2(value_aud2_ * scalar); }

private:
    explicit constexpr AccelerationAUD2(double aud2) noexcept : value_aud2_(aud2) {}
    double value_aud2_;
};

// ============================================================================
// Gravitational Parameter (GM)
// ============================================================================
class GravitationalParameter {
public:
    // Constructors
    constexpr GravitationalParameter() noexcept : value_m3s2_(0.0) {}
    
    // Factory methods
    [[nodiscard]] static constexpr GravitationalParameter from_si(double m3s2)       noexcept { return GravitationalParameter(m3s2); }
    [[nodiscard]] static constexpr GravitationalParameter from_m3_s2(double m3s2)    noexcept { return GravitationalParameter(m3s2); }
    [[nodiscard]] static constexpr GravitationalParameter from_km3_s2(double km3s2)  noexcept { return GravitationalParameter(km3s2 * 1e9); }
    [[nodiscard]] static constexpr GravitationalParameter from_au3_d2(double au3d2)  noexcept { 
        // 1 AU³ / day² = (149597870700)^3 m³ / (86400)² s²
        constexpr double au_m = constants::AU * 1000.0;
        constexpr double factor = (au_m * au_m * au_m) / (86400.0 * 86400.0);
        return GravitationalParameter(au3d2 * factor); 
    }

    // Extraction methods
    [[nodiscard]] constexpr double to_m3_s2()  const noexcept { return value_m3s2_; }
    [[nodiscard]] constexpr double to_km3_s2() const noexcept { return value_m3s2_ / 1e9; }
    [[nodiscard]] constexpr double to_au3_d2() const noexcept { 
        constexpr double au_m = constants::AU * 1000.0;
        constexpr double factor = (au_m * au_m * au_m) / (86400.0 * 86400.0);
        return value_m3s2_ / factor; 
    }

    // Built-in standard GMs
    [[nodiscard]] static constexpr GravitationalParameter sun() noexcept {
        return from_km3_s2(constants::GM_SUN); // 132712440041.939400 km³/s²
    }

private:
    explicit constexpr GravitationalParameter(double m3s2) noexcept : value_m3s2_(m3s2) {}
    double value_m3s2_; // Internal representation is always m³/s² (SI)
};

// ============================================================================
// Constants
// ============================================================================

struct SpeedOfLight {
    // Non-instantiable static class
    SpeedOfLight() = delete;
    
    [[nodiscard]] static constexpr Velocity as_velocity() noexcept { 
        return Velocity::from_km_s(constants::C_LIGHT); 
    }
    [[nodiscard]] static constexpr double to_ms() noexcept { 
        return constants::C_LIGHT * 1000.0; 
    }
    [[nodiscard]] static constexpr double to_au_d() noexcept { 
        return (constants::C_LIGHT * 1000.0) * 86400.0 / (constants::AU * 1000.0);
    }
};

} // namespace astdyn::physics

#endif // ASTDYN_CORE_PHYSICS_TYPES_HPP
