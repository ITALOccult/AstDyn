#ifndef ASTDYN_ASTROMETRY_SKY_TYPES_HPP
#define ASTDYN_ASTROMETRY_SKY_TYPES_HPP

#include "astdyn/core/Constants.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>

namespace astdyn::astrometry {

using namespace astdyn::physics;

// ============================================================================
// Angle Type (Astrometric specific, decoupled from base physics)
// ============================================================================
class Angle {
public:
    // Constructors
    constexpr Angle() noexcept : rad_(0.0) {}
    
    // Factory
    [[nodiscard]] static constexpr Angle from_rad(double rad) noexcept { return Angle(rad); }
    [[nodiscard]] static constexpr Angle from_deg(double deg) noexcept { return Angle(deg * constants::DEG_TO_RAD); }
    [[nodiscard]] static constexpr Angle from_arcsec(double arcsec) noexcept { 
        return Angle(arcsec / 3600.0 * constants::DEG_TO_RAD); 
    }
    [[nodiscard]] static constexpr Angle from_mas(double mas) noexcept { 
        return Angle(mas / (3600.0 * 1000.0) * constants::DEG_TO_RAD); 
    }
    [[nodiscard]] static constexpr Angle zero() noexcept { return Angle(0.0); }

    // Extraction
    [[nodiscard]] constexpr double to_rad()    const noexcept { return rad_; }
    [[nodiscard]] constexpr double to_deg()    const noexcept { return rad_ / constants::DEG_TO_RAD; }
    [[nodiscard]] constexpr double to_arcsec() const noexcept { return to_deg() * 3600.0; }
    [[nodiscard]] constexpr double to_mas()    const noexcept { return to_arcsec() * 1000.0; }

    // Domain normalization [0, 2pi) or [-pi, pi]
    [[nodiscard]] Angle wrap_0_2pi() const noexcept {
        double r = std::fmod(rad_, constants::TWO_PI);
        if (r < 0.0) r += constants::TWO_PI;
        return Angle(r);
    }
    [[nodiscard]] Angle wrap_pi() const noexcept {
        double r = std::fmod(rad_ + constants::PI, constants::TWO_PI);
        if (r < 0.0) r += constants::TWO_PI;
        return Angle(r - constants::PI);
    }

    // Arithmetic
    [[nodiscard]] constexpr Angle operator+(const Angle& o) const noexcept { return Angle(rad_ + o.rad_); }
    [[nodiscard]] constexpr Angle operator-(const Angle& o) const noexcept { return Angle(rad_ - o.rad_); }
    [[nodiscard]] constexpr Angle operator*(double scalar) const noexcept { return Angle(rad_ * scalar); }

protected:
    explicit constexpr Angle(double rad) noexcept : rad_(rad) {}
    double rad_; // internal SI
};

// ============================================================================
// Right Ascension
// ============================================================================
class RightAscension : public Angle {
public:
    constexpr RightAscension() noexcept : Angle(0.0) {}
    explicit constexpr RightAscension(Angle a) noexcept : Angle(a.wrap_0_2pi().to_rad()) {}
    
    [[nodiscard]] static constexpr RightAscension from_rad(double rad) noexcept { 
        return RightAscension(Angle::from_rad(rad)); 
    }
    [[nodiscard]] static constexpr RightAscension from_deg(double deg) noexcept { 
        return RightAscension(Angle::from_deg(deg)); 
    }
    
    // Format: "HH MM SS.ss"
    [[nodiscard]] std::string to_hms() const {
        double total_hours = to_deg() / 15.0;
        int h = static_cast<int>(total_hours);
        double remainder = (total_hours - h) * 60.0;
        int m = static_cast<int>(remainder);
        double s = (remainder - m) * 60.0;
        
        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(2) << h << " "
            << std::setfill('0') << std::setw(2) << m << " "
            << std::fixed << std::setprecision(3) << std::setfill('0') << std::setw(6) << s;
        return oss.str();
    }
};

// ============================================================================
// Declination
// ============================================================================
class Declination : public Angle {
public:
    constexpr Declination() noexcept : Angle(0.0) {}
    explicit constexpr Declination(Angle a) noexcept : Angle(a.wrap_pi().to_rad()) {
        // Enforce [-pi/2, pi/2]
        if (rad_ > constants::PI / 2.0) rad_ = constants::PI - rad_;
        if (rad_ < -constants::PI / 2.0) rad_ = -constants::PI - rad_;
    }

    [[nodiscard]] static constexpr Declination from_rad(double rad) noexcept { 
        return Declination(Angle::from_rad(rad)); 
    }
    [[nodiscard]] static constexpr Declination from_deg(double deg) noexcept { 
        return Declination(Angle::from_deg(deg)); 
    }

    // Format: "+DD MM SS.ss"
    [[nodiscard]] std::string to_dms() const {
        double d_deg = to_deg();
        char sign = (d_deg < 0) ? '-' : '+';
        d_deg = std::abs(d_deg);
        
        int d = static_cast<int>(d_deg);
        double remainder = (d_deg - d) * 60.0;
        int m = static_cast<int>(remainder);
        double s = (remainder - m) * 60.0;
        
        std::ostringstream oss;
        oss << sign << std::setfill('0') << std::setw(2) << d << " "
            << std::setfill('0') << std::setw(2) << m << " "
            << std::fixed << std::setprecision(2) << std::setfill('0') << std::setw(5) << s;
        return oss.str();
    }
};

// ============================================================================
// Parallax
// ============================================================================
class Parallax {
public:
    constexpr Parallax() noexcept : mas_(0.0) {}
    
    [[nodiscard]] static constexpr Parallax from_mas(double mas) noexcept { return Parallax(mas); }
    [[nodiscard]] static constexpr Parallax from_arcsec(double sec) noexcept { return Parallax(sec * 1000.0); }
    
    [[nodiscard]] constexpr double to_mas() const noexcept { return mas_; }
    [[nodiscard]] constexpr double to_arcsec() const noexcept { return mas_ / 1000.0; }
    
    /// Convert to physical Distance (returns 0 if parallax <= 0)
    [[nodiscard]] constexpr physics::Distance to_distance() const noexcept {
        if (mas_ <= 0.0) return physics::Distance::zero();
        // Distance (pc) = 1 / parallax (arcsec)
        // 1 pc in AU is defined as 1 / tan(1 arcsec) ~ 1 / arcsen(1 rad)
        // We use the reciprocal relationship directly.
        double dist_au = 1.0 / (to_arcsec() * constants::ARCSEC_TO_RAD);
        return physics::Distance::from_au(dist_au);
    }

private:
    explicit constexpr Parallax(double mas) noexcept : mas_(mas) {}
    double mas_;
};

// ============================================================================
// Proper Motion
// ============================================================================
class ProperMotion {
public:
    constexpr ProperMotion() noexcept : mas_yr_(0.0) {}
    
    [[nodiscard]] static constexpr ProperMotion from_mas_yr(double val) noexcept { return ProperMotion(val); }
    
    [[nodiscard]] constexpr double to_mas_yr() const noexcept { return mas_yr_; }
    [[nodiscard]] constexpr double to_arcsec_yr() const noexcept { return mas_yr_ / 1000.0; }
    
    /// Compute angular displacement over time
    [[nodiscard]] Angle multiply_time(double years) const noexcept {
        return Angle::from_mas(mas_yr_ * years);
    }

    [[nodiscard]] constexpr double to_rad_yr() const noexcept {
        return (mas_yr_ / 3600000.0) * constants::DEG_TO_RAD;
    }

private:
    explicit constexpr ProperMotion(double val) noexcept : mas_yr_(val) {}
    double mas_yr_;
};

// ============================================================================
// Sky Coordinate (RA/Dec + Distance bound to a specific frame)
// ============================================================================
template <typename TargetFrame>
class SkyCoord {
public:
    [[nodiscard]] static SkyCoord from_rad(double ra_rad, double dec_rad) {
        return SkyCoord(RightAscension(Angle::from_rad(ra_rad)), Declination(Angle::from_rad(dec_rad)));
    }

    [[nodiscard]] static SkyCoord from_deg(double ra_deg, double dec_deg) {
        return SkyCoord(RightAscension(Angle::from_deg(ra_deg)), Declination(Angle::from_deg(dec_deg)));
    }

    // The ONLY way to build a SkyCoord from a distance is from a strongly-typed physical vector
    // The compiler statically verifies that the vector is in TargetFrame
    template <typename PhysUnit>
    [[nodiscard]] static SkyCoord from_vector(const math::Vector3<TargetFrame, PhysUnit>& rho_vec) {
        Eigen::Vector3d rho = rho_vec.to_eigen_si();
        double dist_si = rho.norm();
        if (dist_si < 1e-12) {
            return SkyCoord(RightAscension(Angle::zero()), Declination(Angle::zero()), PhysUnit::from_si(0.0));
        }
        
        Eigen::Vector3d u = rho / dist_si;
        double dec_rad = std::asin(u.z());
        double ra_rad = std::atan2(u.y(), u.x());
        
        return SkyCoord(
            RightAscension(Angle::from_rad(ra_rad)),
            Declination(Angle::from_rad(dec_rad)),
            PhysUnit::from_si(dist_si)
        );
    }

    // Extraction
    [[nodiscard]] const RightAscension& ra() const noexcept { return ra_; }
    [[nodiscard]] const Declination& dec()   const noexcept { return dec_; }
    
    // Separation computation using Haversine or simple dot product
    // Guaranteed to return a valid Angle (e.g. for residuals calculation)
    [[nodiscard]] Angle separation(const SkyCoord& other) const noexcept {
        double d_ra = other.ra().to_rad() - ra().to_rad();
        double d_dec = other.dec().to_rad() - dec().to_rad();
        
        double a = std::sin(d_dec / 2.0) * std::sin(d_dec / 2.0) +
                   std::cos(dec().to_rad()) * std::cos(other.dec().to_rad()) *
                   std::sin(d_ra / 2.0) * std::sin(d_ra / 2.0);
        return Angle::from_rad(2.0 * std::asin(std::sqrt(a)));
    }

private:
    SkyCoord(RightAscension ra, Declination dec) : ra_(ra), dec_(dec) {}

    template <typename PhysUnit>
    SkyCoord(RightAscension ra, Declination dec, PhysUnit dist) 
        : ra_(ra), dec_(dec) {}
    
    RightAscension ra_;
    Declination dec_;
    // We intentionally don't store distance to avoid generic typing complexities.
    // It is evaluated during construction. If distance query is required later,
    // the generic type PhysUnit must be preserved at class template level.
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_SKY_TYPES_HPP
