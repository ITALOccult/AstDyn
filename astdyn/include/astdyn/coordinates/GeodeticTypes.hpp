/**
 * @file GeodeticTypes.hpp
 * @brief Latitudes that know which latitude they are.
 *
 * astrometry::Angle is dimensionally strong but semantically blind: geocentric and geodetic
 * latitude are both "an angle", so
 *
 *     lat = std::atan2(z_ecef, dist_xy);      // geocentric
 *     return GeoPoint{astrometry::Angle::from_rad(lat)};  // consumed as geodetic
 *
 * compiled happily and was wrong by up to 21.4 km at mid-latitudes. These two
 * types make that assignment a compile error, and the only way across is an
 * explicit, named conversion that requires the ellipsoid.
 */
#ifndef ASTDYN_COORDINATES_GEODETIC_TYPES_HPP
#define ASTDYN_COORDINATES_GEODETIC_TYPES_HPP

#include "astdyn/astrometry/sky_types.hpp"
#include <cmath>

namespace astdyn::coordinates {

/// astrometry::Angle at the geocentre between the equatorial plane and the position vector.
class GeocentricLatitude {
public:
    GeocentricLatitude() = default;
    [[nodiscard]] static GeocentricLatitude from_angle(astrometry::Angle a) noexcept { return GeocentricLatitude(a); }
    [[nodiscard]] astrometry::Angle angle() const noexcept { return a_; }
    [[nodiscard]] double to_deg() const noexcept { return a_.to_deg(); }
    [[nodiscard]] double to_rad() const noexcept { return a_.to_rad(); }
private:
    explicit GeocentricLatitude(astrometry::Angle a) noexcept : a_(a) {}
    astrometry::Angle a_{};
};

/// astrometry::Angle between the equatorial plane and the ellipsoid normal. What maps show,
/// what GPS reports, and what Occult4 prints.
class GeodeticLatitude {
public:
    GeodeticLatitude() = default;
    [[nodiscard]] static GeodeticLatitude from_angle(astrometry::Angle a) noexcept { return GeodeticLatitude(a); }
    [[nodiscard]] astrometry::Angle angle() const noexcept { return a_; }
    [[nodiscard]] double to_deg() const noexcept { return a_.to_deg(); }
    [[nodiscard]] double to_rad() const noexcept { return a_.to_rad(); }
private:
    explicit GeodeticLatitude(astrometry::Angle a) noexcept : a_(a) {}
    astrometry::Angle a_{};
};

/// Reference ellipsoid. Centralises the flattening that bug #4 ignored.
class Ellipsoid {
public:
    constexpr Ellipsoid(double a_m, double f) noexcept : a_(a_m), f_(f) {}
    [[nodiscard]] static Ellipsoid wgs84() noexcept { return Ellipsoid(6378137.0, 1.0/298.257223563); }
    [[nodiscard]] constexpr double a_m() const noexcept { return a_; }
    [[nodiscard]] constexpr double f() const noexcept { return f_; }
    [[nodiscard]] constexpr double e2() const noexcept { return 2.0*f_ - f_*f_; }
    [[nodiscard]] constexpr double b_m() const noexcept { return a_ * (1.0 - f_); }

    /// The conversion is only possible with an ellipsoid, so it lives here.
    /// Exact for points lying on the surface.
    [[nodiscard]] GeodeticLatitude to_geodetic(GeocentricLatitude phi) const noexcept {
        return GeodeticLatitude::from_angle(astrometry::Angle::from_rad(
            std::atan2(std::sin(phi.to_rad()), (1.0 - e2()) * std::cos(phi.to_rad()))));
    }
private:
    double a_, f_;
};

} // namespace astdyn::coordinates
#endif
/**
 * @file GeodeticTypes.hpp
 * @brief Latitudes that know which latitude they are.
 *
 * Angle is dimensionally strong but semantically blind: geocentric and geodetic
 * latitude are both "an angle", so
 *
 *     lat = std::atan2(z_ecef, dist_xy);      // geocentric
 *     return GeoPoint{Angle::from_rad(lat)};  // consumed as geodetic
 *
 * compiled happily and was wrong by up to 21.4 km at mid-latitudes. These two
 * types make that assignment a compile error, and the only way across is an
 * explicit, named conversion that requires the ellipsoid.
 */
#ifndef ASTDYN_COORDINATES_GEODETIC_TYPES_HPP
#define ASTDYN_COORDINATES_GEODETIC_TYPES_HPP

#include "astdyn/astrometry/sky_types.hpp"
#include <cmath>

namespace astdyn::coordinates {

/// Angle at the geocentre between the equatorial plane and the position vector.
class GeocentricLatitude {
public:
    GeocentricLatitude() = default;
    [[nodiscard]] static GeocentricLatitude from_angle(Angle a) noexcept { return GeocentricLatitude(a); }
    [[nodiscard]] Angle angle() const noexcept { return a_; }
    [[nodiscard]] double to_deg() const noexcept { return a_.to_deg(); }
    [[nodiscard]] double to_rad() const noexcept { return a_.to_rad(); }
private:
    explicit GeocentricLatitude(Angle a) noexcept : a_(a) {}
    Angle a_{};
};

/// Angle between the equatorial plane and the ellipsoid normal. What maps show,
/// what GPS reports, and what Occult4 prints.
class GeodeticLatitude {
public:
    GeodeticLatitude() = default;
    [[nodiscard]] static GeodeticLatitude from_angle(Angle a) noexcept { return GeodeticLatitude(a); }
    [[nodiscard]] Angle angle() const noexcept { return a_; }
    [[nodiscard]] double to_deg() const noexcept { return a_.to_deg(); }
    [[nodiscard]] double to_rad() const noexcept { return a_.to_rad(); }
private:
    explicit GeodeticLatitude(Angle a) noexcept : a_(a) {}
    Angle a_{};
};

/// Reference ellipsoid. Centralises the flattening that bug #4 ignored.
class Ellipsoid {
public:
    constexpr Ellipsoid(double a_m, double f) noexcept : a_(a_m), f_(f) {}
    [[nodiscard]] static Ellipsoid wgs84() noexcept { return Ellipsoid(6378137.0, 1.0/298.257223563); }
    [[nodiscard]] constexpr double a_m() const noexcept { return a_; }
    [[nodiscard]] constexpr double f() const noexcept { return f_; }
    [[nodiscard]] constexpr double e2() const noexcept { return 2.0*f_ - f_*f_; }
    [[nodiscard]] constexpr double b_m() const noexcept { return a_ * (1.0 - f_); }

    /// The conversion is only possible with an ellipsoid, so it lives here.
    /// Exact for points lying on the surface.
    [[nodiscard]] GeodeticLatitude to_geodetic(GeocentricLatitude phi) const noexcept {
        return GeodeticLatitude::from_angle(Angle::from_rad(
            std::atan2(std::sin(phi.to_rad()), (1.0 - e2()) * std::cos(phi.to_rad()))));
    }
private:
    double a_, f_;
};

} // namespace astdyn::coordinates
#endif
