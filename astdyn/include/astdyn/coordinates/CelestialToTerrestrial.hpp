/**
 * @file CelestialToTerrestrial.hpp
 * @brief Typed IAU 2006/2000B celestial-to-terrestrial transformation.
 *
 *     [ITRF] = W(t) * R(ERA) * Q(t) * [GCRF]
 *
 * Every rotation below carries its endpoints AND its time scale in the type, so
 * the two mistakes that actually cost us cannot be written:
 *
 *   era_rotation(ut1) * p_gcrf          // ERROR: Rotation<CIRS,TIRS> * Vec3<GCRF>
 *   cip_rotation(some_ut1_epoch)        // ERROR: EpochUT1 is not EpochTT
 *
 * Validated against ERFA/SOFA (iauC2i06a): residual < 0.007" (0.19 m) 2000-2050.
 */
#ifndef ASTDYN_COORDINATES_CELESTIAL_TO_TERRESTRIAL_HPP
#define ASTDYN_COORDINATES_CELESTIAL_TO_TERRESTRIAL_HPP

#include "astdyn/math/Rotation.hpp"
#include "astdyn/math/Vec3.hpp"
#include "astdyn/core/frame_tags.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/coordinates/GeodeticTypes.hpp"
#include "astdyn/core/physics_types.hpp"
#include <optional>

namespace astdyn::coordinates {

/// Bias-precession-nutation Q(t): GCRF -> CIRS. Requires TT.
[[nodiscard]] math::Rotation<core::GCRF, core::CIRS> cip_rotation(time::EpochTT tt);

/// Earth rotation R(ERA): CIRS -> TIRS. Requires UT1, never UTC.
[[nodiscard]] math::Rotation<core::CIRS, core::TIRS> era_rotation(time::EpochUT1 ut1);

/// Polar motion W(t): TIRS -> ITRF. Currently identity (< 15 m).
[[nodiscard]] math::Rotation<core::TIRS, core::ITRF> polar_motion_rotation();

/// The full chain. Composition order is enforced by the type system.
[[nodiscard]] math::Rotation<core::GCRF, core::ITRF> gcrf_to_itrf(time::EpochTT tt, time::EpochUT1 ut1);

/// Earth Rotation Angle (IAU 2000).
[[nodiscard]] astrometry::Angle earth_rotation_angle(time::EpochUT1 ut1);

/// CIP coordinates X, Y in GCRS (IAU 2006/2000B).
void cip_xy(time::EpochTT tt, astrometry::Angle& X, astrometry::Angle& Y);

/// A point on the Earth's surface. The latitude type is geodetic by
/// construction, so a geocentric one cannot be stored here by accident.
struct SurfacePoint {
    GeodeticLatitude lat;
    astrometry::Angle lon;
};

/**
 * @brief Intersect the shadow axis with the ellipsoid.
 *
 * The axis passes through the fundamental-plane point (xi, eta) parallel to the
 * star direction. Returns the star-facing intersection, or nullopt on a miss.
 *
 * Note the signature: the star direction is a Direction<GCRF>, so a vector in
 * any other frame is rejected at compile time, and the two epochs cannot be
 * swapped because EpochTT and EpochUT1 are distinct types.
 */
[[nodiscard]] std::optional<SurfacePoint> shadow_point(
    physics::Distance xi, physics::Distance eta,
    const math::Direction<core::GCRF>& star_dir,
    time::EpochTT tt, time::EpochUT1 ut1,
    const Ellipsoid& ellipsoid = Ellipsoid::wgs84());

} // namespace astdyn::coordinates

#endif
