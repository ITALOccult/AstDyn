/**
 * @file CelestialToTerrestrial.hpp
 * @brief GCRF (ICRS) -> ITRF rotation and shadow-axis / ellipsoid intersection.
 *
 * Implements the IAU 2006/2000B CIO-based transformation
 *
 *     [ITRF] = R3(ERA) * Q(t) * [GCRF]
 *
 * where Q(t) is the celestial-to-intermediate (bias-precession-nutation) matrix.
 * Polar motion W(t) is neglected (< 15 m at the Earth's surface).
 *
 * Rationale: the Earth Rotation Angle maps CIRS -> TIRS, NOT GCRF -> ITRF.
 * Applying ERA alone omits Q(t), which by 2026 amounts to a 535" rotation of
 * the CIP, i.e. ~16.6 km of ground error, growing ~0.6 km/yr.
 *
 * Validated against ERFA/SOFA (iauC2i06a, iauGd2gc): residual rotation error
 * < 0.007" (0.19 m at the surface) over 2000-2050.
 */
#ifndef ASTDYN_COORDINATES_CELESTIAL_TO_TERRESTRIAL_HPP
#define ASTDYN_COORDINATES_CELESTIAL_TO_TERRESTRIAL_HPP

#include <Eigen/Dense>

namespace astdyn::coordinates {

/// Fukushima-Williams precession angles (IAU 2006). t = Julian centuries TT from J2000. [rad]
void precession_fw06(double t_jc, double& gamb, double& phib, double& psib, double& epsa);

/// Nutation in longitude/obliquity (IAU 2000B, 20 largest luni-solar terms). [rad]
void nutation_00b(double t_jc, double& dpsi, double& deps);

/// CIP coordinates X, Y in GCRS (IAU 2006/2000B). [rad]
void cip_xy(double jd_tt, double& X, double& Y);

/// Celestial-to-intermediate matrix Q(t): GCRS -> CIRS.
Eigen::Matrix3d celestial_to_intermediate(double jd_tt);

/// Earth Rotation Angle (IAU 2000). [rad]
double earth_rotation_angle(double jd_ut1);

/// Full rotation GCRF -> ITRF (polar motion neglected).
Eigen::Matrix3d gcrf_to_itrf(double jd_tt, double jd_ut1);

/**
 * @brief Intersect the shadow axis with the WGS84 ellipsoid.
 *
 * The shadow axis is the line through the fundamental-plane point (xi, eta)
 * parallel to the star direction. Returns the star-facing intersection.
 *
 * @param xi_m,eta_m  Fundamental-plane coordinates of the axis [m]
 * @param ra_rad,dec_rad  Star direction in GCRF [rad]
 * @param jd_tt, jd_ut1   Epoch [JD]
 * @param lat_rad  GEODETIC latitude [rad] (output)
 * @param lon_rad  Longitude [rad] (output)
 * @return false if the axis misses the Earth.
 */
bool shadow_point_geodetic(double xi_m, double eta_m,
                           double ra_rad, double dec_rad,
                           double jd_tt, double jd_ut1,
                           double& lat_rad, double& lon_rad);

} // namespace astdyn::coordinates

#endif
