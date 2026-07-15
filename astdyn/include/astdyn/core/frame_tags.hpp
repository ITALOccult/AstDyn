#ifndef ASTDYN_CORE_FRAME_TAGS_HPP
#define ASTDYN_CORE_FRAME_TAGS_HPP

namespace astdyn::core {

/**
 * @brief Tag for Geocentric Celestial Reference Frame.
 */
struct GCRF {};

/**
 * @brief Tag for International Terrestrial Reference Frame.
 */
struct ITRF {};

/**
 * @brief Tag for True Equator Mean Equinox frame.
 */
struct TEME {};

/**
 * @brief Tag for the Celestial Intermediate Reference System.
 *
 * The frame the Earth Rotation Angle rotates *from*. Its absence from this list
 * is why the GCRF/CIRS confusion could not even be named, let alone caught by
 * the compiler.
 */
struct CIRS {};

/**
 * @brief Tag for the Terrestrial Intermediate Reference System.
 *
 * ITRF before polar motion is applied. ERA maps CIRS -> TIRS.
 */
struct TIRS {};

/**
 * @brief Tag for Ecliptic J2000 frame.
 */
struct ECLIPJ2000 {};

} // namespace astdyn::core

#endif // ASTDYN_CORE_FRAME_TAGS_HPP
