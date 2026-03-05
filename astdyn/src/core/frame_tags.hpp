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
 * @brief Tag for Ecliptic J2000 frame.
 */
struct ECLIPJ2000 {};

} // namespace astdyn::core

#endif // ASTDYN_CORE_FRAME_TAGS_HPP
