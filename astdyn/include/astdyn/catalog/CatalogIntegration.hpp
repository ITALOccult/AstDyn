/**
 * @file CatalogIntegration.hpp
 * @brief Bridge between AstDyn propagation types and astdyn::catalog queries
 *
 * Provides conversion utilities to build catalog queries from AstDyn's
 * typed orbit states, enabling end-to-end occultation candidate search:
 *
 *   KeplerianStateTyped → propagate → CartesianStateTyped[]
 *                       → make_orbit_query()
 *                       → GaiaDR3Catalog::query_orbit()
 *                       → Star[] (occultation candidates)
 */

#ifndef ASTDYN_CATALOG_INTEGRATION_HPP
#define ASTDYN_CATALOG_INTEGRATION_HPP

#include "astdyn/catalog/CatalogTypes.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/time/epoch.hpp"
#include "src/core/frame_tags.hpp"
#include <vector>

namespace astdyn::catalog {

// ============================================================================
// Star ↔ AstDyn type conversions
// ============================================================================

/**
 * @brief Apply proper motion to propagate a star's position to a target epoch.
 *
 * Uses the linear proper motion model:
 *   RA(t)  = RA(J2016) + pmRA  * (t - J2016) / cos(dec)
 *   Dec(t) = Dec(J2016) + pmDec * (t - J2016)
 *
 * where t is in Julian years from J2016.0 (2457389.0 JD).
 *
 * @param star   Gaia DR3 star (astrometry at J2016.0)
 * @param epoch  Target epoch (TDB)
 * @return Sky position at the target epoch [deg, deg]
 */
[[nodiscard]] SkyPoint propagate_proper_motion(
    const Star& star,
    time::EpochTDB epoch
) noexcept;

/**
 * @brief Convert a Star to an astrometry::SkyCoord at a given epoch.
 *
 * Applies proper motion correction, then constructs a typed SkyCoord
 * in the GCRF frame (= ICRS for stellar objects).
 *
 * @param star   Gaia DR3 star
 * @param epoch  Observation epoch (TDB)
 * @return Typed sky coordinate
 */
[[nodiscard]] astrometry::SkyCoord<core::GCRF> to_sky_coord(
    const Star& star,
    time::EpochTDB epoch
);

// ============================================================================
// Orbit → OrbitQuery bridge
// ============================================================================

/**
 * @brief Build an OrbitQuery from a sequence of propagated Cartesian states.
 *
 * Converts heliocentric Cartesian states (GCRF) to geocentric RA/Dec sky
 * positions, fits Chebyshev polynomials segment-by-segment, and packages
 * the result into an OrbitQuery ready for GaiaDR3Catalog::query_orbit().
 *
 * The states must be evenly sampled in time (or at least monotonically
 * ordered). The Chebyshev approximation uses segments of ~segment_days days,
 * each fitted to degree chebyshev_degree.
 *
 * @param states          Propagated heliocentric Cartesian states (GCRF)
 * @param earth_states    Earth's Cartesian states at the same epochs (GCRF),
 *                        used to compute geocentric direction
 * @param t_start         Start of search window
 * @param t_end           End of search window
 * @param width_arcsec    Corridor half-width [arcsec] (default: 120 = 2')
 * @param max_magnitude   Faint magnitude limit (Gaia G)
 * @param segment_days    Duration of each Chebyshev segment [days]
 * @param chebyshev_degree Degree of the Chebyshev fit per segment
 * @return OrbitQuery ready for GaiaDR3Catalog::query_orbit()
 */
[[nodiscard]] OrbitQuery make_orbit_query(
    const std::vector<physics::CartesianStateTyped<core::GCRF>>& states,
    const std::vector<physics::CartesianStateTyped<core::GCRF>>& earth_states,
    time::EpochTDB t_start,
    time::EpochTDB t_end,
    double width_arcsec    = 120.0,
    double max_magnitude   = 18.0,
    double segment_days    = 10.0,
    int    chebyshev_degree = 12
);

/**
 * @brief Compute geocentric RA/Dec from heliocentric body and Earth positions.
 *
 * body_pos and earth_pos must be in the same frame (GCRF, SI meters).
 * Returns {ra_deg, dec_deg} in the ICRS/J2000 equatorial frame.
 */
[[nodiscard]] SkyPoint heliocentric_to_radec(
    const physics::CartesianStateTyped<core::GCRF>& body,
    const physics::CartesianStateTyped<core::GCRF>& earth
) noexcept;

// ============================================================================
// Convenience: full pipeline in one call
// ============================================================================

/**
 * @brief Search for occultation candidates along a propagated orbit.
 *
 * End-to-end helper that:
 *  1. Converts the trajectory to geocentric RA/Dec
 *  2. Fits Chebyshev polynomials
 *  3. Queries the Gaia DR3 catalog
 *  4. Returns candidate stars with proper motion applied to the mid-epoch
 *
 * @param catalog         Initialized GaiaDR3Catalog instance
 * @param body_states     Propagated heliocentric states of the asteroid (GCRF)
 * @param earth_states    Propagated heliocentric states of the Earth (GCRF)
 * @param t_start         Start of observation window
 * @param t_end           End of observation window
 * @param width_arcsec    Corridor half-width [arcsec]
 * @param max_magnitude   Magnitude limit (Gaia G)
 * @return Candidate stars, proper-motion-corrected to the window mid-epoch
 */
[[nodiscard]] std::vector<Star> find_occultation_candidates(
    const GaiaDR3Catalog& catalog,
    const std::vector<physics::CartesianStateTyped<core::GCRF>>& body_states,
    const std::vector<physics::CartesianStateTyped<core::GCRF>>& earth_states,
    time::EpochTDB t_start,
    time::EpochTDB t_end,
    double width_arcsec  = 120.0,
    double max_magnitude = 18.0
);

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_INTEGRATION_HPP
