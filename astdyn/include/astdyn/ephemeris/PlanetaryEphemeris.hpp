/**
 * @file PlanetaryEphemeris.hpp
 * @brief Planetary ephemerides using analytical approximations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Provides heliocentric positions and velocities for solar system planets.
 * 
 * Implementation uses simplified VSOP87 series expansions suitable for
 * accuracy of ~1 arcsec (Earth) to ~10 arcsec (outer planets) over
 * the period 1800-2050.
 * 
 * For higher accuracy, consider integrating SPICE toolkit (cspice).
 * 
 * References:
 * - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998)
 * - Simon et al. "Numerical expressions for precession formulae..." (A&A 1994)
 * - VSOP87 theory (Bretagnon & Francou, 1988)
 */

#ifndef ASTDYN_EPHEMERIS_PLANETARYEPHEMERIS_HPP
#define ASTDYN_EPHEMERIS_PLANETARYEPHEMERIS_HPP

#include "astdyn/core/Constants.hpp"
#include "src/utils/time_types.hpp"
#include "src/types/vectors.hpp"
#include "src/core/frame_tags.hpp"
#include "src/core/units.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/ephemeris/PlanetaryData.hpp"
#include "astdyn/ephemeris/EphemerisProvider.hpp"
#include <vector>
#include <memory>

namespace astdyn {
namespace ephemeris {

/**
 * @brief Planetary ephemeris calculator using analytical approximations
 * 
 * Provides heliocentric positions and velocities in J2000 ecliptic frame.
 * 
 * Coordinate System:
 * - Reference frame: J2000.0 ecliptic (FK5)
 * - Origin: Solar System Barycenter (approximated as Sun center for planets)
 * - Units: AU for position, AU/day for velocity
 * 
 * Accuracy:
 * - Inner planets (Mercury-Mars): ~1-5 arcsec over 1800-2050
 * - Outer planets (Jupiter-Neptune): ~5-20 arcsec over 1800-2050
 * - Sufficient for asteroid orbit determination, preliminary trajectory design
 * 
 * For sub-arcsecond accuracy, use JPL SPICE kernels (DE440/441).
 */
class PlanetaryEphemeris {
public:
    /**
     * @brief Construct ephemeris calculator
     */
    /**
     * @param body Celestial body
     * @param t Instant
     * @return Position vector [m] in J2000 equatorial frame (ICRF)
     */
    static types::Vector3<core::GCRF, core::Meter> getPosition(CelestialBody body, utils::Instant t);
    
    /**
     * @param body Celestial body
     * @param t Instant
     * @return Velocity vector [m/s] in J2000 equatorial frame (ICRF)
     */
    static types::Vector3<core::GCRF, core::Meter> getVelocity(CelestialBody body, utils::Instant t);
    
    /**
     * @param body Celestial body
     * @param t Instant
     * @return CartesianState in J2000 equatorial frame [m, m/s]
     */
    static coordinates::CartesianState getState(CelestialBody body, utils::Instant t);
    
    /**
     * @param t Instant
     * @return Position vector [m] of Sun relative to SSB
     * 
     * Computed as weighted sum of planetary positions.
     * Useful for high-precision orbit determination.
     */
    static types::Vector3<core::GCRF, core::Meter> getSunBarycentricPosition(utils::Instant t);
    
    /**
     * @brief Convert heliocentric to barycentric state
     * @param heliocentric_state State relative to Sun
     * @param t Time
     * @return State relative to Solar System Barycenter
     */
    static coordinates::CartesianState heliocentricToBarycentric(
        const coordinates::CartesianState& heliocentric_state, 
        utils::Instant t
    );

    /**
     * @brief Set a high-precision ephemeris provider (e.g. DE441)
     * 
     * If set, getPosition/getVelocity will use this provider instead of
     * the internal analytical approximations.
     * 
     * @param provider Shared pointer to provider (can be null to reset to default)
     */
    static void setProvider(std::shared_ptr<EphemerisProvider> provider);

private:
    static std::shared_ptr<EphemerisProvider> global_provider_;

    /**
     * @brief Compute Julian centuries from J2000.0
     * @param jd_tdb Julian Date in TDB
     * @return Julian centuries (T)
     */
    static double julianCenturies(double jd_tdb) {
        return (jd_tdb - constants::JD2000) / 36525.0;
    }
    
    /**
     * @brief Compute planetary orbital elements at epoch
     * @param body Celestial body
     * @param T Julian centuries from J2000.0
     * @param elements Output: [a, e, i, L, omega_bar, Omega] in AU, radians
     * 
     * Uses Simon et al. (1994) low-precision formulae.
     * L = mean longitude, omega_bar = longitude of perihelion
     */
    static void computeOrbitalElements(
        CelestialBody body, 
        double T,
        double elements[6]
    );
    
    /**
     * @brief Convert orbital elements to Cartesian position
     * @param elements [a, e, i, L, omega_bar, Omega]
     * @return Position [AU] in ecliptic frame
     */
    static Eigen::Vector3d elementsToPosition(const double elements[6]);
    
    /**
     * @brief Convert orbital elements to Cartesian velocity
     * @param elements [a, e, i, L, omega_bar, Omega]
     * @param gm Gravitational parameter of central body (Sun) [AU³/day²]
     * @return Velocity [AU/day] in ecliptic frame
     */
    static Eigen::Vector3d elementsToVelocity(const double elements[6], double gm);
    
    /**
     * @brief Compute perturbations for high-precision planets
     * @param body Planet
     * @param T Julian centuries
     * @return Correction to position [AU]
     * 
     * Includes major perturbations (e.g., Jupiter on Saturn).
     * Currently returns zero; can be extended with VSOP87 terms.
     */
    static Eigen::Vector3d computePerturbations(CelestialBody body, double T);
};

} // namespace ephemeris
} // namespace astdyn

#endif // ASTDYN_EPHEMERIS_PLANETARYEPHEMERIS_HPP
