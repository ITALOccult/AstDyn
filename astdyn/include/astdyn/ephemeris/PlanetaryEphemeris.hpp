/**
 * @file PlanetaryEphemeris.hpp
 * @brief Planetary ephemerides using analytical approximations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Provides heliocentric positions and velocities for solar system planets.
 * 
 * @warning Analytical VSOP87/Simon 1994 models are DEPRECATED and DISABLED for 
 * precision reasons. High-precision orbit determination requires 
 * external providers (SPICE/DE441).
 * 
 * References:
 * - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998)
 * - Simon et al. "Numerical expressions for precession formulae..." (A&A 1994)
 * - VSOP87 theory (Bretagnon & Francou, 1988)
 */

#ifndef ASTDYN_EPHEMERIS_PLANETARYEPHEMERIS_HPP
#define ASTDYN_EPHEMERIS_PLANETARYEPHEMERIS_HPP

#include "astdyn/core/Constants.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "src/core/frame_tags.hpp"
#include "src/core/units.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/ephemeris/PlanetaryData.hpp"
#include "astdyn/ephemeris/EphemerisProvider.hpp"
#include <vector>
#include <memory>

namespace astdyn {
namespace ephemeris {

/**
 * @brief Planetary ephemeris gateway (DEPRECATED: internal analytical models)
 * 
 * Provides heliocentric positions and velocities. 
 * 
 * @note This class no longer provides internal analytical approximations (VSOP87).
 * It acts as a gateway for high-precision providers (DE441/Horizons).
 * Call PlanetaryEphemeris::setProvider() before any query.
 */
class PlanetaryEphemeris {
public:
    /**
     * @brief Construct ephemeris calculator
     */
    /**
     * @param body Celestial body
     * @param t Epoch (TDB)
     * @return Position vector in J2000 equatorial frame (ICRF)
     */
    static math::Vector3<core::GCRF, physics::Distance> getPosition(CelestialBody body, time::EpochTDB t);
    
    /**
     * @param body Celestial body
     * @param t Epoch (TDB)
     * @return Velocity vector in J2000 equatorial frame (ICRF)
     */
    static math::Vector3<core::GCRF, physics::Velocity> getVelocity(CelestialBody body, time::EpochTDB t);
    
    /**
     * @param body Celestial body
     * @param t Epoch (TDB)
     * @return CartesianState in J2000 equatorial frame (GCRF)
     */
    static physics::CartesianStateTyped<core::GCRF> getState(CelestialBody body, time::EpochTDB t);
    
    /**
     * @param t Epoch (TDB)
     * @return Position vector of Sun relative to SSB
     * 
     * Computed as weighted sum of planetary positions.
     * Useful for high-precision orbit determination.
     */
    static math::Vector3<core::GCRF, physics::Distance> getSunBarycentricPosition(time::EpochTDB t);
    
    /**
     * @brief Convert heliocentric to barycentric state
     * @param heliocentric_state State relative to Sun
     * @param t Epoch (TDB)
     * @return State relative to Solar System Barycenter
     */
    static physics::CartesianStateTyped<core::GCRF> heliocentricToBarycentric(
        const physics::CartesianStateTyped<core::GCRF>& heliocentric_state, 
        time::EpochTDB t
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
    
    /**
     * @return Current ephemeris provider (usually DE441)
     */
    static std::shared_ptr<EphemerisProvider> getProvider() { return global_provider_; }

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
    [[deprecated("VSOP87 analytical elements are no longer supported. Use a provider.")]]
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
