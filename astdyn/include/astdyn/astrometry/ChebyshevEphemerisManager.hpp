/**
 * @file ChebyshevEphemerisManager.hpp
 * @brief Management of Chebyshev ephemerides for multiple asteroids.
 */

#ifndef ASTDYN_ASTROMETRY_CHEBYSHEV_EPHEMERIS_MANAGER_HPP
#define ASTDYN_ASTROMETRY_CHEBYSHEV_EPHEMERIS_MANAGER_HPP

#include "astdyn/astrometry/AsteroidChebyshevEphemeris.hpp"
#include <map>
#include <string>
#include <memory>
#include <vector>

namespace astdyn::astrometry {

/**
 * @brief High-level manager to handle Chebyshev polynomials for a collection of asteroids.
 *
 * This class coordinates the pre-calculation and retrieval of ephemerides for multiple
 * bodies, optimizing the transition from numerical propagation to analytical interpolation.
 */
class ChebyshevEphemerisManager {
public:
    ChebyshevEphemerisManager(const AstDynConfig& config) : config_(config) {}

    /**
     * @brief Add an asteroid to the managed list and calculate its polynomials.
     * 
     * @param id Asteroid designation or name.
     * @param initial_elements Orbital elements.
     * @param start Start of interval.
     * @param end End of interval.
     * @param degree Polynomial degree (default: 12).
     */
    void add_asteroid(
        const std::string& id,
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
        time::EpochTDB start,
        time::EpochTDB end,
        int degree = 12
    );

    /**
     * @brief Evaluate position and velocity for a specific asteroid.
     * 
     * @param id Asteroid designation.
     * @param t Target epoch.
     * @return pair of tuples: (RA deg, Dec deg, Dist AU), (vRA deg/day, vDec deg/day, vDist AU/day)
     */
    std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>>
    evaluate_full(const std::string& id, time::EpochTDB t) const;

    /**
     * @brief Check if an asteroid is managed.
     */
    bool has_asteroid(const std::string& id) const {
        return ephemerides_.find(id) != ephemerides_.end();
    }

    /**
     * @brief Clear all managed ephemerides.
     */
    void clear() { ephemerides_.clear(); }

    /**
     * @brief Get count of managed asteroids.
     */
    size_t size() const { return ephemerides_.size(); }

    /**
     * @brief Get underlying ephemeris for an asteroid.
     */
    const AsteroidChebyshevEphemeris& get_asteroid(const std::string& id) const;

private:
    const AstDynConfig& config_;
    std::map<std::string, std::unique_ptr<AsteroidChebyshevEphemeris>> ephemerides_;
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_CHEBYSHEV_EPHEMERIS_MANAGER_HPP
