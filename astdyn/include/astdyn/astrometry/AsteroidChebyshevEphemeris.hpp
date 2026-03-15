/**
 * @file AsteroidChebyshevEphemeris.hpp
 * @brief High-performance asteroid ephemeris using daily Chebyshev segments.
 */

#ifndef ASTDYN_ASTROMETRY_ASTEROID_CHEBYSHEV_EPHEMERIS_HPP
#define ASTDYN_ASTROMETRY_ASTEROID_CHEBYSHEV_EPHEMERIS_HPP

#include "astdyn/catalog/CatalogTypes.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/astrometry/IChebyshevEphemeris.hpp"
#include <vector>
#include <tuple>
#include <memory>

namespace astdyn { struct AstDynConfig; }

namespace astdyn::astrometry {

/**
 * @brief Represents an asteroid trajectory as a series of daily Chebyshev polynomials.
 *
 * This class pre-calculates positions (RA, Dec, Distance) for a given asteroid
 * over a specified time interval, fitting one Chebyshev polynomial per 24-hour period.
 * It provides extremely fast position retrieval compared to direct propagation.
 */
class AsteroidChebyshevEphemeris : public IChebyshevEphemeris {
public:
    /**
     * @brief Construct and pre-calculate ephemeris.
     *
     * @param initial_elements Orbit at some epoch.
     * @param start_time       Start of the interval.
     * @param end_time         End of the interval.
     * @param config           AstDyn configuration (propagator settings, etc.).
     * @param degree           Chebyshev polynomial degree (default: 12).
     */
    AsteroidChebyshevEphemeris(
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
        time::EpochTDB start_time,
        time::EpochTDB end_time,
        const astdyn::AstDynConfig& config,
        int degree = 12
    );

    /**
     * @brief Evaluate position at a specific epoch.
     *
     * @param epoch Target time (TDB).
     * @return Tuple containing {RA [deg], Dec [deg], Distance [AU]}.
     * @throw std::out_of_range if epoch is outside the pre-calculated interval.
     */
    std::tuple<double, double, double> evaluate(time::EpochTDB epoch) const override;

    /**
     * @brief Evaluate position and velocity (derivatives) at a specific epoch.
     * 
     * @param epoch Target time (TDB).
     * @return pair of tuples: (RA deg, Dec deg, Dist AU), (vRA deg/day, vDec deg/day, vDist AU/day)
     */
    std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>> 
    evaluate_full(time::EpochTDB epoch) const override;

    /**
     * @brief Get start epoch of the calculated interval.
     */
    time::EpochTDB start_epoch() const override { return start_epoch_; }

    /**
     * @brief Get end epoch of the calculated interval.
     */
    time::EpochTDB end_epoch() const override { return end_epoch_; }

    /**
     * @brief Get number of segments.
     */
    size_t num_segments() const override { return segments_.size(); }

    /**
     * @brief Get segment containing a specific epoch.
     */
    const catalog::ChebyshevSegment& get_segment(time::EpochTDB epoch) const override;

private:
    std::vector<catalog::ChebyshevSegment> segments_;
    time::EpochTDB start_epoch_;
    time::EpochTDB end_epoch_;
    double jd_start_;
    double jd_end_;
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_ASTEROID_CHEBYSHEV_EPHEMERIS_HPP
