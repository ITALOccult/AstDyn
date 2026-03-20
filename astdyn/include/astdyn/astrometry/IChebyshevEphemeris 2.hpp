#ifndef ASTDYN_ASTROMETRY_I_CHEBYSHEV_EPHEMERIS_HPP
#define ASTDYN_ASTROMETRY_I_CHEBYSHEV_EPHEMERIS_HPP

#include "astdyn/catalog/CatalogTypes.hpp"
#include "astdyn/time/epoch.hpp"
#include <tuple>
#include <utility>

namespace astdyn::astrometry {

/**
 * @brief Interface for high-performance asteroid ephemeris using Chebyshev segments.
 */
class IChebyshevEphemeris {
public:
    virtual ~IChebyshevEphemeris() = default;

    /**
     * @brief Evaluate position at a specific epoch.
     * @return {RA [deg], Dec [deg], Distance [AU]}
     */
    virtual std::tuple<double, double, double> evaluate(time::EpochTDB epoch) const = 0;

    /**
     * @brief Evaluate position and analytical velocity at a specific epoch.
     * @return {{RA, Dec, Dist}, {vRA, vDec, vDist}}
     */
    virtual std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>> 
    evaluate_full(time::EpochTDB epoch) const = 0;

    virtual time::EpochTDB start_epoch() const = 0;
    virtual time::EpochTDB end_epoch() const = 0;
    virtual size_t num_segments() const = 0;
    virtual const catalog::ChebyshevSegment& get_segment(time::EpochTDB epoch) const = 0;
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_I_CHEBYSHEV_EPHEMERIS_HPP
