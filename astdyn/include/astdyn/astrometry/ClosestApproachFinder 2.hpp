#ifndef ASTDYN_ASTROMETRY_CLOSEST_APPROACH_FINDER_HPP
#define ASTDYN_ASTROMETRY_CLOSEST_APPROACH_FINDER_HPP

#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/astrometry/IChebyshevEphemeris.hpp"
#include "astdyn/catalog/CatalogTypes.hpp"
#include "astdyn/time/epoch.hpp"
#include <vector>
#include <optional>

namespace astdyn::astrometry {

// ============================================================================
// Result
// ============================================================================

/**
 * @brief Result of a closest approach search between an asteroid and a star.
 *
 * All angular quantities are in the GCRF/ICRS frame, consistent with the
 * Chebyshev ephemeris output (geocentric apparent RA/Dec).
 */
struct ClosestApproachResult {
    time::EpochTDB  t_ca;                  ///< Epoch of closest approach [TDB]
    Angle           separation;            ///< Angular separation at t_ca
    double          position_angle_deg;    ///< PA of star w.r.t. asteroid [deg, N=0, E=90]

    /// Apparent angular velocity of asteroid relative to star at t_ca [arcsec/min]
    double          relative_rate_arcsec_min = 0.0;

    /// True if separation < apparent angular radius of asteroid (requires diameter)
    bool            is_occultation = false;

    /// Apparent diameter of asteroid at t_ca [arcsec]. Zero if not provided.
    double          ast_apparent_diameter_arcsec = 0.0;
};

// ============================================================================
// Finder
// ============================================================================

/**
 * @brief Finds all closest approach minima between an asteroid Chebyshev
 *        ephemeris and a catalog star in a given time window.
 *
 * Algorithm:
 *  1. Coarse uniform sampling of delta(t) = angular_separation(asteroid, star)
 *     over [t1, t2] using n_coarse points.
 *  2. Local minima detection among consecutive triplets.
 *  3. Brent refinement on d(delta)/dt in each bracketing interval, using the
 *     analytic velocity from ChebyshevSegment::evaluate_full().
 *
 * The angular separation derivative is computed analytically from the
 * Chebyshev velocity (vRA, vDec), avoiding numerical differentiation.
 *
 * Proper motion and parallax of the star are applied if the Star record
 * contains valid values (uses Star::predict_at()).
 */
class ClosestApproachFinder {
public:
    /**
     * @brief Find all closest approaches using a full IChebyshevEphemeris.
     *
     * @param ephem        Chebyshev ephemeris covering [t1, t2]
     * @param star         Gaia DR3 star record (proper motion applied if available)
     * @param t1           Search window start
     * @param t2           Search window end
     * @param max_sep      Only return results with separation < max_sep (default: 120 arcsec)
     * @param n_coarse     Number of coarse sampling points (default: 1440, ~1 min for 1-day window)
     * @param ast_diameter_km  Optional: asteroid diameter in km for occultation flag
     * @return             All minima with separation < max_sep, sorted by time
     */
    static std::vector<ClosestApproachResult> find(
        const IChebyshevEphemeris&  ephem,
        const catalog::Star&        star,
        time::EpochTDB              t1,
        time::EpochTDB              t2,
        Angle                       max_sep       = Angle::from_arcsec(120.0),
        int                         n_coarse      = 1440,
        std::optional<double>       ast_diameter_km = std::nullopt);

    /**
     * @brief Find closest approaches using a single ChebyshevSegment.
     *
     * Useful when you already have the relevant segment for the time window.
     * The segment must cover [t1, t2] (in JD TDB).
     */
    static std::vector<ClosestApproachResult> find_in_segment(
        const catalog::ChebyshevSegment& segment,
        const catalog::Star&             star,
        time::EpochTDB                   t1,
        time::EpochTDB                   t2,
        Angle                            max_sep       = Angle::from_arcsec(120.0),
        int                              n_coarse      = 1440,
        std::optional<double>            ast_diameter_km = std::nullopt);

private:
    // -------------------------------------------------------------------------
    // Internal state passed through the refinement lambdas
    // -------------------------------------------------------------------------
    struct EvalContext {
        const catalog::ChebyshevSegment* segment;
        const catalog::Star*             star;
        bool                             use_proper_motion;
        bool                             use_parallax;
    };

    // -------------------------------------------------------------------------
    // Core evaluators
    // -------------------------------------------------------------------------

    /**
     * @brief Compute angular separation delta(jd) and its time derivative.
     *
     * Uses the analytic Chebyshev velocity (vRA [deg/day], vDec [deg/day])
     * to compute d(delta)/dt analytically — no finite differences.
     *
     * @param jd   Julian Date [TDB]
     * @param ctx  Evaluation context
     * @return     {delta [rad], d(delta)/dt [rad/day]}
     */
    static std::pair<double, double> separation_and_derivative(
        double jd, const EvalContext& ctx);

    /**
     * @brief Brent's method to find the root of d(delta)/dt in [ja, jb].
     *
     * Assumes d(delta)/dt changes sign in [ja, jb] (i.e. a minimum exists).
     *
     * @param ja   Left bracket [JD]
     * @param jb   Right bracket [JD]
     * @param ctx  Evaluation context
     * @param tol  Convergence tolerance [days] (default ~0.01 s)
     * @return     JD of closest approach
     */
    static double brent_minimum(
        double ja, double jb,
        const EvalContext& ctx,
        double tol = 1.2e-7);

    /**
     * @brief Build a ClosestApproachResult at a given JD.
     */
    static ClosestApproachResult make_result(
        double jd,
        const EvalContext& ctx,
        std::optional<double> ast_diameter_km);
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_CLOSEST_APPROACH_FINDER_HPP
