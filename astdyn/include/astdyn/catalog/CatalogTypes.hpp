/**
 * @file CatalogTypes.hpp
 * @brief Core types for the astdyn::catalog stellar catalog module
 *
 * Defines the public data types used for stellar catalog queries.
 * The underlying implementation uses the Gaia DR3 catalog via IOC_GaiaLib.
 */

#ifndef ASTDYN_CATALOG_TYPES_HPP
#define ASTDYN_CATALOG_TYPES_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <optional>
#include <functional>

#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "src/core/frame_tags.hpp"

namespace astdyn::catalog {

using namespace astdyn::astrometry;

// ============================================================================
// Module-local constants
// ============================================================================
namespace constants {
    inline constexpr double AU_PER_PARSEC = 206264.806; ///< 1 pc = 206264.806 AU
}

// ============================================================================
// Star — Gaia DR3 stellar record
// ============================================================================

/**
 * @brief A single star record from the Gaia DR3 catalog.
 *
 * Positional data is in the ICRS frame (== GCRF for stellar objects).
 * Proper motions are given at the Gaia reference epoch J2016.0.
 */
struct Star {
    // --- Identification ---
    int64_t     source_id = 0;          ///< Gaia DR3 source identifier
    std::string common_name;            ///< IAU common name (if assigned)
    std::string hd_designation;         ///< Henry Draper catalog number
    std::string hip_designation;        ///< Hipparcos catalog number
    std::string tycho2_designation;     ///< Tycho-2 designation
    std::string sao_designation;        ///< SAO catalog number

    // --- Astrometry (ICRS, epoch J2016.0) ---
    RightAscension ra;                  ///< Right Ascension
    Declination    dec;                 ///< Declination
    Parallax       parallax;            ///< Parallax (mas)
    double         parallax_error_mas = 0.0;
    ProperMotion   pm_ra_cosdec;        ///< Proper motion in RA * cos(dec)
    ProperMotion   pm_dec;              ///< Proper motion in Dec
    double         pmra_error_mas_yr  = 0.0;
    double         pmdec_error_mas_yr = 0.0;

    // --- Photometry ---
    double g_mag  = 99.0;               ///< Gaia G-band magnitude
    double bp_mag = 99.0;               ///< Gaia BP-band magnitude
    double rp_mag = 99.0;               ///< Gaia RP-band magnitude
    double bp_rp  = 0.0;                ///< BP-RP color index

    // --- Quality ---
    double ruwe = 0.0;                  ///< Renormalised unit weight error
    double astrometric_excess_noise = 0.0;

    /// Predict star position at a specific TDB epoch (proper motion + parallax)
    [[nodiscard]] SkyCoord<core::GCRF> predict_at(time::EpochTDB target_time) const;

    [[nodiscard]] bool has_parallax() const noexcept { return parallax.to_mas() > 0.0; }
    [[nodiscard]] bool has_proper_motion() const noexcept {
        return pm_ra_cosdec.to_mas_yr() != 0.0 || pm_dec.to_mas_yr() != 0.0;
    }
};

// ============================================================================
// Query parameter structs
// ============================================================================

/**
 * @brief Parameters for a circular cone search.
 */
struct ConeQuery {
    RightAscension ra;           ///< Center RA
    Declination    dec;          ///< Center Dec
    Angle          radius;       ///< Search radius
    double         max_magnitude = 18.0; ///< Faint limit (Gaia G)
    Parallax       min_parallax; ///< Optional parallax cut
};

/**
 * @brief A point on the sky, used to define a corridor path.
 */
struct SkyPoint {
    double ra_deg  = 0.0; ///< Right Ascension [deg]
    double dec_deg = 0.0; ///< Declination [deg]
};

/**
 * @brief Parameters for a corridor search along a polyline path.
 *
 * Useful for searching stars along the projected path of a solar system body.
 */
struct CorridorQuery {
    std::vector<SkyCoord<core::GCRF>> path; ///< Ordered sky positions
    Angle          width;               ///< Half-width of corridor
    double         max_magnitude = 18.0; ///< Faint magnitude limit (Gaia G)
    Parallax       min_parallax;        ///< Optional parallax cut
    size_t         max_results = 10000;  ///< Maximum number of returned stars
};

/**
 * @brief Chebyshev polynomial representing RA or Dec as a function of time.
 *
 * Covers the interval [t_start, t_end] (Julian Date, TDB).
 * The polynomial is evaluated as:
 *   f(t) = sum_k coeffs[k] * T_k(tau),  tau = 2*(t - t_start)/(t_end - t_start) - 1
 */
struct ChebyshevSegment {
    double t_start = 0.0;               ///< Segment start [JD TDB]
    double t_end   = 0.0;               ///< Segment end   [JD TDB]
    std::vector<double> ra_coeffs;      ///< Chebyshev coefficients for RA [deg]
    std::vector<double> dec_coeffs;     ///< Chebyshev coefficients for Dec [deg]

    /**
     * @brief Evaluate polynomial at specific JD.
     * @param jd Julian Date [TDB]
     * @return pair of RA, Dec in degrees.
     */
    [[nodiscard]] std::pair<double, double> evaluate(double jd) const {
        if (ra_coeffs.empty() || dec_coeffs.empty()) return {0.0, 0.0};
        
        // Normalize JD to [-1, 1]
        double tau = (2.0 * jd - (t_end + t_start)) / (t_end - t_start);
        
        auto eval_at = [](const std::vector<double>& cs, double t) {
            if (cs.empty()) return 0.0;
            if (cs.size() == 1) return cs[0];
            
            double d1 = 0.0, d2 = 0.0;
            for (size_t j = cs.size() - 1; j >= 1; --j) {
                double tmp = d1;
                d1 = cs[j] + 2.0 * t * d1 - d2;
                d2 = tmp;
            }
            return cs[0] + t * d1 - d2;
        };

        return {eval_at(ra_coeffs, tau), eval_at(dec_coeffs, tau)};
    }
};

/**
 * @brief Parameters for searching stars along a propagated orbit.
 *
 * The orbit is represented as a set of Chebyshev segments covering the full
 * time interval. Build this with make_orbit_query() in CatalogIntegration.hpp.
 */
struct OrbitQuery {
    time::EpochTDB t_start;             ///< Search window start
    time::EpochTDB t_end;               ///< Search window end
    std::vector<ChebyshevSegment> segments; ///< Chebyshev approximation of the orbit
    Angle          width;               ///< Corridor half-width
    double         max_magnitude = 18.0; ///< Faint magnitude limit (Gaia G)
    double         step_days     = 1.0;  ///< Discretisation step for orbit sampling [days]
};

// ============================================================================
// Progress callback
// ============================================================================

using ProgressCallback = std::function<void(double fraction, const std::string& message)>;

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_TYPES_HPP
