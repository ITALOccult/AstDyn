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

namespace astdyn::catalog {

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
    double ra_deg  = 0.0;               ///< Right Ascension [deg]
    double dec_deg = 0.0;               ///< Declination [deg]
    double parallax_mas   = 0.0;        ///< Parallax [mas]  (>0 = real star)
    double parallax_error_mas = 0.0;    ///< Parallax uncertainty [mas]
    double pmra_mas_yr    = 0.0;        ///< Proper motion in RA * cos(dec) [mas/yr]
    double pmdec_mas_yr   = 0.0;        ///< Proper motion in Dec [mas/yr]
    double pmra_error_mas_yr  = 0.0;    ///< PM RA uncertainty [mas/yr]
    double pmdec_error_mas_yr = 0.0;    ///< PM Dec uncertainty [mas/yr]

    // --- Photometry ---
    double g_mag  = 99.0;               ///< Gaia G-band magnitude
    double bp_mag = 99.0;               ///< Gaia BP-band magnitude
    double rp_mag = 99.0;               ///< Gaia RP-band magnitude
    double bp_rp  = 0.0;                ///< BP-RP color index

    // --- Quality ---
    double ruwe = 0.0;                  ///< Renormalised unit weight error
    double astrometric_excess_noise = 0.0;

    /// Distance in AU inferred from parallax (returns 0 if parallax <= 0)
    [[nodiscard]] double distance_au() const noexcept {
        if (parallax_mas <= 0.0) return 0.0;
        return 1000.0 / parallax_mas * constants::AU_PER_PARSEC;
    }

    [[nodiscard]] bool has_parallax() const noexcept { return parallax_mas > 0.0; }
    [[nodiscard]] bool has_proper_motion() const noexcept {
        return pmra_mas_yr != 0.0 || pmdec_mas_yr != 0.0;
    }
};

// ============================================================================
// Query parameter structs
// ============================================================================

/**
 * @brief Parameters for a circular cone search.
 */
struct ConeQuery {
    double ra_deg    = 0.0;      ///< Center RA [deg]
    double dec_deg   = 0.0;      ///< Center Dec [deg]
    double radius_arcsec = 60.0; ///< Search radius [arcsec]
    double max_magnitude = 18.0; ///< Faint limit (Gaia G)
    double min_parallax_mas = 0.0; ///< Optional parallax cut [mas]
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
    std::vector<SkyPoint> path;         ///< Ordered sky positions defining the corridor axis
    double width_arcsec  = 120.0;       ///< Half-width of corridor [arcsec]
    double max_magnitude = 18.0;        ///< Faint magnitude limit (Gaia G)
    double min_parallax_mas = 0.0;      ///< Optional parallax cut [mas]
    size_t max_results = 10000;         ///< Maximum number of returned stars
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
};

/**
 * @brief Parameters for searching stars along a propagated orbit.
 *
 * The orbit is represented as a set of Chebyshev segments covering the full
 * time interval. Build this with make_orbit_query() in CatalogIntegration.hpp.
 */
struct OrbitQuery {
    double t_start_jd  = 0.0;              ///< Search window start [JD TDB]
    double t_end_jd    = 0.0;              ///< Search window end   [JD TDB]
    std::vector<ChebyshevSegment> segments; ///< Chebyshev approximation of the orbit
    double width_arcsec  = 120.0;           ///< Corridor half-width [arcsec]
    double max_magnitude = 18.0;            ///< Faint magnitude limit (Gaia G)
    double step_days     = 1.0;             ///< Discretisation step for orbit sampling [days]
};

// ============================================================================
// Progress callback
// ============================================================================

using ProgressCallback = std::function<void(double fraction, const std::string& message)>;

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_TYPES_HPP
