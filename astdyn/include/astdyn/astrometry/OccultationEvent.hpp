/**
 * @file OccultationEvent.hpp
 * @brief One occultation event in Occult4's "occelmnt" XML format.
 *
 * Field names and order follow Occult4's published specification. The previous
 * version of this struct was written by inferring the layout from sample files,
 * and several fields were mis-assigned: the sub-solar point was stored as an
 * altitude and a duration, the mean anomaly and the argument of perihelion were
 * swapped, and the proper motions occupied the slots reserved for the stellar
 * diameter and the double-star code. Keeping the names aligned with the
 * specification is what makes such drift visible.
 *
 *   <Elements> source, duration, year, month, day, UT at closest,
 *              x, y, dX, dY, d2X, d2Y, d3X, d3Y
 *   <Earth>    SubstellarLong, SubstellarLat, SubsolarLong, SubsolarLat, JWST
 *   <Star>     Identifier, RA, Dec, Mb, Mv, Mr, dia, double star code, K2 flag,
 *              Apparent RA, Apparent Dec, MdropV, MdropR,
 *              MagDropsAdjusted, BrightNearbyCount, TotalNearbyCount
 *   <Object>   number, name, magnitude, diameter, distance, #rings, #moons,
 *              dRA, dDec, Taxonomy, DiameterUncertainty,
 *              PlanetMoonInPlanetShadow, MagV, MagR
 *   <Orbit>    equinox, MA, yr, month, day, peri, node, i, e, a, q,
 *              H0, Coeff LogR, G
 *   <Errors>   in Pathwidths, major axis, minor axis, PA, 1 sigma,
 *              Error basis, Reliability, Duplicate, Non-GAIA pm, pm from UCAC4
 *   <ID>       ID, MJD of the prediction calculation
 */

#ifndef ASTDYN_ASTROMETRY_OCCULTATION_EVENT_HPP
#define ASTDYN_ASTROMETRY_OCCULTATION_EVENT_HPP

#include <string>

namespace astdyn::astrometry {

/// A single occultation event, mapped one-to-one onto the occelmnt format.
struct OccultationEvent {
    // ---- <Elements> : Besselian elements ----------------------------------
    std::string elements_source;      ///< Orbit source + prediction date.
    double duration_s = 0.0;          ///< Maximum duration [s].
    int    year = 0, month = 0, day = 0;
    double ut_closest_h = 0.0;        ///< UT of closest approach to the geocentre [h].
    double x = 0.0, y = 0.0;          ///< Shadow axis at closest approach [Earth radii].
    double dx = 0.0, dy = 0.0;        ///< Rate of change [Earth radii/hour].
    double d2x = 0.0, d2y = 0.0;      ///< 2nd-order rate.
    double d3x = 0.0, d3y = 0.0;      ///< 3rd-order rate.

    // ---- <Earth> ----------------------------------------------------------
    /// Sub-star point [deg]. The latitude is "of date", i.e. geocentric, and so
    /// equals the star's apparent declination.
    double substellar_lon_deg = 0.0;
    double substellar_lat_deg = 0.0;
    /// Sub-solar point [deg]; this is what lets the terminator be drawn.
    double subsolar_lon_deg = 0.0;
    double subsolar_lat_deg = 0.0;
    bool   jwst = false;              ///< Prediction is for the JWST location.

    // ---- <Star> -----------------------------------------------------------
    std::string star_id;
    /// BCRS: ICRS frame, epoch of the event (proper motion applied), no parallax.
    double star_ra_h = 0.0;
    double star_dec_deg = 0.0;
    double mag_b = 99.0, mag_v = 99.0, mag_r = 99.0;
    double star_diameter_mas = 0.0;
    int    double_star_code = 0;      ///< 0 none, 1 WDS, 2 other, 4 variable; cumulative.
    std::string k2_flag;              ///< "K" for a Kepler-2 target.
    double star_app_ra_h = 0.0;       ///< Apparent position, of date.
    double star_app_dec_deg = 0.0;
    double mag_drop_v = 0.0, mag_drop_r = 0.0;
    int    mag_drops_adjusted = 0;
    int    bright_nearby_count = -1;  ///< -1 when no check was performed.
    int    total_nearby_count = -1;

    // ---- <Object> ---------------------------------------------------------
    std::string object_number;        ///< Asteroid number, or PxMyy for a moon.
    std::string object_name;
    double object_mag = 0.0;          ///< Apparent magnitude.
    double diameter_km = 0.0;         ///< Augmented by the star's diameter.
    double distance_au = 0.0;         ///< Geocentric distance.
    int    n_rings = 0, n_moons = 0;
    double d_ra_s_hr = 0.0;           ///< Hourly change in RA [s of time/hour].
    double d_dec_as_hr = 0.0;         ///< Hourly change in Dec [arcsec/hour].
    std::string taxonomy;
    double diameter_uncertainty_km = 0.0;
    int    moon_in_planet_shadow = 0;
    double mag_v_asteroid = 0.0, mag_r_asteroid = 0.0;

    // ---- <Orbit> : low-precision, for plotting only -----------------------
    double equinox = 0.0;
    double mean_anomaly_deg = 0.0;
    int    epoch_year = 0, epoch_month = 0, epoch_day = 0;
    double peri_deg = 0.0;
    double node_deg = 0.0;
    double inclination_deg = 0.0;
    double eccentricity = 0.0;
    double semi_major_axis_au = 0.0;
    double perihelion_au = 0.0;
    double h0 = 0.0;
    double coeff_log_r = 0.0;
    double g_param = 0.15;

    // ---- <Errors> ---------------------------------------------------------
    double err_path_widths = 0.0;     ///< Path location error, in path widths.
    double err_major_as = 0.0, err_minor_as = 0.0, err_pa_deg = 0.0;
    double err_1sigma_as = 0.0;
    std::string error_basis;          ///< "Known errors", "Star+PEU", ...
    double reliability = -1.0;        ///< Usually the RUWE; -1 when not set.
    int    duplicate_source = -1;
    int    non_gaia_pm = -1;
    int    pm_added_from_ucac4 = -1;

    // ---- <ID> -------------------------------------------------------------
    std::string event_id;             ///< yyyymmdd_xxxxxx
    /// MJD of the prediction calculation -- NOT the epoch of the event, which
    /// lives in <Elements>.
    double prediction_mjd = 0.0;
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_OCCULTATION_EVENT_HPP
