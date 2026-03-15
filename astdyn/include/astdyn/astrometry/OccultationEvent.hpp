/**
 * @file OccultationEvent.hpp
 * @brief Data structure representing a full occultation event matching external XML formats.
 */

#ifndef ASTDYN_ASTROMETRY_OCCULTATION_EVENT_HPP
#define ASTDYN_ASTROMETRY_OCCULTATION_EVENT_HPP

#include "astdyn/core/physics_state.hpp"
#include "astdyn/astrometry/sky_types.hpp"
#include "astdyn/time/epoch.hpp"
#include <vector>
#include <string>

namespace astdyn::astrometry {

/**
 * @brief Represents a single occultation event with all associated metadata.
 * This structure maps directly to the XML format used by external calculators (like Occult4).
 */
struct OccultationEvent {
    // --- Elements (JPL Integration Data) ---
    std::string elements_source; // e.g., "JPL +INTG:06-Nov-2025"
    std::vector<double> elements_data; // The numeric parameters in the <Elements> tag

    // --- Earth (Shadow Geometry) ---
    double longitude_deg;
    double latitude_deg;
    double alt_or_other; // Third parameter in <Earth>
    double max_duration_sec;
    bool is_daylight;

    // --- Star (Stellar Data) ---
    std::string star_catalog_id;
    double ra_cat_h;      // RA at catalog epoch (hours)
    double dec_cat_deg;   // Dec at catalog epoch (degrees)
    double mag_v, mag_r, mag_k; // Magnitudes
    double pm_ra_as_yr = 0.0;   // Proper motion in RA (arcsec/yr)
    double pm_dec_as_yr = 0.0;  // Proper motion in Dec (arcsec/yr)
    double parallax_as = 0.0;   // Parallax (arcsec)
    double ra_event_h;          // RA at event epoch (hours)
    double dec_event_deg;       // Dec at event epoch (degrees)
    std::vector<double> star_extra_data; // Remaining values in <Star>

    // --- Object (Asteroid Properties) ---
    int object_number;
    std::string object_name;
    double h_mag;
    double diameter_km;
    double apparent_rate_arcsec_hr;
    double g_param = 0.15;      // Slope parameter
    std::string object_type;    // e.g., "Asteroid"
    std::vector<double> object_extra_data; // Remaining values in <Object>

    // --- Orbit (Keplerian Elements) ---
    double orbit_type = 0; // First value (often 0)
    double ra_node_approx; // Second value (sometimes Omega)
    int epoch_year, epoch_month, epoch_day;
    double mean_anomaly_deg;
    double node_deg;
    double inclination_deg;
    double eccentricity;
    double semi_major_axis_au;
    std::vector<double> orbit_extra_data; // Remaining values in <Orbit>

    // --- Errors (Uncertainty) ---
    double combined_error_arcsec;
    double star_error_ra_as = 0.0; // Uncertainty in star RA (arcsec)
    double star_error_dec_as = 0.0; // Uncertainty in star Dec (arcsec)
    double object_error_ra_as = 0.0; // Uncertainty in asteroid RA (arcsec)
    double object_error_dec_as = 0.0; // Uncertainty in asteroid Dec (arcsec)
    std::string uncertainty_method; // e.g., "Star+JplPeakEphemUncert"
    std::vector<double> error_extra_data; // Remaining numeric values in <Errors>

    // --- ID ---
    std::string event_id;
    double mjd;
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_OCCULTATION_EVENT_HPP
