/**
 * @file Observation.hpp
 * @brief Astrometric observation data structures
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Defines structures for various types of observations:
 * - Optical astrometry (RA/Dec, Alt/Az)
 * - Radar observations (range, Doppler)
 * - Photometry (magnitudes)
 * 
 * Supports MPC and ADES formats.
 */

#ifndef ASTDYN_OBSERVATIONS_OBSERVATION_HPP
#define ASTDYN_OBSERVATIONS_OBSERVATION_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/Constants.hpp"
#include <string>
#include <vector>
#include <optional>

namespace astdyn {
namespace observations {

/**
 * @brief Type of observation
 */
enum class ObservationType {
    OPTICAL_RA_DEC,     ///< Optical astrometry in equatorial coordinates
    OPTICAL_AZ_EL,      ///< Optical astrometry in horizontal coordinates
    RADAR_RANGE,        ///< Radar ranging measurement
    RADAR_DOPPLER,      ///< Radar Doppler measurement
    RADAR_RANGE_RATE,   ///< Radar range-rate
    PHOTOMETRY,         ///< Photometric magnitude
    OCCULTATION,        ///< Stellar occultation timing
    SATELLITE_IMAGING,  ///< Space-based imaging
    ROVING_OBSERVER     ///< Roving observer (mobile telescope)
};

/**
 * @brief Observation catalog code
 */
enum class CatalogCode {
    UNKNOWN = 0,
    USNO_A2 = 1,
    USNO_B1 = 2,
    UCAC1 = 3,
    UCAC2 = 4,
    UCAC3 = 5,
    UCAC4 = 6,
    USNO_SA2 = 7,
    TYCHO_1 = 8,
    TYCHO_2 = 9,
    GAIA_DR1 = 10,
    GAIA_DR2 = 11,
    GAIA_DR3 = 12,
    GSC_1_2 = 13,
    GSC_2_3 = 14,
    CMC_14 = 15,
    CMC_15 = 16,
    PPMXL = 17
};

/**
 * @brief Optical astrometric observation
 * 
 * Standard format for optical positional observations.
 * Coordinates can be equatorial (RA/Dec) or horizontal (Az/El).
 */
struct OpticalObservation {
    // Identification
    std::string object_designation;  ///< Object designation (provisional, number, etc.)
    std::string observatory_code;    ///< MPC observatory code (e.g., "568" = Mauna Kea)
    
    // Time
    double mjd_utc;                  ///< Modified Julian Date in UTC
    
    // Astrometry
    double ra;                       ///< Right Ascension [radians]
    double dec;                      ///< Declination [radians]
    double sigma_ra;                 ///< RA uncertainty [radians], σ_α*cos(δ)
    double sigma_dec;                ///< Dec uncertainty [radians]
    
    // Photometry (optional)
    std::optional<double> magnitude; ///< Apparent magnitude
    std::optional<char> mag_band;    ///< Magnitude band (V, R, etc.)
    
    // Metadata
    CatalogCode catalog;             ///< Star catalog used
    std::string note1;               ///< MPC note1 field (discovery, etc.)
    std::string note2;               ///< MPC note2 field (program code)
    
    // Quality flags
    bool is_discovery;               ///< Discovery observation flag
    bool is_offset;                  ///< Offset observation (from nearby reference)
    
    /**
     * @brief Default constructor
     */
    OpticalObservation()
        : mjd_utc(0.0), ra(0.0), dec(0.0), 
          sigma_ra(1e-5), sigma_dec(1e-5),
          catalog(CatalogCode::UNKNOWN),
          is_discovery(false), is_offset(false) {}
};

/**
 * @brief Radar observation
 * 
 * Range and/or Doppler measurements from radar facilities.
 */
struct RadarObservation {
    // Identification
    std::string object_designation;
    std::string transmitter_code;    ///< Transmitting station (e.g., "251" = Arecibo)
    std::string receiver_code;       ///< Receiving station (can differ from TX)
    
    // Time
    double mjd_utc;                  ///< Time of observation (bounce time for range)
    
    // Measurement type
    ObservationType type;            ///< RADAR_RANGE, RADAR_DOPPLER, or RADAR_RANGE_RATE
    
    // Range measurement
    std::optional<double> range;     ///< Range [km] (round-trip time / 2c)
    std::optional<double> sigma_range; ///< Range uncertainty [km]
    
    // Doppler measurement  
    std::optional<double> doppler;   ///< Doppler shift [Hz]
    std::optional<double> sigma_doppler; ///< Doppler uncertainty [Hz]
    
    // Range-rate measurement
    std::optional<double> range_rate; ///< Range rate [km/s]
    std::optional<double> sigma_range_rate; ///< Range rate uncertainty [km/s]
    
    // Frequency
    double frequency_mhz;            ///< Radar frequency [MHz]
    
    /**
     * @brief Default constructor
     */
    RadarObservation()
        : mjd_utc(0.0), type(ObservationType::RADAR_RANGE),
          frequency_mhz(0.0) {}
};

/**
 * @brief Photometric observation
 * 
 * Brightness measurement without positional information.
 */
struct PhotometricObservation {
    // Identification
    std::string object_designation;
    std::string observatory_code;
    
    // Time
    double mjd_utc;
    
    // Photometry
    double magnitude;                ///< Apparent magnitude
    double sigma_magnitude;          ///< Magnitude uncertainty
    char band;                       ///< Photometric band (U, B, V, R, I, etc.)
    
    // Phase information (for lightcurve analysis)
    std::optional<double> phase_angle; ///< Phase angle [degrees]
    std::optional<double> heliocentric_distance; ///< r [AU]
    std::optional<double> geocentric_distance;   ///< Δ [AU]
    
    /**
     * @brief Default constructor
     */
    PhotometricObservation()
        : mjd_utc(0.0), magnitude(0.0), sigma_magnitude(0.1), band('V') {}
};

/**
 * @brief Occultation observation
 * 
 * Timing of stellar occultation events.
 */
struct OccultationObservation {
    // Identification
    std::string object_designation;
    std::string observer_name;
    
    // Location (can be mobile)
    double longitude;                ///< Observer longitude [radians]
    double latitude;                 ///< Observer latitude [radians]
    double altitude;                 ///< Observer altitude [m]
    
    // Timing
    double mjd_utc_start;            ///< Immersion time
    double mjd_utc_end;              ///< Emersion time
    double sigma_time;               ///< Timing uncertainty [seconds]
    
    // Star information
    std::string star_designation;
    double star_ra;                  ///< Star RA [radians]
    double star_dec;                 ///< Star Dec [radians]
    
    /**
     * @brief Default constructor
     */
    OccultationObservation()
        : longitude(0.0), latitude(0.0), altitude(0.0),
          mjd_utc_start(0.0), mjd_utc_end(0.0), sigma_time(0.001),
          star_ra(0.0), star_dec(0.0) {}
};

/**
 * @brief Generic observation container
 * 
 * Unified container that can hold any observation type.
 * Uses std::variant or discriminated union pattern.
 */
struct Observation {
    ObservationType type;
    
    // Storage for different observation types
    std::optional<OpticalObservation> optical;
    std::optional<RadarObservation> radar;
    std::optional<PhotometricObservation> photometric;
    std::optional<OccultationObservation> occultation;
    
    /**
     * @brief Get optical observation
     * @throw std::runtime_error if not optical type
     */
    const OpticalObservation& getOptical() const {
        if (!optical) throw std::runtime_error("Not an optical observation");
        return *optical;
    }
    
    /**
     * @brief Get radar observation
     * @throw std::runtime_error if not radar type
     */
    const RadarObservation& getRadar() const {
        if (!radar) throw std::runtime_error("Not a radar observation");
        return *radar;
    }
    
    /**
     * @brief Get time of observation (works for all types)
     */
    double getMJD() const {
        if (optical) return optical->mjd_utc;
        if (radar) return radar->mjd_utc;
        if (photometric) return photometric->mjd_utc;
        if (occultation) return occultation->mjd_utc_start;
        return 0.0;
    }
    
    /**
     * @brief Get observatory/observer code
     */
    std::string getObservatoryCode() const {
        if (optical) return optical->observatory_code;
        if (radar) return radar->transmitter_code;
        if (photometric) return photometric->observatory_code;
        return "";
    }
};

/**
 * @brief Collection of observations for a single object
 */
struct ObservationSet {
    std::string object_designation;
    std::vector<Observation> observations;
    
    /**
     * @brief Add an observation
     */
    void addObservation(const Observation& obs) {
        observations.push_back(obs);
    }
    
    /**
     * @brief Get number of observations
     */
    size_t size() const {
        return observations.size();
    }
    
    /**
     * @brief Sort observations by time
     */
    void sortByTime();
    
    /**
     * @brief Get time span [days]
     */
    double getTimeSpan() const;
    
    /**
     * @brief Get optical observations only
     */
    std::vector<OpticalObservation> getOpticalObservations() const;
    
    /**
     * @brief Get radar observations only
     */
    std::vector<RadarObservation> getRadarObservations() const;
};

/**
 * @brief Convert MPC catalog code character to enum
 */
CatalogCode parseCatalogCode(char code);

/**
 * @brief Convert catalog enum to string
 */
std::string catalogCodeToString(CatalogCode code);

/**
 * @brief Convert observation type to string
 */
std::string observationTypeToString(ObservationType type);

} // namespace observations
} // namespace astdyn

#endif // ASTDYN_OBSERVATIONS_OBSERVATION_HPP
