#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>
#include <chrono>

namespace ioc {
namespace gaia {

/**
 * Gaia data release versions
 */
enum class GaiaRelease {
    DR2,      ///< Gaia Data Release 2 (2018)
    EDR3,     ///< Gaia Early Data Release 3 (2020)
    DR3       ///< Gaia Data Release 3 (2022)
};

/**
 * Gaia star data structure
 * Contains all essential fields from Gaia DR3
 */
struct GaiaStar {
    int64_t source_id;           ///< Unique Gaia source identifier
    
    // Astrometric data (ICRS, Epoch 2016.0)
    double ra;                   ///< Right Ascension [degrees]
    double dec;                  ///< Declination [degrees]
    double parallax;             ///< Parallax [mas]
    double parallax_error;       ///< Parallax standard error [mas]
    
    // Proper motion
    double pmra;                 ///< Proper motion in RA * cos(dec) [mas/yr]
    double pmdec;                ///< Proper motion in Dec [mas/yr]
    double pmra_error;           ///< PM RA standard error [mas/yr]
    double pmdec_error;          ///< PM Dec standard error [mas/yr]
    
    // Photometry
    double phot_g_mean_mag;      ///< G-band mean magnitude [mag]
    double phot_bp_mean_mag;     ///< BP-band mean magnitude [mag]
    double phot_rp_mean_mag;     ///< RP-band mean magnitude [mag]
    double bp_rp;                ///< BP-RP color [mag]
    
    // Quality flags
    double astrometric_excess_noise;
    double astrometric_chi2_al;
    int visibility_periods_used;
    double ruwe;                 ///< Renormalized Unit Weight Error
    
    // Cross-match identifiers
    std::string sao_designation;     ///< SAO catalog designation
    std::string tycho2_designation;  ///< Tycho-2 catalog designation  
    std::string hd_designation;      ///< Henry Draper catalog designation
    std::string hip_designation;     ///< Hipparcos catalog designation
    std::string common_name;         ///< Common star name (if available)
    
    // Constructor
    GaiaStar();
    
    // Utility methods
    bool isValid() const;
    double getBpRpColor() const;
    std::string getDesignation() const;
    std::string getAllDesignations() const;  ///< Get all available designations
    
    /**
     * @brief Propaga posizione della stella usando il moto proprio
     * @param mjd Epoca target [MJD]
     * @return Stella con coordinate aggiornate all'epoca target
     */
    GaiaStar propagateToEpoch(double mjd) const;
};

/**
 * Equatorial coordinates (J2000.0)
 */
struct EquatorialCoordinates {
    double ra;                   ///< Right Ascension [degrees]
    double dec;                  ///< Declination [degrees]
    
    EquatorialCoordinates() : ra(0.0), dec(0.0) {}
    EquatorialCoordinates(double ra_, double dec_) : ra(ra_), dec(dec_) {}
};

/**
 * Julian Date representation
 */
struct JulianDate {
    double jd;                   ///< Julian Date
    
    JulianDate() : jd(0.0) {}
    explicit JulianDate(double jd_) : jd(jd_) {}
    
    // Convert to/from calendar date
    static JulianDate fromCalendar(int year, int month, int day, double ut);
    void toCalendar(int& year, int& month, int& day, double& ut) const;
    
    // Epoch conversions
    static JulianDate J2000() { return JulianDate(2451545.0); }
    static JulianDate J2016() { return JulianDate(2457389.0); }  // Gaia DR3 epoch
    
    double toYears() const { return (jd - 2451545.0) / 365.25 + 2000.0; }
};

/**
 * Query parameters for Gaia catalog searches
 */
struct QueryParams {
    double ra_center;            ///< Center RA [degrees]
    double dec_center;           ///< Center Dec [degrees]
    double radius;               ///< Search radius [degrees]
    double max_magnitude;        ///< Maximum G magnitude
    double min_parallax;         ///< Minimum parallax [mas], -1 = no limit
    
    QueryParams() 
        : ra_center(0.0), dec_center(0.0), radius(1.0), 
          max_magnitude(20.0), min_parallax(-1.0) {}
};

/**
 * Rectangular search region
 */
struct BoxRegion {
    double ra_min;               ///< Minimum RA [degrees]
    double ra_max;               ///< Maximum RA [degrees]
    double dec_min;              ///< Minimum Dec [degrees]
    double dec_max;              ///< Maximum Dec [degrees]
    
    BoxRegion() : ra_min(0.0), ra_max(0.0), dec_min(0.0), dec_max(0.0) {}
    BoxRegion(double ra_min_, double ra_max_, double dec_min_, double dec_max_)
        : ra_min(ra_min_), ra_max(ra_max_), dec_min(dec_min_), dec_max(dec_max_) {}
};

/**
 * A point on the celestial sphere (for corridor/polyline queries)
 */
struct CelestialPoint {
    double ra;                   ///< Right Ascension [degrees]
    double dec;                  ///< Declination [degrees]
    
    CelestialPoint() : ra(0.0), dec(0.0) {}
    CelestialPoint(double ra_, double dec_) : ra(ra_), dec(dec_) {}
};

/**
 * Parameters for corridor/polyline query
 * Searches for stars within a specified width along a path defined by multiple points
 */
struct CorridorQueryParams {
    std::vector<CelestialPoint> path;  ///< Path defined by ordered points (at least 2)
    double width;                       ///< Corridor half-width [degrees] (total width = 2 * width)
    double max_magnitude;               ///< Maximum G magnitude
    double min_parallax;                ///< Minimum parallax [mas], -1 = no limit
    size_t max_results;                 ///< Maximum results (0 = no limit)
    
    CorridorQueryParams() 
        : width(0.5), max_magnitude(20.0), min_parallax(-1.0), max_results(0) {}
    
    /**
     * Parse corridor query from JSON string
     * Format:
     * {
     *   "path": [
     *     {"ra": 0.0, "dec": 0.0},
     *     {"ra": 10.0, "dec": 5.0},
     *     {"ra": 20.0, "dec": 10.0}
     *   ],
     *   "width": 0.5,
     *   "max_magnitude": 18.0,
     *   "min_parallax": -1,
     *   "max_results": 0
     * }
     */
    static CorridorQueryParams fromJSON(const std::string& json);
    
    /**
     * Serialize to JSON string
     */
    std::string toJSON() const;
    
    /**
     * Validate parameters
     */
    bool isValid() const {
        return path.size() >= 2 && width > 0 && max_magnitude > 0;
    }
    
    /**
     * Get total path length in degrees (approximate great circle)
     */
    double getPathLength() const;
};

/**
 * Chebyshev Polynomial coefficients for orbit interpolation
 * Represents an orbit segment in RA/Dec over a specific time range
 */
struct ChebyshevPolynomial {
    double t_start;             ///< Start time of validity (JD or arbitrary units matching query range)
    double t_end;               ///< End time of validity
    std::vector<double> coeffs_ra;  ///< Chebyshev coefficients for Right Ascension
    std::vector<double> coeffs_dec; ///< Chebyshev coefficients for Declination
    
    ChebyshevPolynomial() : t_start(0), t_end(0) {}
};

/**
 * Parameters for orbit query
 * Searches for stars along an orbit interpolated by Chebyshev polynomials
 */
struct OrbitQueryParams {
    double t_start;                 ///< Query start time
    double t_end;                   ///< Query end time
    std::vector<ChebyshevPolynomial> polynomials; ///< List of polynomials covering the interval
    double width;                   ///< Search width (radius) around the orbit [degrees]
    double max_magnitude;           ///< Maximum G magnitude
    double step_size;               ///< Step size for discretizing the orbit [same units as time] (0 = auto)
    
    OrbitQueryParams() 
        : t_start(0), t_end(0), width(0.1), max_magnitude(20.0), step_size(0) {}
        
    bool isValid() const {
        return t_end > t_start && !polynomials.empty() && width > 0;
    }
};

/**
 * Cache statistics
 */
struct CacheStats {
    size_t total_tiles;          ///< Total HEALPix tiles
    size_t cached_tiles;         ///< Number of cached tiles
    size_t total_stars;          ///< Total cached stars
    size_t disk_size_bytes;      ///< Cache size on disk
    std::chrono::system_clock::time_point last_update;
    
    CacheStats() : total_tiles(0), cached_tiles(0), total_stars(0), disk_size_bytes(0) {}
    
    double getCoveragePercent() const {
        return total_tiles > 0 ? (100.0 * cached_tiles / total_tiles) : 0.0;
    }
};

/**
 * Progress callback for async operations
 */
using ProgressCallback = std::function<void(int percent, const std::string& message)>;

/**
 * Error types for exception handling
 */
enum class ErrorCode {
    SUCCESS = 0,
    NETWORK_ERROR,
    TIMEOUT,
    PARSE_ERROR,
    INVALID_PARAMS,
    CACHE_ERROR,
    RATE_LIMIT_EXCEEDED,
    SERVICE_UNAVAILABLE
};

/**
 * Exception class for Gaia library errors
 */
class GaiaException : public std::exception {
public:
    GaiaException(ErrorCode code, const std::string& message)
        : code_(code), message_(message) {}
    
    const char* what() const noexcept override { return message_.c_str(); }
    ErrorCode code() const noexcept { return code_; }
    
private:
    ErrorCode code_;
    std::string message_;
};

// =============================================================================
// Utility Functions
// =============================================================================

/**
 * Convert string to GaiaRelease enum
 */
GaiaRelease stringToRelease(const std::string& str);

/**
 * Convert GaiaRelease enum to string
 */
std::string releaseToString(GaiaRelease release);

/**
 * Calculate angular distance between two coordinates (haversine formula)
 * @return Distance in degrees
 */
double angularDistance(const EquatorialCoordinates& c1, 
                      const EquatorialCoordinates& c2);

/**
 * Validate coordinate ranges
 */
bool isValidCoordinate(const EquatorialCoordinates& coord);

/**
 * Validate query parameters
 */
bool isValidQueryParams(const QueryParams& params);

/**
 * Validate box region
 */
bool isValidBoxRegion(const BoxRegion& box);

} // namespace gaia
} // namespace ioc
