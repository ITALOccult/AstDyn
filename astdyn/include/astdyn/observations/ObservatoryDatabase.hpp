/**
 * @file ObservatoryDatabase.hpp
 * @brief Observatory location database
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Manages observatory locations using MPC observatory codes.
 * Stores geodetic coordinates (longitude, latitude, altitude) and
 * computes geocentric Cartesian coordinates for observations.
 * 
 * References:
 * - MPC Observatory Codes: https://minorplanetcenter.net/iau/lists/ObsCodes.html
 * - IAU SOFA: Standards of Fundamental Astronomy
 */

#ifndef ASTDYN_OBSERVATIONS_OBSERVATORYDATABASE_HPP
#define ASTDYN_OBSERVATIONS_OBSERVATORYDATABASE_HPP

#include "astdyn/core/Constants.hpp"
#include "astdyn/time/epoch.hpp"
#include "src/types/vectors.hpp"
#include "src/core/frame_tags.hpp"
#include "src/core/units.hpp"
#include <string>
#include <map>
#include <optional>
#include <vector>
#include <Eigen/Dense>

namespace astdyn {
namespace observations {

/**
 * @brief Observatory location data
 */
struct Observatory {
    std::string code;          ///< MPC 3-character code (e.g., "568")
    std::string name;          ///< Observatory name
    
    // Geodetic coordinates
    double longitude;          ///< East longitude [radians]
    double latitude;           ///< Latitude [radians] (geodetic)
    double altitude;           ///< Altitude above sea level [meters]
    
    // Geocentric coordinates (cached)
    types::Vector3<core::ITRF, core::Meter> position;  ///< Geocentric position [m] (ITRF)
    
    // Parallax constants (for speed)
    double rho_cos_phi;        ///< ρ cos φ' (geocentric)
    double rho_sin_phi;        ///< ρ sin φ' (geocentric)
    
    /**
     * @brief Default constructor
     */
    Observatory()
        : longitude(0.0), latitude(0.0), altitude(0.0),
          position(0.0, 0.0, 0.0),
          rho_cos_phi(0.0), rho_sin_phi(0.0) {}
    
    /**
     * @brief Compute geocentric Cartesian coordinates
     * 
     * Converts geodetic (lon, lat, alt) to geocentric ITRF coordinates.
     * Uses WGS84 ellipsoid parameters.
     */
    void computeGeocentricPosition();
    
    /**
     * @brief Get position vector at specific time (includes Earth rotation)
     * @param t Time
     * @return Geocentric position [m] in J2000 frame (GCRF)
     */
    types::Vector3<core::GCRF, core::Meter> getPositionGCRF(time::EpochUTC t) const;
};

/**
 * @brief Database of observatories
 * 
 * Singleton pattern for global access to observatory data.
 */
class ObservatoryDatabase {
public:
    /**
     * @brief Get singleton instance
     */
    static ObservatoryDatabase& getInstance();
    
    /**
     * @brief Load observatories from MPC ObsCodes.txt file
     * @param filepath Path to ObsCodes.txt
     * @return Number of observatories loaded
     */
    size_t loadFromMPCFile(const std::string& filepath);
    
    /**
     * @brief Load observatories from embedded default data
     * 
     * Loads a subset of most common observatories.
     */
    void loadDefaultObservatories();
    
    /**
     * @brief Get observatory by code
     * @param code MPC 3-character code
     * @return Observatory data, or nullopt if not found
     */
    std::optional<Observatory> getObservatory(const std::string& code) const;
    
    /**
     * @brief Check if observatory exists
     */
    bool hasObservatory(const std::string& code) const;
    
    /**
     * @brief Add or update observatory
     */
    void addObservatory(const Observatory& obs);
    
    /**
     * @brief Get all observatory codes
     */
    std::vector<std::string> getAllCodes() const;
    
    /**
     * @brief Get number of observatories in database
     */
    size_t size() const { return observatories_.size(); }
    
    /**
     * @brief Clear all observatories
     */
    void clear() { observatories_.clear(); }

private:
    // Singleton - private constructor
    ObservatoryDatabase() = default;
    ObservatoryDatabase(const ObservatoryDatabase&) = delete;
    ObservatoryDatabase& operator=(const ObservatoryDatabase&) = delete;
    
    std::map<std::string, Observatory> observatories_;
};

/**
 * @brief Convert geodetic to geocentric coordinates
 * 
 * @param lon_geodetic East longitude [radians]
 * @param lat_geodetic Geodetic latitude [radians]
 * @param alt_meters Altitude above WGS84 ellipsoid [meters]
 * @return Geocentric Cartesian position [m] in ITRF frame
 */
types::Vector3<core::ITRF, core::Meter> geodeticToGeocentric(
    double lon_geodetic,
    double lat_geodetic,
    double alt_meters
);

/**
 * @brief Compute parallax constants for observatory
 * 
 * @param lat_geodetic Geodetic latitude [radians]
 * @param alt_meters Altitude [meters]
 * @param rho_cos_phi Output: ρ cos φ'
 * @param rho_sin_phi Output: ρ sin φ'
 */
void computeParallaxConstants(
    double lat_geodetic,
    double alt_meters,
    double& rho_cos_phi,
    double& rho_sin_phi
);

} // namespace observations
} // namespace astdyn

#endif // ASTDYN_OBSERVATIONS_OBSERVATORYDATABASE_HPP
