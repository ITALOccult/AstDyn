/**
 * @file ObservatoryDatabase.cpp
 * @brief Implementation of observatory database
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#include "astdyn/observations/ObservatoryDatabase.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

namespace astdyn {
namespace observations {

// WGS84 ellipsoid parameters
static constexpr double WGS84_A = 6378.137;      // Semi-major axis [km]
static constexpr double WGS84_F = 1.0 / 298.257223563; // Flattening
static constexpr double WGS84_B = WGS84_A * (1.0 - WGS84_F); // Semi-minor axis

// ============================================================================
// Observatory Methods
// ============================================================================

void Observatory::computeGeocentricPosition() {
    position = geodeticToGeocentric(longitude, latitude, altitude);
    computeParallaxConstants(latitude, altitude, rho_cos_phi, rho_sin_phi);
}

Eigen::Vector3d Observatory::getPositionJ2000(double mjd_utc) const {
    // Convert ITRF position to J2000 using Earth rotation
    // This accounts for Earth's rotation during observation
    return coordinates::ReferenceFrame::transform_position(
        position,
        coordinates::FrameType::ITRF,
        coordinates::FrameType::J2000,
        mjd_utc
    );
}

// ============================================================================
// ObservatoryDatabase Methods
// ============================================================================

ObservatoryDatabase& ObservatoryDatabase::getInstance() {
    static ObservatoryDatabase instance;
    return instance;
}

size_t ObservatoryDatabase::loadFromMPCFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open observatory file: " + filepath);
    }
    
    size_t count = 0;
    std::string line;
    
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        // MPC ObsCodes.txt format:
        // Columns 1-3: Code
        // Columns 5-13: Longitude (ddd.ddddd)
        // Columns 14-21: Parallax cos (0.dddddd)
        // Columns 22-30: Parallax sin (±0.dddddd)
        // Columns 31-80: Name
        
        if (line.length() < 30) continue;
        
        Observatory obs;
        
        try {
            // Parse code
            obs.code = line.substr(0, 3);
            
            // Parse longitude (degrees)
            std::string lon_str = line.substr(4, 9);
            double lon_deg = std::stod(lon_str);
            obs.longitude = lon_deg * constants::DEG_TO_RAD;
            
            // Parse parallax constants
            std::string rho_cos_str = line.substr(13, 8);
            std::string rho_sin_str = line.substr(21, 9);
            obs.rho_cos_phi = std::stod(rho_cos_str);
            obs.rho_sin_phi = std::stod(rho_sin_str);
            
            // Parse name (if present)
            if (line.length() > 30) {
                obs.name = line.substr(30);
                // Trim whitespace
                obs.name.erase(0, obs.name.find_first_not_of(" \t"));
                obs.name.erase(obs.name.find_last_not_of(" \t") + 1);
            }
            
            // Compute latitude and altitude from parallax constants
            // ρ sin φ' = (1 - f) sin φ + (h/a)
            // ρ cos φ' = cos φ + (h/a) cos φ (approximately)
            // For simplicity, assume h ≈ 0 for most observatories
            double rho = std::sqrt(obs.rho_cos_phi * obs.rho_cos_phi +
                                  obs.rho_sin_phi * obs.rho_sin_phi);
            obs.latitude = std::atan2(obs.rho_sin_phi, obs.rho_cos_phi);
            obs.altitude = (rho - 1.0) * WGS84_A * 1000.0; // Convert to meters
            
            // Compute geocentric position
            obs.computeGeocentricPosition();
            
            addObservatory(obs);
            count++;
            
        } catch (const std::exception& e) {
            // Skip malformed lines
            continue;
        }
    }
    
    return count;
}

void ObservatoryDatabase::loadDefaultObservatories() {
    // Load a subset of most common observatories
    // Data from MPC ObsCodes.txt
    
    struct DefaultObs {
        const char* code;
        const char* name;
        double lon_deg;    // East longitude
        double rho_cos;
        double rho_sin;
    };
    
    static const DefaultObs defaults[] = {
        {"500", "Geocentric", 0.0, 0.0, 0.0},
        {"568", "Mauna Kea", 204.52833, 0.56497, 0.82400},
        {"G96", "Mt. Lemmon Survey", 249.20892, 0.77071, 0.63098},
        {"F51", "Pan-STARRS 1, Haleakala", 203.74417, 0.56916, 0.81996},
        {"F52", "Pan-STARRS 2, Haleakala", 203.74500, 0.56916, 0.81996},
        {"691", "Steward Observatory, Kitt Peak-Spacewatch", 248.39750, 0.76377, 0.63852},
        {"704", "Lincoln Lab's ETS, New Mexico", 253.39583, 0.73907, 0.66984},
        {"703", "Catalina Sky Survey", 249.20611, 0.77071, 0.63098},
        {"C51", "WISE", 0.0, 0.0, 0.0}, // Space-based
        {"Y28", "Atlas-MLO, Mauna Loa", 204.41600, 0.56630, 0.82236},
        {"T09", "Atlas-HKO, Haleakala", 203.74583, 0.56916, 0.81996},
        {"W84", "Tenerife-ATLAS", 343.44167, 0.86090, 0.49633},
        {"J04", "Cape Town", 18.47750, 0.84537, -0.52923},
        {"474", "Siding Spring Observatory", 149.06611, 0.77014, -0.62868},
        {"A77", "Siding Spring-SSO", 149.07083, 0.77014, -0.62868},
        {"413", "Siding Spring Observatory-DSS", 149.10000, 0.77014, -0.62868}
    };
    
    for (const auto& def : defaults) {
        Observatory obs;
        obs.code = def.code;
        obs.name = def.name;
        obs.longitude = def.lon_deg * constants::DEG_TO_RAD;
        obs.rho_cos_phi = def.rho_cos;
        obs.rho_sin_phi = def.rho_sin;
        
        // Compute latitude from parallax
        if (obs.rho_cos_phi != 0.0 || obs.rho_sin_phi != 0.0) {
            obs.latitude = std::atan2(obs.rho_sin_phi, obs.rho_cos_phi);
            double rho = std::sqrt(obs.rho_cos_phi * obs.rho_cos_phi +
                                  obs.rho_sin_phi * obs.rho_sin_phi);
            obs.altitude = (rho - 1.0) * WGS84_A * 1000.0;
        } else {
            // Geocenter or space-based
            obs.latitude = 0.0;
            obs.altitude = 0.0;
        }
        
        obs.computeGeocentricPosition();
        addObservatory(obs);
    }
}

std::optional<Observatory> ObservatoryDatabase::getObservatory(const std::string& code) const {
    auto it = observatories_.find(code);
    if (it != observatories_.end()) {
        return it->second;
    }
    return std::nullopt;
}

bool ObservatoryDatabase::hasObservatory(const std::string& code) const {
    return observatories_.find(code) != observatories_.end();
}

void ObservatoryDatabase::addObservatory(const Observatory& obs) {
    observatories_[obs.code] = obs;
}

std::vector<std::string> ObservatoryDatabase::getAllCodes() const {
    std::vector<std::string> codes;
    codes.reserve(observatories_.size());
    for (const auto& pair : observatories_) {
        codes.push_back(pair.first);
    }
    return codes;
}

// ============================================================================
// Coordinate Conversion Functions
// ============================================================================

Eigen::Vector3d geodeticToGeocentric(
    double lon_geodetic,
    double lat_geodetic,
    double alt_meters)
{
    // WGS84 ellipsoid
    double a = WGS84_A; // Semi-major axis [km]
    double f = WGS84_F; // Flattening
    double e2 = 2.0 * f - f * f; // Eccentricity squared
    
    // Radius of curvature in prime vertical
    double sin_lat = std::sin(lat_geodetic);
    double N = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
    
    // Geocentric Cartesian coordinates [km]
    double h_km = alt_meters / 1000.0;
    double cos_lat = std::cos(lat_geodetic);
    double cos_lon = std::cos(lon_geodetic);
    double sin_lon = std::sin(lon_geodetic);
    
    double x = (N + h_km) * cos_lat * cos_lon;
    double y = (N + h_km) * cos_lat * sin_lon;
    double z = (N * (1.0 - e2) + h_km) * sin_lat;
    
    return Eigen::Vector3d(x, y, z);
}

void computeParallaxConstants(
    double lat_geodetic,
    double alt_meters,
    double& rho_cos_phi,
    double& rho_sin_phi)
{
    // WGS84 parameters
    double f = WGS84_F;
    double e2 = 2.0 * f - f * f;
    
    double sin_lat = std::sin(lat_geodetic);
    double cos_lat = std::cos(lat_geodetic);
    
    // Radius of curvature
    double N = WGS84_A / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
    
    // Height above ellipsoid in units of Earth radii
    double h_earth_radii = (alt_meters / 1000.0) / WGS84_A;
    
    // Geocentric latitude (phi prime)
    double tan_phi_prime = (1.0 - f) * (1.0 - f) * std::tan(lat_geodetic);
    double phi_prime = std::atan(tan_phi_prime);
    
    // Parallax constants (dimensionless, relative to Earth radius)
    double sin_phi_prime = std::sin(phi_prime);
    double cos_phi_prime = std::cos(phi_prime);
    
    rho_sin_phi = sin_phi_prime + h_earth_radii * sin_lat;
    rho_cos_phi = cos_phi_prime + h_earth_radii * cos_lat;
}

} // namespace observations
} // namespace astdyn
