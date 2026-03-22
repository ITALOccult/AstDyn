#include <astdyn/catalog/gaia/types.h>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

namespace ioc {
namespace gaia {

// =============================================================================
// GaiaStar Implementation
// =============================================================================

GaiaStar::GaiaStar()
    : source_id(0),
      ra(0.0), dec(0.0),
      parallax(0.0), parallax_error(0.0),
      pmra(0.0), pmdec(0.0),
      pmra_error(0.0), pmdec_error(0.0),
      phot_g_mean_mag(0.0), phot_bp_mean_mag(0.0), phot_rp_mean_mag(0.0),
      bp_rp(0.0),
      astrometric_excess_noise(0.0), astrometric_chi2_al(0.0),
      visibility_periods_used(0),
      ruwe(0.0) {}

bool GaiaStar::isValid() const {
    // Check basic validity criteria
    if (source_id <= 0) return false;
    
    // RA range [0, 360)
    if (ra < 0.0 || ra >= 360.0) return false;
    
    // Dec range [-90, 90]
    if (dec < -90.0 || dec > 90.0) return false;
    
    // G magnitude should be reasonable
    if (phot_g_mean_mag < -5.0 || phot_g_mean_mag > 30.0) return false;
    
    // Parallax error should be positive
    if (parallax_error < 0.0) return false;
    
    return true;
}

double GaiaStar::getBpRpColor() const {
    if (std::isnan(phot_bp_mean_mag) || std::isnan(phot_rp_mean_mag)) {
        return std::nan("");
    }
    return phot_bp_mean_mag - phot_rp_mean_mag;
}

std::string GaiaStar::getDesignation() const {
    std::ostringstream ss;
    ss << "Gaia DR3 " << source_id;
    return ss.str();
}

std::string GaiaStar::getAllDesignations() const {
    std::vector<std::string> designations;
    
    // Add Gaia designation (always present)
    designations.push_back("Gaia DR3 " + std::to_string(source_id));
    
    // Add other catalog designations if available
    if (!sao_designation.empty()) {
        designations.push_back("SAO " + sao_designation);
    }
    
    if (!tycho2_designation.empty()) {
        designations.push_back("TYC " + tycho2_designation);
    }
    
    if (!hd_designation.empty()) {
        designations.push_back("HD " + hd_designation);
    }
    
    if (!hip_designation.empty()) {
        designations.push_back("HIP " + hip_designation);
    }
    
    if (!common_name.empty()) {
        designations.push_back(common_name);
    }
    
    // Join all designations with " | "
    std::string result;
    for (size_t i = 0; i < designations.size(); ++i) {
        if (i > 0) result += " | ";
        result += designations[i];
    }
    
    return result;
}

GaiaStar GaiaStar::propagateToEpoch(double mjd_tdb) const {
    GaiaStar result = *this;
    
    // Gaia DR3 epoch is J2016.0 = JD 2457389.0 TCB
    // We approximate to JD 2457389.0 TDB for convenience
    const double MJD_GAIA_DR3 = 2457389.0 - 2400000.5;
    double dt_years = (mjd_tdb - MJD_GAIA_DR3) / 365.25;
    
    if (std::abs(dt_years) < 1e-9) return result;
    
    double cos_dec = std::cos(dec * M_PI / 180.0);
    if (std::abs(cos_dec) > 1e-6) {
        result.ra += (pmra * dt_years / 3600000.0) / cos_dec;
    }
    result.dec += (pmdec * dt_years / 3600000.0);
    
    return result;
}

// =============================================================================
// JulianDate Implementation
// =============================================================================

JulianDate JulianDate::fromCalendar(int year, int month, int day, double ut) {
    // Julian Date calculation (Meeus algorithm)
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    
    double jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    jd += (ut / 24.0) - 0.5;  // Add fractional day from UT (JD starts at noon)
    
    return JulianDate(jd);
}

void JulianDate::toCalendar(int& year, int& month, int& day, double& ut) const {
    // Inverse algorithm (JD starts at noon, so add 0.5 for midnight-based calendar)
    double jd_midnight = jd + 0.5;
    int a = static_cast<int>(jd_midnight);
    int b = a + 1537;
    int c = (b - 122.1) / 365.25;
    int d = 365.25 * c;
    int e = (b - d) / 30.6001;
    
    day = b - d - static_cast<int>(30.6001 * e);
    month = e - (e < 14 ? 1 : 13);
    year = c - (month > 2 ? 4716 : 4715);
    
    // Fractional part is UT
    ut = (jd_midnight - a) * 24.0;
    
}

// =============================================================================
// Utility Functions
// =============================================================================

/**
 * Convert string to GaiaRelease enum
 */
GaiaRelease stringToRelease(const std::string& str) {
    if (str == "DR2") return GaiaRelease::DR2;
    if (str == "EDR3") return GaiaRelease::EDR3;
    if (str == "DR3") return GaiaRelease::DR3;
    
    throw GaiaException(ErrorCode::INVALID_PARAMS, 
                       "Invalid Gaia release: " + str);
}

/**
 * Convert GaiaRelease enum to string
 */
std::string releaseToString(GaiaRelease release) {
    switch (release) {
        case GaiaRelease::DR2:  return "DR2";
        case GaiaRelease::EDR3: return "EDR3";
        case GaiaRelease::DR3:  return "DR3";
        default: return "UNKNOWN";
    }
}

/**
 * Calculate angular distance between two coordinates (haversine formula)
 * @return Distance in degrees
 */
double angularDistance(const EquatorialCoordinates& c1, 
                      const EquatorialCoordinates& c2) {
    constexpr double DEG_TO_RAD = M_PI / 180.0;
    
    double ra1 = c1.ra * DEG_TO_RAD;
    double dec1 = c1.dec * DEG_TO_RAD;
    double ra2 = c2.ra * DEG_TO_RAD;
    double dec2 = c2.dec * DEG_TO_RAD;
    
    double dra = ra2 - ra1;
    double ddec = dec2 - dec1;
    
    double a = std::sin(ddec / 2.0) * std::sin(ddec / 2.0) +
               std::cos(dec1) * std::cos(dec2) *
               std::sin(dra / 2.0) * std::sin(dra / 2.0);
    
    double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    
    return c / DEG_TO_RAD;  // Convert back to degrees
}

/**
 * Validate coordinate ranges
 */
bool isValidCoordinate(const EquatorialCoordinates& coord) {
    return (coord.ra >= 0.0 && coord.ra < 360.0 &&
            coord.dec >= -90.0 && coord.dec <= 90.0);
}

/**
 * Validate query parameters
 */
bool isValidQueryParams(const QueryParams& params) {
    if (params.radius <= 0.0 || params.radius > 180.0) return false;
    if (params.max_magnitude < 0.0 || params.max_magnitude > 30.0) return false;
    
    EquatorialCoordinates center(params.ra_center, params.dec_center);
    return isValidCoordinate(center);
}

/**
 * Validate box region
 */
bool isValidBoxRegion(const BoxRegion& box) {
    if (box.ra_min < 0.0 || box.ra_max > 360.0) return false;
    if (box.dec_min < -90.0 || box.dec_max > 90.0) return false;
    if (box.ra_min >= box.ra_max) return false;
    if (box.dec_min >= box.dec_max) return false;
    return true;
}

// =============================================================================
// CorridorQueryParams Implementation
// =============================================================================

namespace {
    // Simple JSON parser for corridor queries
    std::string extractValue(const std::string& json, const std::string& key) {
        size_t pos = json.find("\"" + key + "\"");
        if (pos == std::string::npos) return "";
        
        pos = json.find(":", pos);
        if (pos == std::string::npos) return "";
        
        // Skip whitespace
        pos++;
        while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t' || json[pos] == '\n')) pos++;
        
        if (pos >= json.size()) return "";
        
        // Handle string value
        if (json[pos] == '"') {
            pos++;
            size_t end = json.find("\"", pos);
            if (end == std::string::npos) return "";
            return json.substr(pos, end - pos);
        }
        
        // Handle number value
        size_t end = pos;
        while (end < json.size() && (std::isdigit(json[end]) || json[end] == '.' || json[end] == '-' || json[end] == '+' || json[end] == 'e' || json[end] == 'E')) {
            end++;
        }
        return json.substr(pos, end - pos);
    }
    
    std::vector<CelestialPoint> extractPath(const std::string& json) {
        std::vector<CelestialPoint> path;
        
        size_t pos = json.find("\"path\"");
        if (pos == std::string::npos) return path;
        
        pos = json.find("[", pos);
        if (pos == std::string::npos) return path;
        
        size_t end = json.find("]", pos);
        if (end == std::string::npos) return path;
        
        std::string path_str = json.substr(pos, end - pos + 1);
        
        // Parse each point object
        size_t obj_start = 0;
        while ((obj_start = path_str.find("{", obj_start)) != std::string::npos) {
            size_t obj_end = path_str.find("}", obj_start);
            if (obj_end == std::string::npos) break;
            
            std::string obj = path_str.substr(obj_start, obj_end - obj_start + 1);
            
            std::string ra_str = extractValue(obj, "ra");
            std::string dec_str = extractValue(obj, "dec");
            
            if (!ra_str.empty() && !dec_str.empty()) {
                CelestialPoint pt;
                pt.ra = std::stod(ra_str);
                pt.dec = std::stod(dec_str);
                path.push_back(pt);
            }
            
            obj_start = obj_end + 1;
        }
        
        return path;
    }
    
    double angularDistancePoints(double ra1, double dec1, double ra2, double dec2) {
        constexpr double DEG_TO_RAD = M_PI / 180.0;
        
        double dra = (ra2 - ra1) * DEG_TO_RAD;
        double ddec = (dec2 - dec1) * DEG_TO_RAD;
        
        double a = std::sin(ddec / 2.0) * std::sin(ddec / 2.0) +
                   std::cos(dec1 * DEG_TO_RAD) * std::cos(dec2 * DEG_TO_RAD) *
                   std::sin(dra / 2.0) * std::sin(dra / 2.0);
        
        double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
        
        return c / DEG_TO_RAD;
    }
}

CorridorQueryParams CorridorQueryParams::fromJSON(const std::string& json) {
    CorridorQueryParams params;
    
    // Parse path
    params.path = extractPath(json);
    
    // Parse width
    std::string width_str = extractValue(json, "width");
    if (!width_str.empty()) {
        params.width = std::stod(width_str);
    }
    
    // Parse max_magnitude
    std::string mag_str = extractValue(json, "max_magnitude");
    if (!mag_str.empty()) {
        params.max_magnitude = std::stod(mag_str);
    }
    
    // Parse min_parallax
    std::string parallax_str = extractValue(json, "min_parallax");
    if (!parallax_str.empty()) {
        params.min_parallax = std::stod(parallax_str);
    }
    
    // Parse max_results
    std::string results_str = extractValue(json, "max_results");
    if (!results_str.empty()) {
        params.max_results = std::stoull(results_str);
    }
    
    return params;
}

std::string CorridorQueryParams::toJSON() const {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6);
    
    ss << "{\n";
    ss << "  \"path\": [\n";
    for (size_t i = 0; i < path.size(); ++i) {
        ss << "    {\"ra\": " << path[i].ra << ", \"dec\": " << path[i].dec << "}";
        if (i < path.size() - 1) ss << ",";
        ss << "\n";
    }
    ss << "  ],\n";
    ss << "  \"width\": " << width << ",\n";
    ss << "  \"max_magnitude\": " << max_magnitude << ",\n";
    ss << "  \"min_parallax\": " << min_parallax << ",\n";
    ss << "  \"max_results\": " << max_results << "\n";
    ss << "}";
    
    return ss.str();
}

double CorridorQueryParams::getPathLength() const {
    double length = 0.0;
    for (size_t i = 1; i < path.size(); ++i) {
        length += angularDistancePoints(path[i-1].ra, path[i-1].dec, 
                                        path[i].ra, path[i].dec);
    }
    return length;
}

} // namespace gaia
} // namespace ioc
