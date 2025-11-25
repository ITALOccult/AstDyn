/**
 * @file MPCReader.cpp
 * @brief Implementation of MPC observation parser
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 */

#include "astdyn/observations/MPCReader.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/utils/StringUtils.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace astdyn {
namespace observations {

using namespace utils;

std::vector<OpticalObservation> MPCReader::readFile(const std::string& filepath) {
    std::vector<OpticalObservation> observations;
    std::ifstream file(filepath);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        auto obs = parseLine(line);
        if (obs) {
            observations.push_back(*obs);
        }
    }
    
    return observations;
}

std::optional<OpticalObservation> MPCReader::parseLine(const std::string& line) {
    // MPC format requires at least 80 characters
    if (line.length() < 80) {
        return std::nullopt;
    }
    
    try {
        OpticalObservation obs;
        
        // Parse designation (columns 1-12)
        std::string packed_des = line.substr(0, 12);
        obs.object_designation = parseDesignation(packed_des);
        
        // Parse discovery/note flags (columns 13-14)
        if (line.length() > 13 && line[12] == '*') {
            obs.is_discovery = true;
        }
        obs.note1 = parseNote1(line[13]);
        obs.note2 = parseNote2(line[14]);
        
        // Parse date (columns 16-32)
        std::string date_str = line.substr(15, 17);
        obs.mjd_utc = parseDate(date_str);
        
        // Parse RA (columns 33-44)
        std::string ra_str = line.substr(32, 12);
        obs.ra = parseRA(ra_str);
        
        // Parse Dec (columns 45-56)
        std::string dec_str = line.substr(44, 12);
        obs.dec = parseDec(dec_str);
        
        // Parse magnitude (columns 66-70, optional)
        if (line.length() >= 70) {
            std::string mag_str = line.substr(65, 5);
            obs.magnitude = parseMagnitude(mag_str);
            
            // Parse band (column 71)
            if (line.length() >= 71 && line[70] != ' ') {
                obs.mag_band = line[70];
            }
        }
        
        // Parse catalog code (column 72, optional)
        if (line.length() >= 72) {
            obs.catalog = parseCatalogCode(line[71]);
        }
        
        // Parse observatory code (columns 78-80)
        obs.observatory_code = trim(line.substr(77, 3));
        
        // Set default uncertainties (can be refined with astrometric catalogs)
        obs.sigma_ra = 0.5 * constants::ARCSEC_TO_RAD;  // 0.5 arcsec
        obs.sigma_dec = 0.5 * constants::ARCSEC_TO_RAD;
        
        return obs;
        
    } catch (const std::exception& e) {
        return std::nullopt;
    }
}

std::map<std::string, ObservationSet> MPCReader::readFileGrouped(const std::string& filepath) {
    std::map<std::string, ObservationSet> grouped;
    
    auto observations = readFile(filepath);
    for (const auto& opt_obs : observations) {
        std::string desig = opt_obs.object_designation;
        
        if (grouped.find(desig) == grouped.end()) {
            grouped[desig].object_designation = desig;
        }
        
        Observation obs;
        obs.type = ObservationType::OPTICAL_RA_DEC;
        obs.optical = opt_obs;
        grouped[desig].addObservation(obs);
    }
    
    // Sort each set by time
    for (auto& pair : grouped) {
        pair.second.sortByTime();
    }
    
    return grouped;
}

// ============================================================================
// Parsing Helper Functions
// ============================================================================

std::string MPCReader::parseDesignation(const std::string& packed) {
    // Check if it's a numbered asteroid (columns 1-5 non-blank)
    std::string num_part = packed.substr(0, 5);
    if (num_part[0] != ' ') {
        int number = unpackNumber(num_part);
        if (number > 0) {
            return "(" + std::to_string(number) + ")";
        }
    }
    
    // Otherwise, it's a provisional designation or comet
    std::string prov_part = trim(packed.substr(5, 7));
    if (!prov_part.empty()) {
        return unpackProvisional(packed);
    }
    
    return trim(packed);
}

double MPCReader::parseDate(const std::string& date_str) {
    // Format: "YYYY MM DD.ddddd"
    std::istringstream iss(date_str);
    int year, month;
    double day;
    
    iss >> year >> month >> day;
    
    // Convert to MJD
    int day_int = static_cast<int>(day);
    double fraction = day - day_int;
    return time::calendar_to_mjd(year, month, day_int, fraction);
}

double MPCReader::parseRA(const std::string& ra_str) {
    // Format: "HH MM SS.ddd"
    std::istringstream iss(ra_str);
    int hours, minutes;
    double seconds;
    
    iss >> hours >> minutes >> seconds;
    
    // Convert to degrees, then radians
    double ra_deg = hours * 15.0 + minutes * 0.25 + seconds * (15.0 / 3600.0);
    return ra_deg * constants::DEG_TO_RAD;
}

double MPCReader::parseDec(const std::string& dec_str) {
    // Format: "sDD MM SS.dd"
    char sign = dec_str[0];
    
    std::istringstream iss(dec_str.substr(1));
    int degrees, minutes;
    double seconds;
    
    iss >> degrees >> minutes >> seconds;
    
    // Convert to degrees
    double dec_deg = degrees + minutes / 60.0 + seconds / 3600.0;
    if (sign == '-') {
        dec_deg = -dec_deg;
    }
    
    return dec_deg * constants::DEG_TO_RAD;
}

std::optional<double> MPCReader::parseMagnitude(const std::string& mag_str) {
    std::string trimmed = trim(mag_str);
    if (trimmed.empty()) {
        return std::nullopt;
    }
    
    try {
        return std::stod(trimmed);
    } catch (...) {
        return std::nullopt;
    }
}

std::string MPCReader::parseNote1(char c) {
    if (c == '*') return "Discovery";
    if (c == 'P') return "Photographic";
    if (c == 'e') return "Encoder";
    if (c == 'C') return "CCD";
    if (c == 'T') return "Meridian/Transit Circle";
    if (c == 'M') return "Micrometer";
    if (c == 'V') return "Roving Observer";
    return "";
}

std::string MPCReader::parseNote2(char c) {
    if (c >= '0' && c <= '9') return std::string(1, c);
    if (c >= 'A' && c <= 'Z') return std::string(1, c);
    return "";
}

int MPCReader::unpackNumber(const std::string& packed) {
    // MPC packed numbers: 00001-99999, then A0000-Z9999, then a0000-z9999
    if (packed.empty()) return 0;
    
    char first = packed[0];
    if (first >= '0' && first <= '9') {
        // Simple number: 00001-99999
        try {
            return std::stoi(packed);
        } catch (...) {
            return 0;
        }
    } else if (first >= 'A' && first <= 'Z') {
        // Extended: A=100000, B=110000, ..., Z=350000
        int base = 100000 + (first - 'A') * 10000;
        try {
            return base + std::stoi(packed.substr(1, 4));
        } catch (...) {
            return 0;
        }
    } else if (first >= 'a' && first <= 'z') {
        // Further extended: a=360000, b=370000, etc.
        int base = 360000 + (first - 'a') * 10000;
        try {
            return base + std::stoi(packed.substr(1, 4));
        } catch (...) {
            return 0;
        }
    }
    
    return 0;
}

std::string MPCReader::unpackProvisional(const std::string& packed) {
    // Simplified unpacking - handles common formats like K14A00A → 2014 AA
    std::string trimmed = trim(packed);
    if (trimmed.length() < 7) return trimmed;
    
    // Century/year encoding
    char century = trimmed[0];
    std::string year_part = trimmed.substr(1, 2);
    
    int year = 0;
    if (century == 'I') year = 1800 + std::stoi(year_part);
    else if (century == 'J') year = 1900 + std::stoi(year_part);
    else if (century == 'K') year = 2000 + std::stoi(year_part);
    else return trimmed;
    
    // Survey/order encoding (simplified)
    std::string suffix = trimmed.substr(3);
    
    return std::to_string(year) + " " + suffix;
}

} // namespace observations
} // namespace astdyn
