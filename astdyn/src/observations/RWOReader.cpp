/**
 * @file RWOReader.cpp
 * @brief Implementation of OrbFit .rwo format parser
 * @author ITALOccult AstDyn Team
 * @date 2025-11-25
 */

#include "astdyn/observations/RWOReader.hpp"
#include "astdyn/observations/MPCReader.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/utils/StringUtils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace astdyn {
namespace observations {

using namespace utils;

std::vector<OpticalObservation> RWOReader::readFile(const std::string& filepath) {
    // RWO format is OrbFit's proprietary format, NOT standard MPC!
    // We must parse it directly using the format specification from OrbFit
    
    std::vector<OpticalObservation> observations;
    std::ifstream file(filepath);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        
        auto obs = parseLine(line);
        if (obs) {
            observations.push_back(*obs);
        }
    }
    
    return observations;
}

std::optional<OpticalObservation> RWOReader::parseLine(const std::string& line) {
    // Skip header lines and ensure line is long enough
    if (line.length() < 150 || line[0] == '!' || line.find("END_OF_HEADER") != std::string::npos) {
        return std::nullopt;
    }
    
    try {
        OpticalObservation obs;
        
        // Fortran FORMAT 110: (1X,A9,1X,A1,1X,A1,1X,A1)
        // Object designation: column 2-10 (1-indexed) = index 1-9 (0-indexed)
        obs.object_designation = trim(line.substr(1, 9));
        
        // Observation type: column 12 (1-indexed) = index 11 (0-indexed)
        if (line.length() < 12 || line[11] != 'O') {
            return std::nullopt;  // Not an optical observation
        }
        
        // Fortran FORMAT 120: (I4,1X,I2,1X,F13.10,1X,E10.3)
        // Date starts at column 18 (1-indexed) = index 17 (0-indexed) in Fortran READ(record(18:),120)
        // But since record starts at column 1, READ(record(18:)) means skip first 17 chars
        // So in our 0-indexed string, date is at position 17
        // Format is: "YYYY MM DD.dddddddddd" (4 + 1 + 2 + 1 + 13 = 21 chars)
        std::string date_str = line.substr(17, 21);
        obs.mjd_utc = parseDate(date_str);
        
        // Fortran FORMAT 521 starts at READ(record(51:),521)
        // In Fortran, record(51:) means from character 51 onwards (1-indexed)
        // In C++, this is index 50 (0-indexed)
        // FORMAT: I2,1X,I2,1X,F6.6,1X,E10.3,...
        size_t ra_start = 50;
        if (line.length() < ra_start + 50) return std::nullopt;
        
        // RA: I2,1X,I2,1X,F6.6 at position 50
        std::string ra_str = trim(line.substr(ra_start, 12));
        obs.ra = parseRA(ra_str);
        
        // After RA: 1X,E10.3 (acc_coord), 1X,A8 (rms), etc.
        // Dec sign is after: 1X,E10.3,1X,A8,1X,L1,1X,F8.3,A9,1X
        // That's about 40 chars after RA start
        // Dec: 1X,A1,I2,1X,I2,1X,F5.2
        size_t dec_start = ra_start + 12 + 1 + 10 + 1 + 8 + 1 + 1 + 1 + 8 + 9 + 1;
        // dec_start ≈ 50 + 54 = 104
        if (line.length() < dec_start + 13) return std::nullopt;
        
        std::string dec_str = line.substr(dec_start, 13);
        obs.dec = parseDec(dec_str);
        
        // Observatory code: parse from end of line
        obs.observatory_code = parseObservatory(line);
        
        // Magnitude
        obs.magnitude = parseMagnitude(line);
        
        // Set default uncertainties (0.5 arcsec)
        obs.sigma_ra = 0.5 * constants::ARCSEC_TO_RAD;
        obs.sigma_dec = 0.5 * constants::ARCSEC_TO_RAD;
        
        return obs;
        
    } catch (const std::exception& e) {
        std::cerr << "Error parsing RWO line: " << e.what() << std::endl;
        return std::nullopt;
    }
}

std::map<std::string, ObservationSet> RWOReader::readFileGrouped(const std::string& filepath) {
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

std::string RWOReader::parseDesignation(const std::string& line) {
    return trim(line.substr(0, 14));
}

double RWOReader::parseDate(const std::string& date_str) {
    // Format: "YYYY MM DD.ddddddddddd" from Fortran FORMAT 120: (I4,1X,I2,1X,F13.10)
    std::istringstream iss(date_str);
    int year, month;
    double day;
    
    iss >> year >> month >> day;
    
    if (year < 1000 || year > 3000 || month < 1 || month > 12 || day < 1 || day > 32) {
        std::cerr << "WARNING: Invalid date parsed: " << year << "/" << month << "/" << day 
                  << " from string: '" << date_str << "'" << std::endl;
    }
    
    // Convert to MJD
    int day_int = static_cast<int>(day);
    double fraction = day - day_int;
    return time::calendar_to_mjd(year, month, day_int, fraction);
}

double RWOReader::parseRA(const std::string& ra_str) {
    // Format: "HH MM SS.ddd"
    std::istringstream iss(ra_str);
    int hours, minutes;
    double seconds;
    
    iss >> hours >> minutes >> seconds;
    
    // Convert to degrees, then radians
    double ra_deg = hours * 15.0 + minutes * 0.25 + seconds * (15.0 / 3600.0);
    return ra_deg * constants::DEG_TO_RAD;
}

double RWOReader::parseDec(const std::string& dec_str) {
    // Format: "sDD MM SS.dd" where s is + or -
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

std::string RWOReader::parseObservatory(const std::string& line) {
    // Observatory code is typically near the end
    // Look for 3-character code after magnitude info
    // In the example: "... V D29 ..." where D29 is observatory
    
    // Search backwards for uppercase letters (observatory codes)
    size_t pos = line.length() - 1;
    while (pos > 100 && std::isspace(line[pos])) {
        pos--;
    }
    
    // Try to find a 3-letter code
    if (pos >= 102) {
        for (size_t i = pos - 20; i < pos - 2; i++) {
            if (std::isalnum(line[i]) && std::isalnum(line[i+1]) && std::isalnum(line[i+2])) {
                // Check if this looks like an observatory code (alphanumeric)
                std::string code = line.substr(i, 3);
                if (code.length() == 3) {
                    return code;
                }
            }
        }
    }
    
    return "500";  // Default: geocentric
}

std::optional<double> RWOReader::parseMagnitude(const std::string& line) {
    // Magnitude is typically around column 120-125
    // Format: "13.67o" or "14.1 R" etc.
    
    // Search for a number followed by optional magnitude band
    for (size_t i = 110; i < line.length() - 5 && i < 140; i++) {
        if (std::isdigit(line[i]) && i + 4 < line.length()) {
            std::string mag_str;
            size_t j = i;
            while (j < line.length() && (std::isdigit(line[j]) || line[j] == '.' || line[j] == '-')) {
                mag_str += line[j];
                j++;
            }
            
            if (!mag_str.empty()) {
                try {
                    double mag = std::stod(mag_str);
                    if (mag >= 0.0 && mag < 30.0) {  // Reasonable range
                        return mag;
                    }
                } catch (...) {
                    // Continue searching
                }
            }
        }
    }
    
    return std::nullopt;
}

} // namespace observations
} // namespace astdyn
