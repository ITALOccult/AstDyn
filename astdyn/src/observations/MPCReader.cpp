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

using namespace astdyn::utils;

std::vector<OpticalObservation> MPCReader::readStream(std::istream& stream) {
    std::vector<OpticalObservation> observations;
    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        auto obs = parseLine(line);
        if (obs) {
            observations.push_back(*obs);
        }
    }
    
    return observations;
}

std::vector<OpticalObservation> MPCReader::readFile(const std::string& filepath) {
    std::ifstream file(filepath);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    return readStream(file);
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
        obs.time = parseDate(date_str);
        
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
        obs.sigma_ra = astrometry::Angle::from_arcsec(0.5);  // 0.5 arcsec
        obs.sigma_dec = astrometry::Angle::from_arcsec(0.5);
        
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

time::EpochUTC MPCReader::parseDate(const std::string& date_str) {
    // Format: "YYYY MM DD.ddddd"
    std::istringstream iss(date_str);
    int year, month;
    double day;
    
    iss >> year >> month >> day;
    
    // Convert to MJD
    int day_int = static_cast<int>(day);
    double fraction = day - day_int;
    double mjd = time::calendar_to_mjd(year, month, day_int, fraction);
    return time::EpochUTC::from_mjd(mjd);
}

astrometry::RightAscension MPCReader::parseRA(const std::string& ra_str) {
    // Format: "HH MM SS.ddd"
    std::istringstream iss(ra_str);
    int hours, minutes;
    double seconds;
    
    iss >> hours >> minutes >> seconds;
    
    // Convert to degrees
    double ra_deg = hours * 15.0 + minutes * 0.25 + seconds * (15.0 / 3600.0);
    return astrometry::RightAscension(astrometry::Angle::from_deg(ra_deg));
}

astrometry::Declination MPCReader::parseDec(const std::string& dec_str) {
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
    
    return astrometry::Declination(astrometry::Angle::from_deg(dec_deg));
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

// ============================================================================
// MPCWriter Implementation
// ============================================================================

void MPCWriter::writeFile(const std::string& filepath, const std::vector<OpticalObservation>& observations) {
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }
    
    for (const auto& obs : observations) {
        file << formatLine(obs) << "\n";
    }
}

std::string MPCWriter::formatLine(const OpticalObservation& obs) {
    std::stringstream ss;
    
    // 1-12: Designation
    ss << std::left << std::setw(12) << packDesignation(obs.object_designation);
    
    // 13: Discovery (*)
    ss << (obs.is_discovery ? "*" : " ");
    
    // 14-15: Notes
    ss << (obs.note1.empty() ? " " : obs.note1.substr(0,1));
    ss << (obs.note2.empty() ? " " : obs.note2.substr(0,1));
    
    // 16-32: Date (YYYY MM DD.ddddd)
    ss << formatDate(obs.time);
    
    // 33-44: RA
    ss << formatRA(obs.ra);
    
    // 45-56: Dec
    ss << formatDec(obs.dec);
    
    // 57-65: Blank
    ss << "         ";
    
    // 66-70: Magnitude
    if (obs.magnitude) {
        ss << std::right << std::fixed << std::setprecision(1) << std::setw(5) << *obs.magnitude;
    } else {
        ss << "     ";
    }
    
    // 71: Band
    ss << (obs.mag_band ? *obs.mag_band : ' ');
    
    // 72-77: Catalog/Offset stuff (blank for now)
    ss << "      ";
    
    // 78-80: Observatory code
    ss << std::right << std::setw(3) << obs.observatory_code;
    
    return ss.str();
}

std::string MPCWriter::packDesignation(const std::string& designation) {
    // Basic implementation: truncate or pad to 12 chars
    if (designation.length() > 12) return designation.substr(0, 12);
    return designation;
}

std::string MPCWriter::formatDate(time::EpochUTC t_utc) {
    auto [year, month, day, fraction] = time::mjd_to_calendar(t_utc.mjd());
    
    std::stringstream ss;
    ss << " " << std::setfill('0') << std::setw(4) << year << " " 
       << std::setw(2) << month << " " << std::fixed << std::setprecision(5) << std::setw(8) << (day + fraction);
    return ss.str();
}

std::string MPCWriter::formatRA(astrometry::RightAscension ra) {
    return " " + ra.to_hms();
}

std::string MPCWriter::formatDec(astrometry::Declination dec) {
    return dec.to_dms();
}

} // namespace observations
} // namespace astdyn
