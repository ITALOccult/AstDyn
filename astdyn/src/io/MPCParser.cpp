/**
 * @file MPCParser.cpp
 * @brief MPC 80-column observation parser.
 */

#include "astdyn/io/MPCParser.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>

namespace astdyn::io {

observations::OpticalObservation MPCParser::parse_line(const std::string& line) {
    observations::OpticalObservation obs = {};
    if (line.size() < 80) return obs;

    // Identification
    obs.object_designation = line.substr(0, 12);
    obs.observatory_code = line.substr(77, 3);
    
    std::string date_part = line.substr(15, 17); // Correctly 17 chars for precision
    obs.time = parse_date(date_part);
    
    // RA (32-44)
    std::string ra_part = line.substr(32, 12);
    obs.ra = parse_ra(ra_part);
    
    // Dec (44-56)
    std::string dec_part = line.substr(44, 12);
    obs.dec = parse_dec(dec_part);
    
    // Magnitude (65-71)
    std::string mag_str = line.substr(65, 6);
    if (!mag_str.empty() && mag_str != "      ") {
        try {
            obs.magnitude = std::stod(mag_str);
        } catch (...) {}
    }
    
    return obs;
}

std::vector<observations::OpticalObservation> MPCParser::parse_file(const std::string& content) {
    std::vector<observations::OpticalObservation> results;
    std::stringstream ss(content);
    std::string line;
    while (std::getline(ss, line)) {
        if (line.size() < 80) continue;
        // Skip auxiliary satellite-offset records (MPC column 15 == 's').
        // These encode the observer spacecraft position in km, not RA/Dec.
        if (line[14] == 's') continue;
        results.push_back(parse_line(line));
    }
    return results;
}

time::EpochUTC MPCParser::parse_date(const std::string& date_str) {
    if (date_str.empty()) return {};

    int y = 0, m = 0;
    double d = 0.0;

    // MPC Dates usually start at column 16. Substring is 16 chars.
    // Examples: "2004 03 15.12" or "K04 03 15.12"
    if (std::isdigit(date_str[0])) {
        // Numeric year (e.g., "2004")
        std::stringstream ss(date_str);
        if (!(ss >> y >> m >> d)) return {};
    } else {
        // Packed year (e.g., "K04")
        char c = date_str[0];
        int century = 0;
        if (c == 'I') century = 1800;
        else if (c == 'J') century = 1900;
        else if (c == 'K') century = 2000;
        else return {};

        try {
            y = century + std::stoi(date_str.substr(1, 2));
            m = std::stoi(date_str.substr(4, 2)); // Month at index 4 (packed format: "KYY MM DD.ddddd")
            d = std::stod(date_str.substr(7));   // Day starts after second space
        } catch(...) { return {}; }
    }
    
    // Julian Date conversion (Meeus formula)
    if (m <= 2) {
        y -= 1;
        m += 12;
    }
    int A = y / 100;
    int B = 2 - A + (A / 4);
    double jd = std::floor(365.25 * (y + 4716)) + std::floor(30.6001 * (m + 1)) + d + B - 1524.5;
    return time::EpochUTC::from_mjd(jd - 2400000.5);
}

astrometry::RightAscension MPCParser::parse_ra(const std::string& ra_str) {
    std::stringstream ss(ra_str);
    double h = 0, m = 0, s = 0;
    if (!(ss >> h >> m >> s)) {
        // Fallback for cases with weird formatting (e.g. trailing characters)
        return astrometry::RightAscension();
    }
    double deg = (h + m / 60.0 + s / 3600.0) * 15.0; // hrs -> deg
    return astrometry::RightAscension(astrometry::Angle::from_deg(deg));
}

astrometry::Declination MPCParser::parse_dec(const std::string& dec_str) {
    // MPC Dec: "+DD MM SS.ss" or "-DD MM SS.ss"
    // Use manual sign check to avoid stream confusion with internal spaces
    double d = 0, m = 0, s = 0;
    double sign = 1.0;
    
    size_t first_non_space = dec_str.find_first_not_of(' ');
    if (first_non_space == std::string::npos) return astrometry::Declination();
    
    if (dec_str[first_non_space] == '-') sign = -1.0;
    
    // Clean string for streaming (replace sign with space)
    std::string clean = dec_str;
    if (clean[first_non_space] == '+' || clean[first_non_space] == '-') {
        clean[first_non_space] = ' ';
    }
    
    std::stringstream ss(clean);
    if (!(ss >> d >> m >> s)) return astrometry::Declination();
    
    double deg = sign * (d + m / 60.0 + s / 3600.0);
    return astrometry::Declination(astrometry::Angle::from_deg(deg));
}

} // namespace astdyn::io
