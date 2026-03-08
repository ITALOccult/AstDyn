/**
 * @file MPCParser.cpp
 * @brief MPC 80-column observation parser.
 */

#include "astdyn/io/MPCParser.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

namespace astdyn::io {

observations::OpticalObservation MPCParser::parse_line(const std::string& line) {
    if (line.size() < 80) return {};

    observations::OpticalObservation obs;
    
    // Identification
    obs.object_designation = line.substr(0, 12);
    obs.observatory_code = line.substr(77, 3);
    
    // Time (15-32)
    obs.time = parse_date(line.substr(15, 17));
    
    // RA (32-44)
    obs.ra = parse_ra(line.substr(32, 12));
    
    // Dec (44-56)
    obs.dec = parse_dec(line.substr(44, 12));
    
    // Magnitude (65-70)
    std::string mag_str = line.substr(65, 5);
    if (!mag_str.empty() && mag_str != "     ") {
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
        if (line.size() >= 80) {
            results.push_back(parse_line(line));
        }
    }
    return results;
}

time::EpochUTC MPCParser::parse_date(const std::string& date_str) {
    int y, m;
    double d;
    std::stringstream ss(date_str);
    ss >> y >> m >> d;
    
    // Convert Y/M/D to MJD
    // Formula from Meeus
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
    double h, m, s;
    if (!(ss >> h >> m >> s)) return astrometry::RightAscension();
    double deg = (h + m / 60.0 + s / 3600.0) * 15.0; // hrs -> deg
    return astrometry::RightAscension(astrometry::Angle::from_deg(deg));
}

astrometry::Declination MPCParser::parse_dec(const std::string& dec_str) {
    std::stringstream ss(dec_str);
    double d, m, s;
    char sign;
    if (!(ss >> std::skipws >> sign >> d >> m >> s)) return astrometry::Declination();
    double deg = d + m / 60.0 + s / 3600.0;
    if (sign == '-') deg = -deg;
    return astrometry::Declination(astrometry::Angle::from_deg(deg));
}

} // namespace astdyn::io
