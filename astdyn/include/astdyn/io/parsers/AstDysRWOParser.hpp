/**
 * @file AstDysRWOParser.hpp
 * @brief Parser for AstDyS .rwo files (MPC extended format)
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * Format: Fixed-width columns from AstDyS database
 * Columns:
 *  - Object name
 *  - Obs type (O=optical)
 *  - Discovery flag
 *  - Note
 *  - Date (YYYY MM DD.ddddd)
 *  - RA (HH MM SS.sss)
 *  - Dec (sDD MM SS.ss)
 *  - Magnitude
 *  - Observatory code
 *  - Residuals and other data
 */

#pragma once

#include "astdyn/io/IParser.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

namespace astdyn {
namespace io {
namespace parsers {

/**
 * @brief Parser for AstDyS .rwo files (MPC extended format with residuals)
 */
class AstDysRWOParser : public IObservationParser {
public:
    std::vector<OpticalObservation> parse(const std::string& filepath, size_t max_count = 0) override {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }

        std::vector<OpticalObservation> observations;
        std::string line;
        size_t count = 0;

        while (std::getline(file, line)) {
            // Skip header and comments
            if (line.empty() || line[0] == '!' || line[0] == '#' || 
                line.find("END_OF_HEADER") != std::string::npos ||
                line.find("version") != std::string::npos ||
                line.find("errmod") != std::string::npos ||
                line.find("RMSast") != std::string::npos ||
                line.find("RMSmag") != std::string::npos) {
                continue;
            }

            // Parse MPC extended format (fixed columns)
            if (line.length() < 80) continue;

            try {
                OpticalObservation obs;

                // Object name (columns 1-6, trimmed)
                obs.object_name = line.substr(0, 6);
                // Trim spaces
                obs.object_name.erase(0, obs.object_name.find_first_not_of(" "));
                obs.object_name.erase(obs.object_name.find_last_not_of(" ") + 1);

                // Date: YYYY MM DD.ddddddddd (columns 18-36)
                int year, month;
                double day;
                std::string date_str = line.substr(17, 19);
                std::istringstream date_iss(date_str);
                date_iss >> year >> month >> day;

                // Convert to MJD
                obs.mjd_utc = date_to_mjd(year, month, day);

                // RA: HH MM SS.sss (columns 55-66)
                std::string ra_str = line.substr(54, 12);
                int ra_h, ra_m;
                double ra_s;
                std::istringstream ra_iss(ra_str);
                ra_iss >> ra_h >> ra_m >> ra_s;
                obs.ra = (ra_h + ra_m / 60.0 + ra_s / 3600.0) * 15.0 * M_PI / 180.0;  // Convert to radians

                // Dec: sDD MM SS.ss (columns 84-95)
                std::string dec_str = line.substr(83, 12);
                char sign;
                int dec_d, dec_m;
                double dec_s;
                std::istringstream dec_iss(dec_str);
                dec_iss >> sign >> dec_d >> dec_m >> dec_s;
                double dec_deg = dec_d + dec_m / 60.0 + dec_s / 3600.0;
                if (sign == '-') dec_deg = -dec_deg;
                obs.dec = dec_deg * M_PI / 180.0;  // Convert to radians

                // Magnitude (columns 113-117, optional)
                if (line.length() > 116) {
                    std::string mag_str = line.substr(112, 5);
                    std::istringstream mag_iss(mag_str);
                    mag_iss >> obs.mag;
                    if (mag_iss.fail()) obs.mag = 0.0;
                } else {
                    obs.mag = 0.0;
                }

                // Observatory code (columns 142-144)
                if (line.length() > 144) {
                    obs.obs_code = line.substr(141, 3);
                    // Trim spaces
                    obs.obs_code.erase(0, obs.obs_code.find_first_not_of(" "));
                    obs.obs_code.erase(obs.obs_code.find_last_not_of(" ") + 1);
                }

                // Default uncertainties (1 arcsec)
                obs.sigma_ra = 1.0;
                obs.sigma_dec = 1.0;

                observations.push_back(obs);

                if (max_count > 0 && ++count >= max_count) {
                    break;
                }

            } catch (const std::exception& e) {
                // Skip malformed lines
                continue;
            }
        }

        if (observations.empty()) {
            throw std::runtime_error("No observations found in file: " + filepath);
        }

        return observations;
    }

    std::string name() const override {
        return "AstDyS RWO Parser (MPC Extended Format)";
    }

    bool can_handle(const std::string& filepath) const override {
        return filepath.find(".rwo") != std::string::npos || 
               filepath.find(".RWO") != std::string::npos;
    }

private:
    /**
     * @brief Convert calendar date to MJD
     */
    double date_to_mjd(int year, int month, double day) const {
        // Julian Day Number calculation
        int a = (14 - month) / 12;
        int y = year + 4800 - a;
        int m = month + 12 * a - 3;
        
        int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
        
        // MJD = JD - 2400000.5
        return jdn - 2400000.5;
    }
};

} // namespace parsers
} // namespace io
} // namespace astdyn
