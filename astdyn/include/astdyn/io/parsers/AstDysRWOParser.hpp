/**
 * @file AstDysRWOParser.hpp  
 * @brief WORKING parser for AstDyS .rwo files with CORRECT Dec parsing
 * @author AstDyn Team
 * @date 2025-12-09
 * 
 * Uses field-based parsing (whitespace-separated) instead of fixed columns
 * 
 * Fields:
 * 0: Object name
 * 1-2: Obs type
 * 3-5: Date (YYYY MM DD.ddddd)
 * 7-9: RA (HH MM SS.sss)
 * 15-17: Dec (sDD MM SS.ss)
 */

#pragma once

#include "astdyn/io/IParser.hpp"
#include <fstream>
#include <sstream>
#include <cmath>

namespace astdyn {
namespace io {
namespace parsers {

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
            // Skip header
            if (line.empty() || line[0] == '!' || line[0] == '#' || 
                line.find("END_OF_HEADER") != std::string::npos ||
                line.find("version") != std::string::npos ||
                line.find("errmod") != std::string::npos ||
                line.find("RMS") != std::string::npos) {
                continue;
            }

            try {
                // Parse fields (whitespace-separated)
                std::istringstream iss(line);
                std::string fields[25];
                int nfields = 0;
                while (nfields < 25 && iss >> fields[nfields]) {
                    nfields++;
                }
                
                if (nfields < 18) continue;  // Need at least 18 fields
                
                OpticalObservation obs;

                // Object name (field 0)
                obs.object_name = fields[0];

                // Date: YYYY MM DD.ddddd (fields 3, 4, 5)
                int year = std::stoi(fields[3]);
                int month = std::stoi(fields[4]);
                double day = std::stod(fields[5]);
                obs.mjd_utc = date_to_mjd(year, month, day);

                // RA: HH MM SS.sss (fields 7, 8, 9)
                int ra_h = std::stoi(fields[7]);
                int ra_m = std::stoi(fields[8]);
                double ra_s = std::stod(fields[9]);
                obs.ra = (ra_h + ra_m / 60.0 + ra_s / 3600.0) * 15.0 * M_PI / 180.0;

                // Dec: sDD MM SS.ss (fields 15, 16, 17)
                std::string dec_d_str = fields[15];
                char sign = dec_d_str[0];
                int dec_d = std::abs(std::stoi(dec_d_str));
                int dec_m = std::stoi(fields[16]);
                double dec_s = std::stod(fields[17]);
                
                double dec_deg = dec_d + dec_m / 60.0 + dec_s / 3600.0;
                if (sign == '-') dec_deg = -dec_deg;
                obs.dec = dec_deg * M_PI / 180.0;

                // Magnitude (field 20, if available)
                if (nfields > 20) {
                    try {
                        obs.mag = std::stod(fields[20]);
                    } catch (...) {
                        obs.mag = 0.0;
                    }
                } else {
                    obs.mag = 0.0;
                }

                // Observatory code (field 23, if available)
                if (nfields > 23) {
                    obs.obs_code = fields[23];
                } else {
                    obs.obs_code = "500";
                }

                // Default uncertainties
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
        return "AstDyS RWO Parser (Field-Based - Correct Dec)";
    }

    bool can_handle(const std::string& filepath) const override {
        return filepath.find(".rwo") != std::string::npos || 
               filepath.find(".RWO") != std::string::npos;
    }

private:
    double date_to_mjd(int year, int month, double day) const {
        int a = (14 - month) / 12;
        int y = year + 4800 - a;
        int m = month + 12 * a - 3;
        
        int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
        return jdn - 2400000.5;
    }
};

} // namespace parsers
} // namespace io
} // namespace astdyn
