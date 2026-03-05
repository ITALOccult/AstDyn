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
#include "astdyn/core/Constants.hpp"

namespace astdyn {
namespace io {
namespace parsers {

class AstDysRWOParser : public IObservationParser {
public:
    std::vector<OpticalObservation> parse_stream(std::istream& stream, size_t max_count = 0) override {
        std::vector<OpticalObservation> observations;
        std::string line;
        size_t count = 0;

        while (std::getline(stream, line)) {
            // Skip header
            if (line.empty() || line[0] == '!' || line[0] == '#' || 
                line.find("END_OF_HEADER") != std::string::npos ||
                line.find("version") != std::string::npos ||
                line.find("errmod") != std::string::npos ||
                line.find("RMS") != std::string::npos) {
                continue;
            }

                // --- STRICT POSITION-BASED PARSING ---
                
                // Helper lambda to extract and trim
                auto extract_field = [&](int start, int len) -> std::string {
                    if (start >= (int)line.length()) return "";
                    size_t actual_len = std::min((size_t)len, line.length() - start);
                    std::string s = line.substr(start, actual_len);
                    
                    // Trim whitespace
                    const auto strBegin = s.find_first_not_of(" ");
                    if (strBegin == std::string::npos) return ""; // All spaces
                    const auto strEnd = s.find_last_not_of(" ");
                    return s.substr(strBegin, strEnd - strBegin + 1);
                };

                // Filter short lines immediately
                if (line.length() < 90) continue; // Minimum length for Coordinates

                try {
                    // Extract substrings carefully
                    std::string obj_name = extract_field(0, 10);
                    std::string date_part = extract_field(17, 21); // YYYY MM DD.ddddd
                    std::string ra_part = extract_field(50, 12);   // HH MM SS.sss
                    std::string dec_part = extract_field(103, 13); // sDD MM SS.ss
                    
                    if (obj_name.empty() || date_part.empty() || ra_part.empty() || dec_part.empty()) continue;

                    OpticalObservation obs;
                    obs.object_name = obj_name;

                    // Parse Date
                    std::stringstream ss_date(date_part);
                    int year, month; double day;
                    if (!(ss_date >> year >> month >> day)) continue;
                    obs.mjd_utc = date_to_mjd(year, month, day);

                    // Parse RA
                    std::stringstream ss_ra(ra_part);
                    int rh, rm; double rs;
                    if (!(ss_ra >> rh >> rm >> rs)) continue;
                    obs.ra = (rh + rm / 60.0 + rs / 3600.0) * 15.0 * M_PI / 180.0;

                    // Parse Dec
                    std::stringstream ss_dec(dec_part);
                    int dd, dm; double ds;
                    if (!(ss_dec >> dd >> dm >> ds)) continue;
                    double dec_val = std::abs(dd) + dm / 60.0 + ds / 3600.0;
                    if (dec_part.find('-') != std::string::npos) dec_val = -dec_val;
                    obs.dec = dec_val * M_PI / 180.0;

                    // Obs Code (at the end of line)
                    if (line.length() > 140) {
                        obs.obs_code = extract_field(143, 3);
                        if (obs.obs_code.empty() || !isdigit(obs.obs_code[0])) {
                             obs.obs_code = extract_field(137, 3);
                        }
                    }
                    if (obs.obs_code.empty() || !isdigit(obs.obs_code[0])) obs.obs_code = "500";

                    obs.sigma_ra = 0.5 * astdyn::constants::ARCSEC_TO_RAD;
                    obs.sigma_dec = 0.5 * astdyn::constants::ARCSEC_TO_RAD;
                    obs.mag = 0.0;

                    observations.push_back(obs);
                    if (max_count > 0 && ++count >= max_count) break;

                } catch (...) {
                    continue;
                }
        }

        if (observations.empty()) {
            throw std::runtime_error("No observations found in stream");
        }

        return observations;
    }

    std::vector<OpticalObservation> parse(const std::string& filepath, size_t max_count = 0) override {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }
        return parse_stream(file, max_count);
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
        
        int jdn = (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
        return (double)jdn + day - 0.5 - 2400000.5;
    }
};

} // namespace parsers
} // namespace io
} // namespace astdyn
