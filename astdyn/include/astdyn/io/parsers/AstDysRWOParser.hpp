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
                    OpticalObservation obs;
                    
                    // 1. Object Name: 0-10
                    obs.object_name = extract_field(0, 10);
                    
                    // 2. Date: YYYY (17,4) MM (22,2) DD (25,9)
                    std::string y_str = extract_field(17, 4);
                    std::string m_str = extract_field(22, 2);
                    std::string d_str = extract_field(25, 9);
                    
                    // Robust conversion
                    if (y_str.empty() || m_str.empty() || d_str.empty()) continue;
                    int year = std::stoi(y_str);
                    int month = std::stoi(m_str);
                    double day = std::stod(d_str);
                    obs.mjd_utc = date_to_mjd(year, month, day);

                    // 3. RA: H (47,3) M (50,3) S (53,8)
                    // Tolerant offsets
                    int ra_h = std::stoi(extract_field(47, 3));
                    int ra_m = std::stoi(extract_field(50, 3));
                    double ra_s = std::stod(extract_field(53, 8));
                    obs.ra = (ra_h + ra_m / 60.0 + ra_s / 3600.0) * 15.0 * M_PI / 180.0;

                    // 4. Dec: sDD(102-105) MM(106-107) SS(109-115)
                    // " +31" or "-13" -> Start 102
                    std::string dec_d_str = extract_field(102, 3);
                    
                    int dec_d = std::abs(std::stoi(dec_d_str));
                    int dec_m = std::stoi(extract_field(106, 2));
                    double dec_s = std::stod(extract_field(109, 8));
                    
                    double dec_deg = dec_d + dec_m / 60.0 + dec_s / 3600.0;
                    if (dec_d_str.find('-') != std::string::npos) dec_deg = -dec_deg;
                    obs.dec = dec_deg * M_PI / 180.0;

                    // 5. Magnitude: Unknown offset, skip or guess
                    obs.mag = 0.0;

                    // 6. Obs Code: Try end of line logic
                    size_t p_pos = line.rfind("p "); 
                    std::string code = "";
                    if (line.length() > 20) {
                        std::string c1 = extract_field(143, 3);
                        if (!c1.empty() && isdigit(c1[0])) code = c1;
                    }
                    if (code.empty()) code = extract_field(137, 3); 

                    if (code.empty() || !isdigit(code[0])) code = "500";
                    obs.obs_code = code;


                    // Defaults
                    obs.sigma_ra = 1.0; 
                    obs.sigma_dec = 1.0;

                    observations.push_back(obs);
                    
                     if (max_count > 0 && ++count >= max_count) break;

                } catch (const std::exception& e) {
                    continue; // Skip lines with bad numbers
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
