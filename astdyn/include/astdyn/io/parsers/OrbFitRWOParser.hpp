#pragma once

#include "astdyn/io/IParser.hpp"
#include <fstream>
#include <sstream>

namespace astdyn {
namespace io {
namespace parsers {

/**
 * @brief Parser for OrbFit .rwo files (optical observations)
 * 
 * Format: objectName MJD_UTC RA DEC mag obscode
 */
class OrbFitRWOParser : public IObservationParser {
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
            // Skip comments and empty lines
            if (line.empty() || line[0] == '!' || line[0] == '#') {
                continue;
            }

            OpticalObservation obs;
            std::istringstream iss(line);

            // Parse: objectName MJD RA DEC mag obscode
            iss >> obs.object_name >> obs.mjd_utc;
            
            // Parse RA (hours:minutes:seconds or decimal degrees)
            double ra_h, ra_m, ra_s;
            char sep1, sep2;
            if (iss.peek() != std::char_traits<char>::eof() && iss >> ra_h) {
                if (iss.peek() == ':') {
                    iss >> sep1 >> ra_m >> sep2 >> ra_s;
                    obs.ra = (ra_h + ra_m / 60.0 + ra_s / 3600.0) * 15.0 * M_PI / 180.0;
                } else {
                    obs.ra = ra_h * M_PI / 180.0;
                }
            }

            // Parse DEC (degrees:arcmin:arcsec or decimal degrees)
            double dec_d, dec_m, dec_s;
            if (iss >> dec_d) {
                if (iss.peek() == ':') {
                    iss >> sep1 >> dec_m >> sep2 >> dec_s;
                    double sign = (dec_d < 0) ? -1.0 : 1.0;
                    obs.dec = sign * (std::abs(dec_d) + dec_m / 60.0 + dec_s / 3600.0) * M_PI / 180.0;
                } else {
                    obs.dec = dec_d * M_PI / 180.0;
                }
            }

            iss >> obs.mag >> obs.obs_code;

            // Default uncertainties (1 arcsec)
            obs.sigma_ra = 1.0;
            obs.sigma_dec = 1.0;

            observations.push_back(obs);
            
            if (max_count > 0 && ++count >= max_count) {
                break;
            }
        }

        if (observations.empty()) {
            throw std::runtime_error("No observations found in file: " + filepath);
        }

        return observations;
    }

    std::string name() const override {
        return "OrbFit RWO Parser";
    }

    bool can_handle(const std::string& filepath) const override {
        return filepath.find(".rwo") != std::string::npos || 
               filepath.find(".RWO") != std::string::npos;
    }
};

} // namespace parsers
} // namespace io
} // namespace astdyn
