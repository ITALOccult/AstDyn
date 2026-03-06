#ifndef ASTDYN_IO_PARSERS_ORBFITEQ1PARSER_HPP
#define ASTDYN_IO_PARSERS_ORBFITEQ1PARSER_HPP

#include "astdyn/io/IParser.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace astdyn::io::parsers {

/**
 * @brief Parser for OrbFit .eq1 files (OEF2.0 format)
 */
class OrbFitEQ1Parser : public IOrbitParser {
public:
    OrbitalElements parse(const std::string& filepath) override {
        std::ifstream file(filepath);
        if (!file.is_open()) throw std::runtime_error("Could not open EQ1 file: " + filepath);

        OrbitalElements elem;
        std::string line;
        
        while (std::getline(file, line)) {
            // Trim leading/trailing whitespace for easier detection
            std::string trimmed = line;
            trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
            trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);

            if (trimmed.starts_with("! Object")) {
                elem.object_name = trimmed.substr(8);
                elem.object_name.erase(0, elem.object_name.find_first_not_of(" "));
            } else if (trimmed.starts_with("EQU")) {
                std::istringstream iss(trimmed.substr(3));
                double a, h, k, p, q, lambda;
                if (iss >> a >> h >> k >> p >> q >> lambda) {
                    elem.semi_major_axis = a;
                    // Equinoctial to Keplerian conversion
                    double e = std::sqrt(h*h + k*k);
                    double i = 2.0 * std::atan(std::sqrt(p*p + q*q));
                    double Omega = std::atan2(p, q);
                    double w_bar = std::atan2(h, k);
                    double omega = w_bar - Omega;
                    double M = lambda - w_bar;

                    // Normalize angles
                    auto norm = [](double ang) {
                        while (ang < 0) ang += 2.0 * M_PI;
                        while (ang >= 2.0 * M_PI) ang -= 2.0 * M_PI;
                        return ang;
                    };

                    elem.eccentricity = e;
                    elem.inclination = i;
                    elem.longitude_asc_node = norm(Omega);
                    elem.argument_perihelion = norm(omega);
                    elem.mean_anomaly = norm(M);
                }
            } else if (trimmed.starts_with("MJD")) {
                std::istringstream iss(trimmed.substr(3));
                iss >> elem.epoch_mjd_tdb;
            } else if (trimmed.starts_with("MAG")) {
                std::istringstream iss(trimmed.substr(3));
                iss >> elem.magnitude >> elem.mag_slope;
            }
        }
        return elem;
    }

    std::string name() const override { return "OrbFit EQ1 Parser (OEF2.0)"; }
    
    bool can_handle(const std::string& filepath) const override {
        return filepath.ends_with(".eq1") || filepath.ends_with(".oel");
    }
};

} // namespace astdyn::io::parsers

#endif // ASTDYN_IO_PARSERS_ORBFITEQ1PARSER_HPP
