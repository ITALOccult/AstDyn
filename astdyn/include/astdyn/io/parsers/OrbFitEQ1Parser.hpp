#pragma once

#include "astdyn/io/IParser.hpp"
#include <cmath>
#include <fstream>
#include <sstream>

namespace astdyn {
namespace io {
namespace parsers {

/**
 * @brief Parser for OrbFit .eq1 files (OEF2.0 format)
 * 
 * Handles equinoctial elements format from AstDyS database:
 *  EQU a h k p q lambda
 *  MJD epoch
 *  MAG H G
 */
class OrbFitEQ1Parser : public IOrbitParser {
public:
    OrbitalElements parse(const std::string& filepath) override {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }

        EquinoctialElements eq_elem;
        bool found_equ = false;
        bool found_mjd = false;

        std::string line;
        while (std::getline(file, line)) {
            // Handle leading spaces in OEF2.0 format
            if (line.find("EQU") == 0 || line.find(" EQU") == 0) {
                size_t start = line.find("EQU") + 3;
                std::istringstream iss(line.substr(start));
                iss >> eq_elem.a >> eq_elem.h >> eq_elem.k 
                    >> eq_elem.p >> eq_elem.q >> eq_elem.lambda;
                found_equ = true;
            }
            else if (line.find("MJD") == 0 || line.find(" MJD") == 0) {
                size_t start = line.find("MJD") + 3;
                std::istringstream iss(line.substr(start));
                iss >> eq_elem.epoch_mjd_tdb;
                found_mjd = true;
            }
            else if (line.find("MAG") == 0 || line.find(" MAG") == 0) {
                size_t start = line.find("MAG") + 3;
                std::istringstream iss(line.substr(start));
                iss >> eq_elem.magnitude >> eq_elem.mag_slope;
            }
            else if (line.find("! Object") == 0 || line.find(" ! Object") == 0) {
                // Extract object name from comment line
                size_t obj_pos = line.find("Object");
                if (obj_pos != std::string::npos) {
                    size_t start = obj_pos + 6;
                    while (start < line.length() && std::isspace(line[start])) ++start;
                    size_t end = start;
                    while (end < line.length() && !std::isspace(line[end])) ++end;
                    eq_elem.object_name = line.substr(start, end - start);
                }
            }
        }

        if (!found_equ || !found_mjd) {
            throw std::runtime_error("Invalid .eq1 file: missing EQU or MJD data");
        }

        // Convert equinoctial to Keplerian
        return equinoctial_to_keplerian(eq_elem);
    }

    std::string name() const override {
        return "OrbFit EQ1 Parser (OEF2.0)";
    }

    bool can_handle(const std::string& filepath) const override {
        return filepath.find(".eq1") != std::string::npos || 
               filepath.find(".EQ1") != std::string::npos;
    }

private:
    struct EquinoctialElements {
        std::string object_name;
        double epoch_mjd_tdb;
        double a, h, k, p, q, lambda;
        double magnitude = 0.0;
        double mag_slope = 0.0;
    };

    OrbitalElements equinoctial_to_keplerian(const EquinoctialElements& eq) {
        OrbitalElements kep;
        kep.object_name = eq.object_name;
        kep.epoch_mjd_tdb = eq.epoch_mjd_tdb;
        kep.magnitude = eq.magnitude;
        kep.mag_slope = eq.mag_slope;

        // Conversions from equinoctial to Keplerian
        kep.semi_major_axis = eq.a;
        
        // Eccentricity: e = sqrt(h^2 + k^2)
        kep.eccentricity = std::sqrt(eq.h * eq.h + eq.k * eq.k);
        
        // Inclination: i = 2 * atan(sqrt(p^2 + q^2))
        kep.inclination = 2.0 * std::atan(std::sqrt(eq.p * eq.p + eq.q * eq.q));
        
        // Longitude of ascending node: Ω = atan2(p, q)  [CORRECTED]
        kep.longitude_asc_node = std::atan2(eq.p, eq.q);
        if (kep.longitude_asc_node < 0) kep.longitude_asc_node += 2.0 * M_PI;
        
        // Argument of perihelion: ω = atan2(h, k) - Ω  [CORRECTED]
        double omega_plus_Omega = std::atan2(eq.h, eq.k);
        if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * M_PI;
        kep.argument_perihelion = omega_plus_Omega - kep.longitude_asc_node;
        if (kep.argument_perihelion < 0) kep.argument_perihelion += 2.0 * M_PI;
        
        // Mean anomaly: M = lambda - atan2(k, h)
        double lambda_rad = eq.lambda * M_PI / 180.0;
        kep.mean_anomaly = lambda_rad - omega_plus_Omega;
        if (kep.mean_anomaly < 0) kep.mean_anomaly += 2.0 * M_PI;

        return kep;
    }
};

} // namespace parsers
} // namespace io
} // namespace astdyn
