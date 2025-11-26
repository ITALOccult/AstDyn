#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>

namespace astdyn {
namespace io {

/**
 * @brief Parser for OrbFit .eq1 equinoctial element files (OEF2.0 format)
 */
class EQ1Parser {
public:
    /**
     * @brief Equinoctial orbital elements
     */
    struct EquinoctialElements {
        std::string object_name;    ///< Object designation
        double epoch_mjd_tdb;       ///< Epoch in MJD TDB
        double a;                   ///< Semi-major axis (AU)
        double h;                   ///< e·sin(ϖ) where ϖ = ω + Ω
        double k;                   ///< e·cos(ϖ)
        double p;                   ///< tan(i/2)·sin(Ω)
        double q;                   ///< tan(i/2)·cos(Ω)
        double lambda;              ///< Mean longitude λ = M + ω + Ω (degrees)
        
        // Optional
        double magnitude = 0.0;     ///< Absolute magnitude
        double mag_slope = 0.0;     ///< Magnitude slope parameter
    };

    /**
     * @brief Parse .eq1 file
     * 
     * Supports OEF2.0 format from AstDyS/OrbFit.
     * 
     * @param filepath Path to .eq1 file
     * @return Equinoctial elements
     * @throws std::runtime_error if file cannot be parsed
     */
    static EquinoctialElements parse(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }

        EquinoctialElements elem;
        std::string line;
        bool found_object = false;
        bool found_equ = false;
        bool found_mjd = false;

        while (std::getline(file, line)) {
            // Skip empty lines and comments
            if (line.empty() || line[0] == '!') continue;

            // Skip header
            if (line.find("format") != std::string::npos ||
                line.find("rectype") != std::string::npos ||
                line.find("refsys") != std::string::npos ||
                line.find("END_OF_HEADER") != std::string::npos) {
                continue;
            }

            // Object name (first non-comment, non-header line)
            if (!found_object) {
                std::istringstream iss(line);
                iss >> elem.object_name;
                found_object = true;
                continue;
            }

            // Equinoctial elements line (may start with space)
            if (line.find("EQU") == 0 || line.find(" EQU") == 0) {
                size_t start = line.find("EQU") + 3;
                std::istringstream iss(line.substr(start));
                if (!(iss >> elem.a >> elem.h >> elem.k >> elem.p >> elem.q >> elem.lambda)) {
                    throw std::runtime_error("Failed to parse EQU line: " + line);
                }
                found_equ = true;
                continue;
            }

            // Epoch line (may start with space)
            if (line.find("MJD") == 0 || line.find(" MJD") == 0) {
                size_t start = line.find("MJD") + 3;
                std::istringstream iss(line.substr(start));
                std::string tdt_str;
                if (!(iss >> elem.epoch_mjd_tdb >> tdt_str)) {
                    throw std::runtime_error("Failed to parse MJD line: " + line);
                }
                found_mjd = true;
                continue;
            }

            // Magnitude line (optional, may start with space)
            if (line.find("MAG") == 0 || line.find(" MAG") == 0) {
                size_t start = line.find("MAG") + 3;
                std::istringstream iss(line.substr(start));
                iss >> elem.magnitude >> elem.mag_slope;
                continue;
            }

            // Stop after getting all required data
            if (found_object && found_equ && found_mjd) {
                break;
            }
        }

        // Validate
        if (!found_object) {
            throw std::runtime_error("No object name found in .eq1 file");
        }
        if (!found_equ) {
            throw std::runtime_error("No EQU line found in .eq1 file");
        }
        if (!found_mjd) {
            throw std::runtime_error("No MJD line found in .eq1 file");
        }

        return elem;
    }

    /**
     * @brief Convert equinoctial elements to Keplerian
     * 
     * @param equ Equinoctial elements
     * @return Keplerian elements (angles in radians)
     */
    static void equinoctial_to_keplerian(
        const EquinoctialElements& equ,
        double& a, double& e, double& i, double& Omega, double& omega, double& M) {
        
        const double PI = 3.14159265358979323846;
        const double DEG_TO_RAD = PI / 180.0;

        a = equ.a;

        // Eccentricity
        e = std::sqrt(equ.h * equ.h + equ.k * equ.k);

        // Inclination
        double tan_half_i = std::sqrt(equ.p * equ.p + equ.q * equ.q);
        i = 2.0 * std::atan(tan_half_i);

        // Longitude of ascending node
        Omega = std::atan2(equ.p, equ.q);
        if (Omega < 0) Omega += 2.0 * PI;

        // Argument of perihelion + longitude of ascending node
        double omega_plus_Omega = std::atan2(equ.h, equ.k);
        if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * PI;

        // Argument of perihelion
        omega = omega_plus_Omega - Omega;
        if (omega < 0) omega += 2.0 * PI;

        // Mean longitude to mean anomaly
        double lambda_rad = equ.lambda * DEG_TO_RAD;
        M = lambda_rad - omega_plus_Omega;
        
        // Normalize to [0, 2π)
        while (M < 0) M += 2.0 * PI;
        while (M >= 2.0 * PI) M -= 2.0 * PI;
    }
};

} // namespace io
} // namespace astdyn
