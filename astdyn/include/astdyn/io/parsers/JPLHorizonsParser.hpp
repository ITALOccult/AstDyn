#pragma once

#include "astdyn/io/IParser.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <optional>

namespace astdyn {
namespace io {
namespace parsers {

/**
 * @brief Parser for JPL HORIZONS orbital elements files or streams.
 * Parses standard tabular text format from JPL Horizons.
 */
class JPLHorizonsParser : public IOrbitParser {
public:
    OrbitalElements parse_stream(std::istream& stream) override {
        OrbitalElements elem;
        
        // Default values
        elem.object_name = "JPL_Object";
        elem.magnitude = 0.0;
        elem.mag_slope = 0.15;
        
        std::string line;
        bool found_epoch = false;
        bool found_a = false, found_e = false, found_i = false;
        bool found_om = false, found_w = false, found_ma = false;

        auto extract_value = [](const std::string& str, const std::string& key) -> std::optional<double> {
            size_t pos = str.find(key);
            if (pos != std::string::npos) {
                pos += key.length();
                // Skip '=' if present specifically (e.g. if key is just "EPOCH")
                while (pos < str.length() && (str[pos] == '=' || str[pos] == ' ')) {
                    pos++;
                }
                if (pos < str.length()) {
                    try {
                        return std::stod(str.substr(pos));
                    } catch (...) {}
                }
            }
            return std::nullopt;
        };

        while (std::getline(stream, line)) {
            // Find Target Body name if present (usually inside "Target body name: 203 Pompeja (A879 UA)")
            if (line.find("Target body name:") != std::string::npos) {
                size_t start = line.find(":") + 1;
                size_t end = line.find("(");
                if (end == std::string::npos) end = line.length();
                if (start < end) {
                    std::string name = line.substr(start, end - start);
                    // trim
                    size_t f = name.find_first_not_of(" \t");
                    size_t l = name.find_last_not_of(" \t");
                    if (f != std::string::npos) elem.object_name = name.substr(f, l - f + 1);
                }
            }

            // Parse EPOCH (JD)
            if (line.find("EPOCH=") != std::string::npos || line.find("Epoch") != std::string::npos) {
                if (auto val = extract_value(line, "EPOCH=")) {
                    elem.epoch_mjd_tdb = *val - 2400000.5; // JD to MJD
                    found_epoch = true;
                } else if (auto val2 = extract_value(line, "Epoch ")) {
                    // sometimes "Epoch 2460000.5"
                    elem.epoch_mjd_tdb = *val2 - 2400000.5;
                    found_epoch = true;
                }
            } else if (!found_epoch && line.find("$$SOE") != std::string::npos) {
                // Next line usually starts with the epoch JD
                std::string next_line;
                if (std::getline(stream, next_line)) {
                    try {
                        size_t end_num = next_line.find_first_of(" =");
                        double jd = std::stod(next_line.substr(0, end_num));
                        if (jd > 2000000) {
                            elem.epoch_mjd_tdb = jd - 2400000.5;
                            found_epoch = true;
                        }
                        line = next_line; // Might contain elements on same line 
                    } catch(...) {}
                }
            } else if (!found_epoch && std::isdigit(line[0]) && line.find(" = A.D.") != std::string::npos) {
                try {
                    double jd = std::stod(line);
                    if (jd > 2000000) {
                        elem.epoch_mjd_tdb = jd - 2400000.5;
                        found_epoch = true;
                    }
                } catch(...) {}
            }

            // Parse EC (Eccentricity)
            if (auto val = extract_value(line, "EC=")) { elem.eccentricity = *val; found_e = true; }
            
            // Parse IN (Inclination) in degrees -> radians
            if (auto val = extract_value(line, "IN=")) { elem.inclination = *val * M_PI / 180.0; found_i = true; }
            
            // Parse OM (Longitude of Ascending Node) -> radians
            if (auto val = extract_value(line, "OM=")) { elem.longitude_asc_node = *val * M_PI / 180.0; found_om = true; }
            
            // Parse W (Argument of Perihelion) -> radians
            if (auto val = extract_value(line, "W =")) { elem.argument_perihelion = *val * M_PI / 180.0; found_w = true; }
            else if (auto val = extract_value(line, "W=")) { elem.argument_perihelion = *val * M_PI / 180.0; found_w = true; }
            
            // Parse MA (Mean Anomaly) -> radians
            if (auto val = extract_value(line, "MA=")) { elem.mean_anomaly = *val * M_PI / 180.0; found_ma = true; }
            
            // Parse A (Semi-major axis) in AU
            if (auto val = extract_value(line, "A =")) { elem.semi_major_axis = *val; found_a = true; }
            else if (auto val = extract_value(line, "A=")) { elem.semi_major_axis = *val; found_a = true; }
        }

        if (!found_epoch || !found_a || !found_e || !found_i || !found_om || !found_w || !found_ma) {
            throw std::runtime_error("Invalid JPL Horizons stream: missing essential orbital parameters.");
        }

        return elem;
    }

    OrbitalElements parse(const std::string& filepath) override {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }
        return parse_stream(file);
    }

    std::string name() const override {
        return "JPL HORIZONS Parser";
    }

    bool can_handle(const std::string& filepath) const override {
        std::string lower = filepath;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        return lower.find("horizons") != std::string::npos || lower.find(".txt") != std::string::npos;
    }
};

} // namespace parsers
} // namespace io
} // namespace astdyn
