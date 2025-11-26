#include "astdyn/io/IParser.hpp"
#include "astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "astdyn/io/parsers/OrbFitRWOParser.hpp"
#include <algorithm>
#include <cctype>

namespace astdyn {
namespace io {

std::string ParserFactory::detect_format(const std::string& filepath) {
    // Convert to lowercase for case-insensitive comparison
    std::string lower_path = filepath;
    std::transform(lower_path.begin(), lower_path.end(), lower_path.begin(), 
                   [](unsigned char c) { return std::tolower(c); });

    // Detect by extension
    if (lower_path.find(".eq1") != std::string::npos) return "eq1";
    if (lower_path.find(".rwo") != std::string::npos) return "rwo";
    if (lower_path.find(".oel") != std::string::npos) return "oel";
    if (lower_path.find(".mpc") != std::string::npos) return "mpc";

    return "unknown";
}

std::unique_ptr<IOrbitParser> ParserFactory::create_orbit_parser(
    const std::string& filepath,
    const std::string& format_hint) {
    
    std::string format = (format_hint == "auto") ? detect_format(filepath) : format_hint;

    // Convert to lowercase
    std::string lower_format = format;
    std::transform(lower_format.begin(), lower_format.end(), lower_format.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (lower_format == "eq1") {
        return std::make_unique<parsers::OrbFitEQ1Parser>();
    }
    // Future: add more parsers here
    // else if (lower_format == "mpc") {
    //     return std::make_unique<parsers::MPCOrbitParser>();
    // }
    // else if (lower_format == "oel") {
    //     return std::make_unique<parsers::OELParser>();
    // }

    throw std::runtime_error("Unsupported orbit format: " + format + 
                           " for file: " + filepath);
}

std::unique_ptr<IObservationParser> ParserFactory::create_observation_parser(
    const std::string& filepath,
    const std::string& format_hint) {
    
    std::string format = (format_hint == "auto") ? detect_format(filepath) : format_hint;

    // Convert to lowercase
    std::string lower_format = format;
    std::transform(lower_format.begin(), lower_format.end(), lower_format.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (lower_format == "rwo") {
        return std::make_unique<parsers::OrbFitRWOParser>();
    }
    // Future: add more parsers here
    // else if (lower_format == "mpc") {
    //     return std::make_unique<parsers::MPCObservationParser>();
    // }

    throw std::runtime_error("Unsupported observation format: " + format + 
                           " for file: " + filepath);
}

} // namespace io
} // namespace astdyn
