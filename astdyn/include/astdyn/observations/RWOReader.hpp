/**
 * @file RWOReader.hpp
 * @brief Reader for OrbFit .rwo (residual/observation) format
 * @author ITALOccult AstDyn Team
 * @date 2025-11-25
 */

#pragma once

#include "astdyn/observations/Observation.hpp"
#include <string>
#include <vector>
#include <optional>
#include <map>

namespace astdyn {
namespace observations {

/**
 * @brief Parser for OrbFit .rwo format
 * 
 * The .rwo format contains observations with residuals and weights.
 * Format (columns approximate):
 * - Object name (1-14)
 * - Obs type 'O' (15-16) 
 * - Note/status 'C' (17-18)
 * - Date YYYY MM DD.ddddd (19-36)
 * - RA error (38-47)
 * - RA HH MM SS.ddd (48-60)
 * - RA weight (62-71)
 * - [more fields...]
 * - Dec +DD MM SS.dd (86-99)
 * - Dec weight (101-110)
 * - [magnitude, observatory, etc.]
 */
class RWOReader {
public:
    /**
     * @brief Read observations from .rwo file
     * @param filepath Path to .rwo file
     * @return Vector of optical observations
     */
    static std::vector<OpticalObservation> readFile(const std::string& filepath);
    
    /**
     * @brief Parse single line from .rwo file
     * @param line Line text
     * @return Observation if successfully parsed
     */
    static std::optional<OpticalObservation> parseLine(const std::string& line);
    
    /**
     * @brief Read and group observations by object designation
     * @param filepath Path to .rwo file
     * @return Map of object designation to observation sets
     */
    static std::map<std::string, ObservationSet> readFileGrouped(const std::string& filepath);

private:
    /**
     * @brief Parse object designation
     */
    static std::string parseDesignation(const std::string& line);
    
    /**
     * @brief Parse observation date
     */
    static double parseDate(const std::string& date_str);
    
    /**
     * @brief Parse RA from "HH MM SS.ddd" format
     */
    static double parseRA(const std::string& ra_str);
    
    /**
     * @brief Parse Dec from "sDD MM SS.dd" format
     */
    static double parseDec(const std::string& dec_str);
    
    /**
     * @brief Parse observatory code (usually 3 chars near end)
     */
    static std::string parseObservatory(const std::string& line);
    
    /**
     * @brief Parse magnitude if present
     */
    static std::optional<double> parseMagnitude(const std::string& line);
};

} // namespace observations
} // namespace astdyn
