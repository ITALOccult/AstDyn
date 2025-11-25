/**
 * @file MPCReader.hpp
 * @brief Parser for MPC 80-column observation format
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Reads astrometric observations in MPC's classic 80-column format.
 * 
 * Format specification:
 * https://minorplanetcenter.net/iau/info/OpticalObs.html
 * 
 * Example line:
 * 00433         C2023 01 15.41667 10 34 23.45 +19 40 25.8          17.5 V      568
 * 
 * Columns:
 *  1-5:   Packed provisional designation or number
 *  6-12:  Discovery (*), note1
 *  13:    Note2 (program code)
 *  16-32: Date (YYYY MM DD.ddddd)
 *  33-44: RA (HH MM SS.ddd)
 *  45-56: Dec (sDD MM SS.dd)
 *  66-70: Magnitude
 *  71:    Band
 *  78-80: Observatory code
 */

#ifndef ASTDYN_OBSERVATIONS_MPCREADER_HPP
#define ASTDYN_OBSERVATIONS_MPCREADER_HPP

#include "astdyn/observations/Observation.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <map>

namespace astdyn {
namespace observations {

/**
 * @brief MPC 80-column observation reader
 */
class MPCReader {
public:
    /**
     * @brief Read observations from MPC file
     * @param filepath Path to MPC observations file
     * @return Vector of observations
     */
    static std::vector<OpticalObservation> readFile(const std::string& filepath);
    
    /**
     * @brief Parse a single MPC observation line
     * @param line 80-column MPC observation line
     * @return OpticalObservation, or nullopt if parsing fails
     */
    static std::optional<OpticalObservation> parseLine(const std::string& line);
    
    /**
     * @brief Read observations into ObservationSet grouped by object
     * @param filepath Path to MPC file
     * @return Map of object designation to observation set
     */
    static std::map<std::string, ObservationSet> readFileGrouped(const std::string& filepath);

private:
    /**
     * @brief Parse packed designation (columns 1-12)
     * 
     * MPC uses packed format for both numbers and provisional designations.
     * Examples:
     * - "00433       " → "(433)" Eros
     * - "K14A00A     " → "2014 AA"
     * - "    PLS1193 " → "P/1993 A1"
     */
    static std::string parseDesignation(const std::string& packed);
    
    /**
     * @brief Parse date (columns 16-32)
     * Format: "YYYY MM DD.ddddd"
     * @return MJD in UTC
     */
    static double parseDate(const std::string& date_str);
    
    /**
     * @brief Parse RA (columns 33-44)
     * Format: "HH MM SS.ddd"
     * @return RA in radians
     */
    static double parseRA(const std::string& ra_str);
    
    /**
     * @brief Parse Dec (columns 45-56)
     * Format: "sDD MM SS.dd"
     * @return Dec in radians
     */
    static double parseDec(const std::string& dec_str);
    
    /**
     * @brief Parse magnitude (columns 65-69)
     */
    static std::optional<double> parseMagnitude(const std::string& mag_str);
    
    /**
     * @brief Parse note1 field (column 14)
     * Discovery asterisk, etc.
     */
    static std::string parseNote1(char c);
    
    /**
     * @brief Parse note2 field (column 15)
     * Program code
     */
    static std::string parseNote2(char c);
    
    /**
     * @brief Unpack minor planet number
     * Examples: "00433" → 433, "A0001" → 100001
     */
    static int unpackNumber(const std::string& packed);
    
    /**
     * @brief Unpack provisional designation
     * Examples: "K14A00A" → "2014 AA"
     */
    static std::string unpackProvisional(const std::string& packed);
};

/**
 * @brief Write observations to MPC 80-column format
 */
class MPCWriter {
public:
    /**
     * @brief Write observations to file
     * @param filepath Output file path
     * @param observations Observations to write
     */
    static void writeFile(
        const std::string& filepath,
        const std::vector<OpticalObservation>& observations
    );
    
    /**
     * @brief Format single observation as MPC 80-column line
     */
    static std::string formatLine(const OpticalObservation& obs);

private:
    /**
     * @brief Pack designation to MPC format
     */
    static std::string packDesignation(const std::string& designation);
    
    /**
     * @brief Format date to MPC format
     */
    static std::string formatDate(double mjd_utc);
    
    /**
     * @brief Format RA to MPC format
     */
    static std::string formatRA(double ra_rad);
    
    /**
     * @brief Format Dec to MPC format
     */
    static std::string formatDec(double dec_rad);
};

} // namespace observations
} // namespace astdyn

#endif // ASTDYN_OBSERVATIONS_MPCREADER_HPP
