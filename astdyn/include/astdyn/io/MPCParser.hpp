/**
 * @file MPCParser.hpp
 * @brief MPC 80-column observation parser.
 */

#ifndef ASTDYN_IO_MPC_PARSER_HPP
#define ASTDYN_IO_MPC_PARSER_HPP

#include "astdyn/observations/Observation.hpp"
#include <string>
#include <vector>

namespace astdyn::io {

/**
 * @brief Utility class to parse MPC (Minor Planet Center) 80-column format.
 */
class MPCParser {
public:
    /**
     * @brief Parse a single MPC line.
     * @param line 80-char line.
     * @return OpticalObservation.
     */
    static observations::OpticalObservation parse_line(const std::string& line);

    /**
     * @brief Parse a multi-line MPC file.
     * @param content File content.
     * @return Vector of observations.
     */
    static std::vector<observations::OpticalObservation> parse_file(const std::string& content);

private:
    /** @brief Helper to parse MPC date (YYYY MM DD.ddddd). */
    static time::EpochUTC parse_date(const std::string& date_str);
    /** @brief Helper to parse MPC RA (HH MM SS.sss). */
    static double parse_ra(const std::string& ra_str);
    /** @brief Helper to parse MPC Dec (DD MM SS.ss). */
    static double parse_dec(const std::string& dec_str);
};

} // namespace astdyn::io

#endif // ASTDYN_IO_MPC_PARSER_HPP
