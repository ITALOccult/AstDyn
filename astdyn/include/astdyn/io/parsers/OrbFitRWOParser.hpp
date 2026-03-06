#ifndef ASTDYN_IO_PARSERS_ORBFITRWOPARSER_HPP
#define ASTDYN_IO_PARSERS_ORBFITRWOPARSER_HPP

#include "astdyn/io/IParser.hpp"
#include <stdexcept>

namespace astdyn::io::parsers {

/**
 * @brief Parser for OrbFit .rwo files (Stub)
 */
class OrbFitRWOParser : public IObservationParser {
public:
    std::vector<OpticalObservation> parse(const std::string& filepath, size_t max_count = 0) override {
        // Implementation would normally use RWOReader
        throw std::runtime_error("OrbFitRWOParser::parse not fully implemented yet. Please use legacy RWOReader for now.");
        return {};
    }

    std::string name() const override { return "OrbFit RWO Parser"; }
    
    bool can_handle(const std::string& filepath) const override {
        return filepath.ends_with(".rwo");
    }
};

} // namespace astdyn::io::parsers

#endif // ASTDYN_IO_PARSERS_ORBFITRWOPARSER_HPP
