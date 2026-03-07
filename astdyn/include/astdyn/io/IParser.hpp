#ifndef ASTDYN_IO_IPARSER_HPP
#define ASTDYN_IO_IPARSER_HPP

#include "astdyn/time/epoch.hpp"
#include <string>
#include <vector>
#include <memory>

namespace astdyn {
namespace io {

/**
 * @brief Interface for orbital element parsers
 */
class IOrbitParser {
public:
    struct OrbitalElements {
        std::string object_name;
        time::EpochTDB epoch;
        double semi_major_axis;      // AU
        double eccentricity;
        double inclination;          // radianti
        double longitude_asc_node;   // Ω (radianti)
        double argument_perihelion;  // ω (radianti)
        double mean_anomaly;         // M (radianti)
        double magnitude = 0;
        double mag_slope = 0.15;
    };

    virtual ~IOrbitParser() = default;
    virtual OrbitalElements parse(const std::string& filepath) = 0;
    virtual std::string name() const = 0;
    virtual bool can_handle(const std::string& filepath) const = 0;
};

/**
 * @brief Interface for observation parsers
 */
class IObservationParser {
public:
    struct OpticalObservation {
        std::string object_name;
        time::EpochUTC time;
        double ra, dec;          // radianti
        double mag;
        std::string obs_code;
        double sigma_ra, sigma_dec;  // arcsec
    };

    virtual ~IObservationParser() = default;
    virtual std::vector<OpticalObservation> parse(const std::string& filepath, size_t max_count = 0) = 0;
    virtual std::string name() const = 0;
    virtual bool can_handle(const std::string& filepath) const = 0;
};

/**
 * @brief Factory for creating appropriate parsers
 */
class ParserFactory {
public:
    static std::string detect_format(const std::string& filepath);
    
    static std::unique_ptr<IOrbitParser> create_orbit_parser(
        const std::string& filepath,
        const std::string& format_hint = "auto");

    static std::unique_ptr<IObservationParser> create_observation_parser(
        const std::string& filepath,
        const std::string& format_hint = "auto");
};

} // namespace io
} // namespace astdyn

#endif // ASTDYN_IO_IPARSER_HPP
