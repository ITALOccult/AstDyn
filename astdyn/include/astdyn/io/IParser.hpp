#pragma once

#include <string>
#include <memory>
#include <vector>

namespace astdyn {
namespace io {

/**
 * @brief Base interface for orbital elements parsers
 * 
 * This allows switching between different input formats (OrbFit .eq1, MPC, etc.)
 */
class IOrbitParser {
public:
    /**
     * @brief Orbital elements structure (format-agnostic)
     */
    struct OrbitalElements {
        std::string object_name;
        double epoch_mjd_tdb;
        
        // Keplerian elements (in SI/standard units)
        double semi_major_axis;        ///< AU
        double eccentricity;           ///< dimensionless
        double inclination;            ///< radians
        double longitude_asc_node;     ///< radians (Ω)
        double argument_perihelion;    ///< radians (ω)
        double mean_anomaly;           ///< radians (M)
        
        // Optional
        double magnitude = 0.0;
        double mag_slope = 0.0;
    };

    virtual ~IOrbitParser() = default;

    /**
     * @brief Parse orbital elements from file
     * 
     * @param filepath Path to input file
     * @return Parsed orbital elements
     * @throws std::runtime_error if parsing fails
     */
    virtual OrbitalElements parse(const std::string& filepath) = 0;

    /**
     * @brief Get parser name/description
     */
    virtual std::string name() const = 0;

    /**
     * @brief Check if this parser can handle the file (based on extension/format)
     */
    virtual bool can_handle(const std::string& filepath) const = 0;
};

/**
 * @brief Base interface for observation parsers
 */
class IObservationParser {
public:
    /**
     * @brief Optical observation structure
     */
    struct OpticalObservation {
        std::string object_name;
        double mjd_utc;
        double ra;              ///< Right ascension (radians)
        double dec;             ///< Declination (radians)
        double mag;             ///< Apparent magnitude
        std::string obs_code;   ///< Observatory code
        
        // Uncertainties
        double sigma_ra = 1.0;  ///< arcsec
        double sigma_dec = 1.0; ///< arcsec
    };

    virtual ~IObservationParser() = default;

    /**
     * @brief Parse observations from file
     * 
     * @param filepath Path to input file
     * @param max_count Maximum observations to read (0 = all)
     * @return Vector of observations
     * @throws std::runtime_error if parsing fails
     */
    virtual std::vector<OpticalObservation> parse(
        const std::string& filepath, 
        size_t max_count = 0) = 0;

    /**
     * @brief Get parser name/description
     */
    virtual std::string name() const = 0;

    /**
     * @brief Check if this parser can handle the file
     */
    virtual bool can_handle(const std::string& filepath) const = 0;
};

/**
 * @brief Factory for creating parsers
 */
class ParserFactory {
public:
    /**
     * @brief Create appropriate orbit parser based on file extension/format
     * 
     * @param filepath Path to file (used to detect format)
     * @param format_hint Optional format hint ("eq1", "mpc", "oel", etc.)
     * @return Unique pointer to parser
     */
    static std::unique_ptr<IOrbitParser> create_orbit_parser(
        const std::string& filepath,
        const std::string& format_hint = "auto");

    /**
     * @brief Create appropriate observation parser based on file extension/format
     * 
     * @param filepath Path to file
     * @param format_hint Optional format hint ("rwo", "mpc", etc.)
     * @return Unique pointer to parser
     */
    static std::unique_ptr<IObservationParser> create_observation_parser(
        const std::string& filepath,
        const std::string& format_hint = "auto");

private:
    static std::string detect_format(const std::string& filepath);
};

} // namespace io
} // namespace astdyn
