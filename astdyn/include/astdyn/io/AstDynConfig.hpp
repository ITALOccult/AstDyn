/**
 * @file AstDynConfig.hpp
 * @brief OrbFit configuration file parser and manager
 * 
 * Handles configuration compatible with Fortran orbit fitting software:
 * - .opt files (options)
 * - .key files (keywords)
 * - .rwo files (observations with weights)
 * - .oef files (orbital element files)
 * - Mean → Osculating element conversion
 * 
 * Reference: OrbFit Software (Milani et al.)
 */

#ifndef ASTDYN_CONFIG_PARSER_HPP
#define ASTDYN_CONFIG_PARSER_HPP

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <fstream>
#include "astdyn/core/Types.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/observations/Observation.hpp"

namespace astdyn::config {

/**
 * @brief Orbital element subtype (MEAN/OSCULATING/PROPER)
 * 
 * Distinct from astdyn::ElementType which indicates coordinate system.
 * This enum specifies whether Keplerian elements are instantaneous,
 * averaged, or proper elements.
 */
enum class OrbitalElementSubType {
    MEAN,          ///< Mean elements (with J2 secular terms removed)
    OSCULATING,    ///< Osculating (instantaneous) elements
    PROPER         ///< Proper elements (long-term averaged)
};

/**
 * @brief OEF (Orbital Element File) format
 * 
 * Compatible with OrbFit Fortran format:
 * Line 1: Object name
 * Line 2: Element type (KEP/EQU/CAR), Epoch (MJD), Reference system
 * Line 3-8: Elements and covariance (if available)
 */
struct OrbitalElementFile {
    std::string object_name;
    std::string element_format;      ///< KEP, EQU, CAR, COM
    OrbitalElementSubType element_type;  ///< MEAN, OSCULATING, PROPER
    double epoch_mjd;
    std::string time_scale;          ///< TDT, UTC, TDB
    std::string reference_frame;     ///< ECLM J2000, EQUA J2000
    
    // Elements (format-dependent)
    propagation::KeplerianElements keplerian;
    
    // Optional covariance
    bool has_covariance;
    Eigen::Matrix<double, 6, 6> covariance;
    
    // Optional proper elements info
    std::map<std::string, double> metadata;
};

/**
 * @brief RWO (Residuals and Weights for Observations) format
 * 
 * Extended MPC format with weights and residuals:
 * - Standard MPC 80-column format
 * - Additional fields: weights, residuals, flags
 */
struct RWOObservation {
    observations::OpticalObservation observation;
    
    // Weights
    double weight_ra;   ///< Weight for RA (0=rejected)
    double weight_dec;  ///< Weight for Dec (0=rejected)
    
    // Residuals (if computed)
    double residual_ra;  ///< O-C in RA [arcsec]
    double residual_dec; ///< O-C in Dec [arcsec]
    
    // Flags
    bool outlier;
    char selection_flag; ///< S=selected, R=rejected, F=forced
    
    RWOObservation() 
        : weight_ra(1.0), weight_dec(1.0),
          residual_ra(0.0), residual_dec(0.0),
          outlier(false), selection_flag('S') {}
};

// OptionFileParser removed - usage replaced by JSON


/**
 * @brief OEF file reader/writer
 */
class OEFFileHandler {
public:
    /**
     * @brief Read orbital elements from OEF file
     */
    static OrbitalElementFile read(const std::string& filename);
    
    /**
     * @brief Write orbital elements to OEF file
     */
    static void write(const std::string& filename, const OrbitalElementFile& oef);
    
    /**
     * @brief Convert mean elements to osculating
     * 
     * Uses J2 theory to convert mean elements (with secular terms removed)
     * back to osculating elements at the same epoch.
     * 
     * @param mean_elements Mean Keplerian elements
     * @return Osculating Keplerian elements
     */
    static propagation::KeplerianElements meanToOsculating(
        const propagation::KeplerianElements& mean_elements);
    
    /**
     * @brief Convert osculating elements to mean
     */
    static propagation::KeplerianElements osculatingToMean(
        const propagation::KeplerianElements& osc_elements);
    
private:
    static OrbitalElementSubType parseElementType(const std::string& type_str);
    static std::string elementTypeToString(OrbitalElementSubType type);
};

/**
 * @brief RWO file reader/writer
 */
class RWOFileHandler {
public:
    /**
     * @brief Read observations from RWO file
     */
    static std::vector<RWOObservation> read(const std::string& filename);
    
    /**
     * @brief Write observations to RWO file
     */
    static void write(const std::string& filename, 
                     const std::vector<RWOObservation>& observations);
    
    /**
     * @brief Parse single RWO line (80-column MPC + extensions)
     */
    static RWOObservation parseLine(const std::string& line);
    
    /**
     * @brief Format RWO line for output
     */
    static std::string formatLine(const RWOObservation& obs);
};

/**
 * @brief Complete OrbFit configuration manager
 * 
 * Manages all configuration aspects:
 * - Option files (.opt)
 * - Orbital elements (.oef)
 * - Observations (.rwo or .obs)
 * - Mean/Osculating conversion
 */
class AstDynConfigManager {
public:
    AstDynConfigManager() = default;
    
    /**
     * @brief Load complete configuration from directory
     * 
     * Looks for:
     * - <name>.opt (options)
     * - <name>.oef (orbital elements)
     * - <name>.rwo or <name>.obs (observations)
     */
    bool loadConfiguration(const std::string& base_path, const std::string& object_name);
    
    /**
     * @brief Get orbital elements (automatically converts MEAN → OSCULATING)
     */
    propagation::KeplerianElements getOsculatingElements() const;
    
    /**
     * @brief Get original elements (as stored in file)
     */
    propagation::KeplerianElements getOriginalElements() const;
    
    /**
     * @brief Get element type from file
     */
    OrbitalElementSubType getElementType() const { return oef_data_.element_type; }
    
    /**
     * @brief Get observations with weights
     */
    const std::vector<RWOObservation>& getObservations() const { return observations_; }
    
    /**
     * @brief Get filtered observations (weight > 0)
     */
    std::vector<observations::OpticalObservation> getValidObservations() const;
    
    /**
     * @brief Save current configuration
     */
    bool saveConfiguration(const std::string& base_path, const std::string& object_name) const;

    /**
     * @brief Get option value
     */
    bool exportForLegacyFortran(const std::string& output_dir, 
                                const std::string& object_name) const;
    
private:
    // OptionFileParser options_; // Removed
    OrbitalElementFile oef_data_;
    std::vector<RWOObservation> observations_;
    std::string object_name_;
    
    bool loadOEF(const std::string& filename);
    bool loadRWO(const std::string& filename);
    bool loadOBS(const std::string& filename);
};

} // namespace astdyn::config

#endif // ASTDYN_CONFIG_PARSER_HPP
