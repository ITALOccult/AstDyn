#ifndef ASTDYN_IO_ASTDYNCONFIG_HPP
#define ASTDYN_IO_ASTDYNCONFIG_HPP

#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/observations/OpticalObservation.hpp"
#include <string>
#include <vector>

namespace astdyn::config {

using namespace astdyn::propagation;
using namespace astdyn::observations;

enum class OrbitalElementSubType {
    OSCULATING,
    MEAN,
    PROPER
};

struct OrbitalElementFile {
    std::string object_name;
    std::string element_format; // e.g. "KEP", "EQ1"
    double epoch_mjd;
    std::string time_scale;     // TDB, UTC
    std::string reference_frame; // ECLM J2000, EQUA J2000
    OrbitalElementSubType element_type;
    KeplerianElements keplerian;
    bool has_covariance = false;
};

struct RWOObservation {
    OpticalObservation observation;
    double weight_ra = 1.0;
    double weight_dec = 1.0;
    double residual_ra = 0.0;
    double residual_dec = 0.0;
    int selection_flag = 1;
    bool outlier = false;
};

class OEFFileHandler {
public:
    static OrbitalElementFile read(const std::string& filename);
    static void write(const std::string& filename, const OrbitalElementFile& oef);
    
    static KeplerianElements meanToOsculating(const KeplerianElements& mean_elements);
    static KeplerianElements osculatingToMean(const KeplerianElements& osc_elements);
    
private:
    static OrbitalElementSubType parseElementType(const std::string& type_str);
    static std::string elementTypeToString(OrbitalElementSubType type);
};

class RWOFileHandler {
public:
    static std::vector<RWOObservation> read(const std::string& filename);
    static void write(const std::string& filename, const std::vector<RWOObservation>& observations);
    
private:
    static RWOObservation parseLine(const std::string& line);
    static std::string formatLine(const RWOObservation& obs);
};

class AstDynConfigManager {
public:
    AstDynConfigManager() = default;
    
    bool loadConfiguration(const std::string& base_path, const std::string& object_name);
    bool saveConfiguration(const std::string& base_path, const std::string& object_name) const;
    
    KeplerianElements getOsculatingElements() const;
    KeplerianElements getOriginalElements() const;
    std::vector<OpticalObservation> getValidObservations() const;
    
    bool exportForLegacyFortran(const std::string& output_dir, const std::string& object_name) const;

private:
    std::string object_name_;
    OrbitalElementFile oef_data_;
    std::vector<RWOObservation> observations_;
    
    bool loadOEF(const std::string& filename);
    bool loadRWO(const std::string& filename);
    bool loadOBS(const std::string& filename);
};

} // namespace astdyn::config

#endif // ASTDYN_IO_ASTDYNCONFIG_HPP
