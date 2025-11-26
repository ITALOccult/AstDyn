/**
 * @file AstDynConfig.cpp
 * @brief Implementation of AstDyn configuration parsers
 */

#include "astdyn/io/AstDynConfig.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/observations/MPCReader.hpp"
#include "astdyn/core/Constants.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace astdyn::config {

using namespace astdyn::constants;
using namespace astdyn::propagation;
using namespace astdyn::observations;

// ============================================================================
// OptionFileParser Implementation
// ============================================================================

bool OptionFileParser::load(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        line = trim(line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#' || line[0] == '!') {
            continue;
        }
        
        // Section header: [section]
        if (line[0] == '[' && line.back() == ']') {
            current_section_ = line.substr(1, line.length() - 2);
            continue;
        }
        
        // Parse KEY = VALUE
        auto [key, value] = parseLine(line);
        if (!key.empty()) {
            std::string full_key = current_section_.empty() ? key : 
                                   current_section_ + "." + key;
            options_[full_key] = value;
        }
    }
    
    return true;
}

std::string OptionFileParser::getString(const std::string& key, 
                                       const std::string& default_val) const {
    auto it = options_.find(key);
    return (it != options_.end()) ? it->second : default_val;
}

double OptionFileParser::getDouble(const std::string& key, double default_val) const {
    auto it = options_.find(key);
    if (it == options_.end()) return default_val;
    try {
        return std::stod(it->second);
    } catch (...) {
        return default_val;
    }
}

int OptionFileParser::getInt(const std::string& key, int default_val) const {
    auto it = options_.find(key);
    if (it == options_.end()) return default_val;
    try {
        return std::stoi(it->second);
    } catch (...) {
        return default_val;
    }
}

bool OptionFileParser::getBool(const std::string& key, bool default_val) const {
    auto it = options_.find(key);
    if (it == options_.end()) return default_val;
    
    std::string val = it->second;
    std::transform(val.begin(), val.end(), val.begin(), ::tolower);
    
    return (val == "true" || val == "yes" || val == "1" || val == "on");
}

bool OptionFileParser::hasKey(const std::string& key) const {
    return options_.find(key) != options_.end();
}

void OptionFileParser::set(const std::string& key, const std::string& value) {
    options_[key] = value;
}

bool OptionFileParser::save(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    for (const auto& [key, value] : options_) {
        file << key << " = " << value << "\n";
    }
    
    return true;
}

std::string OptionFileParser::trim(const std::string& str) const {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::pair<std::string, std::string> OptionFileParser::parseLine(const std::string& line) const {
    size_t eq_pos = line.find('=');
    if (eq_pos == std::string::npos) {
        return {"", ""};
    }
    
    std::string key = trim(line.substr(0, eq_pos));
    std::string value_part = line.substr(eq_pos + 1);
    
    // Remove comments (starting with ! or #)
    size_t comment_pos = value_part.find_first_of("!#");
    if (comment_pos != std::string::npos) {
        value_part = value_part.substr(0, comment_pos);
    }
    
    std::string value = trim(value_part);
    
    return {key, value};
}

// ============================================================================
// OEFFileHandler Implementation
// ============================================================================

OrbitalElementFile OEFFileHandler::read(const std::string& filename) {
    OrbitalElementFile oef;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open OEF file: " + filename);
    }
    
    std::string line;
    
    // Line 1: Object name
    std::getline(file, line);
    oef.object_name = line;
    
    // Line 2: Format, Epoch, Type
    // Example: KEP 60310.0 TDB ECLM J2000 MEAN
    std::getline(file, line);
    std::istringstream iss(line);
    std::string ref2, type_str;  // J2000 is split into "ECLM" and "J2000" or "EQUA" and "J2000"
    iss >> oef.element_format >> oef.epoch_mjd >> oef.time_scale 
        >> oef.reference_frame >> ref2 >> type_str;
    
    oef.reference_frame = oef.reference_frame + " " + ref2;  // Reconstruct "ECLM J2000"
    oef.element_type = parseElementType(type_str);
    
    // Lines 3-8: Elements (KEP format: a, e, i, Omega, omega, M)
    if (oef.element_format == "KEP") {
        double a, e, i, Omega, omega, M;
        file >> a >> e >> i >> Omega >> omega >> M;
        
        oef.keplerian.epoch_mjd_tdb = oef.epoch_mjd;
        oef.keplerian.semi_major_axis = a;
        oef.keplerian.eccentricity = e;
        oef.keplerian.inclination = i * DEG_TO_RAD;
        oef.keplerian.longitude_ascending_node = Omega * DEG_TO_RAD;
        oef.keplerian.argument_perihelion = omega * DEG_TO_RAD;
        oef.keplerian.mean_anomaly = M * DEG_TO_RAD;
        oef.keplerian.gravitational_parameter = GMS;
    }
    
    oef.has_covariance = false;
    
    return oef;
}

void OEFFileHandler::write(const std::string& filename, const OrbitalElementFile& oef) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot write OEF file: " + filename);
    }
    
    // Line 1: Object name
    file << oef.object_name << "\n";
    
    // Line 2: Format info
    file << oef.element_format << " "
         << std::fixed << std::setprecision(6) << oef.epoch_mjd << " "
         << oef.time_scale << " "
         << oef.reference_frame << " "
         << elementTypeToString(oef.element_type) << "\n";
    
    // Elements
    if (oef.element_format == "KEP") {
        file << std::setprecision(9) << oef.keplerian.semi_major_axis << "\n"
             << std::setprecision(9) << oef.keplerian.eccentricity << "\n"
             << std::setprecision(6) << (oef.keplerian.inclination * RAD_TO_DEG) << "\n"
             << std::setprecision(6) << (oef.keplerian.longitude_ascending_node * RAD_TO_DEG) << "\n"
             << std::setprecision(6) << (oef.keplerian.argument_perihelion * RAD_TO_DEG) << "\n"
             << std::setprecision(6) << (oef.keplerian.mean_anomaly * RAD_TO_DEG) << "\n";
    }
}

KeplerianElements OEFFileHandler::meanToOsculating(const KeplerianElements& mean_elements) {
    // Use the mean-to-osculating converter with Sun's J2
    // For heliocentric orbits, J2_sun ≈ 2e-7 is negligible but we apply it for consistency
    
    // Convert mean → osculating
    constexpr double j2_sun = 2e-7;
    return mean_to_osculating(mean_elements, j2_sun);
}

KeplerianElements OEFFileHandler::osculatingToMean(const KeplerianElements& osc_elements) {
    constexpr double j2_sun = 2e-7;
    return osculating_to_mean(osc_elements, j2_sun);
}

OrbitalElementSubType OEFFileHandler::parseElementType(const std::string& type_str) {
    if (type_str == "MEAN" || type_str == "MEA") {
        return OrbitalElementSubType::MEAN;
    } else if (type_str == "PROPER" || type_str == "PROP") {
        return OrbitalElementSubType::PROPER;
    }
    return OrbitalElementSubType::OSCULATING;
}

std::string OEFFileHandler::elementTypeToString(OrbitalElementSubType type) {
    switch (type) {
        case OrbitalElementSubType::MEAN: return "MEAN";
        case OrbitalElementSubType::PROPER: return "PROPER";
        default: return "OSC";
    }
}

// ============================================================================
// RWOFileHandler Implementation
// ============================================================================

std::vector<RWOObservation> RWOFileHandler::read(const std::string& filename) {
    std::vector<RWOObservation> observations;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open RWO file: " + filename);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.length() < 80) continue;
        
        try {
            RWOObservation rwo = parseLine(line);
            observations.push_back(rwo);
        } catch (...) {
            // Skip malformed lines
            continue;
        }
    }
    
    return observations;
}

void RWOFileHandler::write(const std::string& filename, 
                           const std::vector<RWOObservation>& observations) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot write RWO file: " + filename);
    }
    
    for (const auto& obs : observations) {
        file << formatLine(obs) << "\n";
    }
}

RWOObservation RWOFileHandler::parseLine(const std::string& line) {
    RWOObservation rwo;
    
    // First 80 columns: standard MPC format
    // Use MPCReader to parse the base observation
    std::vector<OpticalObservation> obs_vec;
    try {
        // Create temporary file with single line
        std::string temp = line.substr(0, 80);
        // Parse using MPC format (simplified - would need full MPCReader integration)
        
        // For now, basic parsing:
        rwo.observation.object_designation = line.substr(0, 12);
        
    } catch (...) {
        throw std::runtime_error("Failed to parse RWO line");
    }
    
    // Extended columns (if present): weights, residuals
    if (line.length() > 80) {
        std::istringstream iss(line.substr(80));
        iss >> rwo.weight_ra >> rwo.weight_dec 
            >> rwo.residual_ra >> rwo.residual_dec
            >> rwo.selection_flag;
        
        rwo.outlier = (rwo.weight_ra == 0.0 || rwo.weight_dec == 0.0);
    }
    
    return rwo;
}

std::string RWOFileHandler::formatLine(const RWOObservation& obs) {
    // Format standard MPC 80-column line + extensions
    std::ostringstream oss;
    
    // MPC format (simplified)
    oss << std::left << std::setw(12) << obs.observation.object_designation;
    // ... (full MPC formatting would go here)
    
    // Extended fields
    oss << std::fixed << std::setprecision(2)
        << " " << std::setw(6) << obs.weight_ra
        << " " << std::setw(6) << obs.weight_dec
        << " " << std::setw(8) << obs.residual_ra
        << " " << std::setw(8) << obs.residual_dec
        << " " << obs.selection_flag;
    
    return oss.str();
}

// ============================================================================
// AstDynConfigManager Implementation
// ============================================================================

bool AstDynConfigManager::loadConfiguration(const std::string& base_path, 
                                            const std::string& object_name) {
    object_name_ = object_name;
    
    // Try to load .opt file
    std::string opt_file = base_path + "/" + object_name + ".opt";
    options_.load(opt_file);  // OK if fails
    
    // Load .oef file (required)
    std::string oef_file = base_path + "/" + object_name + ".oef";
    if (!loadOEF(oef_file)) {
        return false;
    }
    
    // Try .rwo first, then .obs
    std::string rwo_file = base_path + "/" + object_name + ".rwo";
    if (!loadRWO(rwo_file)) {
        std::string obs_file = base_path + "/" + object_name + ".obs";
        if (!loadOBS(obs_file)) {
            std::cerr << "Warning: No observation file found\n";
        }
    }
    
    return true;
}

KeplerianElements AstDynConfigManager::getOsculatingElements() const {
    if (oef_data_.element_type == OrbitalElementSubType::MEAN) {
        return OEFFileHandler::meanToOsculating(oef_data_.keplerian);
    }
    return oef_data_.keplerian;
}

KeplerianElements AstDynConfigManager::getOriginalElements() const {
    return oef_data_.keplerian;
}

std::vector<OpticalObservation> AstDynConfigManager::getValidObservations() const {
    std::vector<OpticalObservation> valid_obs;
    for (const auto& rwo : observations_) {
        if (rwo.weight_ra > 0.0 && rwo.weight_dec > 0.0) {
            valid_obs.push_back(rwo.observation);
        }
    }
    return valid_obs;
}

bool AstDynConfigManager::saveConfiguration(const std::string& base_path, 
                                            const std::string& object_name) const {
    try {
        // Save .opt
        options_.save(base_path + "/" + object_name + ".opt");
        
        // Save .oef
        OEFFileHandler::write(base_path + "/" + object_name + ".oef", oef_data_);
        
        // Save .rwo
        RWOFileHandler::write(base_path + "/" + object_name + ".rwo", observations_);
        
        return true;
    } catch (...) {
        return false;
    }
}

bool AstDynConfigManager::exportForLegacyFortran(const std::string& output_dir, 
                                                 const std::string& object_name) const {
    return saveConfiguration(output_dir, object_name);
}

bool AstDynConfigManager::loadOEF(const std::string& filename) {
    try {
        oef_data_ = OEFFileHandler::read(filename);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Failed to load OEF: " << e.what() << "\n";
        return false;
    }
}

bool AstDynConfigManager::loadRWO(const std::string& filename) {
    try {
        observations_ = RWOFileHandler::read(filename);
        return true;
    } catch (...) {
        return false;
    }
}

bool AstDynConfigManager::loadOBS(const std::string& filename) {
    try {
        // Load standard MPC observations
        auto mpc_obs = MPCReader::readFile(filename);
        
        // Convert to RWO format
        observations_.clear();
        for (const auto& obs : mpc_obs) {
            RWOObservation rwo;
            rwo.observation = obs;
            rwo.weight_ra = 1.0;
            rwo.weight_dec = 1.0;
            observations_.push_back(rwo);
        }
        return true;
    } catch (...) {
        return false;
    }
}

} // namespace astdyn::config
