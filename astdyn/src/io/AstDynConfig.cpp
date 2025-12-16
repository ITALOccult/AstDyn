/**
 * @file AstDynConfig.cpp
 * @brief Implementation of AstDyn configuration parsers
 */

#include "astdyn/io/AstDynConfig.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/observations/MPCReader.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/time/TimeScale.hpp"
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
        // Skip short lines
        if (line.length() < 80) continue;
        
        // Skip comment lines (starting with !)
        if (!line.empty() && line[0] == '!') continue;
        
        // Skip header lines (starting with % or #)
        if (!line.empty() && (line[0] == '%' || line[0] == '#')) continue;
        
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
    
    // RWO format (AstDyS format - similar to MPC but with different column positions)
    // Columns 1-12: Object designation
    // Column 14: Observation type (O=optical, S=satellite, etc.)
    // Column 16: Note/discovery asterisk
    // Columns 18-37: Date (YYYY MM DD.dddddddddd)
    // Columns 55-66: RA (HH MM SS.sss)
    // Columns 68-79: Dec (sDD MM SS.ss)
    
    try {
        if (line.length() < 114) {
            throw std::runtime_error("RWO line too short (need at least 114 chars for Dec)");
        }
        
        // Parse object designation (columns 1-12)
        rwo.observation.object_designation = line.substr(0, 12);
        // Trim trailing spaces
        size_t end = rwo.observation.object_designation.find_last_not_of(" ");
        if (end != std::string::npos) {
            rwo.observation.object_designation = rwo.observation.object_designation.substr(0, end + 1);
        }
        
        // Parse discovery asterisk (column 16)
        if (line.length() > 15 && line[15] != ' ') {
            rwo.observation.is_discovery = (line[15] == '*' || line[15] == 'A');
        }
        
        // Parse date (columns 18-34, starting at index 17): "YYYY MM DD.dddddddddd"  
        std::string date_str = line.substr(17, 17);
        
        std::istringstream date_iss(date_str);
        int year, month;
        double day;
        date_iss >> year >> month >> day;
        
        // Convert to MJD
        int day_int = static_cast<int>(day);
        double fraction = day - day_int;
        rwo.observation.mjd_utc = astdyn::time::calendar_to_mjd(year, month, day_int, fraction);
        
        // Parse RA (columns 51-62, starting at index 50): "HH MM SS.sss"
        std::string ra_str = line.substr(50, 12);
        std::istringstream ra_iss(ra_str);
        int ra_h, ra_m;
        double ra_s;
        ra_iss >> ra_h >> ra_m >> ra_s;
        double ra_deg = ra_h * 15.0 + ra_m * 0.25 + ra_s * (15.0 / 3600.0);
        rwo.observation.ra = ra_deg * constants::DEG_TO_RAD;
        
        // Parse Dec (columns 104-115, starting at index 103): "sDD MM SS.ss"
        std::string dec_str = line.substr(103, 12);
        char sign = dec_str[0];
        std::istringstream dec_iss(dec_str.substr(1));
        int dec_d, dec_m;
        double dec_s;
        dec_iss >> dec_d >> dec_m >> dec_s;
        double dec_deg = dec_d + dec_m / 60.0 + dec_s / 3600.0;
        if (sign == '-') dec_deg = -dec_deg;
        rwo.observation.dec = dec_deg * constants::DEG_TO_RAD;
        
        // Set default uncertainties
        rwo.observation.sigma_ra = 0.5 * constants::ARCSEC_TO_RAD;
        rwo.observation.sigma_dec = 0.5 * constants::ARCSEC_TO_RAD;
        
    } catch (...) {
        throw std::runtime_error("Failed to parse RWO line");
    }
    
    // Extended columns (if present): weights, residuals
    if (line.length() > 80) {
        try {
            std::istringstream iss(line.substr(80));
            iss >> rwo.weight_ra >> rwo.weight_dec 
                >> rwo.residual_ra >> rwo.residual_dec
                >> rwo.selection_flag;
            
            rwo.outlier = (rwo.weight_ra == 0.0 || rwo.weight_dec == 0.0);
        } catch (...) {
            // If extended columns fail to parse, use defaults
            rwo.weight_ra = 1.0;
            rwo.weight_dec = 1.0;
            rwo.residual_ra = 0.0;
            rwo.residual_dec = 0.0;
            rwo.selection_flag = 1;
            rwo.outlier = false;
        }
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
    
    // Option parsing (.opt) removed in favor of JSON
    // Only loading elements and observations here
    
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
        // Option saving removed
        
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
