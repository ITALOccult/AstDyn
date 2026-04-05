#pragma once

#include <string>
#include <unordered_map>
#include <optional>
#include <vector>
#include "types.h"

namespace ioc::gaia {

/**
 * @brief Cross-match information for a star
 */
struct CrossMatchInfo {
    std::string common_name;
    std::string sao_designation;
    std::string hd_designation;
    std::string hip_designation;
    std::string tycho2_designation;
    std::string constellation;
    std::string notes;
    
    // IAU designation fields
    std::string flamsteed_designation;  ///< Flamsteed designation (e.g., "61 Cyg")
    std::string bayer_designation;      ///< Bayer designation (e.g., "α CMa")
    bool has_flamsteed = false;         ///< True if star has Flamsteed designation
    bool has_bayer = false;             ///< True if star has Bayer designation
};

/**
 * @brief Common star names and cross-match database
 * 
 * Provides lookup functionality for:
 * - Common star names (Sirius, Vega, etc.)
 * - Cross-catalog identifiers (SAO, HD, HIP, TYC)
 * - Constellation information
 * - Historical notes
 */
class CommonStarNames {
public:
    /**
     * @brief Load star names database from file
     * @param database_path Path to CSV file with star names
     * @return true if loaded successfully
     */
    bool loadDatabase(const std::string& database_path);
    
    /**
     * @brief Load default embedded database
     * @return true if loaded successfully
     */
    bool loadDefaultDatabase();
    
    /**
     * @brief Get cross-match info by Gaia source ID
     * @param source_id Gaia DR3 source identifier
     * @return Cross-match information if available
     */
    std::optional<CrossMatchInfo> getCrossMatch(uint64_t source_id) const;
    
    /**
     * @brief Find star by common name
     * @param name Common star name (case-insensitive)
     * @return Gaia source ID if found
     */
    std::optional<uint64_t> findByName(const std::string& name) const;
    
    /**
     * @brief Find star by SAO designation
     * @param sao_number SAO catalog number
     * @return Gaia source ID if found
     */
    std::optional<uint64_t> findBySAO(const std::string& sao_number) const;
    
    /**
     * @brief Find star by HD designation
     * @param hd_number HD catalog number
     * @return Gaia source ID if found
     */
    std::optional<uint64_t> findByHD(const std::string& hd_number) const;
    
    /**
     * @brief Find star by Hipparcos designation
     * @param hip_number Hipparcos catalog number
     * @return Gaia source ID if found
     */
    std::optional<uint64_t> findByHipparcos(const std::string& hip_number) const;
    
    /**
     * @brief Find star by Tycho-2 designation
     * @param tycho2_id Tycho-2 designation
     * @return Gaia source ID if found
     */
    std::optional<uint64_t> findByTycho2(const std::string& tycho2_id) const;
    
    /**
     * @brief Get all stars in a constellation
     * @param constellation Constellation name
     * @return Vector of Gaia source IDs
     */
    std::vector<uint64_t> getConstellationStars(const std::string& constellation) const;
    
    /**
     * @brief Get number of stars in database
     */
    size_t size() const { return source_id_to_info_.size(); }
    
    /**
     * @brief Check if database is loaded
     */
    bool isLoaded() const { return !source_id_to_info_.empty(); }
    
    /**
     * @brief Get all available common names
     * @return Vector of all star names in database
     */
    std::vector<std::string> getAllNames() const;

    /**
     * @brief Get access to internal source_id map (for IAU parser)
     * @return Const reference to internal map
     */
    const std::unordered_map<uint64_t, CrossMatchInfo>& getSourceIdMap() const {
        return source_id_to_info_;
    }
    
    /**
     * @brief Add star to internal indices (for IAU parser)
     * @param source_id Gaia source ID
     * @param info Cross-match information
     */
    void addToIndices(uint64_t source_id, const CrossMatchInfo& info);

private:
    // Primary storage: source_id -> cross-match info
    std::unordered_map<uint64_t, CrossMatchInfo> source_id_to_info_;
    
    // Reverse lookup indices
    std::unordered_map<std::string, uint64_t> name_to_source_id_;
    std::unordered_map<std::string, uint64_t> sao_to_source_id_;
    std::unordered_map<std::string, uint64_t> hd_to_source_id_;
    std::unordered_map<std::string, uint64_t> hip_to_source_id_;
    std::unordered_map<std::string, uint64_t> tycho2_to_source_id_;
    std::unordered_map<std::string, std::vector<uint64_t>> constellation_to_source_ids_;
    
    /**
     * @brief Parse CSV line into cross-match info
     * @param line CSV line
     * @param info Output cross-match info
     * @param source_id Output source ID
     * @return true if parsed successfully
     */
    bool parseLine(const std::string& line, CrossMatchInfo& info, uint64_t& source_id);
    
    /**
     * @brief Convert string to lowercase for case-insensitive lookup
     */
    std::string toLower(const std::string& str) const;
    
    /**
     * @brief Load embedded star database
     */
    void loadEmbeddedDatabase();
};

} // namespace ioc::gaia