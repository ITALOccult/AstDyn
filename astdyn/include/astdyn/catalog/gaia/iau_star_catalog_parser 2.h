#pragma once

#include <astdyn/catalog/gaia/common_star_names.h>
#include <string>
#include <memory>

namespace ioc::gaia {

/**
 * @brief Parser for official IAU Catalog of Star Names (IAU-CSN)
 * 
 * This class loads and parses the official IAU star names catalog
 * to provide comprehensive cross-match information including
 * Flamsteed numbers, Bayer designations, and proper names.
 */
class IAUStarCatalogParser {
public:
    /**
     * @brief Load IAU catalog from JSON file
     * @param json_file_path Path to IAU-CSN.json file
     * @return True if loading was successful
     */
    static bool loadFromJSON(const std::string& json_file_path, CommonStarNames& star_names);
    
    /**
     * @brief Get number of stars loaded from IAU catalog
     */
    static size_t getLoadedStarsCount() { return loaded_stars_count_; }
    
    /**
     * @brief Get statistics about loaded catalog
     */
    struct CatalogStats {
        size_t total_stars = 0;
        size_t stars_with_hr = 0;
        size_t stars_with_hd = 0;
        size_t stars_with_hip = 0;
        size_t stars_with_flamsteed = 0;
        size_t stars_with_bayer = 0;
        size_t stars_with_proper_names = 0;
        size_t exoplanet_host_stars = 0;
    };
    
    static CatalogStats getStatistics() { return catalog_stats_; }

private:
    static size_t loaded_stars_count_;
    static CatalogStats catalog_stats_;
    
    /**
     * @brief Parse a single star entry from JSON
     * @param json_star JSON object representing one star
     * @param cross_match Output cross-match info
     * @param gaia_source_id Output Gaia source ID (if available)
     * @return True if parsing was successful
     */
    static bool parseStar(const std::string& json_star, CrossMatchInfo& cross_match, uint64_t& gaia_source_id);
    
    /**
     * @brief Try to resolve Gaia source ID from HR/HD/HIP numbers
     * This is a placeholder - in a real implementation you'd have
     * a cross-match table or use SIMBAD/CDS services
     */
    static uint64_t resolveGaiaSourceId(const std::string& hr, const std::string& hd, const std::string& hip);
    
    /**
     * @brief Parse Flamsteed designation from ID field
     * @param id_field The ID field from IAU catalog
     * @param constellation The constellation abbreviation
     * @return Flamsteed designation (e.g., "61 Cyg") or empty string
     */
    static std::string parseFlammsteed(const std::string& id_field, const std::string& constellation);
    
    /**
     * @brief Parse Bayer designation from ID field
     * @param id_field The ID field from IAU catalog
     * @param constellation The constellation abbreviation
     * @return Bayer designation (e.g., "α CMa") or empty string
     */
    static std::string parseBayer(const std::string& id_field, const std::string& constellation);
};

} // namespace ioc::gaia