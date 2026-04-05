/**
 * @file gaia_mag18_catalog.h
 * @brief Fast access to compressed Gaia catalog subset (mag ≤ 18)
 * 
 * Provides efficient access to a pre-filtered subset of Gaia EDR3 containing
 * only stars with G magnitude ≤ 18. The catalog is stored in a compressed
 * binary format optimized for quick queries.
 * 
 * Features:
 * - Compressed storage (~50-70% size reduction)
 * - Binary search by source_id (O(log N))
 * - Cone search with spatial indexing
 * - All essential Gaia parameters
 * - Memory-efficient streaming
 */

#ifndef IOC_GAIALIB_GAIA_MAG18_CATALOG_H
#define IOC_GAIALIB_GAIA_MAG18_CATALOG_H

#include "types.h"
#include <string>
#include <vector>
#include <optional>
#include <cstdint>
#include <zlib.h>

namespace ioc_gaialib {

/**
 * @brief Header for compressed catalog file
 */
#pragma pack(push, 1)
struct Mag18CatalogHeader {
    char magic[8];           // "GAIA18\0\0" (8 bytes)
    uint32_t version;        // Format version (1) (4 bytes)
    uint64_t total_stars;    // Total number of stars (8 bytes)
    double mag_limit;        // Magnitude limit (18.0) (8 bytes)
    uint64_t data_offset;    // Offset to compressed data (8 bytes)
    uint64_t data_size;      // Size of uncompressed data (8 bytes)
    uint8_t reserved[20];    // Reserved for future use (20 bytes) = 64 total
};
#pragma pack(pop)

static_assert(sizeof(Mag18CatalogHeader) == 64, "Header must be 64 bytes");

/**
 * @brief Binary record format (52 bytes per star)
 * 
 * Optimized format preserving essential Gaia parameters:
 * - source_id: 8 bytes (uint64)
 * - RA: 8 bytes (double, degrees)
 * - Dec: 8 bytes (double, degrees)
 * - G mag: 8 bytes (double)
 * - BP mag: 8 bytes (double)
 * - RP mag: 8 bytes (double)
 * - Parallax: 4 bytes (float, mas)
 */
#pragma pack(push, 1)
struct Mag18Record {
    uint64_t source_id;
    double ra;
    double dec;
    double g_mag;
    double bp_mag;
    double rp_mag;
    float parallax;
};
#pragma pack(pop)

static_assert(sizeof(Mag18Record) == 52, "Record must be 52 bytes");

/**
 * @brief Reader for compressed Gaia mag≤18 catalog
 * 
 * Efficient access to pre-filtered Gaia subset with stars up to magnitude 18.
 * The catalog is sorted by source_id for fast binary search lookups.
 * 
 * Usage:
 * @code
 * GaiaMag18Catalog catalog("gaia_mag18.cat.gz");
 * 
 * // Query by source_id (fast binary search)
 * auto star = catalog.queryBySourceId(12345678901234);
 * 
 * // Cone search
 * auto stars = catalog.queryCone(180.0, 0.0, 1.0, 100);
 * 
 * // Get statistics
 * auto stats = catalog.getStatistics();
 * std::cout << "Total stars: " << stats.total_stars << std::endl;
 * @endcode
 */
class GaiaMag18Catalog {
public:
    /**
     * @brief Catalog statistics
     */
    struct Statistics {
        uint64_t total_stars;
        double mag_limit;
        size_t file_size;
        size_t uncompressed_size;
        double compression_ratio;
    };
    
    /**
     * @brief Constructor
     * @param catalog_file Path to compressed catalog file (.gz)
     */
    explicit GaiaMag18Catalog(const std::string& catalog_file);
    
    /**
     * @brief Destructor
     */
    ~GaiaMag18Catalog();
    
    /**
     * @brief Check if catalog is loaded
     * @return true if catalog is ready to use
     */
    bool isLoaded() const { return loaded_; }
    
    /**
     * @brief Get catalog statistics
     * @return Statistics structure with catalog info
     */
    Statistics getStatistics() const;
    
    /**
     * @brief Query star by Gaia source_id
     * @param source_id Gaia DR3 source identifier
     * @return GaiaStar if found, std::nullopt otherwise
     * 
     * Uses binary search on sorted catalog - O(log N) complexity.
     * Very fast even for large catalogs.
     */
    std::optional<ioc::gaia::GaiaStar> queryBySourceId(uint64_t source_id) const;
    
    /**
     * @brief Cone search around position
     * @param ra Right Ascension (degrees, J2016.0)
     * @param dec Declination (degrees, J2016.0)
     * @param radius Search radius (degrees)
     * @param max_results Maximum number of results (0 = unlimited)
     * @return Vector of stars within cone
     * 
     * Linear scan of all stars - suitable for radius < 5°.
     * For larger regions, consider using spatial indexing.
     */
    std::vector<ioc::gaia::GaiaStar> queryCone(double ra, double dec, double radius,
                                                size_t max_results = 0) const;
    
    /**
     * @brief Cone search with magnitude filter
     * @param ra Right Ascension (degrees)
     * @param dec Declination (degrees)
     * @param radius Search radius (degrees)
     * @param mag_min Minimum G magnitude
     * @param mag_max Maximum G magnitude
     * @param max_results Maximum number of results (0 = unlimited)
     * @return Vector of stars matching criteria
     */
    std::vector<ioc::gaia::GaiaStar> queryConeWithMagnitude(
        double ra, double dec, double radius,
        double mag_min, double mag_max,
        size_t max_results = 0) const;
    
    /**
     * @brief Find N brightest stars in cone
     * @param ra Right Ascension (degrees)
     * @param dec Declination (degrees)
     * @param radius Search radius (degrees)
     * @param n_brightest Number of brightest stars to return
     * @return Vector of N brightest stars, sorted by magnitude
     */
    std::vector<ioc::gaia::GaiaStar> queryBrightest(double ra, double dec, double radius,
                                                     size_t n_brightest) const;
    
    /**
     * @brief Count stars in cone
     * @param ra Right Ascension (degrees)
     * @param dec Declination (degrees)
     * @param radius Search radius (degrees)
     * @return Number of stars in cone
     */
    size_t countInCone(double ra, double dec, double radius) const;
    
    /**
     * @brief Get catalog file path
     * @return Path to catalog file
     */
    std::string getCatalogPath() const { return catalog_file_; }

private:
    /**
     * @brief Load catalog header and prepare for queries
     * @return true if successful
     */
    bool loadHeader();
    
    /**
     * @brief Read record at specific index
     * @param index Record index (0-based)
     * @return Record data if successful
     */
    std::optional<Mag18Record> readRecord(uint64_t index) const;
    
    /**
     * @brief Convert Mag18Record to GaiaStar
     * @param record Binary record
     * @return GaiaStar structure
     */
    ioc::gaia::GaiaStar recordToStar(const Mag18Record& record) const;
    
    /**
     * @brief Calculate angular distance between two positions
     * @param ra1 Right Ascension 1 (degrees)
     * @param dec1 Declination 1 (degrees)
     * @param ra2 Right Ascension 2 (degrees)
     * @param dec2 Declination 2 (degrees)
     * @return Angular distance (degrees)
     */
    double angularDistance(double ra1, double dec1, 
                          double ra2, double dec2) const;
    
    std::string catalog_file_;
    bool loaded_;
    Mag18CatalogHeader header_;
    mutable gzFile gz_file_;
};

} // namespace ioc_gaialib

#endif // IOC_GAIALIB_GAIA_MAG18_CATALOG_H
