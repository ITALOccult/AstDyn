#pragma once

#ifndef IOC_GAIALIB_UNIFIED_CATALOG_H
#define IOC_GAIALIB_UNIFIED_CATALOG_H

#include <string>
#include <vector>
#include <memory>
#include <optional>
#include <functional>
#include <future>
#include <filesystem>
#include <map>
#include "types.h"

namespace ioc::gaia {

// Forward declarations
struct GaiaStar;
// Use existing QueryParams from types.h

/**
 * @brief Configuration for Gaia catalog access
 */
struct GaiaCatalogConfig {
    enum class CatalogType {
        MULTIFILE_V2,      // Local multi-file V2 catalog
        COMPRESSED_V2,     // Local compressed V2 catalog
        SQLITE_DR3,        // SQLite with R*Tree and cross-references (recommended)
        ONLINE_ESA,        // ESA Gaia Archive online
        ONLINE_VIZIER      // VizieR online service
    };
    
    enum class CacheStrategy {
        DISABLED,          // No caching
        MEMORY_ONLY,       // In-memory cache only
        DISK_CACHE,        // Persistent disk cache
        HYBRID             // Memory + disk cache (default)
    };

    CatalogType catalog_type = CatalogType::MULTIFILE_V2;
    CacheStrategy cache_strategy = CacheStrategy::HYBRID;
    
    // Local catalog paths
        std::string multifile_directory;
        // Path to the SQLite database file
        std::string sqlite_file_path;
      // Required for MULTIFILE_V2
    std::string compressed_file_path;     // Required for COMPRESSED_V2
    
    // Online configuration
    std::string server_url;               // Custom server URL (optional)
    int timeout_seconds = 30;             // Network timeout
    std::string api_key;                  // API key if required
    
    // Performance settings
    size_t max_cached_chunks = 50;        // Memory cache size (~4GB)
    size_t max_concurrent_requests = 8;   // Max parallel requests
    bool enable_compression = true;       // Enable response compression
    
    // Cache directories
    std::string cache_directory = "~/.cache/gaia_catalog";
    size_t max_cache_size_gb = 10;        // Max disk cache size
    int cache_expiry_days = 30;           // Cache expiration
    
    // Logging
    enum class LogLevel { SILENT, ERROR, WARNING, INFO, DEBUG };
    LogLevel log_level = LogLevel::WARNING;
    std::string log_file;                 // Optional log file path
};

/**
 * @brief Statistics about catalog operations
 */
struct CatalogStats {
    size_t total_queries = 0;
    double average_query_time_ms = 0.0;
    size_t total_stars_returned = 0;
    double cache_hit_rate = 0.0;
    size_t memory_used_mb = 0;
    size_t disk_cache_used_mb = 0;
};

/**
 * @brief Information about loaded catalog
 */
struct CatalogInfo {
    std::string catalog_name;
    std::string version;
    size_t total_stars = 0;
    double magnitude_limit = 0.0;
    bool is_online = false;
    std::vector<std::string> available_fields;
};

/**
 * @brief Progress callback for async operations
 */
using ProgressCallback = std::function<void(int percent_complete, const std::string& status)>;

/**
 * @brief Unified Gaia catalog interface
 * 
 * This is the ONLY API that should be used to access Gaia catalogs.
 * It replaces all previous catalog access methods and provides:
 * - JSON-based configuration
 * - Automatic optimization for concurrent access
 * - Unified interface for all catalog types
 * - Built-in caching and performance monitoring
 */
class UnifiedGaiaCatalog {
public:
    /**
     * @brief Initialize the catalog with JSON configuration
     * @param json_config JSON string with catalog configuration
     * @return true if initialization successful
     * 
     * Example JSON configurations:
     * 
     * Multi-file V2 (recommended for performance):
     * {
     *   "catalog_type": "multifile_v2",
     *   "multifile_directory": "/path/to/multifile/catalog",
     *   "max_cached_chunks": 100,
     *   "log_level": "info"
     * }
     * 
     * Compressed V2:
     * {
     *   "catalog_type": "compressed_v2",
     *   "compressed_file_path": "/path/to/catalog.mag18v2",
     *   "log_level": "info"  
     * }
     * 
     * Online ESA:
     * {
     *   "catalog_type": "online_esa",
     *   "timeout_seconds": 30,
     *   "log_level": "info"
     * }
     */
    static bool initialize(const std::string& json_config);
    
    /**
     * @brief Get singleton instance (must call initialize first)
     */
    static UnifiedGaiaCatalog& getInstance();
    
    /**
     * @brief Shutdown and cleanup catalog
     */
    static void shutdown();
    
    /**
     * @brief Check if catalog is properly initialized
     */
    bool isInitialized() const;
    
    /**
     * @brief Get catalog information
     */
    CatalogInfo getCatalogInfo() const;
    
    /**
     * @brief Perform cone search query
     * @param params Query parameters
     * @return Vector of found stars
     */
    std::vector<GaiaStar> queryCone(const QueryParams& params) const;
    
    /**
     * @brief Asynchronous cone search
     * @param params Query parameters
     * @param progress_callback Optional progress callback
     * @return Future containing results
     */
    std::future<std::vector<GaiaStar>> queryAsync(
        const QueryParams& params,
        ProgressCallback progress_callback = nullptr
    ) const;
    
    /**
     * @brief Batch query multiple regions
     * @param param_list List of query parameters
     * @return Vector of result vectors (one per query)
     */
    std::vector<std::vector<GaiaStar>> batchQuery(
        const std::vector<QueryParams>& param_list
    ) const;
    
    /**
     * @brief Query star by source ID
     * @param source_id Gaia source identifier
     * @return Star data if found
     */
    std::optional<GaiaStar> queryBySourceId(uint64_t source_id) const;
    
    /**
     * @brief Query star by SAO designation
     * @param sao_number SAO catalog number (e.g., "123456")
     * @return Star data if found
     */
    std::optional<GaiaStar> queryBySAO(const std::string& sao_number) const;
    
    /**
     * @brief Query star by Tycho-2 designation
     * @param tycho2_designation Tycho-2 designation (e.g., "1234-567-1")
     * @return Star data if found
     */
    std::optional<GaiaStar> queryByTycho2(const std::string& tycho2_designation) const;
    
    /**
     * @brief Query star by HD designation
     * @param hd_number Henry Draper catalog number (e.g., "12345")
     * @return Star data if found
     */
    std::optional<GaiaStar> queryByHD(const std::string& hd_number) const;
    
    /**
     * @brief Query star by Hipparcos designation
     * @param hip_number Hipparcos catalog number (e.g., "12345")
     * @return Star data if found
     */
    std::optional<GaiaStar> queryByHipparcos(const std::string& hip_number) const;
    
    /**
     * @brief Query star by common name
     * @param common_name Star common name (e.g., "Sirius", "Polaris")
     * @return Star data if found
     */
    std::optional<GaiaStar> queryByName(const std::string& common_name) const;
    
    /**
     * @brief Query stars along a corridor/path
     * 
     * Searches for stars within a specified width along a path defined by
     * multiple waypoints. The corridor is formed by connecting the waypoints
     * with great circle arcs and including all stars within 'width' degrees
     * of the path.
     * 
     * @param params Corridor query parameters including path points and width
     * @return Vector of found stars
     * 
     * Example JSON input:
     * {
     *   "path": [
     *     {"ra": 0.0, "dec": 0.0},
     *     {"ra": 10.0, "dec": 5.0},
     *     {"ra": 20.0, "dec": 10.0}
     *   ],
     *   "width": 0.5,
     *   "max_magnitude": 16.0
     * }
     */
    std::vector<GaiaStar> queryCorridor(const CorridorQueryParams& params) const;
    
    /**
     * @brief Query stars along a corridor using JSON configuration
     * @param json_config JSON string with corridor parameters
     * @return Vector of found stars
     */
    std::vector<GaiaStar> queryCorridorJSON(const std::string& json_config) const;
    
    /**
     * @brief Asynchronous corridor query
     * @param params Corridor query parameters
     * @param progress_callback Optional progress callback
     * @return Future containing results
     */
    std::future<std::vector<GaiaStar>> queryCorridorAsync(
        const CorridorQueryParams& params,
        ProgressCallback progress_callback = nullptr
    ) const;
    
    /**
     * @brief Query stars along an orbit interpolated by Chebyshev polynomials
     * 
     * This method discretizes the orbit defined by Chebyshev polynomials over the
     * specified time interval and searches for stars within 'width' degrees of
     * the path.
     * 
     * @param params Orbit query parameters including polynomials, time range, and width
     * @return Vector of found stars
     */
    std::vector<GaiaStar> queryOrbit(const OrbitQueryParams& params) const;
    
    /**
     * @brief Get performance statistics
     */
    CatalogStats getStatistics() const;
    
    /**
     * @brief Clear all caches
     */
    void clearCache();
    
    /**
     * @brief Reconfigure catalog (requires reinitialization)
     * @param json_config New configuration
     * @return true if reconfiguration successful
     */
    bool reconfigure(const std::string& json_config);

private:
    // Singleton pattern
    UnifiedGaiaCatalog() = default;
    UnifiedGaiaCatalog(const UnifiedGaiaCatalog&) = delete;
    UnifiedGaiaCatalog& operator=(const UnifiedGaiaCatalog&) = delete;
    
    // Private implementation pattern
    class Impl;
    std::unique_ptr<Impl> pimpl_;
    
    // Static members for singleton
    static std::unique_ptr<UnifiedGaiaCatalog> instance_;
    static bool initialized_;
    
    friend class std::default_delete<UnifiedGaiaCatalog>; // Allow std::unique_ptr to delete
};

} // namespace ioc::gaia

#endif // IOC_GAIALIB_UNIFIED_CATALOG_H