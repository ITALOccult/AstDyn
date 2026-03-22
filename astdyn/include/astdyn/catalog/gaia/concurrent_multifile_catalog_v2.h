#pragma once

#ifndef IOC_GAIALIB_CONCURRENT_MULTIFILE_CATALOG_V2_H
#define IOC_GAIALIB_CONCURRENT_MULTIFILE_CATALOG_V2_H

#include <string>
#include <vector>
#include <memory>
#include <shared_mutex>
#include <unordered_map>
#include <chrono>
#include <optional>
#include <set>
#include "gaia_mag18_catalog_v2.h"

namespace ioc {
namespace gaia {

/**
 * @brief Thread-safe multi-file catalog optimized for concurrent access
 */
class ConcurrentMultiFileCatalogV2 {
public:
    /**
     * @brief Load catalog from multi-file directory
     * @param catalog_dir Directory containing metadata.dat and chunks/
     * @param max_cached_chunks Maximum chunks to keep in memory (default: 50)
     */
    explicit ConcurrentMultiFileCatalogV2(const std::string& catalog_dir, 
                                          size_t max_cached_chunks = 50);
    
    /**
     * @brief Destructor
     */
    ~ConcurrentMultiFileCatalogV2();
    
    // Non-copyable, non-movable 
    ConcurrentMultiFileCatalogV2(const ConcurrentMultiFileCatalogV2&) = delete;
    ConcurrentMultiFileCatalogV2& operator=(const ConcurrentMultiFileCatalogV2&) = delete;
    ConcurrentMultiFileCatalogV2(ConcurrentMultiFileCatalogV2&&) = delete;
    ConcurrentMultiFileCatalogV2& operator=(ConcurrentMultiFileCatalogV2&&) = delete;
    
    /**
     * @brief Thread-safe cone search
     */
    std::vector<GaiaStar> queryCone(double ra, double dec, double radius, 
                                    size_t max_results = 0);
    
    /**
     * @brief Thread-safe search by source_id
     */
    std::optional<GaiaStar> queryBySourceId(uint64_t source_id);
    
    /**
     * @brief Get basic catalog info (thread-safe, no locking needed)
     */
    uint64_t getTotalStars() const { return header_.total_stars; }
    uint64_t getNumChunks() const { return header_.total_chunks; }
    uint32_t getNumPixels() const { return header_.num_healpix_pixels; }
    double getMagLimit() const { return header_.mag_limit; }
    
    /**
     * @brief Get performance statistics
     */
    struct ConcurrencyStats {
        size_t cache_hits;
        size_t cache_misses;
        size_t chunks_loaded;
        size_t memory_used_mb;
        size_t active_readers;
        double hit_rate;
    };
    ConcurrencyStats getStats() const;
    
    /**
     * @brief Preload chunks for better concurrent performance
     */
    void preloadChunks(const std::vector<uint64_t>& chunk_ids);
    
    /**
     * @brief Clear cache and free memory
     */
    void clearCache();

private:
    struct ChunkData {
        uint64_t chunk_id;
        std::vector<Mag18RecordV2> records;
        std::chrono::steady_clock::time_point last_access;
        mutable std::shared_mutex access_mutex;  // Allow multiple readers
        
        ChunkData(uint64_t id, std::vector<Mag18RecordV2> data) 
            : chunk_id(id), records(std::move(data)), 
              last_access(std::chrono::steady_clock::now()) {}
    };
    
    // Core data
    std::string catalog_dir_;
    Mag18CatalogHeaderV2 header_;
    
    // NEW HEALPix index format: maps pixel -> list of chunk IDs
    std::vector<PixelChunkEntry> pixel_index_;
    std::vector<uint32_t> chunk_lists_;  // Flattened array of chunk IDs
    
    std::vector<ChunkInfo> chunk_index_;
    
    // Thread-safe cache with read-write locks
    mutable std::shared_mutex cache_mutex_;  // Protects cache structure
    mutable std::unordered_map<uint64_t, std::shared_ptr<ChunkData>> chunk_cache_;
    size_t max_cached_chunks_;
    
    // Statistics
    mutable std::atomic<size_t> cache_hits_{0};
    mutable std::atomic<size_t> cache_misses_{0};
    mutable std::atomic<size_t> active_readers_{0};
    
    // Internal methods
    bool loadMetadata();
    std::shared_ptr<ChunkData> loadChunk(uint64_t chunk_id);
    std::shared_ptr<ChunkData> getOrLoadChunk(uint64_t chunk_id);
    std::set<uint32_t> getChunksForCone(double ra, double dec, double radius) const;
    std::vector<uint32_t> getPixelsInCone(double ra, double dec, double radius) const;
    uint32_t getHEALPixPixel(double ra, double dec) const;
    uint32_t ang2pix_nest(double theta, double phi) const;
    GaiaStar recordToStar(const Mag18RecordV2& record) const;
    double angularDistance(double ra1, double dec1, double ra2, double dec2) const;
    void evictLRUChunks();
    std::string getChunkPath(uint64_t chunk_id) const;
};

} // namespace gaia
} // namespace ioc

#endif // IOC_GAIALIB_CONCURRENT_MULTIFILE_CATALOG_V2_H