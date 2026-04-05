#ifndef GAIA_MAG18_CATALOG_V2_H
#define GAIA_MAG18_CATALOG_V2_H

#include "types.h"
#include <string>
#include <vector>
#include <memory>
#include <zlib.h>
#include <cstdint>
#include <mutex>
#include <shared_mutex>
#include <atomic>

namespace ioc {
namespace gaia {

/**
 * @brief Gaia Mag18 Catalog V2 - Optimized version with HEALPix spatial index
 * @deprecated Use UnifiedGaiaCatalog with "compressed_v2" configuration instead.
 * 
 * Major improvements over V1:
 * - HEALPix NSIDE=64 spatial index for 100-300x faster cone searches
 * - Chunk-based compression (1M records/chunk) for 5x faster random access
 * - Extended 80-byte record format with proper motion and quality flags
 * 
 * Performance targets:
 * - Cone search 0.5°: 50 ms (vs 15 sec in V1)
 * - Cone search 5°: 500 ms (vs 48 sec in V1)
 * - Query source_id: <1 ms (unchanged)
 * 
 * File size: ~14 GB (vs 9 GB in V1, +55% for extended data)
 */

/**
 * @brief Extended record format with proper motion and quality parameters
 * Total: 84 bytes per star (with 4-byte alignment)
 */
#pragma pack(push, 4)
struct Mag18RecordV2 {
    uint64_t source_id;        // 8 bytes - Gaia DR3 source_id
    double ra;                 // 8 bytes - Right Ascension (degrees)
    double dec;                // 8 bytes - Declination (degrees)
    
    // Photometry (28 bytes)
    float g_mag;              // 4 bytes - G magnitude
    float bp_mag;             // 4 bytes - BP magnitude
    float rp_mag;             // 4 bytes - RP magnitude
    float g_mag_error;        // 4 bytes - G magnitude error
    float bp_mag_error;       // 4 bytes - BP magnitude error
    float rp_mag_error;       // 4 bytes - RP magnitude error
    float bp_rp;              // 4 bytes - BP-RP color
    
    // Astrometry (20 bytes)
    float parallax;           // 4 bytes - Parallax (mas)
    float parallax_error;     // 4 bytes - Parallax error (mas)
    float pmra;               // 4 bytes - Proper motion RA (mas/yr)
    float pmdec;              // 4 bytes - Proper motion Dec (mas/yr)
    float pmra_error;         // 4 bytes - PM RA error (mas/yr)
    
    // Quality (8 bytes)
    float ruwe;               // 4 bytes - Renormalized Unit Weight Error
    uint16_t phot_bp_n_obs;   // 2 bytes - Number of BP observations
    uint16_t phot_rp_n_obs;   // 2 bytes - Number of RP observations
    
    // HEALPix pixel (4 bytes)
    uint32_t healpix_pixel;   // NSIDE=64, NESTED scheme
};
#pragma pack(pop)

// Static assertion to verify correct structure size
static_assert(sizeof(Mag18RecordV2) == 84, "Mag18RecordV2 must be exactly 84 bytes");

/**
 * @brief HEALPix index entry - OLD format (deprecated, kept for compatibility)
 */
#pragma pack(push, 4)
struct HEALPixIndexEntry {
    uint32_t pixel_id;        // HEALPix pixel number (NSIDE=64, NESTED)
    uint64_t first_star_idx;  // Index of first star in this pixel
    uint32_t num_stars;       // Number of stars in this pixel
    uint32_t reserved;        // Reserved for future use
};
#pragma pack(pop)

// Static assertion for HEALPixIndexEntry size
static_assert(sizeof(HEALPixIndexEntry) == 20, "HEALPixIndexEntry must be exactly 20 bytes");

/**
 * @brief NEW HEALPix index entry - maps pixel to chunks (16 bytes)
 * This format is more efficient as stars are stored sorted by source_id, not by pixel
 */
#pragma pack(push, 4)
struct PixelChunkEntry {
    uint32_t pixel_id;           // HEALPix pixel number (NSIDE=64, NESTED)
    uint32_t num_chunks;         // Number of chunks containing stars in this pixel
    uint64_t chunk_list_offset;  // Offset into chunk list array (in elements, not bytes)
};
#pragma pack(pop)

static_assert(sizeof(PixelChunkEntry) == 16, "PixelChunkEntry must be exactly 16 bytes");

/**
 * @brief Chunk compression info - one per chunk (40 bytes)
 */
struct ChunkInfo {
    uint64_t chunk_id;              // Chunk number (0-based)
    uint64_t first_star_idx;        // Global index of first star in chunk
    uint32_t num_stars;             // Number of stars in this chunk
    uint32_t compressed_size;       // Size of compressed data (bytes)
    uint64_t file_offset;           // Offset in file where compressed chunk starts
    uint32_t uncompressed_size;     // Original size before compression
    uint32_t reserved;              // Reserved for future use
};

/**
 * @brief Catalog header with spatial index
 * Total: 268 bytes (with 4-byte alignment)
 */
#pragma pack(push, 4)
struct Mag18CatalogHeaderV2 {
    char magic[8];                  // "GAIA18V2"
    uint32_t version;               // Version number (2)
    uint32_t format_flags;          // Format options bitfield
    
    uint64_t total_stars;           // Total number of stars in catalog
    uint64_t total_chunks;          // Number of compressed chunks
    uint32_t stars_per_chunk;       // Stars per chunk (1,000,000)
    uint32_t healpix_nside;         // HEALPix NSIDE (64)
    
    double mag_limit;               // Maximum G magnitude (18.0)
    double ra_min, ra_max;          // RA bounds
    double dec_min, dec_max;        // Dec bounds
    
    uint64_t header_size;           // Size of this header (268)
    uint64_t healpix_index_offset;  // Offset to HEALPix index
    uint64_t healpix_index_size;    // Size of HEALPix index (bytes)
    uint32_t num_healpix_pixels;    // Number of pixels with data
    
    uint64_t chunk_index_offset;    // Offset to chunk index
    uint64_t chunk_index_size;      // Size of chunk index (bytes)
    
    uint64_t data_offset;           // Offset to compressed star data
    uint64_t data_size;             // Total size of compressed data
    
    char creation_date[32];         // ISO 8601 date
    char source_catalog[32];        // Source catalog name (GRAPPA3E)
    char reserved[64];              // Reserved for future use
};
#pragma pack(pop)

// Static assertion for header size
static_assert(sizeof(Mag18CatalogHeaderV2) == 268, "Mag18CatalogHeaderV2 must be exactly 268 bytes");

/**
 * @brief Gaia Mag18 Catalog V2 Reader with spatial indexing
 * @deprecated Use UnifiedGaiaCatalog with "compressed_v2" configuration instead.
 */
class [[deprecated("Use UnifiedGaiaCatalog instead")]] Mag18CatalogV2 {
public:
    Mag18CatalogV2();
    explicit Mag18CatalogV2(const std::string& catalog_path);
    ~Mag18CatalogV2();
    
    /**
     * @brief Open and load catalog with indices
     */
    bool open(const std::string& catalog_path);
    
    /**
     * @brief Close catalog and free resources
     */
    void close();
    
    /**
     * @brief Get catalog statistics
     */
    uint64_t getTotalStars() const { return header_.total_stars; }
    double getMagLimit() const { return header_.mag_limit; }
    uint32_t getNumPixels() const { return header_.num_healpix_pixels; }
    uint64_t getNumChunks() const { return header_.total_chunks; }
    
    /**
     * @brief Query by Gaia source_id (binary search, <1 ms)
     */
    std::optional<GaiaStar> queryBySourceId(uint64_t source_id);
    
    /**
     * @brief Cone search using HEALPix index (target: 50 ms for 0.5°)
     * @param ra Right Ascension (degrees)
     * @param dec Declination (degrees)
     * @param radius Search radius (degrees)
     * @param max_results Maximum number of results (0 = no limit)
     * @return Vector of stars within cone
     */
    std::vector<GaiaStar> queryCone(double ra, double dec, double radius, 
                                     size_t max_results = 0);
    
    /**
     * @brief Cone search with magnitude filter (optimized)
     */
    std::vector<GaiaStar> queryConeWithMagnitude(double ra, double dec, double radius,
                                                   double mag_min, double mag_max,
                                                   size_t max_results = 0);
    
    /**
     * @brief Get N brightest stars in cone (optimized)
     */
    std::vector<GaiaStar> queryBrightest(double ra, double dec, double radius,
                                          size_t num_stars);
    
    /**
     * @brief Count stars in cone (fast, no star loading)
     */
    size_t countInCone(double ra, double dec, double radius);
    
    /**
     * @brief Get HEALPix pixel for position
     */
    uint32_t getHEALPixPixel(double ra, double dec) const;
    
    /**
     * @brief Get all pixels that intersect with cone
     */
    std::vector<uint32_t> getPixelsInCone(double ra, double dec, double radius) const;
    
    /**
     * Get extended star information (V2 format)
     */
    std::optional<Mag18RecordV2> getExtendedRecord(uint64_t source_id);
    
    /**
     * @brief Enable/disable parallel processing
     * @param enable True to enable parallel queries
     * @param num_threads Number of threads (0 = auto-detect)
     */
    void setParallelProcessing(bool enable, size_t num_threads = 0);
    
    /**
     * @brief Check if parallel processing is enabled
     */
    bool isParallelEnabled() const { return enable_parallel_.load(); }
    
    /**
     * @brief Get number of threads used for parallel processing
     */
    size_t getNumThreads() const { return num_threads_; }
    
private:
    std::string catalog_path_;
    FILE* file_;
    mutable std::mutex file_mutex_;  // Protects file operations
    Mag18CatalogHeaderV2 header_;
    
    // HEALPix spatial index (loaded in memory for fast queries) - READ ONLY after load
    std::vector<HEALPixIndexEntry> healpix_index_;
    
    // Chunk index (loaded in memory) - READ ONLY after load
    std::vector<ChunkInfo> chunk_index_;
    
    // Parallelization settings
    std::atomic<bool> enable_parallel_;
    size_t num_threads_;
    
    // Chunk cache (LRU, keep last N decompressed chunks) - THREAD SAFE
    struct ChunkCache {
        uint64_t chunk_id;
        std::vector<Mag18RecordV2> records;
        std::atomic<size_t> access_count;
        
        ChunkCache() : chunk_id(0), access_count(0) {}
        ChunkCache(uint64_t id, std::vector<Mag18RecordV2> recs, size_t count)
            : chunk_id(id), records(std::move(recs)), access_count(count) {}
        
        // Move constructor (atomic requires explicit move)
        ChunkCache(ChunkCache&& other) noexcept
            : chunk_id(other.chunk_id),
              records(std::move(other.records)),
              access_count(other.access_count.load()) {}
        
        // Move assignment (atomic requires explicit move)
        ChunkCache& operator=(ChunkCache&& other) noexcept {
            if (this != &other) {
                chunk_id = other.chunk_id;
                records = std::move(other.records);
                access_count.store(other.access_count.load());
            }
            return *this;
        }
        
        // Delete copy operations (atomics can't be copied)
        ChunkCache(const ChunkCache&) = delete;
        ChunkCache& operator=(const ChunkCache&) = delete;
    };
    static constexpr size_t MAX_CACHED_CHUNKS = 10;
    mutable std::vector<ChunkCache> chunk_cache_;
    mutable std::shared_mutex cache_mutex_;  // Reader-writer lock for cache
    
    // Internal methods
    bool loadHeader();
    bool loadHEALPixIndex();
    bool loadChunkIndex();
    
    std::optional<Mag18RecordV2> readRecord(uint64_t index);
    std::vector<Mag18RecordV2> readChunk(uint64_t chunk_id);
    
    GaiaStar recordToStar(const Mag18RecordV2& record) const;
    
    double angularDistance(double ra1, double dec1, double ra2, double dec2) const;
    
    // HEALPix helpers
    uint32_t ang2pix_nest(double theta, double phi) const;
    void pix2ang_nest(uint32_t pixel, double& theta, double& phi) const;
    void query_disc(double theta, double phi, double radius,
                    std::vector<uint32_t>& pixels) const;
};

} // namespace gaia
} // namespace ioc

#endif // GAIA_MAG18_CATALOG_V2_H
