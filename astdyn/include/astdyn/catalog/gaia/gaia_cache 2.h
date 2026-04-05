#pragma once

#include "types.h"
#include <string>
#include <vector>
#include <memory>
#include <unordered_set>

namespace ioc {
namespace gaia {

// Forward declaration
class GaiaClient;

/**
 * GaiaCache - HEALPix-based local cache for Gaia catalog data
 * 
 * Implements a disk-based cache using HEALPix tessellation (NSIDE=32) to store
 * Gaia stars locally. Significantly speeds up repeated queries by avoiding
 * network requests.
 * 
 * Cache structure:
 *   cache_dir/
 *     ├── index.json           (tile coverage and metadata)
 *     ├── tiles/
 *     │   ├── tile_00000.json
 *     │   ├── tile_00001.json
 *     │   └── ...
 *     └── metadata.json        (cache statistics)
 * 
 * Example usage:
 * @code
 *   GaiaCache cache("~/.ioccultcalc/gaia_cache");
 *   
 *   // Download region once
 *   cache.downloadRegion(67.0, 15.0, 5.0, 12.0);
 *   
 *   // Fast cached queries
 *   auto stars = cache.queryRegion(67.0, 15.0, 0.5, 12.0);  // < 100ms
 * @endcode
 */
class GaiaCache {
public:
    /**
     * Construct a cache manager
     * 
     * @param cache_dir Directory to store cache files (created if doesn't exist)
     * @param nside HEALPix NSIDE parameter (default: 32, tiles ~13.4° square)
     * @param release Gaia data release (default: DR3)
     */
    explicit GaiaCache(
        const std::string& cache_dir,
        int nside = 32,
        GaiaRelease release = GaiaRelease::DR3
    );
    
    /**
     * Destructor - saves index if modified
     */
    ~GaiaCache();
    
    // No copy, allow move
    GaiaCache(const GaiaCache&) = delete;
    GaiaCache& operator=(const GaiaCache&) = delete;
    GaiaCache(GaiaCache&&) noexcept;
    GaiaCache& operator=(GaiaCache&&) noexcept;
    
    /**
     * Set the GaiaClient to use for downloading data
     * Must be called before downloadRegion() or downloadBox()
     * 
     * @param client Pointer to GaiaClient instance (must remain valid)
     */
    void setClient(GaiaClient* client);
    
    /**
     * Download and cache stars in a circular region
     * 
     * Downloads all stars brighter than max_magnitude within radius of the center
     * and stores them in HEALPix tiles. Skips tiles already cached.
     * 
     * @param ra_center Center RA [degrees]
     * @param dec_center Center Dec [degrees]
     * @param radius Search radius [degrees]
     * @param max_magnitude Maximum G-band magnitude
     * @param callback Optional progress callback
     * @return Number of stars downloaded
     * @throws GaiaException on network or I/O errors
     */
    size_t downloadRegion(
        double ra_center,
        double dec_center,
        double radius,
        double max_magnitude,
        ProgressCallback callback = nullptr
    );
    
    /**
     * Download and cache stars in a rectangular region
     * 
     * @param ra_min Minimum RA [degrees]
     * @param ra_max Maximum RA [degrees]
     * @param dec_min Minimum Dec [degrees]
     * @param dec_max Maximum Dec [degrees]
     * @param max_magnitude Maximum G-band magnitude
     * @param callback Optional progress callback
     * @return Number of stars downloaded
     * @throws GaiaException on network or I/O errors
     */
    size_t downloadBox(
        double ra_min,
        double ra_max,
        double dec_min,
        double dec_max,
        double max_magnitude,
        ProgressCallback callback = nullptr
    );
    
    /**
     * Query cached stars in a circular region
     * 
     * Returns stars from local cache only. If region is not fully cached,
     * throws GaiaException with ErrorCode::CACHE_ERROR.
     * 
     * @param ra_center Center RA [degrees]
     * @param dec_center Center Dec [degrees]
     * @param radius Search radius [degrees]
     * @param max_magnitude Maximum G-band magnitude
     * @return Vector of GaiaStar objects from cache
     * @throws GaiaException if region not cached
     */
    std::vector<GaiaStar> queryRegion(
        double ra_center,
        double dec_center,
        double radius,
        double max_magnitude
    );
    
    /**
     * Query cached stars in a rectangular region
     * 
     * @param ra_min Minimum RA [degrees]
     * @param ra_max Maximum RA [degrees]
     * @param dec_min Minimum Dec [degrees]
     * @param dec_max Maximum Dec [degrees]
     * @param max_magnitude Maximum G-band magnitude
     * @return Vector of GaiaStar objects from cache
     * @throws GaiaException if region not cached
     */
    std::vector<GaiaStar> queryBox(
        double ra_min,
        double ra_max,
        double dec_min,
        double dec_max,
        double max_magnitude
    );
    
    /**
     * Check if a point is covered by the cache
     * 
     * @param ra Right Ascension [degrees]
     * @param dec Declination [degrees]
     * @return true if the HEALPix tile containing this point is cached
     */
    bool isCovered(double ra, double dec) const;
    
    /**
     * Check if a circular region is fully covered by the cache
     * 
     * @param ra_center Center RA [degrees]
     * @param dec_center Center Dec [degrees]
     * @param radius Search radius [degrees]
     * @return true if all required HEALPix tiles are cached
     */
    bool isRegionCovered(double ra_center, double dec_center, double radius) const;
    
    /**
     * Clear all cached data
     * Removes all tile files and resets the index
     */
    void clear();
    
    /**
     * Remove specific tiles from cache
     * 
     * @param tile_ids Vector of HEALPix tile IDs to remove
     * @return Number of tiles removed
     */
    size_t removeTiles(const std::vector<int>& tile_ids);
    
    /**
     * Verify cache integrity
     * Checks that all index entries have corresponding files
     * 
     * @return true if cache is consistent
     */
    bool verify();
    
    /**
     * Get cache statistics
     * 
     * @return CacheStats structure with tile counts, sizes, etc.
     */
    CacheStats getStatistics() const;
    
    /**
     * Get list of all cached tile IDs
     * 
     * @return Vector of HEALPix tile IDs
     */
    std::vector<int> getCachedTiles() const;
    
    /**
     * Get cache directory path
     * 
     * @return Absolute path to cache directory
     */
    std::string getCacheDir() const;
    
    /**
     * Get HEALPix NSIDE parameter
     * 
     * @return NSIDE value
     */
    int getNside() const;

private:
    class Impl;
    std::unique_ptr<Impl> pImpl_;
    
    // Internal methods
    std::vector<int> getTileIdsForRegion(double ra, double dec, double radius) const;
    std::string getTilePath(int tile_id) const;
    void loadIndex();
    void saveIndex();
    std::vector<GaiaStar> loadTile(int tile_id);
    void saveTile(int tile_id, const std::vector<GaiaStar>& stars);
};

} // namespace gaia
} // namespace ioc
