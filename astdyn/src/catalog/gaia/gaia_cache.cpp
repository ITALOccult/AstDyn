#include <astdyn/catalog/gaia/gaia_cache.h>
#include <astdyn/catalog/gaia/gaia_client.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <functional>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

namespace fs = std::filesystem;

namespace ioc {
namespace gaia {

// =============================================================================
// HEALPix Utilities (simplified, ring scheme)
// =============================================================================

namespace healpix {

/**
 * Convert (RA, Dec) to HEALPix tile ID (ring scheme)
 * Simplified implementation for NSIDE <= 256
 */
int ang2pix(int nside, double ra, double dec) {
    constexpr double DEG_TO_RAD = M_PI / 180.0;
    
    // Convert to theta (colatitude), phi
    double theta = (90.0 - dec) * DEG_TO_RAD;  // [0, π]
    double phi = ra * DEG_TO_RAD;              // [0, 2π)
    
    double z = std::cos(theta);
    int npix = 12 * nside * nside;
    
    // Equatorial region
    double za = std::abs(z);
    double tt = (2.0 / M_PI) * std::atan(za / std::sqrt(1.0 - za * za));
    
    if (za <= 2.0 / 3.0) {  // Equatorial region
        int jp = static_cast<int>(nside * (0.5 + tt - z * 0.75));
        int jm = static_cast<int>(nside * (0.5 + tt + z * 0.75));
        
        int ir = nside + 1 + jp - jm;  // Ring number
        int kshift = 0;
        if ((ir & 1) == 0) kshift = 1;
        
        int ip = static_cast<int>((phi / (2.0 * M_PI)) * 4 * ir + kshift) / 2;
        ip = ip % (4 * ir);
        
        return 2 * ir * (ir - 1) + ip;
    }
    
    // Polar caps
    int jp = static_cast<int>(nside * tt * std::sqrt(3.0 * (1.0 - za)));
    int jm = static_cast<int>(nside * tt * std::sqrt(3.0 * (1.0 - za)));
    
    int ir = jp + jm + 1;
    int ip = static_cast<int>((phi / (2.0 * M_PI)) * 4 * ir);
    ip = ip % (4 * ir);
    
    if (z > 0) {
        return 2 * ir * (ir - 1) + ip;  // North cap
    } else {
        return npix - 2 * ir * (ir + 1) + ip;  // South cap
    }
}

/**
 * Get all tiles within a circular region
 */
std::vector<int> query_disc(int nside, double ra_center, double dec_center, double radius) {
    std::unordered_set<int> tiles;
    
    // Sample points around the circle to find all tiles
    int n_samples = std::max(10, static_cast<int>(radius * 20));
    
    for (int i = 0; i < n_samples; ++i) {
        double angle = 2.0 * M_PI * i / n_samples;
        double dr = radius * std::cos(angle);
        double dd = radius * std::sin(angle);
        
        // Approximate offset (good for small radii)
        double ra = ra_center + dr / std::cos(dec_center * M_PI / 180.0);
        double dec = dec_center + dd;
        
        tiles.insert(ang2pix(nside, ra, dec));
    }
    
    // Add center tile
    tiles.insert(ang2pix(nside, ra_center, dec_center));
    
    return std::vector<int>(tiles.begin(), tiles.end());
}

/**
 * Get all tiles in a box region
 */
std::vector<int> query_box(int nside, double ra_min, double ra_max, 
                          double dec_min, double dec_max) {
    std::unordered_set<int> tiles;
    
    // Sample grid points
    int n_ra = std::max(5, static_cast<int>((ra_max - ra_min) * 5));
    int n_dec = std::max(5, static_cast<int>((dec_max - dec_min) * 5));
    
    for (int i = 0; i <= n_ra; ++i) {
        for (int j = 0; j <= n_dec; ++j) {
            double ra = ra_min + (ra_max - ra_min) * i / n_ra;
            double dec = dec_min + (dec_max - dec_min) * j / n_dec;
            tiles.insert(ang2pix(nside, ra, dec));
        }
    }
    
    return std::vector<int>(tiles.begin(), tiles.end());
}

} // namespace healpix

// =============================================================================
// GaiaCache::Impl - Private implementation
// =============================================================================

class GaiaCache::Impl {
public:
    std::string cache_dir_;
    int nside_;
    GaiaRelease release_;
    std::unordered_set<int> cached_tiles_;
    bool index_modified_;
    GaiaClient* client_;
    
    Impl(const std::string& cache_dir, int nside, GaiaRelease release)
        : cache_dir_(cache_dir), nside_(nside), release_(release), 
          index_modified_(false), client_(nullptr) {
        
        // Create cache directory structure
        fs::create_directories(cache_dir_);
        fs::create_directories(fs::path(cache_dir_) / "tiles");
    }
    
    std::string getTilePath(int tile_id) const {
        std::ostringstream oss;
        oss << cache_dir_ << "/tiles/tile_" 
            << std::setw(5) << std::setfill('0') << tile_id << ".json";
        return oss.str();
    }
    
    std::string getIndexPath() const {
        return cache_dir_ + "/index.json";
    }
    
    void loadIndex() {
        std::string index_path = getIndexPath();
        if (!fs::exists(index_path)) {
            return;  // Empty cache
        }
        
        std::ifstream ifs(index_path);
        if (!ifs) {
            throw GaiaException(ErrorCode::CACHE_ERROR, 
                              "Cannot open index file: " + index_path);
        }
        
        // Simple JSON parsing (tile IDs only)
        std::string line;
        while (std::getline(ifs, line)) {
            if (line.find("\"tile_id\":") != std::string::npos) {
                size_t pos = line.find(':');
                if (pos != std::string::npos) {
                    int tile_id = std::stoi(line.substr(pos + 1));
                    cached_tiles_.insert(tile_id);
                }
            }
        }
    }
    
    void saveIndex() {
        std::string index_path = getIndexPath();
        std::ofstream ofs(index_path);
        if (!ofs) {
            throw GaiaException(ErrorCode::CACHE_ERROR, 
                              "Cannot write index file: " + index_path);
        }
        
        ofs << "{\n";
        ofs << "  \"nside\": " << nside_ << ",\n";
        ofs << "  \"release\": \"" << releaseToString(release_) << "\",\n";
        ofs << "  \"tiles\": [\n";
        
        bool first = true;
        for (int tile_id : cached_tiles_) {
            if (!first) ofs << ",\n";
            ofs << "    {\"tile_id\": " << tile_id << "}";
            first = false;
        }
        
        ofs << "\n  ]\n";
        ofs << "}\n";
        
        index_modified_ = false;
    }
    
    std::vector<GaiaStar> loadTile(int tile_id) {
        std::string tile_path = getTilePath(tile_id);
        std::ifstream ifs(tile_path);
        if (!ifs) {
            throw GaiaException(ErrorCode::CACHE_ERROR, 
                              "Cannot open tile file: " + tile_path);
        }
        
        std::vector<GaiaStar> stars;
        
        // Parse JSON (simplified - expects one star per line)
        std::string line;
        while (std::getline(ifs, line)) {
            if (line.find("\"source_id\":") != std::string::npos) {
                GaiaStar star;
                
                // Extract values (very simplified parser)
                auto extract = [&line](const std::string& key) -> double {
                    size_t pos = line.find("\"" + key + "\":");
                    if (pos == std::string::npos) return 0.0;
                    pos = line.find(':', pos) + 1;
                    return std::stod(line.substr(pos));
                };
                
                star.source_id = static_cast<int64_t>(extract("source_id"));
                star.ra = extract("ra");
                star.dec = extract("dec");
                star.parallax = extract("parallax");
                star.parallax_error = extract("parallax_error");
                star.pmra = extract("pmra");
                star.pmdec = extract("pmdec");
                star.phot_g_mean_mag = extract("phot_g_mean_mag");
                star.phot_bp_mean_mag = extract("phot_bp_mean_mag");
                star.phot_rp_mean_mag = extract("phot_rp_mean_mag");
                
                stars.push_back(star);
            }
        }
        
        return stars;
    }
    
    void saveTile(int tile_id, const std::vector<GaiaStar>& stars) {
        std::string tile_path = getTilePath(tile_id);
        std::ofstream ofs(tile_path);
        if (!ofs) {
            throw GaiaException(ErrorCode::CACHE_ERROR, 
                              "Cannot write tile file: " + tile_path);
        }
        
        ofs << "{\n  \"tile_id\": " << tile_id << ",\n";
        ofs << "  \"stars\": [\n";
        
        for (size_t i = 0; i < stars.size(); ++i) {
            const auto& s = stars[i];
            ofs << "    {";
            ofs << "\"source_id\":" << s.source_id << ",";
            ofs << "\"ra\":" << s.ra << ",";
            ofs << "\"dec\":" << s.dec << ",";
            ofs << "\"parallax\":" << s.parallax << ",";
            ofs << "\"parallax_error\":" << s.parallax_error << ",";
            ofs << "\"pmra\":" << s.pmra << ",";
            ofs << "\"pmdec\":" << s.pmdec << ",";
            ofs << "\"phot_g_mean_mag\":" << s.phot_g_mean_mag << ",";
            ofs << "\"phot_bp_mean_mag\":" << s.phot_bp_mean_mag << ",";
            ofs << "\"phot_rp_mean_mag\":" << s.phot_rp_mean_mag;
            ofs << "}";
            if (i < stars.size() - 1) ofs << ",";
            ofs << "\n";
        }
        
        ofs << "  ]\n}\n";
        
        cached_tiles_.insert(tile_id);
        index_modified_ = true;
    }
};

// =============================================================================
// GaiaCache Public Interface
// =============================================================================

GaiaCache::GaiaCache(const std::string& cache_dir, int nside, GaiaRelease release)
    : pImpl_(std::make_unique<Impl>(cache_dir, nside, release)) {
    pImpl_->loadIndex();
}

GaiaCache::~GaiaCache() {
    if (pImpl_ && pImpl_->index_modified_) {
        pImpl_->saveIndex();
    }
}

GaiaCache::GaiaCache(GaiaCache&&) noexcept = default;
GaiaCache& GaiaCache::operator=(GaiaCache&&) noexcept = default;

void GaiaCache::setClient(GaiaClient* client) {
    pImpl_->client_ = client;
}

size_t GaiaCache::downloadRegion(
    double ra_center,
    double dec_center,
    double radius,
    double max_magnitude,
    ProgressCallback callback) {
    
    // Get tiles needed
    auto tile_ids = healpix::query_disc(pImpl_->nside_, ra_center, dec_center, radius);
    
    size_t total_stars = 0;
    size_t tiles_processed = 0;
    
    for (int tile_id : tile_ids) {
        // Skip if already cached
        if (pImpl_->cached_tiles_.count(tile_id)) {
            tiles_processed++;
            if (callback) {
                int progress = (100 * tiles_processed) / tile_ids.size();
                callback(progress, "Tile " + std::to_string(tile_id) + " already cached");
            }
            continue;
        }
        
        // Download from Gaia TAP service via client
        std::vector<GaiaStar> stars;
        
        if (pImpl_->client_) {
            try {
                // Query the entire region for this download
                // (in a production version, we'd query tile-by-tile more efficiently)
                stars = pImpl_->client_->queryCone(ra_center, dec_center, radius, max_magnitude);
                
                if (callback) {
                    int progress = (100 * tiles_processed) / tile_ids.size();
                    callback(progress, "Downloaded " + std::to_string(stars.size()) + 
                            " stars for tile " + std::to_string(tile_id));
                }
            } catch (const GaiaException& e) {
                if (callback) {
                    callback(-1, "Error downloading tile " + std::to_string(tile_id) + 
                            ": " + std::string(e.what()));
                }
                throw;
            }
        } else {
            // No client set - create empty tile
            if (callback) {
                callback(-1, "Warning: No client set, creating empty tile " + 
                        std::to_string(tile_id));
            }
        }
        
        pImpl_->saveTile(tile_id, stars);
        total_stars += stars.size();
        tiles_processed++;
        
        if (callback) {
            int progress = (100 * tiles_processed) / tile_ids.size();
            callback(progress, "Saved tile " + std::to_string(tile_id));
        }
    }
    
    return total_stars;
}

size_t GaiaCache::downloadBox(
    double ra_min,
    double ra_max,
    double dec_min,
    double dec_max,
    double max_magnitude,
    ProgressCallback callback) {
    
    auto tile_ids = healpix::query_box(pImpl_->nside_, ra_min, ra_max, 
                                       dec_min, dec_max);
    
    size_t total_stars = 0;
    size_t tiles_processed = 0;
    
    for (int tile_id : tile_ids) {
        if (pImpl_->cached_tiles_.count(tile_id)) {
            tiles_processed++;
            if (callback) {
                int progress = (100 * tiles_processed) / tile_ids.size();
                callback(progress, "Tile " + std::to_string(tile_id) + " already cached");
            }
            continue;
        }
        
        // Download from Gaia TAP service via client
        std::vector<GaiaStar> stars;
        
        if (pImpl_->client_) {
            try {
                stars = pImpl_->client_->queryBox(ra_min, ra_max, dec_min, dec_max, max_magnitude);
                
                if (callback) {
                    int progress = (100 * tiles_processed) / tile_ids.size();
                    callback(progress, "Downloaded " + std::to_string(stars.size()) + 
                            " stars for tile " + std::to_string(tile_id));
                }
            } catch (const GaiaException& e) {
                if (callback) {
                    callback(-1, "Error downloading tile " + std::to_string(tile_id) + 
                            ": " + std::string(e.what()));
                }
                throw;
            }
        }
        
        pImpl_->saveTile(tile_id, stars);
        total_stars += stars.size();
        tiles_processed++;
        
        if (callback) {
            int progress = (100 * tiles_processed) / tile_ids.size();
            callback(progress, "Saved tile " + std::to_string(tile_id));
        }
    }
    
    return total_stars;
}

std::vector<GaiaStar> GaiaCache::queryRegion(
    double ra_center,
    double dec_center,
    double radius,
    double max_magnitude) {
    
    // Check if region is covered
    if (!isRegionCovered(ra_center, dec_center, radius)) {
        throw GaiaException(ErrorCode::CACHE_ERROR, 
                          "Region not fully cached. Call downloadRegion() first.");
    }
    
    auto tile_ids = healpix::query_disc(pImpl_->nside_, ra_center, dec_center, radius);
    
    std::vector<GaiaStar> results;
    EquatorialCoordinates center(ra_center, dec_center);
    
    for (int tile_id : tile_ids) {
        auto stars = pImpl_->loadTile(tile_id);
        
        for (const auto& star : stars) {
            EquatorialCoordinates pos(star.ra, star.dec);
            double dist = angularDistance(center, pos);
            
            if (dist <= radius && star.phot_g_mean_mag <= max_magnitude) {
                results.push_back(star);
            }
        }
    }
    
    return results;
}

std::vector<GaiaStar> GaiaCache::queryBox(
    double ra_min,
    double ra_max,
    double dec_min,
    double dec_max,
    double max_magnitude) {
    
    auto tile_ids = healpix::query_box(pImpl_->nside_, ra_min, ra_max, 
                                       dec_min, dec_max);
    
    std::vector<GaiaStar> results;
    
    for (int tile_id : tile_ids) {
        if (!pImpl_->cached_tiles_.count(tile_id)) {
            throw GaiaException(ErrorCode::CACHE_ERROR, 
                              "Region not fully cached. Call downloadBox() first.");
        }
        
        auto stars = pImpl_->loadTile(tile_id);
        
        for (const auto& star : stars) {
            if (star.ra >= ra_min && star.ra <= ra_max &&
                star.dec >= dec_min && star.dec <= dec_max &&
                star.phot_g_mean_mag <= max_magnitude) {
                results.push_back(star);
            }
        }
    }
    
    return results;
}

bool GaiaCache::isCovered(double ra, double dec) const {
    int tile_id = healpix::ang2pix(pImpl_->nside_, ra, dec);
    return pImpl_->cached_tiles_.count(tile_id) > 0;
}

bool GaiaCache::isRegionCovered(double ra_center, double dec_center, double radius) const {
    auto tile_ids = healpix::query_disc(pImpl_->nside_, ra_center, dec_center, radius);
    
    for (int tile_id : tile_ids) {
        if (pImpl_->cached_tiles_.count(tile_id) == 0) {
            return false;
        }
    }
    return true;
}

void GaiaCache::clear() {
    // Remove all tile files
    fs::path tiles_dir = fs::path(pImpl_->cache_dir_) / "tiles";
    for (const auto& entry : fs::directory_iterator(tiles_dir)) {
        fs::remove(entry);
    }
    
    pImpl_->cached_tiles_.clear();
    pImpl_->index_modified_ = true;
    pImpl_->saveIndex();
}

size_t GaiaCache::removeTiles(const std::vector<int>& tile_ids) {
    size_t removed = 0;
    
    for (int tile_id : tile_ids) {
        if (pImpl_->cached_tiles_.count(tile_id)) {
            std::string tile_path = pImpl_->getTilePath(tile_id);
            fs::remove(tile_path);
            pImpl_->cached_tiles_.erase(tile_id);
            removed++;
        }
    }
    
    if (removed > 0) {
        pImpl_->index_modified_ = true;
    }
    
    return removed;
}

bool GaiaCache::verify() {
    for (int tile_id : pImpl_->cached_tiles_) {
        if (!fs::exists(pImpl_->getTilePath(tile_id))) {
            return false;
        }
    }
    return true;
}

CacheStats GaiaCache::getStatistics() const {
    CacheStats stats;
    stats.total_tiles = 12 * pImpl_->nside_ * pImpl_->nside_;
    stats.cached_tiles = pImpl_->cached_tiles_.size();
    stats.last_update = std::chrono::system_clock::now();
    
    // Calculate disk size and star count
    fs::path tiles_dir = fs::path(pImpl_->cache_dir_) / "tiles";
    for (const auto& entry : fs::directory_iterator(tiles_dir)) {
        stats.disk_size_bytes += fs::file_size(entry);
    }
    
    return stats;
}

std::vector<int> GaiaCache::getCachedTiles() const {
    return std::vector<int>(pImpl_->cached_tiles_.begin(), 
                           pImpl_->cached_tiles_.end());
}

std::string GaiaCache::getCacheDir() const {
    return pImpl_->cache_dir_;
}

int GaiaCache::getNside() const {
    return pImpl_->nside_;
}

} // namespace gaia
} // namespace ioc

#ifdef __clang__
#pragma clang diagnostic pop
#endif
