#include <astdyn/catalog/gaia/gaia_mag18_catalog_v2.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <ctime>
#include <iostream>
#include <omp.h>

namespace ioc {
namespace gaia {

// HEALPix constants
static constexpr double PI = 3.14159265358979323846;
static constexpr double TWOPI = 2.0 * PI;
static constexpr double HALFPI = 0.5 * PI;
static constexpr double DEG2RAD = PI / 180.0;
static constexpr double RAD2DEG = 180.0 / PI;

Mag18CatalogV2::Mag18CatalogV2() 
    : file_(nullptr), 
      enable_parallel_(true),
      num_threads_(omp_get_max_threads()) {
    chunk_cache_.reserve(MAX_CACHED_CHUNKS);
    if (num_threads_ == 0) num_threads_ = 4;  // Fallback
}

Mag18CatalogV2::Mag18CatalogV2(const std::string& catalog_path)
    : Mag18CatalogV2() {
    if (!open(catalog_path)) {
        throw std::runtime_error("Failed to open catalog: " + catalog_path);
    }
}

Mag18CatalogV2::~Mag18CatalogV2() {
    close();
}

void Mag18CatalogV2::setParallelProcessing(bool enable, size_t num_threads) {
    enable_parallel_.store(enable);
    if (num_threads > 0) {
        num_threads_ = num_threads;
        omp_set_num_threads(num_threads);
    } else {
        num_threads_ = omp_get_max_threads();
        if (num_threads_ == 0) num_threads_ = 4;
    }
}

bool Mag18CatalogV2::open(const std::string& catalog_path) {
    catalog_path_ = catalog_path;
    
    file_ = fopen(catalog_path.c_str(), "rb");
    if (!file_) {
        std::cerr << "Failed to open catalog: " << catalog_path << std::endl;
        return false;
    }
    
    if (!loadHeader()) {
        fclose(file_);
        file_ = nullptr;
        return false;
    }
    
    if (!loadHEALPixIndex()) {
        fclose(file_);
        file_ = nullptr;
        return false;
    }
    
    if (!loadChunkIndex()) {
        fclose(file_);
        file_ = nullptr;
        return false;
    }
    
    std::cout << "✅ Opened Mag18 V2 catalog: " << header_.total_stars << " stars\n";
    std::cout << "   HEALPix: NSIDE=" << header_.healpix_nside 
              << ", " << header_.num_healpix_pixels << " pixels with data\n";
    std::cout << "   Chunks: " << header_.total_chunks 
              << " (" << header_.stars_per_chunk << " stars/chunk)\n";
    
    return true;
}

void Mag18CatalogV2::close() {
    if (file_) {
        fclose(file_);
        file_ = nullptr;
    }
    healpix_index_.clear();
    chunk_index_.clear();
    chunk_cache_.clear();
}

bool Mag18CatalogV2::loadHeader() {
    if (fread(&header_, sizeof(Mag18CatalogHeaderV2), 1, file_) != 1) {
        std::cerr << "Failed to read catalog header\n";
        return false;
    }
    
    if (strncmp(header_.magic, "GAIA18V2", 8) != 0) {
        std::cerr << "Invalid catalog magic. Expected GAIA18V2\n";
        return false;
    }
    
    if (header_.version != 2) {
        std::cerr << "Unsupported catalog version: " << header_.version << "\n";
        return false;
    }
    
    return true;
}

bool Mag18CatalogV2::loadHEALPixIndex() {
    if (header_.num_healpix_pixels == 0) {
        std::cerr << "No HEALPix pixels in catalog\n";
        return false;
    }
    
    healpix_index_.resize(header_.num_healpix_pixels);
    
    if (fseek(file_, header_.healpix_index_offset, SEEK_SET) != 0) {
        std::cerr << "Failed to seek to HEALPix index\n";
        return false;
    }
    
    if (fread(healpix_index_.data(), sizeof(HEALPixIndexEntry), 
              header_.num_healpix_pixels, file_) != header_.num_healpix_pixels) {
        std::cerr << "Failed to read HEALPix index\n";
        return false;
    }
    
    return true;
}

bool Mag18CatalogV2::loadChunkIndex() {
    if (header_.total_chunks == 0) {
        std::cerr << "No chunks in catalog\n";
        return false;
    }
    
    chunk_index_.resize(header_.total_chunks);
    
    if (fseek(file_, header_.chunk_index_offset, SEEK_SET) != 0) {
        std::cerr << "Failed to seek to chunk index\n";
        return false;
    }
    
    if (fread(chunk_index_.data(), sizeof(ChunkInfo), 
              header_.total_chunks, file_) != header_.total_chunks) {
        std::cerr << "Failed to read chunk index\n";
        return false;
    }
    
    return true;
}

uint32_t Mag18CatalogV2::ang2pix_nest(double theta, double phi) const {
    // Simplified HEALPix ang2pix_nest for NSIDE=64
    // theta: colatitude [0, pi], phi: longitude [0, 2*pi]
    
    const uint32_t nside = header_.healpix_nside;
    const uint32_t nl2 = 2 * nside;
    const double z = cos(theta);
    const double za = fabs(z);
    
    // Equatorial region
    if (za <= 2.0/3.0) {
        const double temp1 = nside * (0.5 + phi / TWOPI);
        const double temp2 = nside * z * 0.75;
        const uint32_t jp = static_cast<uint32_t>(temp1 - temp2);
        const uint32_t jm = static_cast<uint32_t>(temp1 + temp2);
        const uint32_t ir = nside + 1 + jp - jm;
        const uint32_t kshift = 1 - (ir & 1);
        const uint32_t ip = (jp + jm - nside + kshift + 1) / 2;
        
        const uint32_t face_num = ((jp < nside) ? 0 : 
                                   (jm < nside) ? 2 : 
                                   (jp >= nl2) ? 4 : 
                                   (jm >= nl2) ? 6 : 8) + (ir - nside);
        
        return face_num * nside * nside + (ip * nside + (ir - nside + 1));
    }
    
    // Polar caps
    const uint32_t ntt = static_cast<uint32_t>(phi * 4.0 / TWOPI);
    const double tp = phi - ntt * TWOPI / 4.0;
    const double tmp = sqrt(3.0 * (1.0 - za));
    const uint32_t jp = static_cast<uint32_t>(nside * tp / HALFPI * tmp);
    const uint32_t jm = static_cast<uint32_t>(nside * (1.0 - tp / HALFPI) * tmp);
    const uint32_t ir = jp + jm + 1;
    const uint32_t ip = ntt + 1;
    
    if (z > 0) {
        return 2 * ir * (ir - 1) + ip - 1;
    } else {
        const uint32_t npix = 12 * nside * nside;
        return npix - 2 * ir * (ir + 1) + ip - 1;
    }
}

uint32_t Mag18CatalogV2::getHEALPixPixel(double ra, double dec) const {
    // Convert RA/Dec to HEALPix pixel
    double theta = (90.0 - dec) * DEG2RAD;  // Colatitude
    double phi = ra * DEG2RAD;               // Longitude
    
    return ang2pix_nest(theta, phi);
}

void Mag18CatalogV2::query_disc(double theta, double phi, double radius,
                                 std::vector<uint32_t>& pixels) const {
    // Simplified disc query - returns pixels that intersect with cone
    // For production, use proper HEALPix query_disc from chealpix
    
    pixels.clear();
    const uint32_t nside = header_.healpix_nside;
    const uint32_t npix = 12 * nside * nside;
    
    // Brute force check all pixels (for NSIDE=64, this is 49,152 pixels)
    // In production, use proper hierarchical algorithm
    for (uint32_t pix = 0; pix < npix; ++pix) {
        double pix_theta, pix_phi;
        pix2ang_nest(pix, pix_theta, pix_phi);
        
        // Angular distance between disc center and pixel center
        double cos_dist = sin(theta) * sin(pix_theta) * cos(phi - pix_phi) +
                         cos(theta) * cos(pix_theta);
        double dist = acos(std::max(-1.0, std::min(1.0, cos_dist)));
        
        // Pixel angular size ~ sqrt(4*pi / npix)
        double pixel_size = sqrt(4.0 * PI / npix) * 1.5; // 1.5x for safety margin
        
        if (dist <= radius * DEG2RAD + pixel_size) {
            pixels.push_back(pix);
        }
    }
}

void Mag18CatalogV2::pix2ang_nest(uint32_t pixel, double& theta, double& phi) const {
    // Simplified HEALPix pix2ang_nest for NSIDE=64
    const uint32_t nside = header_.healpix_nside;
    const uint32_t npface = nside * nside;
    const uint32_t ncap = 2 * nside * (nside - 1);
    
    if (pixel < ncap) {
        // North polar cap
        const uint32_t iring = static_cast<uint32_t>(0.5 * (1 + sqrt(1 + 2 * pixel)));
        const uint32_t iphi = pixel - 2 * iring * (iring - 1);
        theta = acos(1.0 - iring * iring / (3.0 * npface));
        phi = (iphi + 0.5) * PI / (2.0 * iring);
    } else if (pixel < 12 * npface - ncap) {
        // Equatorial region
        const uint32_t ip = pixel - ncap;
        const uint32_t iring = ip / (4 * nside) + nside;
        const uint32_t iphi = ip % (4 * nside);
        theta = acos((2 * nside - iring) * 2.0 / (3.0 * nside));
        phi = (iphi + 0.5) * PI / (2.0 * nside);
    } else {
        // South polar cap
        const uint32_t ip = 12 * npface - pixel;
        const uint32_t iring = static_cast<uint32_t>(0.5 * (1 + sqrt(2 * ip - 1)));
        const uint32_t iphi = 2 * iring * (iring + 1) - ip - 1;
        theta = acos(-1.0 + iring * iring / (3.0 * npface));
        phi = (iphi + 0.5) * PI / (2.0 * iring);
    }
}

std::vector<uint32_t> Mag18CatalogV2::getPixelsInCone(double ra, double dec, 
                                                        double radius) const {
    double theta = (90.0 - dec) * DEG2RAD;
    double phi = ra * DEG2RAD;
    
    std::vector<uint32_t> pixels;
    query_disc(theta, phi, radius * DEG2RAD, pixels);  // BUG FIX: Convert radius to radians
    
    return pixels;
}

std::vector<Mag18RecordV2> Mag18CatalogV2::readChunk(uint64_t chunk_id) {
    // Check cache first (shared read lock)
    {
        std::shared_lock<std::shared_mutex> lock(cache_mutex_);
        for (auto& cached : chunk_cache_) {
            if (cached.chunk_id == chunk_id) {
                cached.access_count.fetch_add(1, std::memory_order_relaxed);
                return cached.records;
            }
        }
    }
    
    // Not in cache, load from file
    if (chunk_id >= header_.total_chunks) {
        return {};
    }
    
    const ChunkInfo& chunk = chunk_index_[chunk_id];
    
    // Read compressed data (exclusive file lock)
    std::vector<uint8_t> compressed(chunk.compressed_size);
    {
        std::lock_guard<std::mutex> file_lock(file_mutex_);
        if (fseek(file_, chunk.file_offset, SEEK_SET) != 0) {
            return {};
        }
        if (fread(compressed.data(), 1, chunk.compressed_size, file_) != chunk.compressed_size) {
            return {};
        }
    }
    
    // Decompress
    std::vector<uint8_t> uncompressed(chunk.uncompressed_size);
    uLongf dest_len = chunk.uncompressed_size;
    if (uncompress(uncompressed.data(), &dest_len, compressed.data(), chunk.compressed_size) != Z_OK) {
        return {};
    }
    
    // Parse records
    std::vector<Mag18RecordV2> records(chunk.num_stars);
    memcpy(records.data(), uncompressed.data(), chunk.num_stars * sizeof(Mag18RecordV2));
    
    // Add to cache (exclusive write lock)
    {
        std::unique_lock<std::shared_mutex> lock(cache_mutex_);
        if (chunk_cache_.size() >= MAX_CACHED_CHUNKS) {
            // Evict least recently used
            auto lru = std::min_element(chunk_cache_.begin(), chunk_cache_.end(),
                [](const ChunkCache& a, const ChunkCache& b) {
                    return a.access_count.load() < b.access_count.load();
                });
            // Move assign (explicit for atomic members)
            lru->chunk_id = chunk_id;
            lru->records = std::move(records);
            lru->access_count.store(1);
        } else {
            chunk_cache_.emplace_back(chunk_id, std::move(records), 1);
        }
    }
    
    return records;
}

std::optional<Mag18RecordV2> Mag18CatalogV2::readRecord(uint64_t index) {
    if (index >= header_.total_stars) {
        return std::nullopt;
    }
    
    // Find chunk containing this record
    uint64_t chunk_id = index / header_.stars_per_chunk;
    uint64_t offset_in_chunk = index % header_.stars_per_chunk;
    
    auto records = readChunk(chunk_id);
    if (records.empty() || offset_in_chunk >= records.size()) {
        return std::nullopt;
    }
    
    return records[offset_in_chunk];
}

std::optional<GaiaStar> Mag18CatalogV2::queryBySourceId(uint64_t source_id) {
    // Binary search (catalog must be sorted by source_id)
    uint64_t left = 0;
    uint64_t right = header_.total_stars - 1;
    
    while (left <= right) {
        uint64_t mid = left + (right - left) / 2;
        auto record = readRecord(mid);
        
        if (!record) {
            return std::nullopt;
        }
        
        if (record->source_id == source_id) {
            return recordToStar(*record);
        } else if (record->source_id < source_id) {
            left = mid + 1;
        } else {
            if (mid == 0) break;
            right = mid - 1;
        }
    }
    
    return std::nullopt;
}

std::vector<GaiaStar> Mag18CatalogV2::queryCone(double ra, double dec, double radius,
                                                  size_t max_results) {
    // Get HEALPix pixels intersecting cone
    auto pixels = getPixelsInCone(ra, dec, radius);
    
    // If parallel processing is disabled or few pixels, use sequential
    if (!enable_parallel_.load() || pixels.size() < 4) {
        std::vector<GaiaStar> results;
        for (uint32_t pixel : pixels) {
            // Find pixel in index
            auto it = std::lower_bound(healpix_index_.begin(), healpix_index_.end(), pixel,
                [](const HEALPixIndexEntry& entry, uint32_t pix) {
                    return entry.pixel_id < pix;
                });
            
            if (it == healpix_index_.end() || it->pixel_id != pixel) {
                continue;
            }
            
            // Scan all stars in pixel
            for (uint32_t i = 0; i < it->num_stars; ++i) {
                auto record = readRecord(it->first_star_idx + i);
                if (!record) continue;
                
                double dist = angularDistance(ra, dec, record->ra, record->dec);
                if (dist <= radius) {
                    results.push_back(recordToStar(*record));
                    
                    if (max_results > 0 && results.size() >= max_results) {
                        return results;
                    }
                }
            }
        }
        return results;
    }
    
    // PARALLEL PROCESSING with OpenMP
    std::vector<GaiaStar> results;
    std::mutex results_mutex;
    std::atomic<bool> limit_reached(false);
    
    #pragma omp parallel if(enable_parallel_.load())
    {
        std::vector<GaiaStar> thread_results;
        
        #pragma omp for schedule(dynamic)
        for (size_t p = 0; p < pixels.size(); ++p) {
            if (limit_reached.load()) continue;
            
            uint32_t pixel = pixels[p];
            
            // Find pixel in index
            auto it = std::lower_bound(healpix_index_.begin(), healpix_index_.end(), pixel,
                [](const HEALPixIndexEntry& entry, uint32_t pix) {
                    return entry.pixel_id < pix;
                });
            
            if (it == healpix_index_.end() || it->pixel_id != pixel) {
                continue;
            }
            
            // Scan all stars in pixel
            for (uint32_t i = 0; i < it->num_stars; ++i) {
                if (limit_reached.load()) break;
                
                auto record = readRecord(it->first_star_idx + i);
                if (!record) continue;
                
                double dist = angularDistance(ra, dec, record->ra, record->dec);
                if (dist <= radius) {
                    thread_results.push_back(recordToStar(*record));
                }
            }
        }
        
        // Merge thread results (critical section)
        if (!thread_results.empty()) {
            #pragma omp critical
            {
                results.insert(results.end(), thread_results.begin(), thread_results.end());
                
                if (max_results > 0 && results.size() >= max_results) {
                    limit_reached.store(true);
                }
            }
        }
    }
    
    // Truncate if needed
    if (max_results > 0 && results.size() > max_results) {
        results.resize(max_results);
    }
    
    return results;
}

std::vector<GaiaStar> Mag18CatalogV2::queryConeWithMagnitude(double ra, double dec, double radius,
                                                               double mag_min, double mag_max,
                                                               size_t max_results) {
    auto pixels = getPixelsInCone(ra, dec, radius);
    
    // Sequential for small queries
    if (!enable_parallel_.load() || pixels.size() < 4) {
        std::vector<GaiaStar> results;
        
        for (uint32_t pixel : pixels) {
            auto it = std::lower_bound(healpix_index_.begin(), healpix_index_.end(), pixel,
                [](const HEALPixIndexEntry& entry, uint32_t pix) {
                    return entry.pixel_id < pix;
                });
            
            if (it == healpix_index_.end() || it->pixel_id != pixel) continue;
            
            for (uint32_t i = 0; i < it->num_stars; ++i) {
                auto record = readRecord(it->first_star_idx + i);
                if (!record) continue;
                
                if (record->g_mag < mag_min || record->g_mag > mag_max) continue;
                
                double dist = angularDistance(ra, dec, record->ra, record->dec);
                if (dist <= radius) {
                    results.push_back(recordToStar(*record));
                    
                    if (max_results > 0 && results.size() >= max_results) {
                        return results;
                    }
                }
            }
        }
        return results;
    }
    
    // PARALLEL PROCESSING with OpenMP
    std::vector<GaiaStar> results;
    std::atomic<bool> limit_reached(false);
    
    #pragma omp parallel if(enable_parallel_.load())
    {
        std::vector<GaiaStar> thread_results;
        
        #pragma omp for schedule(dynamic)
        for (size_t p = 0; p < pixels.size(); ++p) {
            if (limit_reached.load()) continue;
            
            uint32_t pixel = pixels[p];
            
            auto it = std::lower_bound(healpix_index_.begin(), healpix_index_.end(), pixel,
                [](const HEALPixIndexEntry& entry, uint32_t pix) {
                    return entry.pixel_id < pix;
                });
            
            if (it == healpix_index_.end() || it->pixel_id != pixel) continue;
            
            for (uint32_t i = 0; i < it->num_stars; ++i) {
                if (limit_reached.load()) break;
                
                auto record = readRecord(it->first_star_idx + i);
                if (!record) continue;
                
                if (record->g_mag < mag_min || record->g_mag > mag_max) continue;
                
                double dist = angularDistance(ra, dec, record->ra, record->dec);
                if (dist <= radius) {
                    thread_results.push_back(recordToStar(*record));
                }
            }
        }
        
        // Merge results (critical section)
        if (!thread_results.empty()) {
            #pragma omp critical
            {
                results.insert(results.end(), thread_results.begin(), thread_results.end());
                
                if (max_results > 0 && results.size() >= max_results) {
                    limit_reached.store(true);
                }
            }
        }
    }
    
    if (max_results > 0 && results.size() > max_results) {
        results.resize(max_results);
    }
    
    return results;
}

std::vector<GaiaStar> Mag18CatalogV2::queryBrightest(double ra, double dec, double radius,
                                                       size_t num_stars) {
    auto all_stars = queryCone(ra, dec, radius, 0);
    
    // Partial sort to get N brightest
    if (all_stars.size() > num_stars) {
        std::partial_sort(all_stars.begin(), all_stars.begin() + num_stars, all_stars.end(),
            [](const GaiaStar& a, const GaiaStar& b) {
                return a.phot_g_mean_mag < b.phot_g_mean_mag;
            });
        all_stars.resize(num_stars);
    } else {
        std::sort(all_stars.begin(), all_stars.end(),
            [](const GaiaStar& a, const GaiaStar& b) {
                return a.phot_g_mean_mag < b.phot_g_mean_mag;
            });
    }
    
    return all_stars;
}

size_t Mag18CatalogV2::countInCone(double ra, double dec, double radius) {
    size_t count = 0;
    auto pixels = getPixelsInCone(ra, dec, radius);
    
    for (uint32_t pixel : pixels) {
        auto it = std::lower_bound(healpix_index_.begin(), healpix_index_.end(), pixel,
            [](const HEALPixIndexEntry& entry, uint32_t pix) {
                return entry.pixel_id < pix;
            });
        
        if (it == healpix_index_.end() || it->pixel_id != pixel) continue;
        
        for (uint32_t i = 0; i < it->num_stars; ++i) {
            auto record = readRecord(it->first_star_idx + i);
            if (!record) continue;
            
            double dist = angularDistance(ra, dec, record->ra, record->dec);
            if (dist <= radius) {
                count++;
            }
        }
    }
    
    return count;
}

std::optional<Mag18RecordV2> Mag18CatalogV2::getExtendedRecord(uint64_t source_id) {
    // Binary search
    uint64_t left = 0;
    uint64_t right = header_.total_stars - 1;
    
    while (left <= right) {
        uint64_t mid = left + (right - left) / 2;
        auto record = readRecord(mid);
        
        if (!record) {
            return std::nullopt;
        }
        
        if (record->source_id == source_id) {
            return record;
        } else if (record->source_id < source_id) {
            left = mid + 1;
        } else {
            if (mid == 0) break;
            right = mid - 1;
        }
    }
    
    return std::nullopt;
}

GaiaStar Mag18CatalogV2::recordToStar(const Mag18RecordV2& record) const {
    GaiaStar star;
    star.source_id = record.source_id;
    star.ra = record.ra;
    star.dec = record.dec;
    star.parallax = record.parallax;
    star.pmra = record.pmra;
    star.pmdec = record.pmdec;
    star.phot_g_mean_mag = record.g_mag;
    star.phot_bp_mean_mag = record.bp_mag;
    star.phot_rp_mean_mag = record.rp_mag;
    star.bp_rp = record.bp_rp;
    star.ruwe = record.ruwe;
    return star;
}

double Mag18CatalogV2::angularDistance(double ra1, double dec1, double ra2, double dec2) const {
    // Haversine formula
    double dra = (ra2 - ra1) * DEG2RAD;
    double ddec = (dec2 - dec1) * DEG2RAD;
    double a = sin(ddec/2) * sin(ddec/2) +
               cos(dec1 * DEG2RAD) * cos(dec2 * DEG2RAD) *
               sin(dra/2) * sin(dra/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return c * RAD2DEG;
}

} // namespace gaia
} // namespace ioc
