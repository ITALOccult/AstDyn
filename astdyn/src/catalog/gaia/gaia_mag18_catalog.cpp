#include <astdyn/catalog/gaia/gaia_mag18_catalog.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <cstring>

namespace fs = std::filesystem;

namespace ioc_gaialib {

GaiaMag18Catalog::GaiaMag18Catalog(const std::string& catalog_file)
    : catalog_file_(catalog_file), loaded_(false), gz_file_(nullptr) {
    loaded_ = loadHeader();
}

GaiaMag18Catalog::~GaiaMag18Catalog() {
    if (gz_file_) {
        gzclose(gz_file_);
        gz_file_ = nullptr;
    }
}

bool GaiaMag18Catalog::loadHeader() {
    // Open gzip file
    gz_file_ = gzopen(catalog_file_.c_str(), "rb");
    if (!gz_file_) {
        return false;
    }
    
    // Read header
    int bytes_read = gzread(gz_file_, &header_, sizeof(header_));
    if (bytes_read != sizeof(header_)) {
        gzclose(gz_file_);
        gz_file_ = nullptr;
        return false;
    }
    
    // Verify magic
    if (std::strncmp(header_.magic, "GAIA18", 6) != 0) {
        gzclose(gz_file_);
        gz_file_ = nullptr;
        return false;
    }
    
    // Verify version
    if (header_.version != 1) {
        gzclose(gz_file_);
        gz_file_ = nullptr;
        return false;
    }
    
    return true;
}

GaiaMag18Catalog::Statistics GaiaMag18Catalog::getStatistics() const {
    Statistics stats;
    stats.total_stars = header_.total_stars;
    stats.mag_limit = header_.mag_limit;
    stats.uncompressed_size = header_.data_size;
    
    if (fs::exists(catalog_file_)) {
        stats.file_size = fs::file_size(catalog_file_);
        if (stats.uncompressed_size > 0) {
            stats.compression_ratio = 100.0 * 
                (1.0 - (double)stats.file_size / stats.uncompressed_size);
        } else {
            stats.compression_ratio = 0.0;
        }
    } else {
        stats.file_size = 0;
        stats.compression_ratio = 0.0;
    }
    
    return stats;
}

std::optional<Mag18Record> GaiaMag18Catalog::readRecord(uint64_t index) const {
    if (!loaded_ || !gz_file_) {
        return std::nullopt;
    }
    
    if (index >= header_.total_stars) {
        return std::nullopt;
    }
    
    // Calculate position in file
    z_off_t pos = header_.data_offset + (index * sizeof(Mag18Record));
    
    // Seek to position
    if (gzseek(gz_file_, pos, SEEK_SET) < 0) {
        return std::nullopt;
    }
    
    // Read record
    Mag18Record record;
    int bytes_read = gzread(gz_file_, &record, sizeof(record));
    if (bytes_read != sizeof(record)) {
        return std::nullopt;
    }
    
    return record;
}

ioc::gaia::GaiaStar GaiaMag18Catalog::recordToStar(const Mag18Record& record) const {
    ioc::gaia::GaiaStar star;
    star.source_id = record.source_id;
    star.ra = record.ra;
    star.dec = record.dec;
    star.phot_g_mean_mag = record.g_mag;
    star.phot_bp_mean_mag = record.bp_mag;
    star.phot_rp_mean_mag = record.rp_mag;
    star.parallax = record.parallax;
    
    // Set other fields to default/unknown
    star.pmra = 0.0;
    star.pmdec = 0.0;
    star.parallax_error = 0.0;
    star.pmra_error = 0.0;
    star.pmdec_error = 0.0;
    
    return star;
}

std::optional<ioc::gaia::GaiaStar> GaiaMag18Catalog::queryBySourceId(uint64_t source_id) const {
    if (!loaded_) {
        return std::nullopt;
    }
    
    // Binary search (catalog is sorted by source_id)
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

double GaiaMag18Catalog::angularDistance(double ra1, double dec1,
                                         double ra2, double dec2) const {
    // Convert to radians
    const double deg2rad = M_PI / 180.0;
    ra1 *= deg2rad;
    dec1 *= deg2rad;
    ra2 *= deg2rad;
    dec2 *= deg2rad;
    
    // Haversine formula
    double dra = ra2 - ra1;
    double ddec = dec2 - dec1;
    
    double a = std::sin(ddec/2) * std::sin(ddec/2) +
               std::cos(dec1) * std::cos(dec2) *
               std::sin(dra/2) * std::sin(dra/2);
    
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1-a));
    
    return c / deg2rad;  // Convert back to degrees
}

std::vector<ioc::gaia::GaiaStar> GaiaMag18Catalog::queryCone(double ra, double dec,
                                                              double radius,
                                                              size_t max_results) const {
    std::vector<ioc::gaia::GaiaStar> results;
    
    if (!loaded_) {
        return results;
    }
    
    // Reserve space
    results.reserve(max_results > 0 ? max_results : 1000);
    
    std::cout << "[GaiaMag18Catalog] DEBUG: queryCone RA=" << ra << " Dec=" << dec 
              << " Radius=" << radius << " TotalStars=" << header_.total_stars << "\n";

    // Linear scan through all stars
    for (uint64_t i = 0; i < header_.total_stars; ++i) {
        auto record = readRecord(i);
        if (!record) {
            continue;
        }
        
        // Check distance
        double dist = angularDistance(ra, dec, record->ra, record->dec);
        if (dist <= radius) {
            results.push_back(recordToStar(*record));
            
            // Check if we've reached max results
            if (max_results > 0 && results.size() >= max_results) {
                break;
            }
        }
    }
    
    return results;
}

std::vector<ioc::gaia::GaiaStar> GaiaMag18Catalog::queryConeWithMagnitude(
    double ra, double dec, double radius,
    double mag_min, double mag_max,
    size_t max_results) const {
    
    std::vector<ioc::gaia::GaiaStar> results;
    
    if (!loaded_) {
        return results;
    }
    
    results.reserve(max_results > 0 ? max_results : 1000);
    
    for (uint64_t i = 0; i < header_.total_stars; ++i) {
        auto record = readRecord(i);
        if (!record) {
            continue;
        }
        
        // Check magnitude first (faster filter)
        if (record->g_mag < mag_min || record->g_mag > mag_max) {
            continue;
        }
        
        // Check distance
        double dist = angularDistance(ra, dec, record->ra, record->dec);
        if (dist <= radius) {
            results.push_back(recordToStar(*record));
            
            if (max_results > 0 && results.size() >= max_results) {
                break;
            }
        }
    }
    
    return results;
}

std::vector<ioc::gaia::GaiaStar> GaiaMag18Catalog::queryBrightest(double ra, double dec,
                                                                   double radius,
                                                                   size_t n_brightest) const {
    // Get all stars in cone
    auto all_stars = queryCone(ra, dec, radius, 0);
    
    // Sort by magnitude
    std::partial_sort(all_stars.begin(),
                     all_stars.begin() + std::min(n_brightest, all_stars.size()),
                     all_stars.end(),
                     [](const ioc::gaia::GaiaStar& a, const ioc::gaia::GaiaStar& b) {
                         return a.phot_g_mean_mag < b.phot_g_mean_mag;
                     });
    
    // Return only N brightest
    if (all_stars.size() > n_brightest) {
        all_stars.resize(n_brightest);
    }
    
    return all_stars;
}

size_t GaiaMag18Catalog::countInCone(double ra, double dec, double radius) const {
    if (!loaded_) {
        return 0;
    }
    
    size_t count = 0;
    
    for (uint64_t i = 0; i < header_.total_stars; ++i) {
        auto record = readRecord(i);
        if (!record) {
            continue;
        }
        
        double dist = angularDistance(ra, dec, record->ra, record->dec);
        if (dist <= radius) {
            count++;
        }
    }
    
    return count;
}

} // namespace ioc_gaialib
