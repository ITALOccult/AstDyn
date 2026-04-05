#include <astdyn/catalog/gaia/unified_gaia_catalog.h>
#include <astdyn/catalog/gaia/concurrent_multifile_catalog_v2.h>
#include <astdyn/catalog/gaia/gaia_mag18_catalog.h>
#include <astdyn/catalog/gaia/gaia_mag18_catalog_v2.h>
#include <astdyn/catalog/gaia/gaia_client.h>
#include <astdyn/catalog/gaia/common_star_names.h>
#include <astdyn/catalog/gaia/iau_star_catalog_parser.h>
#include <astdyn/catalog/gaia/gaia_sqlite_catalog.h>
#include <iostream>
#include <sstream>
#include <thread>
#include <chrono>
#include <algorithm>
#include <atomic>
#include <set>
#include <cmath>
#include <limits>
#include <filesystem>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

namespace ioc::gaia {

namespace {

std::string find_existing_path(const std::initializer_list<const char*>& candidates) {
    for (const char* candidate : candidates) {
        if (std::filesystem::exists(candidate)) {
            return candidate;
        }
    }
    return {};
}

} // namespace

// Simple JSON parser to avoid external dependencies
class SimpleJSON {
public:
    static std::map<std::string, std::string> parse(const std::string& json_str) {
        std::map<std::string, std::string> result;
        
        std::string clean_str = json_str;
        
        // Remove all whitespace except within quoted strings
        std::string cleaned;
        bool in_quotes = false;
        for (char c : clean_str) {
            if (c == '"') {
                in_quotes = !in_quotes;
                cleaned += c;
            } else if (in_quotes || (c != ' ' && c != '\t' && c != '\n' && c != '\r')) {
                cleaned += c;
            }
        }
        
        // Remove { } 
        size_t start = cleaned.find('{');
        size_t end = cleaned.find_last_of('}');
        if (start != std::string::npos && end != std::string::npos) {
            cleaned = cleaned.substr(start + 1, end - start - 1);
        }
        
        // Split by comma
        std::stringstream ss(cleaned);
        std::string item;
        
        while (std::getline(ss, item, ',')) {
            size_t colon = item.find(':');
            if (colon != std::string::npos) {
                std::string key = item.substr(0, colon);
                std::string value = item.substr(colon + 1);
                
                // Remove quotes
                key = removeQuotes(key);
                value = removeQuotes(value);
                
                result[key] = value;
            }
        }
        
        return result;
    }
    
private:
    static std::string trim(const std::string& str) {
        size_t start = str.find_first_not_of(" \t\n\r");
        if (start == std::string::npos) return "";
        size_t end = str.find_last_not_of(" \t\n\r");
        return str.substr(start, end - start + 1);
    }
    
    static std::string removeQuotes(const std::string& str) {
        std::string result = str;
        if (result.front() == '"' && result.back() == '"') {
            result = result.substr(1, result.length() - 2);
        }
        return result;
    }
};

// Private implementation class
class UnifiedGaiaCatalog::Impl {
public:
    GaiaCatalogConfig config_;
    
    // Catalog implementations
    std::unique_ptr<ConcurrentMultiFileCatalogV2> multifile_catalog_;
    std::unique_ptr<ioc_gaialib::GaiaMag18Catalog> compressed_catalog_;
    std::unique_ptr<ioc::gaia::Mag18CatalogV2> compressed_catalog_v2_;
    std::unique_ptr<GaiaSqliteCatalog> sqlite_catalog_;
    std::unique_ptr<GaiaClient> online_client_;
    
    // Star names cross-match system
    CommonStarNames star_names_;
    bool star_names_loaded_{false};
    
    // Statistics
    std::atomic<size_t> total_queries_{0};
    std::atomic<size_t> total_stars_returned_{0};
    std::atomic<double> total_query_time_{0.0};
    
    Impl() {
        // Try to load star names database on construction
        try {
            initializeStarNames();
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to initialize star names database: " << e.what() << std::endl;
        }
    }
    
    bool initializeMultiFile(const std::string& directory) {
        try {
            multifile_catalog_ = std::make_unique<ConcurrentMultiFileCatalogV2>(
                directory, config_.max_cached_chunks
            );
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Failed to initialize multi-file catalog: " << e.what() << std::endl;
            return false;
        }
    }
    
    bool initializeCompressed(const std::string& file_path) {
        try {
            // Try to open as V2 first
            auto cat_v2 = std::make_unique<ioc::gaia::Mag18CatalogV2>(file_path);
            if (cat_v2->getTotalStars() > 0) { // Simple check if loaded successfully
                 compressed_catalog_v2_ = std::move(cat_v2);
                 std::cout << "Loaded Gaia Compressed V2 catalog: " << file_path << std::endl;
                 return true;
            }
            
            // Fallback to V1
            compressed_catalog_ = std::make_unique<ioc_gaialib::GaiaMag18Catalog>(file_path);
            if (compressed_catalog_->isLoaded()) {
                std::cout << "Loaded Gaia Compressed V1 catalog (no proper motions)" << std::endl;
                return true;
            }
            
            return false;
        } catch (const std::exception& e) {
            std::cerr << "Failed to initialize compressed catalog: " << e.what() << std::endl;
            return false;
        }
    }
    
    bool initializeSqlite(const std::string& db_path) {
        try {
            sqlite_catalog_ = std::make_unique<GaiaSqliteCatalog>(db_path);
            if (sqlite_catalog_->isOpen()) {
                std::cout << "Loaded Gaia SQLite DR3 catalog: " << db_path << std::endl;
                return true;
            }
            return false;
        } catch (const std::exception& e) {
            std::cerr << "Failed to initialize SQLite catalog: " << e.what() << std::endl;
            return false;
        }
    }
    
    bool initializeOnline() {
        try {
            online_client_ = std::make_unique<GaiaClient>();
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Failed to initialize online client: " << e.what() << std::endl;
            return false;
        }
    }
    
    std::vector<GaiaStar> performQuery(const QueryParams& params) {
        auto start_time = std::chrono::high_resolution_clock::now();
        total_queries_++;
        
        std::vector<GaiaStar> results;
        
        try {
            switch (config_.catalog_type) {
                case GaiaCatalogConfig::CatalogType::MULTIFILE_V2:
                    if (multifile_catalog_) {
                        results = multifile_catalog_->queryCone(
                            params.ra_center, params.dec_center, params.radius
                        );
                    }
                    break;
                    
                case GaiaCatalogConfig::CatalogType::COMPRESSED_V2:
                    if (compressed_catalog_v2_) {
                         results = compressed_catalog_v2_->queryCone(
                            params.ra_center, params.dec_center, params.radius,
                             // Map limiting magnitude if set (primitive optimization)
                            0 
                        );
                        // Apply magnitude filter manually if needed, or V2 might support it
                        if (params.max_magnitude < 21.0) {
                             results.erase(
                                std::remove_if(results.begin(), results.end(),
                                    [&](const GaiaStar& star) {
                                        return star.phot_g_mean_mag > params.max_magnitude;
                                    }),
                                results.end()
                            );
                        }
                    } else if (compressed_catalog_) {
                        results = compressed_catalog_->queryCone(
                            params.ra_center, params.dec_center, params.radius
                        );
                    }
                    break;
                    
                case GaiaCatalogConfig::CatalogType::ONLINE_ESA:
                    if (online_client_) {
                        results = online_client_->queryCone(
                            params.ra_center, params.dec_center, params.radius,
                            params.max_magnitude
                        );
                    }
                    break;
                    
                case GaiaCatalogConfig::CatalogType::SQLITE_DR3:
                    if (sqlite_catalog_) {
                        results = sqlite_catalog_->queryCone(
                            params.ra_center, params.dec_center, params.radius,
                            params.max_magnitude
                        );
                    }
                    break;

                case GaiaCatalogConfig::CatalogType::ONLINE_VIZIER:
                    // TODO: Implement VizieR support
                    throw std::runtime_error("VizieR support not yet implemented");
                    break;
            }
            
            // Apply additional filters
            if (params.max_magnitude > 0) {
                results.erase(
                    std::remove_if(results.begin(), results.end(),
                        [&](const GaiaStar& star) {
                            return star.phot_g_mean_mag > params.max_magnitude;
                        }),
                    results.end()
                );
            }
            
            if (params.min_parallax >= 0) {
                results.erase(
                    std::remove_if(results.begin(), results.end(),
                        [&](const GaiaStar& star) {
                            return star.parallax < params.min_parallax;
                        }),
                    results.end()
                );
            }
            
        } catch (const std::exception& e) {
            std::cerr << "Query failed: " << e.what() << std::endl;
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double duration_ms = duration.count();
        
        // Atomic operations for statistics
        double old_time = total_query_time_.load();
        while (!total_query_time_.compare_exchange_weak(old_time, old_time + duration_ms)) {
            // Retry if another thread modified the value
        }
        total_stars_returned_ += results.size();
        
        return results;
    }
    
    std::optional<GaiaStar> queryBySourceId(uint64_t source_id) {
        std::optional<GaiaStar> result;
        if (config_.catalog_type == GaiaCatalogConfig::CatalogType::SQLITE_DR3 && sqlite_catalog_) {
            result = sqlite_catalog_->queryBySourceId(source_id);
        } else {
            result = performSourceIdQuery(source_id);
        }

        if (result.has_value()) {
            // Fallback to internal star names if the catalog didn't provide them
            if (star_names_loaded_ && result->common_name.empty()) {
                auto cross_match = star_names_.getCrossMatch(source_id);
                if (cross_match.has_value()) {
                    if (result->sao_designation.empty()) result->sao_designation = cross_match->sao_designation;
                    if (result->hd_designation.empty()) result->hd_designation = cross_match->hd_designation;
                    if (result->hip_designation.empty()) result->hip_designation = cross_match->hip_designation;
                    if (result->tycho2_designation.empty()) result->tycho2_designation = cross_match->tycho2_designation;
                    if (result->common_name.empty()) result->common_name = cross_match->common_name;
                }
            }
        }
        return result;
    }
    
    std::optional<GaiaStar> performSourceIdQuery(uint64_t source_id) {
        switch (config_.catalog_type) {
            case GaiaCatalogConfig::CatalogType::MULTIFILE_V2:
                if (multifile_catalog_) {
                    return multifile_catalog_->queryBySourceId(source_id);
                }
                break;
                
            case GaiaCatalogConfig::CatalogType::COMPRESSED_V2:
                if (compressed_catalog_v2_) {
                    return compressed_catalog_v2_->queryBySourceId(source_id);
                } else if (compressed_catalog_) {
                    return compressed_catalog_->queryBySourceId(source_id);
                }
                break;
                
            case GaiaCatalogConfig::CatalogType::ONLINE_ESA:
                if (online_client_) {
                    auto results = online_client_->queryBySourceIds({static_cast<int64_t>(source_id)});
                    if (!results.empty()) return results[0];
                }
                break;
                
            case GaiaCatalogConfig::CatalogType::SQLITE_DR3:
                if (sqlite_catalog_) {
                    return sqlite_catalog_->queryBySourceId(source_id);
                }
                break;

            case GaiaCatalogConfig::CatalogType::ONLINE_VIZIER:
                break;
        }
        return std::nullopt;
    }
    
    std::optional<GaiaStar> queryByCatalogDesignation(const std::string& catalog_type, const std::string& designation) {
        if (config_.catalog_type == GaiaCatalogConfig::CatalogType::SQLITE_DR3 && sqlite_catalog_) {
            auto result = sqlite_catalog_->queryByDesignation(catalog_type, designation);
            if (result) return result;
            // Fallback to internal CSN if not found in SQLite
        }

        if (!star_names_loaded_) {
            return std::nullopt;
        }
        
        std::optional<uint64_t> source_id;
        
        // Find source_id based on catalog type
        if (catalog_type == "SAO") {
            source_id = star_names_.findBySAO(designation);
        } else if (catalog_type == "HD") {
            source_id = star_names_.findByHD(designation);
        } else if (catalog_type == "HIP") {
            source_id = star_names_.findByHipparcos(designation);
        } else if (catalog_type == "TYC") {
            source_id = star_names_.findByTycho2(designation);
        } else if (catalog_type == "NAME") {
            source_id = star_names_.findByName(designation);
        }
        
        if (!source_id.has_value()) {
            return std::nullopt;
        }
        
        // Query by source_id and populate cross-match information
        return queryBySourceId(source_id.value());
    }

private:
    void initializeStarNames() {
        try {
            const std::string iau_catalog_path = find_existing_path({
                "data/IAU-CSN.json",
                "external/ioc_gaialib/data/IAU-CSN.json",
                "astdyn/external/ioc_gaialib/data/IAU-CSN.json",
                "astdyn/data/IAU-CSN.json"
            });

            // Try to load from IAU catalog first (if exists)
            if (!iau_catalog_path.empty() && IAUStarCatalogParser::loadFromJSON(iau_catalog_path, star_names_)) {
                star_names_loaded_ = true;
                std::cout << "Loaded star names from official IAU catalog" << std::endl;
                return;
            }

            const std::string common_names_path = find_existing_path({
                "data/common_star_names.csv",
                "external/ioc_gaialib/data/common_star_names.csv",
                "astdyn/external/ioc_gaialib/data/common_star_names.csv",
                "astdyn/data/common_star_names.csv"
            });
            
            // Try to load from custom CSV file
            if (!common_names_path.empty() && star_names_.loadDatabase(common_names_path)) {
                star_names_loaded_ = true;
                std::cout << "Loaded star names from custom CSV file" << std::endl;
                return;
            } 
            
            // Fall back to embedded database
            if (star_names_.loadDefaultDatabase()) {
                star_names_loaded_ = true;
                std::cout << "Loaded star names from embedded database" << std::endl;
            } else {
                std::cerr << "Warning: Could not load star names database" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error initializing star names: " << e.what() << std::endl;
        }
    }
};

// Static members
std::unique_ptr<UnifiedGaiaCatalog> UnifiedGaiaCatalog::instance_{nullptr};
bool UnifiedGaiaCatalog::initialized_{false};

bool UnifiedGaiaCatalog::initialize(const std::string& json_config) {
    try {
        auto config_map = SimpleJSON::parse(json_config);
        
        // Create instance if needed
        if (!instance_) {
            instance_ = std::unique_ptr<UnifiedGaiaCatalog>(new UnifiedGaiaCatalog());
            instance_->pimpl_ = std::make_unique<Impl>();
        }
        
        auto& impl = *instance_->pimpl_;
        
        // Parse configuration
        std::string catalog_type_str = config_map["catalog_type"];
        if (catalog_type_str == "multifile_v2") {
            impl.config_.catalog_type = GaiaCatalogConfig::CatalogType::MULTIFILE_V2;
            impl.config_.multifile_directory = config_map["multifile_directory"];
            
            if (config_map.find("max_cached_chunks") != config_map.end()) {
                impl.config_.max_cached_chunks = std::stoull(config_map["max_cached_chunks"]);
            }
            
            if (!impl.initializeMultiFile(impl.config_.multifile_directory)) {
                return false;
            }
            
        } else if (catalog_type_str == "compressed_v2") {
            impl.config_.catalog_type = GaiaCatalogConfig::CatalogType::COMPRESSED_V2;
            impl.config_.compressed_file_path = config_map["compressed_file_path"];
            
            if (!impl.initializeCompressed(impl.config_.compressed_file_path)) {
                return false;
            }
            
        } else if (catalog_type_str == "online_esa") {
            impl.config_.catalog_type = GaiaCatalogConfig::CatalogType::ONLINE_ESA;
            
            if (config_map.find("timeout_seconds") != config_map.end()) {
                impl.config_.timeout_seconds = std::stoi(config_map["timeout_seconds"]);
            }
            
            if (!impl.initializeOnline()) {
                return false;
            }
            
        } else if (catalog_type_str == "sqlite_dr3") {
            impl.config_.catalog_type = GaiaCatalogConfig::CatalogType::SQLITE_DR3;
            impl.config_.sqlite_file_path = config_map["sqlite_file_path"];
            
            if (!impl.initializeSqlite(impl.config_.sqlite_file_path)) {
                return false;
            }

        } else {
            std::cerr << "Unknown catalog type: " << catalog_type_str << std::endl;
            return false;
        }
        
        initialized_ = true;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Initialization failed: " << e.what() << std::endl;
        return false;
    }
}

UnifiedGaiaCatalog& UnifiedGaiaCatalog::getInstance() {
    if (!initialized_ || !instance_) {
        throw std::runtime_error("UnifiedGaiaCatalog not initialized. Call initialize() first.");
    }
    return *instance_;
}

void UnifiedGaiaCatalog::shutdown() {
    if (instance_) {
        instance_.reset();
        initialized_ = false;
    }
}

bool UnifiedGaiaCatalog::isInitialized() const {
    return pimpl_ != nullptr;
}

CatalogInfo UnifiedGaiaCatalog::getCatalogInfo() const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    
    CatalogInfo info;
    
    switch (pimpl_->config_.catalog_type) {
        case GaiaCatalogConfig::CatalogType::MULTIFILE_V2:
            info.catalog_name = "Multi-File Gaia DR3 V2";
            info.version = "2.0";
            info.magnitude_limit = 21.0;
            info.is_online = false;
            break;
            
        case GaiaCatalogConfig::CatalogType::COMPRESSED_V2:
            info.catalog_name = "Compressed Gaia DR3 V2";
            info.version = "2.0";
            info.magnitude_limit = 18.0;
            info.is_online = false;
            break;
            
        case GaiaCatalogConfig::CatalogType::ONLINE_ESA:
            info.catalog_name = "ESA Gaia Archive";
            info.version = "DR3";
            info.magnitude_limit = 21.0;
            info.is_online = true;
            break;
            
        case GaiaCatalogConfig::CatalogType::SQLITE_DR3:
            info.catalog_name = "Gaia DR3 SQLite (Occult Pro)";
            info.version = "DR3";
            info.magnitude_limit = 21.0;
            info.is_online = false;
            if (pimpl_->sqlite_catalog_) {
                info.total_stars = pimpl_->sqlite_catalog_->getTotalStars();
            }
            break;

        case GaiaCatalogConfig::CatalogType::ONLINE_VIZIER:
            info.catalog_name = "VizieR Gaia DR3";
            info.version = "DR3";
            info.magnitude_limit = 21.0;
            info.is_online = true;
            break;
    }
    
    info.available_fields = {
        "source_id", "ra", "dec", "phot_g_mean_mag", 
        "phot_bp_mean_mag", "phot_rp_mean_mag", "parallax"
    };
    
    return info;
}

std::vector<GaiaStar> UnifiedGaiaCatalog::queryCone(const QueryParams& params) const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    return pimpl_->performQuery(params);
}

std::future<std::vector<GaiaStar>> UnifiedGaiaCatalog::queryAsync(
    const QueryParams& params,
    ProgressCallback progress_callback
) const {
    return std::async(std::launch::async, [this, params, progress_callback]() {
        if (progress_callback) {
            progress_callback(0, "Starting query...");
        }
        
        auto result = queryCone(params);
        
        if (progress_callback) {
            progress_callback(100, "Query completed");
        }
        
        return result;
    });
}

std::vector<std::vector<GaiaStar>> UnifiedGaiaCatalog::batchQuery(
    const std::vector<QueryParams>& param_list
) const {
    std::vector<std::vector<GaiaStar>> results;
    results.reserve(param_list.size());
    
    for (const auto& params : param_list) {
        results.push_back(queryCone(params));
    }
    
    return results;
}

std::optional<GaiaStar> UnifiedGaiaCatalog::queryBySourceId(uint64_t source_id) const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    return pimpl_->queryBySourceId(source_id);
}

std::optional<GaiaStar> UnifiedGaiaCatalog::queryBySAO(const std::string& sao_number) const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    return pimpl_->queryByCatalogDesignation("SAO", sao_number);
}

std::optional<GaiaStar> UnifiedGaiaCatalog::queryByTycho2(const std::string& tycho2_designation) const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    return pimpl_->queryByCatalogDesignation("TYC", tycho2_designation);
}

std::optional<GaiaStar> UnifiedGaiaCatalog::queryByHD(const std::string& hd_number) const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    return pimpl_->queryByCatalogDesignation("HD", hd_number);
}

std::optional<GaiaStar> UnifiedGaiaCatalog::queryByHipparcos(const std::string& hip_number) const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    return pimpl_->queryByCatalogDesignation("HIP", hip_number);
}

std::optional<GaiaStar> UnifiedGaiaCatalog::queryByName(const std::string& common_name) const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    return pimpl_->queryByCatalogDesignation("NAME", common_name);
}

// =============================================================================
// Corridor Query Implementation
// =============================================================================

namespace {
    // Calculate distance from a point to a line segment on a sphere
    // Uses the cross-track distance formula
    double pointToSegmentDistance(double point_ra, double point_dec,
                                   double seg_start_ra, double seg_start_dec,
                                   double seg_end_ra, double seg_end_dec) {
        constexpr double DEG_TO_RAD = M_PI / 180.0;
        
        // Convert to radians
        double p_ra = point_ra * DEG_TO_RAD;
        double p_dec = point_dec * DEG_TO_RAD;
        double s_ra = seg_start_ra * DEG_TO_RAD;
        double s_dec = seg_start_dec * DEG_TO_RAD;
        double e_ra = seg_end_ra * DEG_TO_RAD;
        double e_dec = seg_end_dec * DEG_TO_RAD;
        
        // Angular distance from start to point
        double d_sp = std::acos(
            std::sin(s_dec) * std::sin(p_dec) +
            std::cos(s_dec) * std::cos(p_dec) * std::cos(p_ra - s_ra)
        );
        
        // Angular distance from start to end
        double d_se = std::acos(
            std::sin(s_dec) * std::sin(e_dec) +
            std::cos(s_dec) * std::cos(e_dec) * std::cos(e_ra - s_ra)
        );
        
        // If segment has zero length, return distance to start point
        if (d_se < 1e-10) {
            return d_sp / DEG_TO_RAD;
        }
        
        // Bearing from start to end
        double theta_se = std::atan2(
            std::sin(e_ra - s_ra) * std::cos(e_dec),
            std::cos(s_dec) * std::sin(e_dec) - std::sin(s_dec) * std::cos(e_dec) * std::cos(e_ra - s_ra)
        );
        
        // Bearing from start to point
        double theta_sp = std::atan2(
            std::sin(p_ra - s_ra) * std::cos(p_dec),
            std::cos(s_dec) * std::sin(p_dec) - std::sin(s_dec) * std::cos(p_dec) * std::cos(p_ra - s_ra)
        );
        
        // Cross-track distance (perpendicular distance to great circle)
        double dxt = std::asin(std::sin(d_sp) * std::sin(theta_sp - theta_se));
        
        // Along-track distance (how far along the segment)
        double dat = std::acos(std::cos(d_sp) / std::cos(dxt));
        
        // Check if the closest point is within the segment
        if (dat > d_se) {
            // Closest point is past the end - return distance to end
            double d_pe = std::acos(
                std::sin(e_dec) * std::sin(p_dec) +
                std::cos(e_dec) * std::cos(p_dec) * std::cos(p_ra - e_ra)
            );
            return d_pe / DEG_TO_RAD;
        } else if (dat < 0 || std::cos(theta_sp - theta_se) < 0) {
            // Closest point is before the start - return distance to start
            return d_sp / DEG_TO_RAD;
        }
        
        // Return cross-track distance
        return std::abs(dxt) / DEG_TO_RAD;
    }
    
    // Calculate minimum distance from a point to a polyline
    double pointToPolylineDistance(double ra, double dec, 
                                    const std::vector<CelestialPoint>& path) {
        double min_dist = std::numeric_limits<double>::max();
        
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            double dist = pointToSegmentDistance(
                ra, dec,
                path[i].ra, path[i].dec,
                path[i+1].ra, path[i+1].dec
            );
            min_dist = std::min(min_dist, dist);
        }
        
        return min_dist;
    }
}

std::vector<GaiaStar> UnifiedGaiaCatalog::queryCorridor(const CorridorQueryParams& params) const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    
    if (!params.isValid()) {
        throw std::runtime_error("Invalid corridor query parameters");
    }
    
    std::vector<GaiaStar> results;
    std::set<uint64_t> seen_source_ids;  // Avoid duplicates
    
    // Optimized strategy: Adaptive Cone Covering
    // Instead of querying every point, we generate a minimal set of query cones
    // that covers the entire corridor tube.
    
    struct SearchCone {
        double ra;
        double dec;
        double radius;
    };
    std::vector<SearchCone> search_cones;
    
    // Heuristic: Use a query radius that is at least a reasonable chunk size (e.g. 0.5 deg)
    // but scales with width. A larger radius means fewer queries but more filtering.
    // 0.5 degrees is generally efficient for file-based partial reads.
    double chunk_radius = std::max(params.width * 2.0, 0.5);
    
    // Start first cone
    if (!params.path.empty()) {
        search_cones.push_back({params.path[0].ra, params.path[0].dec, chunk_radius});
    }
    
    // Iterate through points and add cones when we step out of current coverage
    // We also check segment midpoints for long segments
    for (size_t i = 0; i < params.path.size(); ++i) {
        // Check current point
        SearchCone& current_cone = search_cones.back();
        double dist = angularDistance({params.path[i].ra, params.path[i].dec}, 
                                      {current_cone.ra, current_cone.dec});
        
        // If the tube section at this point sticks out of the cone, start a new one
        // Tube sticks out if: center_dist + width > cone_radius
        if (dist + params.width > current_cone.radius) {
            search_cones.push_back({params.path[i].ra, params.path[i].dec, chunk_radius});
        }
        
        // Handle long segments between i and i+1
        if (i + 1 < params.path.size()) {
            double seg_len = angularDistance(
                {params.path[i].ra, params.path[i].dec},
                {params.path[i+1].ra, params.path[i+1].dec}
            );
            
            // If segment is long, walk along it
            if (seg_len > chunk_radius) {
                int steps = static_cast<int>(seg_len / (chunk_radius * 0.5)); // 50% overlap approx
                for (int k = 1; k < steps; ++k) {
                    double t = static_cast<double>(k) / steps;
                    double mid_ra = params.path[i].ra + t * (params.path[i+1].ra - params.path[i].ra);
                    double mid_dec = params.path[i].dec + t * (params.path[i+1].dec - params.path[i].dec);
                    
                    SearchCone& active = search_cones.back();
                    double mid_dist = angularDistance({mid_ra, mid_dec}, {active.ra, active.dec});
                    
                    if (mid_dist + params.width > active.radius) {
                        search_cones.push_back({mid_ra, mid_dec, chunk_radius});
                    }
                }
            }
        }
    }
    
    // Execute batched queries
    // We collect ALL candidates first, then filter.
    // This allows using the batchQuery API if available, or just loop.
    std::cout << "Optimized Corridor: " << params.path.size() << " points -> " 
              << search_cones.size() << " query cones." << std::endl;
              
    for (const auto& cone : search_cones) {
        QueryParams cone_params;
        cone_params.ra_center = cone.ra;
        cone_params.dec_center = cone.dec;
        cone_params.radius = cone.radius;
        cone_params.max_magnitude = params.max_magnitude;
        cone_params.min_parallax = params.min_parallax;
        
        // Perform query
        // Note: In a real high-perf scenario, we would use batchQuery/async here
        auto cone_results = pimpl_->performQuery(cone_params);
        
        // Accumulate and Filter
        for (const auto& star : cone_results) {
            if (seen_source_ids.count(star.source_id) > 0) continue;
            
            // Precise geometric check
            // Only calculate expensive distance if star is roughly close
            // (performQuery already did a coarse check via cone radius)
            
            double dist = pointToPolylineDistance(star.ra, star.dec, params.path);
            if (dist <= params.width) {
                seen_source_ids.insert(star.source_id);
                results.push_back(star);
                
                if (params.max_results > 0 && results.size() >= params.max_results) {
                     return results;
                }
            }
        }
    }
    
    return results;
}



namespace {
    double evaluateChebyshev(double t, double t_start, double t_end, const std::vector<double>& coeffs) {
        if (coeffs.empty()) return 0.0;
        
        // Normalize time to [-1, 1]
        // Handle singularity if t_start == t_end
        if (std::abs(t_end - t_start) < 1e-9) return coeffs[0];

        double u = 2.0 * (t - t_start) / (t_end - t_start) - 1.0;
        
        // Clamping to handle potential floating point epsilon issues at boundaries
        if (u < -1.0) u = -1.0;
        if (u > 1.0) u = 1.0;
        
        double T_prev = 1.0;
        double T_curr = u;
        
        double result = coeffs[0] * T_prev;
        if (coeffs.size() > 1) {
            result += coeffs[1] * T_curr;
        }
        
        for (size_t k = 2; k < coeffs.size(); ++k) {
            double T_next = 2.0 * u * T_curr - T_prev;
            result += coeffs[k] * T_next;
            T_prev = T_curr;
            T_curr = T_next;
        }
        
        return result;
    }
}

std::vector<GaiaStar> UnifiedGaiaCatalog::queryOrbit(const OrbitQueryParams& params) const {
    if (!params.isValid()) {
        throw std::runtime_error("Invalid orbit query parameters");
    }

    // Convert orbit to corridor path
    CorridorQueryParams corridor_params;
    corridor_params.width = params.width;
    corridor_params.max_magnitude = params.max_magnitude;
    corridor_params.min_parallax = -1.0; // Default
    corridor_params.max_results = 0; // Default

    double step = params.step_size;
    if (step <= 0) {
        // Default sampling: 100 points across the interval
        step = (params.t_end - params.t_start) / 100.0;
        if (step <= 0) step = 1.0; // Should not happen given isValid check
    }

    // Sort polynomials by start time to optimize search (optional but good)
    // We make a copy to sort
    std::vector<ChebyshevPolynomial> sorted_polys = params.polynomials;
    std::sort(sorted_polys.begin(), sorted_polys.end(), 
        [](const ChebyshevPolynomial& a, const ChebyshevPolynomial& b) {
            return a.t_start < b.t_start;
        });

    for (double t = params.t_start; t <= params.t_end; t += step) {
        // Find valid polynomial
        const ChebyshevPolynomial* active_poly = nullptr;
        for (const auto& poly : sorted_polys) {
            if (t >= poly.t_start && t <= poly.t_end) {
                active_poly = &poly;
                break;
            }
            // Optimization: if we passed the time, stop (since sorted)
            if (poly.t_start > t) break; 
        }

        if (active_poly) {
            // Evaluate RA/Dec
            double ra = evaluateChebyshev(t, active_poly->t_start, active_poly->t_end, active_poly->coeffs_ra);
            double dec = evaluateChebyshev(t, active_poly->t_start, active_poly->t_end, active_poly->coeffs_dec);
            
            corridor_params.path.emplace_back(ra, dec);
        }
    }
    
    // Ensure we process the exact end time if loop missed it
    if (std::abs(params.t_end - params.t_start) > 1e-9) { // If not a point query
         double t = params.t_end;
         const ChebyshevPolynomial* active_poly = nullptr;
         for (const auto& poly : sorted_polys) {
             if (t >= poly.t_start && t <= poly.t_end) {
                 active_poly = &poly;
                 break;
             }
         }
         
         if (active_poly) {
             double ra = evaluateChebyshev(t, active_poly->t_start, active_poly->t_end, active_poly->coeffs_ra);
             double dec = evaluateChebyshev(t, active_poly->t_start, active_poly->t_end, active_poly->coeffs_dec);
             
             // Avoid duplicate if step landed exactly on end
             if (corridor_params.path.empty() || 
                 std::abs(corridor_params.path.back().ra - ra) > 1e-6 || 
                 std::abs(corridor_params.path.back().dec - dec) > 1e-6) {
                 corridor_params.path.emplace_back(ra, dec);
             }
         }
    }

    if (corridor_params.path.size() < 2) {
         // Need at least 2 points for a corridor.
         // If path has 1 point, treat as a cone search?
         // queryCorridor might handle 1 point if modified, but let's check.
         // Effectively, if we have 1 point, it's a cone search at that point.
         if (corridor_params.path.size() == 1) {
             QueryParams cone_params;
             cone_params.ra_center = corridor_params.path[0].ra;
             cone_params.dec_center = corridor_params.path[0].dec;
             cone_params.radius = params.width;
             cone_params.max_magnitude = params.max_magnitude;
             return queryCone(cone_params);
         }
         
         // If 0 points, return empty
         return {};
    }

    return queryCorridor(corridor_params);
}

std::vector<GaiaStar> UnifiedGaiaCatalog::queryCorridorJSON(const std::string& json_config) const {
    CorridorQueryParams params = CorridorQueryParams::fromJSON(json_config);
    return queryCorridor(params);
}

std::future<std::vector<GaiaStar>> UnifiedGaiaCatalog::queryCorridorAsync(
    const CorridorQueryParams& params,
    ProgressCallback progress_callback
) const {
    return std::async(std::launch::async, [this, params, progress_callback]() {
        if (progress_callback) {
            progress_callback(0, "Starting corridor query...");
        }
        
        auto result = queryCorridor(params);
        
        if (progress_callback) {
            progress_callback(100, "Corridor query completed");
        }
        
        return result;
    });
}

CatalogStats UnifiedGaiaCatalog::getStatistics() const {
    if (!pimpl_) {
        throw std::runtime_error("Catalog not initialized");
    }
    
    CatalogStats stats;
    stats.total_queries = pimpl_->total_queries_.load();
    stats.total_stars_returned = pimpl_->total_stars_returned_.load();
    
    if (stats.total_queries > 0) {
        stats.average_query_time_ms = pimpl_->total_query_time_.load() / stats.total_queries;
    }
    
    // Cache hit rate calculation would depend on specific catalog implementation
    stats.cache_hit_rate = 85.0; // Placeholder
    stats.memory_used_mb = 512;   // Placeholder
    stats.disk_cache_used_mb = 0; // Placeholder
    
    return stats;
}

void UnifiedGaiaCatalog::clearCache() {
    if (!pimpl_) {
        return;
    }
    
    // Clear caches based on catalog type
    switch (pimpl_->config_.catalog_type) {
        case GaiaCatalogConfig::CatalogType::MULTIFILE_V2:
            if (pimpl_->multifile_catalog_) {
                pimpl_->multifile_catalog_->clearCache();
            }
            break;
        case GaiaCatalogConfig::CatalogType::SQLITE_DR3:
            // SQLite handles its own caching via the OS and its internal page cache.
            // In a more complex scenario, we might call sqlite3_db_release_memory().
            break;
        // Other catalog types would implement clearCache differently
        case GaiaCatalogConfig::CatalogType::COMPRESSED_V2:
        case GaiaCatalogConfig::CatalogType::ONLINE_ESA:
        case GaiaCatalogConfig::CatalogType::ONLINE_VIZIER:
            break;
    }
}

bool UnifiedGaiaCatalog::reconfigure(const std::string& json_config) {
    shutdown();
    return initialize(json_config);
}

} // namespace ioc::gaia

#ifdef __clang__
#pragma clang diagnostic pop
#endif