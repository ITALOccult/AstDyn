#include <astdyn/catalog/gaia/gaia_client.h>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <cmath>
#include <algorithm>
#include <curl/curl.h>

namespace ioc {
namespace gaia {

// =============================================================================
// Rate Limiter - Simple token bucket implementation
// =============================================================================

class RateLimiter {
public:
    explicit RateLimiter(int queries_per_minute)
        : queries_per_minute_(queries_per_minute),
          last_query_time_(std::chrono::steady_clock::now()) {}
    
    void waitIfNeeded() {
        if (queries_per_minute_ <= 0) return;  // No limit
        
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            now - last_query_time_).count();
        
        int min_interval_ms = 60000 / queries_per_minute_;
        
        if (elapsed < min_interval_ms) {
            std::this_thread::sleep_for(
                std::chrono::milliseconds(min_interval_ms - elapsed));
        }
        
        last_query_time_ = std::chrono::steady_clock::now();
    }
    
    void setRate(int queries_per_minute) {
        queries_per_minute_ = queries_per_minute;
    }

private:
    int queries_per_minute_;
    std::chrono::steady_clock::time_point last_query_time_;
};

// =============================================================================
// Query Builder - Generate ADQL queries
// =============================================================================

class QueryBuilder {
public:
    explicit QueryBuilder(GaiaRelease release) : release_(release) {}
    
    std::string buildConeQuery(double ra, double dec, double radius, double max_mag) {
        std::ostringstream adql;
        adql << std::fixed << std::setprecision(6);
        
        adql << "SELECT TOP 100000 ";
        adql << "source_id, ra, dec, parallax, parallax_error, ";
        adql << "pmra, pmdec, pmra_error, pmdec_error, ";
        adql << "phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ";
        adql << "astrometric_excess_noise, astrometric_chi2_al, visibility_periods_used ";
        adql << "FROM " << getTableName() << " ";
        adql << "WHERE 1=CONTAINS(POINT('ICRS', ra, dec), ";
        adql << "CIRCLE('ICRS', " << ra << ", " << dec << ", " << radius << ")) ";
        adql << "AND phot_g_mean_mag < " << max_mag;
        
        return adql.str();
    }
    
    std::string buildBoxQuery(double ra_min, double ra_max, 
                             double dec_min, double dec_max, double max_mag) {
        std::ostringstream adql;
        adql << std::fixed << std::setprecision(6);
        
        adql << "SELECT TOP 100000 ";
        adql << "source_id, ra, dec, parallax, parallax_error, ";
        adql << "pmra, pmdec, pmra_error, pmdec_error, ";
        adql << "phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ";
        adql << "astrometric_excess_noise, astrometric_chi2_al, visibility_periods_used ";
        adql << "FROM " << getTableName() << " ";
        adql << "WHERE ra >= " << ra_min << " AND ra <= " << ra_max << " ";
        adql << "AND dec >= " << dec_min << " AND dec <= " << dec_max << " ";
        adql << "AND phot_g_mean_mag < " << max_mag;
        
        return adql.str();
    }
    
    std::string buildSourceIdQuery(const std::vector<int64_t>& source_ids) {
        std::ostringstream adql;
        
        adql << "SELECT ";
        adql << "source_id, ra, dec, parallax, parallax_error, ";
        adql << "pmra, pmdec, pmra_error, pmdec_error, ";
        adql << "phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ";
        adql << "astrometric_excess_noise, astrometric_chi2_al, visibility_periods_used ";
        adql << "FROM " << getTableName() << " ";
        adql << "WHERE source_id IN (";
        
        for (size_t i = 0; i < source_ids.size(); ++i) {
            if (i > 0) adql << ", ";
            adql << source_ids[i];
        }
        
        adql << ")";
        
        return adql.str();
    }
    
    std::string getTableName() const {
        switch (release_) {
            case GaiaRelease::DR2:  return "gaiadr2.gaia_source";
            case GaiaRelease::EDR3: return "gaiaedr3.gaia_source";
            case GaiaRelease::DR3:  return "gaiadr3.gaia_source";
            default: return "gaiadr3.gaia_source";
        }
    }

private:
    GaiaRelease release_;
};

// =============================================================================
// CURL Helper Functions
// =============================================================================

namespace {

size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

std::string curlPost(const std::string& url, const std::string& data) {
    CURL* curl = curl_easy_init();
    if (!curl) {
        throw GaiaException(ErrorCode::NETWORK_ERROR, "Failed to initialize CURL");
    }
    
    std::string response;
    struct curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "Content-Type: application/x-www-form-urlencoded");
    
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, data.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 60L);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    
    CURLcode res = curl_easy_perform(curl);
    
    if (res != CURLE_OK) {
        curl_slist_free_all(headers);
        curl_easy_cleanup(curl);
        throw GaiaException(ErrorCode::NETWORK_ERROR, 
                          std::string("CURL error: ") + curl_easy_strerror(res));
    }
    
    long http_code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
    
    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
    
    if (http_code != 200) {
        throw GaiaException(ErrorCode::NETWORK_ERROR, 
                          "HTTP error " + std::to_string(http_code));
    }
    
    return response;
}

} // anonymous namespace

// =============================================================================
// CSV Parser for VOTable/CSV responses
// =============================================================================

class CSVParser {
public:
    static std::vector<GaiaStar> parseCSV(const std::string& csv_data) {
        std::vector<GaiaStar> stars;
        std::istringstream stream(csv_data);
        std::string line;
        
        // Skip header line
        if (!std::getline(stream, line)) {
            return stars;
        }
        
        // Parse data lines
        while (std::getline(stream, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            GaiaStar star = parseLine(line);
            if (star.isValid()) {
                stars.push_back(star);
            }
        }
        
        return stars;
    }

private:
    static GaiaStar parseLine(const std::string& line) {
        GaiaStar star;
        std::istringstream ss(line);
        std::string field;
        int col = 0;
        
        while (std::getline(ss, field, ',')) {
            try {
                switch (col) {
                    case 0: star.source_id = std::stoll(field); break;
                    case 1: star.ra = std::stod(field); break;
                    case 2: star.dec = std::stod(field); break;
                    case 3: star.parallax = std::stod(field); break;
                    case 4: star.parallax_error = std::stod(field); break;
                    case 5: star.pmra = std::stod(field); break;
                    case 6: star.pmdec = std::stod(field); break;
                    case 7: star.pmra_error = std::stod(field); break;
                    case 8: star.pmdec_error = std::stod(field); break;
                    case 9: star.phot_g_mean_mag = std::stod(field); break;
                    case 10: star.phot_bp_mean_mag = std::stod(field); break;
                    case 11: star.phot_rp_mean_mag = std::stod(field); break;
                    case 12: star.astrometric_excess_noise = std::stod(field); break;
                    case 13: star.astrometric_chi2_al = std::stod(field); break;
                    case 14: star.visibility_periods_used = std::stoi(field); break;
                }
            } catch (...) {
                // Handle missing/invalid values
            }
            col++;
        }
        
        return star;
    }
};

// =============================================================================
// GaiaClient::Impl - Private implementation
// =============================================================================

class GaiaClient::Impl {
public:
    GaiaRelease release_;
    std::unique_ptr<RateLimiter> rate_limiter_;
    std::unique_ptr<QueryBuilder> query_builder_;
    int timeout_seconds_;
    int max_retries_;
    std::string tap_url_;
    
    explicit Impl(GaiaRelease release)
        : release_(release),
          rate_limiter_(std::make_unique<RateLimiter>(10)),  // 10 queries/min default
          query_builder_(std::make_unique<QueryBuilder>(release)),
          timeout_seconds_(60),
          max_retries_(3) {
        
        // Set TAP URL based on release
        tap_url_ = "https://gea.esac.esa.int/tap-server/tap/sync";
    }
    
    std::vector<GaiaStar> executeQuery(const std::string& adql, int retry_count = 0) {
        try {
            // Rate limiting
            rate_limiter_->waitIfNeeded();
            
            // Prepare POST data
            std::ostringstream post_data;
            post_data << "REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=";
            
            // URL-encode the query (simplified)
            std::string encoded_query = urlEncode(adql);
            post_data << encoded_query;
            
            // Execute HTTP POST
            std::string response = curlPost(tap_url_, post_data.str());
            
            // Parse CSV response
            return CSVParser::parseCSV(response);
            
        } catch (const GaiaException& e) {
            // Retry on network errors
            if (e.code() == ErrorCode::NETWORK_ERROR && retry_count < max_retries_) {
                std::this_thread::sleep_for(std::chrono::seconds(2 * (retry_count + 1)));
                return executeQuery(adql, retry_count + 1);
            }
            throw;
        }
    }
    
    static std::string urlEncode(const std::string& str) {
        std::ostringstream escaped;
        escaped.fill('0');
        escaped << std::hex;
        
        for (char c : str) {
            if (std::isalnum(c) || c == '-' || c == '_' || c == '.' || c == '~') {
                escaped << c;
            } else if (c == ' ') {
                escaped << '+';
            } else {
                escaped << '%' << std::setw(2) << int((unsigned char)c);
            }
        }
        
        return escaped.str();
    }
};

// =============================================================================
// GaiaClient Public Interface
// =============================================================================

GaiaClient::GaiaClient(GaiaRelease release)
    : pImpl_(std::make_unique<Impl>(release)) {
    // Initialize CURL globally
    static bool curl_initialized = false;
    if (!curl_initialized) {
        curl_global_init(CURL_GLOBAL_DEFAULT);
        curl_initialized = true;
    }
}

GaiaClient::~GaiaClient() = default;

GaiaClient::GaiaClient(GaiaClient&&) noexcept = default;
GaiaClient& GaiaClient::operator=(GaiaClient&&) noexcept = default;

std::vector<GaiaStar> GaiaClient::queryCone(
    double ra_center, 
    double dec_center, 
    double radius, 
    double max_magnitude) {
    
    std::string adql = pImpl_->query_builder_->buildConeQuery(
        ra_center, dec_center, radius, max_magnitude);
    
    return pImpl_->executeQuery(adql);
}

std::vector<GaiaStar> GaiaClient::queryBox(
    double ra_min,
    double ra_max,
    double dec_min,
    double dec_max,
    double max_magnitude) {
    
    std::string adql = pImpl_->query_builder_->buildBoxQuery(
        ra_min, ra_max, dec_min, dec_max, max_magnitude);
    
    return pImpl_->executeQuery(adql);
}

std::vector<GaiaStar> GaiaClient::queryADQL(const std::string& adql) {
    return pImpl_->executeQuery(adql);
}

std::vector<GaiaStar> GaiaClient::queryBySourceIds(
    const std::vector<int64_t>& source_ids) {
    
    if (source_ids.empty()) {
        return {};
    }
    
    std::string adql = pImpl_->query_builder_->buildSourceIdQuery(source_ids);
    return pImpl_->executeQuery(adql);
}

std::vector<GaiaStar> GaiaClient::queryAsync(
    const QueryParams& params,
    ProgressCallback callback) {
    
    if (callback) {
        callback(0, "Starting query...");
    }
    
    auto stars = queryCone(params.ra_center, params.dec_center, 
                          params.radius, params.max_magnitude);
    
    if (callback) {
        callback(100, "Query complete");
    }
    
    return stars;
}

void GaiaClient::setRateLimit(int queries_per_minute) {
    pImpl_->rate_limiter_->setRate(queries_per_minute);
}

void GaiaClient::setTimeout(int seconds) {
    pImpl_->timeout_seconds_ = seconds;
}

void GaiaClient::setMaxRetries(int max_retries) {
    pImpl_->max_retries_ = max_retries;
}

std::string GaiaClient::getTapUrl() const {
    return pImpl_->tap_url_;
}

GaiaRelease GaiaClient::getRelease() const {
    return pImpl_->release_;
}

std::string GaiaClient::getTableName() const {
    return pImpl_->query_builder_->getTableName();
}

} // namespace gaia
} // namespace ioc
