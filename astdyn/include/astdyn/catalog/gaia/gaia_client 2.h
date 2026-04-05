#pragma once

#include "types.h"
#include <string>
#include <vector>
#include <memory>
#include <functional>

namespace ioc {
namespace gaia {

// Forward declarations
class RateLimiter;
class QueryBuilder;

/**
 * GaiaClient - Main interface for querying the Gaia catalog via TAP/ADQL
 * @deprecated Use UnifiedGaiaCatalog with "online_esa" configuration instead.
 * 
 * Provides high-level methods for common queries (cone, box, ADQL) with
 * automatic rate limiting, error handling, and retry logic.
 * 
 * Example usage:
 * @code
 *   GaiaClient client(GaiaRelease::DR3);
 *   auto stars = client.queryCone(67.0, 15.0, 0.5, 12.0);
 *   for (const auto& star : stars) {
 *       std::cout << "Gaia DR3 " << star.source_id << "\n";
 *   }
 * @endcode
 */
class [[deprecated("Use UnifiedGaiaCatalog instead")]] GaiaClient {
public:
    /**
     * Construct a Gaia client for the specified data release
     * @param release Which Gaia data release to query (default: DR3)
     */
    explicit GaiaClient(GaiaRelease release = GaiaRelease::DR3);
    
    /**
     * Destructor
     */
    ~GaiaClient();
    
    // No copy, allow move
    GaiaClient(const GaiaClient&) = delete;
    GaiaClient& operator=(const GaiaClient&) = delete;
    GaiaClient(GaiaClient&&) noexcept;
    GaiaClient& operator=(GaiaClient&&) noexcept;
    
    /**
     * Query stars in a circular region (cone search)
     * 
     * @param ra_center Center Right Ascension [degrees]
     * @param dec_center Center Declination [degrees]
     * @param radius Search radius [degrees]
     * @param max_magnitude Maximum G-band magnitude (default: 20.0)
     * @return Vector of GaiaStar objects
     * @throws GaiaException on network or parse errors
     */
    std::vector<GaiaStar> queryCone(
        double ra_center, 
        double dec_center, 
        double radius, 
        double max_magnitude = 20.0
    );
    
    /**
     * Query stars in a rectangular region (box search)
     * 
     * @param ra_min Minimum RA [degrees]
     * @param ra_max Maximum RA [degrees]
     * @param dec_min Minimum Dec [degrees]
     * @param dec_max Maximum Dec [degrees]
     * @param max_magnitude Maximum G-band magnitude (default: 20.0)
     * @return Vector of GaiaStar objects
     * @throws GaiaException on network or parse errors
     */
    std::vector<GaiaStar> queryBox(
        double ra_min,
        double ra_max,
        double dec_min,
        double dec_max,
        double max_magnitude = 20.0
    );
    
    /**
     * Execute a custom ADQL query
     * 
     * Allows full control over the query for advanced use cases.
     * See: https://www.cosmos.esa.int/web/gaia-users/archive/writing-queries
     * 
     * @param adql ADQL query string
     * @return Vector of GaiaStar objects
     * @throws GaiaException on network or parse errors
     * 
     * Example:
     * @code
     *   std::string adql = "SELECT TOP 100 * FROM gaiadr3.gaia_source "
     *                      "WHERE parallax > 10 AND phot_g_mean_mag < 10";
     *   auto stars = client.queryADQL(adql);
     * @endcode
     */
    std::vector<GaiaStar> queryADQL(const std::string& adql);
    
    /**
     * Query stars by Gaia source IDs
     * 
     * @param source_ids Vector of Gaia DR3 source_id values
     * @return Vector of GaiaStar objects (may be smaller if some IDs not found)
     * @throws GaiaException on network or parse errors
     */
    std::vector<GaiaStar> queryBySourceIds(const std::vector<int64_t>& source_ids);
    
    /**
     * Asynchronous cone query with progress callback
     * 
     * @param params Query parameters
     * @param callback Progress callback function (percent, message)
     * @return Vector of GaiaStar objects
     * @throws GaiaException on network or parse errors
     */
    std::vector<GaiaStar> queryAsync(
        const QueryParams& params,
        ProgressCallback callback
    );
    
    /**
     * Set rate limiting (queries per minute)
     * Default: 10 queries/minute (Gaia Archive limit)
     * 
     * @param queries_per_minute Maximum queries per minute (0 = unlimited)
     */
    void setRateLimit(int queries_per_minute);
    
    /**
     * Set network timeout for queries
     * Default: 60 seconds
     * 
     * @param seconds Timeout in seconds
     */
    void setTimeout(int seconds);
    
    /**
     * Set maximum number of retry attempts
     * Default: 3 retries
     * 
     * @param max_retries Number of retries (0 = no retries)
     */
    void setMaxRetries(int max_retries);
    
    /**
     * Get the TAP service URL being used
     * @return TAP endpoint URL
     */
    std::string getTapUrl() const;
    
    /**
     * Get the Gaia release being queried
     * @return GaiaRelease enum value
     */
    GaiaRelease getRelease() const;
    
    /**
     * Get table name for current release
     * @return Table name (e.g., "gaiadr3.gaia_source")
     */
    std::string getTableName() const;

private:
    class Impl;
    std::unique_ptr<Impl> pImpl_;
    
    // Internal query execution
    std::vector<GaiaStar> executeQuery(const std::string& adql);
    std::string httpPost(const std::string& url, const std::string& data);
};

} // namespace gaia
} // namespace ioc
