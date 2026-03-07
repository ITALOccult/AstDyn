/**
 * @file GaiaDR3Catalog.hpp
 * @brief Singleton interface to the Gaia DR3 local stellar catalog
 *
 * Wraps IOC_GaiaLib's UnifiedGaiaCatalog, exposing its functionality under
 * the astdyn::catalog namespace with AstDyn-compatible types.
 *
 * @note Requires local Gaia DR3 catalog data (~19 GB).
 *       Build with -DASTDYN_BUILD_CATALOG=ON.
 *
 * Usage:
 * @code
 *   // Initialize once at application startup
 *   GaiaDR3Catalog::initialize("/path/to/gaia_config.json");
 *
 *   auto& cat = GaiaDR3Catalog::instance();
 *   ConeQuery q;
 *   q.ra_deg = 123.45; q.dec_deg = -10.0; q.radius_arcsec = 60.0;
 *   auto stars = cat.query_cone(q);
 * @endcode
 */

#ifndef ASTDYN_CATALOG_GAIADR3CATALOG_HPP
#define ASTDYN_CATALOG_GAIADR3CATALOG_HPP

#include "astdyn/catalog/CatalogTypes.hpp"
#include <string>
#include <vector>
#include <optional>
#include <memory>
#include <stdexcept>

// Forward-declare the upstream type to avoid polluting the public header
namespace ioc::gaia { class UnifiedGaiaCatalog; }

namespace astdyn::catalog {

// ============================================================================
// Exceptions
// ============================================================================

struct CatalogError : std::runtime_error {
    explicit CatalogError(const std::string& msg) : std::runtime_error(msg) {}
};

struct CatalogNotInitialized : CatalogError {
    CatalogNotInitialized()
        : CatalogError("GaiaDR3Catalog not initialized — call GaiaDR3Catalog::initialize() first") {}
};

// ============================================================================
// Statistics
// ============================================================================

struct CatalogStats {
    uint64_t total_queries   = 0;
    uint64_t stars_returned  = 0;
    double   avg_query_ms    = 0.0;
    double   cache_hit_rate  = 0.0;
    uint64_t cache_size_mb   = 0;
};

// ============================================================================
// GaiaDR3Catalog — main API
// ============================================================================

/**
 * @brief Thread-safe singleton interface to the Gaia DR3 local catalog.
 *
 * Backed by IOC_GaiaLib's UnifiedGaiaCatalog. Supports multiple catalog
 * backends (multifile HEALPix, SQLite, online ESA/VizieR) configured via JSON.
 */
class GaiaDR3Catalog {
public:
    // -----------------------------------------------------------------------
    // Lifecycle
    // -----------------------------------------------------------------------

    /**
     * @brief Initialize the catalog from a JSON configuration file.
     *
     * Must be called once before any query. Subsequent calls reconfigure.
     *
     * @param json_config_path Path to the IOC_GaiaLib JSON config file.
     * @throws CatalogError on failure to load or parse the catalog.
     */
    static void initialize(const std::string& json_config_path);

    /**
     * @brief Access the singleton instance.
     * @throws CatalogNotInitialized if initialize() was not called.
     */
    static GaiaDR3Catalog& instance();

    /**
     * @brief Release catalog resources. Safe to call multiple times.
     */
    static void shutdown() noexcept;

    [[nodiscard]] bool is_initialized() const noexcept;

    // -----------------------------------------------------------------------
    // Spatial queries
    // -----------------------------------------------------------------------

    /**
     * @brief Circular cone search around a sky position.
     * @param params  Center, radius, magnitude and parallax cuts.
     * @return Stars within the cone, sorted by G magnitude.
     */
    [[nodiscard]] std::vector<Star> query_cone(const ConeQuery& params) const;

    /**
     * @brief Search along a polyline corridor on the sky.
     *
     * Suitable for searching stars near the apparent path of an asteroid
     * when the corridor is pre-computed externally.
     *
     * @param params  Path, half-width, magnitude limit.
     * @return Stars within the corridor.
     */
    [[nodiscard]] std::vector<Star> query_corridor(const CorridorQuery& params) const;

    /**
     * @brief Search stars along a propagated orbit (Chebyshev representation).
     *
     * Use make_orbit_query() from CatalogIntegration.hpp to build an OrbitQuery
     * from AstDyn propagation output.
     *
     * @param params  Orbit segments + corridor width + magnitude limit.
     * @param progress Optional callback for progress reporting (async).
     * @return Stars within the orbit corridor over the time range.
     */
    [[nodiscard]] std::vector<Star> query_orbit(
        const OrbitQuery& params,
        ProgressCallback progress = nullptr) const;

    // -----------------------------------------------------------------------
    // Direct lookups
    // -----------------------------------------------------------------------

    [[nodiscard]] std::optional<Star> by_source_id(int64_t source_id) const;
    [[nodiscard]] std::optional<Star> by_name(const std::string& common_name) const;
    [[nodiscard]] std::optional<Star> by_hd(const std::string& hd_number) const;
    [[nodiscard]] std::optional<Star> by_hipparcos(const std::string& hip_number) const;
    [[nodiscard]] std::optional<Star> by_tycho2(const std::string& tycho2_id) const;

    // -----------------------------------------------------------------------
    // Diagnostics
    // -----------------------------------------------------------------------

    [[nodiscard]] CatalogStats statistics() const noexcept;
    void clear_cache() noexcept;

    // Non-copyable singleton
    GaiaDR3Catalog(const GaiaDR3Catalog&) = delete;
    GaiaDR3Catalog& operator=(const GaiaDR3Catalog&) = delete;

private:
    GaiaDR3Catalog();
    ~GaiaDR3Catalog();

    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace astdyn::catalog

#endif // ASTDYN_CATALOG_GAIADR3CATALOG_HPP
