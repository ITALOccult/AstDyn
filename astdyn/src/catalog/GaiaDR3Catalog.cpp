/**
 * @file GaiaDR3Catalog.cpp
 * @brief Implementation of GaiaDR3Catalog — wrapper over IOC_GaiaLib's UnifiedGaiaCatalog
 */

#include "astdyn/catalog/GaiaDR3Catalog.hpp"

// Gaia catalog — now internalized under astdyn/catalog/gaia/
#include <astdyn/catalog/gaia/unified_gaia_catalog.h>
#include <astdyn/catalog/gaia/types.h>

// GaiaClient is deprecated in favour of UnifiedGaiaCatalog, but it is the
// direct route to a per-source_id online query for the error model; the
// unified layer does not expose one. Use is contained to this file.
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
#include <astdyn/catalog/gaia/gaia_client.h>
#ifdef __clang__
#  pragma clang diagnostic pop
#endif

#include <sqlite3.h>
#include <atomic>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <string>

namespace astdyn::catalog {

// ============================================================================
// Internal helpers: convert between ioc_gaialib types and astdyn::catalog types
// ============================================================================

namespace {

Star from_upstream(const ioc::gaia::GaiaStar& s) {
    Star out;
    out.source_id            = static_cast<int64_t>(s.source_id);
    out.common_name          = s.common_name;
    out.hd_designation       = s.hd_designation;
    out.hip_designation      = s.hip_designation;
    out.tycho2_designation   = s.tycho2_designation;
    out.sao_designation      = s.sao_designation;
    
    out.ra                   = RightAscension(Angle::from_deg(s.ra));
    out.dec                  = Declination(Angle::from_deg(s.dec));
    out.parallax             = Parallax::from_mas(s.parallax);
    out.parallax_error_mas   = s.parallax_error;
    out.pm_ra_cosdec         = ProperMotion::from_mas_yr(s.pmra);
    out.pm_dec               = ProperMotion::from_mas_yr(s.pmdec);
    out.pmra_error_mas_yr    = s.pmra_error;
    out.pmdec_error_mas_yr   = s.pmdec_error;

    // Position errors (ra_error is alpha*, Gaia convention) + full correlation
    // model for C_star. Meaningful only when astrometric_params_solved >= 5.
    out.ra_error_mas               = s.ra_error;
    out.dec_error_mas              = s.dec_error;
    out.astrometric_params_solved  = s.astrometric_params_solved;
    out.ra_dec_corr          = s.ra_dec_corr;
    out.ra_parallax_corr     = s.ra_parallax_corr;
    out.ra_pmra_corr         = s.ra_pmra_corr;
    out.ra_pmdec_corr        = s.ra_pmdec_corr;
    out.dec_parallax_corr    = s.dec_parallax_corr;
    out.dec_pmra_corr        = s.dec_pmra_corr;
    out.dec_pmdec_corr       = s.dec_pmdec_corr;
    out.parallax_pmra_corr   = s.parallax_pmra_corr;
    out.parallax_pmdec_corr  = s.parallax_pmdec_corr;
    out.pmra_pmdec_corr      = s.pmra_pmdec_corr;
    
    out.g_mag                = s.phot_g_mean_mag;
    out.bp_mag               = s.phot_bp_mean_mag;
    out.rp_mag               = s.phot_rp_mean_mag;
    out.bp_rp                = s.bp_rp;
    out.ruwe                 = s.ruwe;
    out.astrometric_excess_noise = s.astrometric_excess_noise;
    return out;
}

std::vector<Star> from_upstream(const std::vector<ioc::gaia::GaiaStar>& src) {
    std::vector<Star> out;
    out.reserve(src.size());
    for (const auto& s : src) out.push_back(from_upstream(s));
    return out;
}

ioc::gaia::QueryParams to_upstream(const ConeQuery& q) {
    ioc::gaia::QueryParams p;
    p.ra_center      = q.ra.to_deg();
    p.dec_center     = q.dec.to_deg();
    p.radius         = q.radius.to_deg();
    p.max_magnitude  = q.max_magnitude;
    p.min_parallax   = q.min_parallax.to_mas();
    return p;
}

ioc::gaia::CorridorQueryParams to_upstream(const CorridorQuery& q) {
    ioc::gaia::CorridorQueryParams p;
    for (const auto& pt : q.path) {
        p.path.push_back(ioc::gaia::CelestialPoint(pt.ra().to_deg(), pt.dec().to_deg()));
    }
    p.width         = q.width.to_deg();
    p.max_magnitude = q.max_magnitude;
    p.min_parallax  = q.min_parallax.to_mas();
    p.max_results   = q.max_results;
    return p;
}

ioc::gaia::OrbitQueryParams to_upstream(const OrbitQuery& q) {
    ioc::gaia::OrbitQueryParams p;
    p.t_start       = q.t_start.jd();
    p.t_end         = q.t_end.jd();
    p.width         = q.width.to_deg();
    p.max_magnitude = q.max_magnitude;
    p.step_size     = q.step_days;
    for (const auto& seg : q.segments) {
        ioc::gaia::ChebyshevPolynomial poly;
        poly.t_start    = seg.t_start;
        poly.t_end      = seg.t_end;
        poly.coeffs_ra  = seg.ra_coeffs;
        poly.coeffs_dec = seg.dec_coeffs;
        p.polynomials.push_back(std::move(poly));
    }
    return p;
}

} // namespace

// ============================================================================
// Star::predict_at implementation
// ============================================================================

SkyCoord<core::GCRF> Star::predict_at(time::EpochTDB target_time, const std::optional<Eigen::Vector3d>& observer_pos_ssb_m) const {
    // 1. Setup Time
    double dt_years = (target_time.jd() - ::astdyn::constants::GAIA_EPOCH_JD) / ::astdyn::constants::DAYS_PER_YEAR;
    
    // 2. Setup Distance
    double dist_m = ::astdyn::constants::PARSEC_TO_M * 1000.0; // Default: 1kpc
    if (parallax.to_mas() > ::astdyn::constants::EPSILON) {
        dist_m = ::astdyn::constants::PARSEC_TO_M / (parallax.to_mas() / 1000.0);
    }

    // 3. Compute 3D Position at Gaia Epoch
    double r_a = ra.to_rad();
    double r_d = dec.to_rad();
    Eigen::Vector3d pos_epoch(
        dist_m * std::cos(r_a) * std::cos(r_d),
        dist_m * std::sin(r_a) * std::cos(r_d),
        dist_m * std::sin(r_d)
    );

    // 4. Proper Motion Velocity (m/yr)
    // Basis vectors in spherical coordinates
    Eigen::Vector3d e_ra(-std::sin(r_a), std::cos(r_a), 0.0);
    Eigen::Vector3d e_dec(-std::sin(r_d) * std::cos(r_a), -std::sin(r_d) * std::sin(r_a), std::cos(r_d));
    
    double pm_ra_rad_yr  = pm_ra_cosdec.to_rad_yr();
    double pm_dec_rad_yr = pm_dec.to_rad_yr();
    Eigen::Vector3d vel_m_yr = dist_m * (pm_ra_rad_yr * e_ra + pm_dec_rad_yr * e_dec);
    
    // 5. Propagate and Apply Annual Parallax
    Eigen::Vector3d pos_star_ssb = pos_epoch + vel_m_yr * dt_years;
    Eigen::Vector3d rho_vec = observer_pos_ssb_m ? (pos_star_ssb - *observer_pos_ssb_m) : pos_star_ssb;
    
    return SkyCoord<core::GCRF>::from_vector(math::Vector3<core::GCRF, Distance>::from_si(rho_vec.x(), rho_vec.y(), rho_vec.z()));
}

// ============================================================================
// GaiaDR3Catalog::Impl — PIMPL hiding UnifiedGaiaCatalog
// ============================================================================

struct GaiaDR3Catalog::Impl {
    // UnifiedGaiaCatalog is a singleton itself; we just hold a reference to it
    // after initialization.
    bool initialized = false;
};

// ============================================================================
// Static singleton machinery
// ============================================================================

namespace {
    std::mutex g_mutex;
    std::atomic<bool> g_initialized{false};
    std::unique_ptr<GaiaDR3Catalog> g_instance;
}

void GaiaDR3Catalog::initialize(const std::string& json_config_path) {
    std::lock_guard<std::mutex> lock(g_mutex);
    bool ok = ioc::gaia::UnifiedGaiaCatalog::initialize(json_config_path);
    if (!ok) {
        throw CatalogError("GaiaDR3Catalog: failed to initialize from config: " + json_config_path);
    }
    if (!g_instance) {
        g_instance.reset(new GaiaDR3Catalog());
    }
    g_instance->impl_->initialized = true;
    g_initialized.store(true, std::memory_order_release);
}

GaiaDR3Catalog& GaiaDR3Catalog::instance() {
    if (!g_initialized.load(std::memory_order_acquire)) {
        throw CatalogNotInitialized{};
    }
    std::lock_guard<std::mutex> lock(g_mutex);
    if (!g_instance) {
        throw CatalogNotInitialized{};
    }
    return *g_instance;
}

void GaiaDR3Catalog::shutdown() noexcept {
    std::lock_guard<std::mutex> lock(g_mutex);
    if (g_instance) {
        // Mark as uninitialized before destruction to prevent concurrent access.
        g_initialized.store(false, std::memory_order_release);
        ioc::gaia::UnifiedGaiaCatalog::shutdown();
        g_instance.reset();
    }
}

// ============================================================================
// Constructor / Destructor
// ============================================================================

GaiaDR3Catalog::GaiaDR3Catalog() : impl_(std::make_unique<Impl>()) {}
GaiaDR3Catalog::~GaiaDR3Catalog() = default;

bool GaiaDR3Catalog::is_initialized() const noexcept {
    return impl_->initialized;
}

// ============================================================================
// Spatial queries
// ============================================================================

std::vector<Star> GaiaDR3Catalog::query_cone(const ConeQuery& params) const {
    auto& upstream = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    return from_upstream(upstream.queryCone(to_upstream(params)));
}

std::vector<Star> GaiaDR3Catalog::query_corridor(const CorridorQuery& params) const {
    auto& upstream = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    return from_upstream(upstream.queryCorridor(to_upstream(params)));
}

std::vector<Star> GaiaDR3Catalog::query_orbit(
    const OrbitQuery& params,
    ProgressCallback /*progress*/) const
{
    auto& upstream = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    return from_upstream(upstream.queryOrbit(to_upstream(params)));
}

// ============================================================================
// Direct lookups
// ============================================================================

std::optional<Star> GaiaDR3Catalog::by_source_id(int64_t source_id) const {
    auto& upstream = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    auto result = upstream.queryBySourceId(static_cast<uint64_t>(source_id));
    if (!result) return std::nullopt;
    return from_upstream(*result);
}

std::optional<Star> GaiaDR3Catalog::by_name(const std::string& name) const {
    auto& upstream = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    auto result = upstream.queryByName(name);
    if (!result) return std::nullopt;
    return from_upstream(*result);
}

std::optional<Star> GaiaDR3Catalog::by_hd(const std::string& hd) const {
    auto& upstream = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    auto result = upstream.queryByHD(hd);
    if (!result) return std::nullopt;
    return from_upstream(*result);
}

std::optional<Star> GaiaDR3Catalog::by_hipparcos(const std::string& hip) const {
    auto& upstream = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    auto result = upstream.queryByHipparcos(hip);
    if (!result) return std::nullopt;
    return from_upstream(*result);
}

std::optional<Star> GaiaDR3Catalog::by_tycho2(const std::string& tyc) const {
    auto& upstream = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    auto result = upstream.queryByTycho2(tyc);
    if (!result) return std::nullopt;
    return from_upstream(*result);
}

// ============================================================================
// Diagnostics
// ============================================================================

CatalogStats GaiaDR3Catalog::statistics() const noexcept {
    // IOC_GaiaLib exposes stats via getStatistics() returning its own struct.
    // We do a best-effort mapping; unmapped fields remain zero-initialised.
    auto& upstream = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    auto s = upstream.getStatistics();
    CatalogStats out;
    out.total_queries  = s.total_queries;
    out.stars_returned = s.total_stars_returned;
    out.avg_query_ms   = s.average_query_time_ms;
    out.cache_hit_rate = s.cache_hit_rate;
    return out;
}

void GaiaDR3Catalog::clear_cache() noexcept {
    ioc::gaia::UnifiedGaiaCatalog::getInstance().clearCache();
}


// ============================================================================
// Astrometric errors: local SQLite cache + online per-source_id fetch
// ============================================================================

namespace {

// Cache path: alongside the catalogue config, in ~/.ioccultcalc/.
std::string error_cache_path() {
    const char* home = std::getenv("HOME");
    if (!home) return "star_errors_cache.db";
    return std::string(home) + "/.ioccultcalc/star_errors_cache.db";
}

// Open (creating if needed) the cache DB and ensure the schema exists.
// Returns nullptr on failure; callers then degrade to online-only.
sqlite3* open_error_cache() {
    sqlite3* db = nullptr;
    if (sqlite3_open_v2(error_cache_path().c_str(), &db,
            SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, nullptr) != SQLITE_OK) {
        if (db) sqlite3_close(db);
        return nullptr;
    }
    const char* schema =
        "CREATE TABLE IF NOT EXISTS star_errors ("
        " sid INTEGER PRIMARY KEY,"
        " ra_err REAL, dec_err REAL, plx_err REAL, pmra_err REAL, pmdec_err REAL,"
        " params_solved INTEGER,"
        " c_ra_dec REAL, c_ra_plx REAL, c_ra_pmra REAL, c_ra_pmdec REAL,"
        " c_dec_plx REAL, c_dec_pmra REAL, c_dec_pmdec REAL,"
        " c_plx_pmra REAL, c_plx_pmdec REAL, c_pmra_pmdec REAL);";
    sqlite3_exec(db, schema, nullptr, nullptr, nullptr);
    return db;
}

// Read one row from the cache. Returns nullopt on miss or error.
std::optional<StarAstrometricErrors> cache_read(sqlite3* db, int64_t sid) {
    if (!db) return std::nullopt;
    const char* sql = "SELECT ra_err,dec_err,plx_err,pmra_err,pmdec_err,params_solved,"
        "c_ra_dec,c_ra_plx,c_ra_pmra,c_ra_pmdec,c_dec_plx,c_dec_pmra,c_dec_pmdec,"
        "c_plx_pmra,c_plx_pmdec,c_pmra_pmdec FROM star_errors WHERE sid=?;";
    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &st, nullptr) != SQLITE_OK) return std::nullopt;
    sqlite3_bind_int64(st, 1, sid);
    std::optional<StarAstrometricErrors> out;
    if (sqlite3_step(st) == SQLITE_ROW) {
        StarAstrometricErrors e;
        e.ra_error_mas       = sqlite3_column_double(st, 0);
        e.dec_error_mas      = sqlite3_column_double(st, 1);
        e.parallax_error_mas = sqlite3_column_double(st, 2);
        e.pmra_error_mas_yr  = sqlite3_column_double(st, 3);
        e.pmdec_error_mas_yr = sqlite3_column_double(st, 4);
        e.astrometric_params_solved = sqlite3_column_int(st, 5);
        e.ra_dec_corr        = sqlite3_column_double(st, 6);
        e.ra_parallax_corr   = sqlite3_column_double(st, 7);
        e.ra_pmra_corr       = sqlite3_column_double(st, 8);
        e.ra_pmdec_corr      = sqlite3_column_double(st, 9);
        e.dec_parallax_corr  = sqlite3_column_double(st, 10);
        e.dec_pmra_corr      = sqlite3_column_double(st, 11);
        e.dec_pmdec_corr     = sqlite3_column_double(st, 12);
        e.parallax_pmra_corr = sqlite3_column_double(st, 13);
        e.parallax_pmdec_corr= sqlite3_column_double(st, 14);
        e.pmra_pmdec_corr    = sqlite3_column_double(st, 15);
        out = e;
    }
    sqlite3_finalize(st);
    return out;
}

void cache_write(sqlite3* db, int64_t sid, const StarAstrometricErrors& e) {
    if (!db) return;
    const char* sql = "INSERT OR REPLACE INTO star_errors VALUES"
        "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &st, nullptr) != SQLITE_OK) return;
    sqlite3_bind_int64(st, 1, sid);
    sqlite3_bind_double(st, 2, e.ra_error_mas);
    sqlite3_bind_double(st, 3, e.dec_error_mas);
    sqlite3_bind_double(st, 4, e.parallax_error_mas);
    sqlite3_bind_double(st, 5, e.pmra_error_mas_yr);
    sqlite3_bind_double(st, 6, e.pmdec_error_mas_yr);
    sqlite3_bind_int   (st, 7, e.astrometric_params_solved);
    sqlite3_bind_double(st, 8,  e.ra_dec_corr);
    sqlite3_bind_double(st, 9,  e.ra_parallax_corr);
    sqlite3_bind_double(st, 10, e.ra_pmra_corr);
    sqlite3_bind_double(st, 11, e.ra_pmdec_corr);
    sqlite3_bind_double(st, 12, e.dec_parallax_corr);
    sqlite3_bind_double(st, 13, e.dec_pmra_corr);
    sqlite3_bind_double(st, 14, e.dec_pmdec_corr);
    sqlite3_bind_double(st, 15, e.parallax_pmra_corr);
    sqlite3_bind_double(st, 16, e.parallax_pmdec_corr);
    sqlite3_bind_double(st, 17, e.pmra_pmdec_corr);
    sqlite3_step(st);
    sqlite3_finalize(st);
}

} // namespace

std::optional<StarAstrometricErrors>
GaiaDR3Catalog::fetch_astrometric_errors(int64_t source_id) const {
    // 1. Local cache first (offline path).
    sqlite3* cache = open_error_cache();
    if (auto hit = cache_read(cache, source_id)) {
        if (cache) sqlite3_close(cache);
        return hit;
    }

    // 2. Miss: one online query by source_id for the full error model.
    StarAstrometricErrors e;
    try {
        ioc::gaia::GaiaClient client(ioc::gaia::GaiaRelease::DR3);
        auto rows = client.queryBySourceIds({static_cast<int64_t>(source_id)});
        if (rows.empty()) { if (cache) sqlite3_close(cache); return std::nullopt; }
        const auto& s = rows.front();
        e.ra_error_mas       = s.ra_error;
        e.dec_error_mas      = s.dec_error;
        e.parallax_error_mas = s.parallax_error;
        e.pmra_error_mas_yr  = s.pmra_error;
        e.pmdec_error_mas_yr = s.pmdec_error;
        e.astrometric_params_solved = s.astrometric_params_solved;
        e.ra_dec_corr        = s.ra_dec_corr;
        e.ra_parallax_corr   = s.ra_parallax_corr;
        e.ra_pmra_corr       = s.ra_pmra_corr;
        e.ra_pmdec_corr      = s.ra_pmdec_corr;
        e.dec_parallax_corr  = s.dec_parallax_corr;
        e.dec_pmra_corr      = s.dec_pmra_corr;
        e.dec_pmdec_corr     = s.dec_pmdec_corr;
        e.parallax_pmra_corr = s.parallax_pmra_corr;
        e.parallax_pmdec_corr= s.parallax_pmdec_corr;
        e.pmra_pmdec_corr    = s.pmra_pmdec_corr;
    } catch (const std::exception&) {
        if (cache) sqlite3_close(cache);
        return std::nullopt;  // offline and not cached: honest failure
    }

    // 3. Cache for next time, then return.
    cache_write(cache, source_id, e);
    if (cache) sqlite3_close(cache);
    return e;
}

} // namespace astdyn::catalog
