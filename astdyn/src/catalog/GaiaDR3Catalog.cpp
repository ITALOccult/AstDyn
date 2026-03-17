/**
 * @file GaiaDR3Catalog.cpp
 * @brief Implementation of GaiaDR3Catalog — wrapper over IOC_GaiaLib's UnifiedGaiaCatalog
 */

#include "astdyn/catalog/GaiaDR3Catalog.hpp"

// IOC_GaiaLib upstream headers
#include <ioc_gaialib/unified_gaia_catalog.h>
#include <ioc_gaialib/types.h>

#include <atomic>

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
    std::atomic<bool> g_initialized{false};
    GaiaDR3Catalog* g_instance = nullptr;
}

void GaiaDR3Catalog::initialize(const std::string& json_config_path) {
    bool ok = ioc::gaia::UnifiedGaiaCatalog::initialize(json_config_path);
    if (!ok) {
        throw CatalogError("GaiaDR3Catalog: failed to initialize from config: " + json_config_path);
    }
    if (!g_instance) {
        g_instance = new GaiaDR3Catalog();
    }
    g_instance->impl_->initialized = true;
    g_initialized.store(true);
}

GaiaDR3Catalog& GaiaDR3Catalog::instance() {
    if (!g_initialized.load()) {
        throw CatalogNotInitialized{};
    }
    return *g_instance;
}

void GaiaDR3Catalog::shutdown() noexcept {
    if (g_instance) {
        ioc::gaia::UnifiedGaiaCatalog::shutdown();
        delete g_instance;
        g_instance = nullptr;
        g_initialized.store(false);
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

} // namespace astdyn::catalog
