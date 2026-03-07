/**
 * @file CatalogIntegration.cpp
 * @brief Bridge between AstDyn propagation output and astdyn::catalog queries
 */

#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/core/Constants.hpp"

#include <Eigen/Dense>

// Bring astdyn math constants into scope without shadowing astdyn::catalog::constants
namespace c = astdyn::constants;
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace astdyn::catalog {

// ============================================================================
// Constants
// ============================================================================

namespace {
    // Gaia DR3 reference epoch: J2016.0 = JD 2457389.0 (TDB)
    inline constexpr double GAIA_EPOCH_JD = 2457389.0;
    // Julian year in days
    inline constexpr double JY = 365.25;
    // arcsec -> deg
    inline constexpr double ARCSEC_TO_DEG = 1.0 / 3600.0;
} // namespace

// ============================================================================
// propagate_proper_motion
// ============================================================================

SkyPoint propagate_proper_motion(const Star& star, time::EpochTDB epoch) noexcept {
    if (!star.has_proper_motion()) {
        return SkyPoint{star.ra_deg, star.dec_deg};
    }

    const double dt_yr = (epoch.jd() - GAIA_EPOCH_JD) / JY;
    const double dec_rad = star.dec_deg * c::DEG_TO_RAD;
    const double cos_dec = std::cos(dec_rad);

    // Linear model: pmRA is already pmRA*cos(dec) in Gaia convention
    double ra_new  = star.ra_deg  + (star.pmra_mas_yr  * ARCSEC_TO_DEG / 1000.0 / std::max(cos_dec, 1e-6)) * dt_yr;
    double dec_new = star.dec_deg + (star.pmdec_mas_yr * ARCSEC_TO_DEG / 1000.0) * dt_yr;

    // Wrap RA to [0, 360)
    ra_new = std::fmod(ra_new, 360.0);
    if (ra_new < 0.0) ra_new += 360.0;

    return SkyPoint{ra_new, dec_new};
}

// ============================================================================
// to_sky_coord
// ============================================================================

astrometry::SkyCoord<core::GCRF> to_sky_coord(const Star& star, time::EpochTDB epoch) {
    SkyPoint pos = propagate_proper_motion(star, epoch);
    // Build a unit vector in GCRF and wrap it in SkyCoord via from_vector()
    const double ra_rad  = pos.ra_deg  * c::DEG_TO_RAD;
    const double dec_rad = pos.dec_deg * c::DEG_TO_RAD;
    const double cos_dec = std::cos(dec_rad);

    Eigen::Vector3d unit_vec(
        cos_dec * std::cos(ra_rad),
        cos_dec * std::sin(ra_rad),
        std::sin(dec_rad)
    );

    // Scale to 1 AU (arbitrary — SkyCoord only uses direction)
    constexpr double AU_M = 1.495978707e11; // 1 AU in meters
    auto rho = math::Vector3<core::GCRF, physics::Distance>::from_si(
        unit_vec.x() * AU_M, unit_vec.y() * AU_M, unit_vec.z() * AU_M);

    return astrometry::SkyCoord<core::GCRF>::from_vector(rho);
}

// ============================================================================
// heliocentric_to_radec
// ============================================================================

SkyPoint heliocentric_to_radec(
    const physics::CartesianStateTyped<core::GCRF>& body,
    const physics::CartesianStateTyped<core::GCRF>& earth) noexcept
{
    // Geocentric direction vector (SI meters)
    const double rx = body.position.x_si() - earth.position.x_si();
    const double ry = body.position.y_si() - earth.position.y_si();
    const double rz = body.position.z_si() - earth.position.z_si();
    const double rho = std::sqrt(rx*rx + ry*ry + rz*rz);

    if (rho < 1.0) return SkyPoint{0.0, 0.0}; // degenerate

    const double ux = rx / rho;
    const double uy = ry / rho;
    const double uz = rz / rho;

    double dec_rad = std::asin(std::clamp(uz, -1.0, 1.0));
    double ra_rad  = std::atan2(uy, ux);
    if (ra_rad < 0.0) ra_rad += 2.0 * c::PI;

    return SkyPoint{ra_rad / c::DEG_TO_RAD, dec_rad / c::DEG_TO_RAD};
}

// ============================================================================
// Chebyshev fitting helpers
// ============================================================================

namespace {

// Fit Chebyshev coefficients to a set of (tau_k, f_k) samples via DCT-I
std::vector<double> cheby_fit(
    const std::vector<double>& t_samples, // normalised tau in [-1,1]
    const std::vector<double>& f_samples,
    int degree)
{
    const int N = static_cast<int>(t_samples.size());
    const int M = std::min(degree + 1, N);

    // Vandermonde-style least-squares fit
    // Build matrix T where T[i][j] = T_j(tau_i)
    Eigen::MatrixXd A(N, M);
    Eigen::VectorXd b(N);
    for (int i = 0; i < N; ++i) {
        double tau = t_samples[i];
        b(i) = f_samples[i];
        A(i, 0) = 1.0;
        if (M > 1) A(i, 1) = tau;
        for (int j = 2; j < M; ++j) {
            A(i, j) = 2.0 * tau * A(i, j-1) - A(i, j-2);
        }
    }
    // Solve via normal equations (small M, this is fine)
    Eigen::VectorXd x = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    std::vector<double> coeffs(M);
    for (int j = 0; j < M; ++j) coeffs[j] = x(j);
    return coeffs;
}

} // namespace

// ============================================================================
// make_orbit_query
// ============================================================================

OrbitQuery make_orbit_query(
    const std::vector<physics::CartesianStateTyped<core::GCRF>>& states,
    const std::vector<physics::CartesianStateTyped<core::GCRF>>& earth_states,
    time::EpochTDB t_start,
    time::EpochTDB t_end,
    double width_arcsec,
    double max_magnitude,
    double segment_days,
    int    chebyshev_degree)
{
    if (states.size() != earth_states.size() || states.empty()) {
        throw std::invalid_argument("make_orbit_query: states and earth_states must be non-empty and same size");
    }

    // Compute geocentric RA/Dec at each sample epoch
    const std::size_t N = states.size();
    std::vector<double> jd_samples(N), ra_samples(N), dec_samples(N);
    for (std::size_t i = 0; i < N; ++i) {
        jd_samples[i] = states[i].epoch.jd();
        auto sky = heliocentric_to_radec(states[i], earth_states[i]);
        ra_samples[i]  = sky.ra_deg;
        dec_samples[i] = sky.dec_deg;
    }

    const double jd0 = t_start.jd();
    const double jd1 = t_end.jd();

    // Partition into segments
    std::vector<ChebyshevSegment> segments;
    double seg_start = jd0;
    while (seg_start < jd1) {
        double seg_end = std::min(seg_start + segment_days, jd1);

        // Collect sample points within this segment
        std::vector<double> tau_pts, ra_pts, dec_pts;
        for (std::size_t i = 0; i < N; ++i) {
            if (jd_samples[i] >= seg_start && jd_samples[i] <= seg_end) {
                // Normalise to [-1, 1]
                double tau = 2.0 * (jd_samples[i] - seg_start) / (seg_end - seg_start) - 1.0;
                tau_pts.push_back(tau);
                ra_pts.push_back(ra_samples[i]);
                dec_pts.push_back(dec_samples[i]);
            }
        }

        if (tau_pts.size() >= 2) {
            ChebyshevSegment seg;
            seg.t_start  = seg_start;
            seg.t_end    = seg_end;
            seg.ra_coeffs  = cheby_fit(tau_pts, ra_pts,  chebyshev_degree);
            seg.dec_coeffs = cheby_fit(tau_pts, dec_pts, chebyshev_degree);
            segments.push_back(std::move(seg));
        }

        seg_start = seg_end;
    }

    OrbitQuery q;
    q.t_start_jd  = jd0;
    q.t_end_jd    = jd1;
    q.segments    = std::move(segments);
    q.width_arcsec  = width_arcsec;
    q.max_magnitude = max_magnitude;
    q.step_days     = segment_days / std::max(chebyshev_degree, 1);
    return q;
}

// ============================================================================
// find_occultation_candidates
// ============================================================================

std::vector<Star> find_occultation_candidates(
    const GaiaDR3Catalog& catalog,
    const std::vector<physics::CartesianStateTyped<core::GCRF>>& body_states,
    const std::vector<physics::CartesianStateTyped<core::GCRF>>& earth_states,
    time::EpochTDB t_start,
    time::EpochTDB t_end,
    double width_arcsec,
    double max_magnitude)
{
    OrbitQuery oq = make_orbit_query(
        body_states, earth_states,
        t_start, t_end,
        width_arcsec, max_magnitude);

    std::vector<Star> candidates = catalog.query_orbit(oq);

    // Apply proper motion to mid-epoch for each candidate
    time::EpochTDB mid_epoch = time::EpochTDB::from_jd(
        0.5 * (t_start.jd() + t_end.jd()));

    for (auto& star : candidates) {
        auto pos = propagate_proper_motion(star, mid_epoch);
        star.ra_deg  = pos.ra_deg;
        star.dec_deg = pos.dec_deg;
    }

    return candidates;
}

} // namespace astdyn::catalog
