#include "astdyn/catalog/CatalogIntegration.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/astrometry/Astrometry.hpp"

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

// Local constants removed in favor of astdyn::constants

// ============================================================================
// to_sky_coord
// ============================================================================

astrometry::SkyCoord<core::GCRF> to_sky_coord(const Star& star, time::EpochTDB epoch) {
    return star.predict_at(epoch);
}

SkyCoord<core::GCRF> heliocentric_to_skycoord(
    const physics::CartesianStateTyped<core::GCRF>& body,
    const physics::CartesianStateTyped<core::GCRF>& earth) noexcept
{
    // Geocentric direction vector (SI meters)
    Eigen::Vector3d rho_vec = body.position.to_eigen_si() - earth.position.to_eigen_si();
    Distance rho = Distance::from_si(rho_vec.norm());

    if (rho.to_m() < 1.0) {
        return SkyCoord<core::GCRF>::from_vector(math::Vector3<core::GCRF, Distance>());
    }

    return SkyCoord<core::GCRF>::from_vector(math::Vector3<core::GCRF, Distance>::from_si(rho_vec.x(), rho_vec.y(), rho_vec.z()));
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
    Angle          width,
    double         max_magnitude,
    double         segment_days,
    int            chebyshev_degree)
{
    if (states.size() != earth_states.size() || states.empty()) {
        throw std::invalid_argument("make_orbit_query: states and earth_states must be non-empty and same size");
    }

    // Compute geocentric SkyCoord at each sample epoch
    const std::size_t N = states.size();
    std::vector<double> jd_samples(N), ra_samples(N), dec_samples(N);
    for (std::size_t i = 0; i < N; ++i) {
        jd_samples[i] = states[i].epoch.jd();
        auto sky = heliocentric_to_skycoord(states[i], earth_states[i]);
        ra_samples[i]  = sky.ra().to_deg();
        dec_samples[i] = sky.dec().to_deg();
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
    q.t_start     = t_start;
    q.t_end       = t_end;
    q.segments    = std::move(segments);
    q.width       = width;
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
    Angle          width,
    double         max_magnitude)
{
    OrbitQuery oq = make_orbit_query(
        body_states, earth_states,
        t_start, t_end,
        width, max_magnitude);

    std::vector<Star> candidates = catalog.query_orbit(oq);

    // Apply proper motion correction to the center of the window
    time::EpochTDB mid_epoch = time::EpochTDB::from_jd(
        0.5 * (t_start.jd() + t_end.jd()));

    for (auto& star : candidates) {
        auto sky = star.predict_at(mid_epoch);
        star.ra  = sky.ra();
        star.dec = sky.dec();
    }

    return candidates;
}

ChebyshevSegment fit_chebyshev(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_elements,
    time::EpochTDB center_epoch,
    double         duration_days,
    const AstDynConfig& cfg,
    int            degree)
{
    const double t_start = center_epoch.jd() - duration_days / 2.0;
    const double t_end   = center_epoch.jd() + duration_days / 2.0;

    // Sample the orbit (using +2 points to ensure enough data for degree D)
    const int N = degree + 5;
    std::vector<double> jd_samples(N), ra_samples(N), dec_samples(N), dist_samples(N), tau_samples(N);

    astrometry::AstrometricSettings a_settings;
    a_settings.light_time_correction = true;
    a_settings.aberrazione_differenziale = false;
    a_settings.frame_conversion_to_equatorial = true;

    for (int i = 0; i < N; ++i) {
        double tau = -1.0 + 2.0 * i / (N - 1); // Normalized time [-1, 1]
        double jd  = t_start + (tau + 1.0) * (t_end - t_start) / 2.0;
        time::EpochTDB t_obs = time::EpochTDB::from_jd(jd);

        auto obs = astrometry::AstrometryReducer::compute_observation(
            initial_elements, initial_elements.epoch, t_obs, cfg, a_settings
        );

        if (obs) {
            jd_samples[i]  = jd;
            tau_samples[i] = tau;
            ra_samples[i]  = obs->ra.to_deg();
            dec_samples[i] = obs->dec.to_deg();
            dist_samples[i] = obs->distance.to_au();
        } else {
            // Degenerate or failure, should not happen in nominal conditions
            throw std::runtime_error("fit_chebyshev: propagation failed at JD " + std::to_string(jd));
        }
    }

    ChebyshevSegment seg;
    seg.t_start = t_start;
    seg.t_end   = t_end;
    seg.ra_coeffs  = cheby_fit(tau_samples, ra_samples, degree);
    seg.dec_coeffs = cheby_fit(tau_samples, dec_samples, degree);
    seg.dist_coeffs = cheby_fit(tau_samples, dist_samples, degree);
    return seg;
}

ChebyshevSegment fit_chebyshev_spk(
    astdyn::io::SPKReader& reader,
    int            target_id,
    time::EpochTDB center_epoch,
    double         duration_days,
    int            degree,
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem)
{
    const double t_start = center_epoch.jd() - duration_days / 2.0;
    const double t_end   = center_epoch.jd() + duration_days / 2.0;

    // Sample the orbit
    const int N = degree + 5;
    std::vector<double> jd_samples(N), ra_samples(N), dec_samples(N), dist_samples(N), tau_samples(N);

    // Settings for the reduction
    astrometry::AstrometricSettings a_settings;
    a_settings.light_time_correction = true;
    a_settings.aberrazione_differenziale = false; 
    a_settings.frame_conversion_to_equatorial = true;

    for (int i = 0; i < N; ++i) {
        double tau = -1.0 + 2.0 * i / (N - 1);
        double jd  = t_start + (tau + 1.0) * (t_end - t_start) / 2.0;
        time::EpochTDB t_obs = time::EpochTDB::from_jd(jd);
        double et_obs = (jd - c::JD2000) * c::SECONDS_PER_DAY;

        // 1. Get Earth position at observation time
        auto actual_ephem = ephem ? ephem : std::make_shared<ephemeris::PlanetaryEphemeris>();
        auto earth_state = actual_ephem->getState(ephemeris::CelestialBody::EARTH, t_obs);
        Eigen::Vector3d r_earth = earth_state.position.to_eigen_si() / 1000.0; // km

        // 2. Iterative light-time correction
        double lt = 0.0;
        Eigen::Vector3d rho_vec;
        for (int iter = 0; iter < 5; ++iter) {
            double et_emit = et_obs - lt;
            auto s_emit = reader.getState(target_id, et_emit);
            rho_vec = Eigen::Vector3d(s_emit[0], s_emit[1], s_emit[2]) - r_earth;
            lt = rho_vec.norm() / (c::C_LIGHT / 1000.0);
        }

        // 3. Convert to Equatorial J2000 if needed (AstrometryReducer helper)
        // Note: SPK data is already in J2000 Ecliptic (native for AstDyn) or ICRF?
        // SPKReader.hpp says Ecliptic J2000 frame.
        // We need to rotate to Equatorial for catalog matching.
        
        // Use the rotation from AstrometryReducer if private access is problematic, 
        // we can do it manually here.
        double obliquity = c::OBLIQUITY_J2000;
        Eigen::Matrix3d rot;
        rot << 1, 0, 0,
               0, std::cos(-obliquity), -std::sin(-obliquity),
               0, std::sin(-obliquity),  std::cos(-obliquity);
        
        Eigen::Vector3d rho_eq = rot * rho_vec;
        
        // 4. Final conversion to RA/Dec
        double r = rho_eq.norm();
        double ra = std::atan2(rho_eq.y(), rho_eq.x()) * c::RAD_TO_DEG;
        double dec = std::asin(rho_eq.z() / r) * c::RAD_TO_DEG;
        if (ra < 0) ra += 2.0 * c::PI * c::RAD_TO_DEG;

        jd_samples[i]  = jd;
        tau_samples[i] = tau;
        ra_samples[i]  = ra;
        dec_samples[i] = dec;
        dist_samples[i] = r / c::AU;
    }

    ChebyshevSegment seg;
    seg.t_start = t_start;
    seg.t_end   = t_end;
    seg.ra_coeffs  = cheby_fit(tau_samples, ra_samples, degree);
    seg.dec_coeffs = cheby_fit(tau_samples, dec_samples, degree);
    seg.dist_coeffs = cheby_fit(tau_samples, dist_samples, degree);
    return seg;
}

std::vector<Star> find_stars_near_segment(
    const GaiaDR3Catalog& catalog,
    const ChebyshevSegment& segment,
    Angle          width,
    double         max_magnitude)
{
    OrbitQuery oq;
    oq.t_start     = time::EpochTDB::from_jd(segment.t_start);
    oq.t_end       = time::EpochTDB::from_jd(segment.t_end);
    oq.segments    = { segment };
    oq.width       = width;
    oq.max_magnitude = max_magnitude;
    oq.step_days     = (segment.t_end - segment.t_start) / std::max((int)segment.ra_coeffs.size(), 5);
    
    return catalog.query_orbit(oq);
}

} // namespace astdyn::catalog
