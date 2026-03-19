#include "astdyn/astrometry/ClosestApproachFinder.hpp"
#include "astdyn/core/Constants.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace astdyn::astrometry {

// ============================================================================
// Constants
// ============================================================================
namespace {
    constexpr double DEG_TO_RAD = constants::DEG_TO_RAD;
    constexpr double RAD_TO_DEG = 1.0 / constants::DEG_TO_RAD;
    constexpr double AU_TO_KM   = 1.495978707e8;   // 1 AU in km
}

// ============================================================================
// separation_and_derivative
// ============================================================================
// 
// Let:
//   alpha(t), delta(t)  = asteroid RA, Dec [rad] from Chebyshev
//   alpha_s, delta_s    = star RA, Dec [rad] (fixed or with proper motion)
//
// Versors:
//   p = (cos(delta)*cos(alpha), cos(delta)*sin(alpha), sin(delta))   asteroid
//   s = (cos(delta_s)*cos(alpha_s), cos(delta_s)*sin(alpha_s), sin(delta_s))  star
//
// Separation:
//   theta = atan2( |p x s|, p . s )          (numerically stable for all angles)
//
// Derivative:
//   d(theta)/dt = d/dt [ atan2(|p x s|, p.s) ]
//
// Let  c = p.s,  q = |p x s|  (= sin(theta)),  so theta = atan2(q, c).
//
//   d(theta)/dt = ( c * dq/dt - q * dc/dt ) / (q^2 + c^2)
//              = ( c * dq/dt - q * dc/dt )      [since q^2+c^2 = 1 when |p|=|s|=1]
//
// With  p' = dp/dt  (asteroid angular velocity vector, from Chebyshev vRA, vDec):
//
//   dc/dt = p'.s
//   cross = p x s
//   dq/dt = (cross . (p' x s)) / q      [if q != 0]
//         = d/dt |p x s|
//
// For q -> 0 (very close approach) we use the small-angle limit:
//   theta ≈ |p x s|   =>   d(theta)/dt ≈ (p x s).(p' x s) / |p x s|
// which is the same formula — so the expression is valid throughout.
//
// p' in Cartesian from (vRA [deg/day], vDec [deg/day]):
//   p'_x = -sin(delta)*cos(alpha)*vDec_rad - cos(delta)*sin(alpha)*vRA_rad
//   p'_y = -sin(delta)*sin(alpha)*vDec_rad + cos(delta)*cos(alpha)*vRA_rad
//   p'_z =  cos(delta)*vDec_rad
//
std::pair<double, double> ClosestApproachFinder::separation_and_derivative(
    double jd, const EvalContext& ctx)
{
    // --- Asteroid position and velocity from Chebyshev ---
    auto [pos, vel] = ctx.segment->evaluate_full(jd);
    auto [ra_deg, dec_deg, dist_au] = pos;
    auto [vra_degd, vdec_degd, vdist]   = vel;   // deg/day

    const double ra  = ra_deg  * DEG_TO_RAD;
    const double dec = dec_deg * DEG_TO_RAD;
    const double vra  = vra_degd  * DEG_TO_RAD;  // rad/day
    const double vdec = vdec_degd * DEG_TO_RAD;

    // --- Star position (with proper motion if available) ---
    double star_ra_rad, star_dec_rad;

    if (ctx.use_proper_motion && ctx.star->has_proper_motion()) {
        // Use Star::predict_at() — handles proper motion + parallax
        time::EpochTDB t = time::EpochTDB::from_jd(jd);
        auto coord = ctx.star->predict_at(t);        // SkyCoord<GCRF>
        star_ra_rad  = coord.ra().to_rad();
        star_dec_rad = coord.dec().to_rad();
    } else {
        star_ra_rad  = ctx.star->ra.to_rad();
        star_dec_rad = ctx.star->dec.to_rad();
    }

    // --- Unit versors ---
    // Asteroid versor p
    const double cos_dec = std::cos(dec);
    const double sin_dec = std::sin(dec);
    const double cos_ra  = std::cos(ra);
    const double sin_ra  = std::sin(ra);

    const double px = cos_dec * cos_ra;
    const double py = cos_dec * sin_ra;
    const double pz = sin_dec;

    // Star versor s
    const double cos_ds = std::cos(star_dec_rad);
    const double sin_ds = std::sin(star_dec_rad);
    const double cos_rs = std::cos(star_ra_rad);
    const double sin_rs = std::sin(star_ra_rad);

    const double sx = cos_ds * cos_rs;
    const double sy = cos_ds * sin_rs;
    const double sz = sin_ds;

    // --- Separation: theta = atan2(|p x s|, p.s) ---
    const double dot = px*sx + py*sy + pz*sz;

    // cross = p x s
    const double cx = py*sz - pz*sy;
    const double cy = pz*sx - px*sz;
    const double cz = px*sy - py*sx;
    const double cross_norm = std::sqrt(cx*cx + cy*cy + cz*cz);

    const double theta = std::atan2(cross_norm, dot);   // [rad], always >= 0

    // --- Derivative: d(theta)/dt ---
    // p' = dp/dt in Cartesian (vra, vdec in rad/day)
    const double dpx = -sin_dec*cos_ra*vdec - cos_dec*sin_ra*vra;
    const double dpy = -sin_dec*sin_ra*vdec + cos_dec*cos_ra*vra;
    const double dpz =  cos_dec*vdec;

    // dc/dt = p'.s
    const double dc_dt = dpx*sx + dpy*sy + dpz*sz;

    // d|cross|/dt = (cross . (p' x s)) / |cross|      if |cross| > eps
    //             = 0                                   otherwise (exact minimum)
    double dtheta_dt;
    if (cross_norm > 1e-15) {
        // p' x s
        const double pcx = dpy*sz - dpz*sy;
        const double pcy = dpz*sx - dpx*sz;
        const double pcz = dpx*sy - dpy*sx;
        const double dcross_dt = (cx*pcx + cy*pcy + cz*pcz) / cross_norm;
        // d(theta)/dt = (dot * dcross_dt - cross_norm * dc_dt) / 1  [|p|=|s|=1]
        dtheta_dt = dot * dcross_dt - cross_norm * dc_dt;
    } else {
        // Very close to zero separation — derivative is dominated by p'.s term
        // Use L'Hopital limit: dtheta/dt -> p'_perp magnitude
        // (sufficient precision for Brent bracketing near zero)
        dtheta_dt = -dc_dt;  // sign: when dot -> 1, atan2 -> 0, derivative sign flips
    }

    return {theta, dtheta_dt};
}

// ============================================================================
// brent_minimum
// ============================================================================
double ClosestApproachFinder::brent_minimum(
    double ja, double jb,
    const EvalContext& ctx,
    double tol)
{
    // Brent's method for root of f(t) = d(theta)/dt in [ja, jb].
    // Assumes f(ja) and f(jb) have opposite signs.
    constexpr int MAX_ITER = 100;

    auto f = [&](double jd) {
        return separation_and_derivative(jd, ctx).second;
    };

    double fa = f(ja);
    double fb = f(jb);

    // Sanity: if no sign change, return the lower-separation endpoint
    if (fa * fb > 0.0) {
        auto [ta, _a] = separation_and_derivative(ja, ctx);
        auto [tb, _b] = separation_and_derivative(jb, ctx);
        return (ta < tb) ? ja : jb;
    }

    double c = ja, fc = fa;
    double d = jb - ja, e = d;

    for (int i = 0; i < MAX_ITER; ++i) {
        if (fb * fc > 0.0) { c = ja; fc = fa; d = e = jb - ja; }
        if (std::abs(fc) < std::abs(fb)) {
            ja = jb; fa = fb;
            jb = c;  fb = fc;
            c  = ja; fc = fa;
        }

        double tol1 = 2.0 * 2.2e-16 * std::abs(jb) + 0.5 * tol;
        double xm   = 0.5 * (c - jb);

        if (std::abs(xm) <= tol1 || fb == 0.0) return jb;

        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
            double s = fb / fa;
            double p, q;
            if (ja == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                double r = fb / fc;
                q = fa / fc;
                s = fb / fa;
                p = s * (2.0 * xm * q * (q - r) - (jb - ja) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q; else p = -p;
            if (2.0 * p < std::min(3.0 * xm * q - std::abs(tol1 * q), std::abs(e * q))) {
                e = d; d = p / q;
            } else {
                d = xm; e = d;
            }
        } else {
            d = xm; e = d;
        }

        ja = jb; fa = fb;
        jb += (std::abs(d) > tol1) ? d : (xm > 0 ? tol1 : -tol1);
        fb = f(jb);
    }
    return jb;
}

// ============================================================================
// make_result
// ============================================================================
ClosestApproachResult ClosestApproachFinder::make_result(
    double jd,
    const EvalContext& ctx,
    std::optional<double> ast_diameter_km)
{
    auto [pos, vel] = ctx.segment->evaluate_full(jd);
    auto [ra_deg, dec_deg, dist_au] = pos;
    auto [vra_degd, vdec_degd, vdist] = vel;

    // Star position at t_ca
    double star_ra_rad, star_dec_rad;
    if (ctx.use_proper_motion && ctx.star->has_proper_motion()) {
        time::EpochTDB t = time::EpochTDB::from_jd(jd);
        auto coord = ctx.star->predict_at(t);
        star_ra_rad  = coord.ra().to_rad();
        star_dec_rad = coord.dec().to_rad();
    } else {
        star_ra_rad  = ctx.star->ra.to_rad();
        star_dec_rad = ctx.star->dec.to_rad();
    }

    // Separation (reuse separation_and_derivative)
    auto [theta, _] = separation_and_derivative(jd, ctx);

    // Position angle of star w.r.t. asteroid
    // PA = atan2( sin(dRA)*cos(dec_s), cos(dec_a)*sin(dec_s) - sin(dec_a)*cos(dec_s)*cos(dRA) )
    // measured N=0, E=90 (standard astrometric convention)
    const double dra      = star_ra_rad - ra_deg * DEG_TO_RAD;
    const double cos_ds   = std::cos(star_dec_rad);
    const double sin_ds   = std::sin(star_dec_rad);
    const double cos_da   = std::cos(dec_deg * DEG_TO_RAD);
    const double sin_da   = std::sin(dec_deg * DEG_TO_RAD);

    const double pa_rad = std::atan2(
        std::sin(dra) * cos_ds,
        cos_da * sin_ds - sin_da * cos_ds * std::cos(dra)
    );
    const double pa_deg = std::fmod(pa_rad * RAD_TO_DEG + 360.0, 360.0);

    // Relative angular rate [arcsec/min]
    // omega = sqrt( (vRA * cos(dec))^2 + vDec^2 )  [arcsec/min from deg/day]
    const double cos_dec    = std::cos(dec_deg * DEG_TO_RAD);
    const double rate_arcsec_min = std::sqrt(
        std::pow(vra_degd * cos_dec * 3600.0, 2) +
        std::pow(vdec_degd        * 3600.0, 2)
    ) / (24.0 * 60.0);

    // Apparent angular diameter of asteroid [arcsec]
    double app_diam_arcsec = 0.0;
    if (ast_diameter_km.has_value() && dist_au > 0.0) {
        const double diam_au   = ast_diameter_km.value() / AU_TO_KM;
        app_diam_arcsec = (diam_au / dist_au) * RAD_TO_DEG * 3600.0;
    }

    ClosestApproachResult r;
    r.t_ca                        = time::EpochTDB::from_jd(jd);
    r.separation                  = Angle::from_rad(theta);
    r.position_angle_deg          = pa_deg;
    r.relative_rate_arcsec_min    = rate_arcsec_min;
    r.ast_apparent_diameter_arcsec = app_diam_arcsec;
    r.is_occultation              = (app_diam_arcsec > 0.0) &&
                                    (theta * RAD_TO_DEG * 3600.0 < app_diam_arcsec / 2.0);
    return r;
}

// ============================================================================
// find_in_segment  (core implementation)
// ============================================================================
std::vector<ClosestApproachResult> ClosestApproachFinder::find_in_segment(
    const catalog::ChebyshevSegment& segment,
    const catalog::Star&             star,
    time::EpochTDB                   t1,
    time::EpochTDB                   t2,
    Angle                            max_sep,
    int                              n_coarse,
    std::optional<double>            ast_diameter_km)
{
    if (n_coarse < 3) n_coarse = 3;

    const double jd1 = t1.jd();
    const double jd2 = t2.jd();
    const double dt  = (jd2 - jd1) / static_cast<double>(n_coarse - 1);

    EvalContext ctx { &segment, &star,
                      /*use_proper_motion=*/true,
                      /*use_parallax=*/true };

    // -------------------------------------------------------------------------
    // 1. Coarse sampling — collect (jd, theta, dtheta_dt)
    // -------------------------------------------------------------------------
    struct Sample { double jd, theta, dtheta; };
    std::vector<Sample> samples;
    samples.reserve(n_coarse);

    for (int i = 0; i < n_coarse; ++i) {
        double jd = jd1 + i * dt;
        auto [theta, dtheta] = separation_and_derivative(jd, ctx);
        samples.push_back({jd, theta, dtheta});
    }

    // -------------------------------------------------------------------------
    // 2. Detect sign changes in dtheta (= local minima AND maxima of theta)
    //    Keep only intervals where theta is decreasing then increasing (minimum)
    // -------------------------------------------------------------------------
    std::vector<ClosestApproachResult> results;

    for (int i = 0; i + 1 < static_cast<int>(samples.size()); ++i) {
        const auto& s0 = samples[i];
        const auto& s1 = samples[i + 1];

        // Sign change in derivative => extremum in [s0.jd, s1.jd]
        if (s0.dtheta * s1.dtheta >= 0.0) continue;

        // Minimum: derivative goes from negative to positive
        if (s0.dtheta >= 0.0) continue;   // maximum, skip

        // Quick pre-filter: if both endpoints are above max_sep * 2, skip
        if (s0.theta > max_sep.to_rad() * 3.0 &&
            s1.theta > max_sep.to_rad() * 3.0) continue;

        // -----------------------------------------------------------------------
        // 3. Brent refinement
        // -----------------------------------------------------------------------
        double jd_min = brent_minimum(s0.jd, s1.jd, ctx);

        auto [theta_min, _] = separation_and_derivative(jd_min, ctx);

        if (theta_min > max_sep.to_rad()) continue;

        results.push_back(make_result(jd_min, ctx, ast_diameter_km));
    }

    // Sort by time (should already be sorted, but enforce)
    std::sort(results.begin(), results.end(),
        [](const ClosestApproachResult& a, const ClosestApproachResult& b) {
            return a.t_ca.jd() < b.t_ca.jd();
        });

    return results;
}

// ============================================================================
// find  (multi-segment dispatch)
// ============================================================================
std::vector<ClosestApproachResult> ClosestApproachFinder::find(
    const IChebyshevEphemeris& ephem,
    const catalog::Star&       star,
    time::EpochTDB             t1,
    time::EpochTDB             t2,
    Angle                      max_sep,
    int                        n_coarse,
    std::optional<double>      ast_diameter_km)
{
    if (t1.jd() >= t2.jd())
        throw std::invalid_argument("ClosestApproachFinder::find: t1 must be < t2");

    if (t1.jd() < ephem.start_epoch().jd() || t2.jd() > ephem.end_epoch().jd())
        throw std::out_of_range("ClosestApproachFinder::find: time window outside ephemeris coverage");

    std::vector<ClosestApproachResult> all_results;

    // Iterate over segments that overlap [t1, t2]
    // We walk by probing get_segment at the midpoint of each known segment.
    // Since segments are contiguous, we sweep t from t1 to t2.
    double jd = t1.jd();

    while (jd < t2.jd()) {
        const auto& seg = ephem.get_segment(time::EpochTDB::from_jd(jd + 1e-9));

        // Clamp to our window
        double seg_start = std::max(seg.t_start, t1.jd());
        double seg_end   = std::min(seg.t_end,   t2.jd());

        if (seg_end <= seg_start) {
            jd = seg.t_end + 1e-9;
            continue;
        }

        // Distribute n_coarse proportionally to segment fraction
        double frac = (seg_end - seg_start) / (t2.jd() - t1.jd());
        int n_seg   = std::max(3, static_cast<int>(n_coarse * frac));

        auto seg_results = find_in_segment(
            seg, star,
            time::EpochTDB::from_jd(seg_start),
            time::EpochTDB::from_jd(seg_end),
            max_sep, n_seg, ast_diameter_km);

        all_results.insert(all_results.end(),
                           seg_results.begin(), seg_results.end());

        jd = seg.t_end + 1e-9;
    }

    // Final sort and dedup (results near segment boundaries may appear twice)
    std::sort(all_results.begin(), all_results.end(),
        [](const ClosestApproachResult& a, const ClosestApproachResult& b) {
            return a.t_ca.jd() < b.t_ca.jd();
        });

    // Remove duplicates closer than 1 minute
    constexpr double ONE_MIN_JD = 1.0 / (24.0 * 60.0);
    all_results.erase(
        std::unique(all_results.begin(), all_results.end(),
            [](const ClosestApproachResult& a, const ClosestApproachResult& b) {
                return std::abs(a.t_ca.jd() - b.t_ca.jd()) < ONE_MIN_JD;
            }),
        all_results.end());

    return all_results;
}

} // namespace astdyn::astrometry
