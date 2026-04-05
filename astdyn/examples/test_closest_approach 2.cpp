/**
 * @file test_closest_approach.cpp
 * @brief Minimal verification of ClosestApproachFinder::find_in_segment().
 *
 * Geometry:
 *   Asteroid: linear RA motion at 0.3 deg/day, constant Dec 10°, dist 2.5 AU.
 *   Star: 0.5 arcsec to the North at the midpoint of the interval.
 *   Expected: 1 minimum at t_mid, separation ≈ 0.5 arcsec, PA ≈ 0° (North).
 */

#include "astdyn/astrometry/ClosestApproachFinder.hpp"
#include "astdyn/catalog/CatalogTypes.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/time/epoch.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace astdyn::astrometry;
using namespace astdyn::catalog;
namespace ac = astdyn::constants;
using astdyn::time::EpochTDB;

// RA(t) = RA_mid + rate*(t-t_mid):  c0=RA_mid, c1=rate*dt/2  (Chebyshev linear)
// Dec = 10° (constant), Dist = 2.5 AU (constant)
static ChebyshevSegment make_synthetic_segment() {
    constexpr double rate_deg_day  = 0.3;    // MBA-like motion
    constexpr double dt_days       = 2.0;
    constexpr double ra_mid_deg    = 100.3;  // RA at tau=0 (t_mid)
    constexpr double dec_const_deg = 10.0;
    constexpr double dist_const_au = 2.5;

    ChebyshevSegment seg;
    seg.t_start     = ac::JD2000;
    seg.t_end       = ac::JD2000 + dt_days;
    seg.ra_coeffs   = {ra_mid_deg, rate_deg_day * dt_days / 2.0};
    seg.dec_coeffs  = {dec_const_deg};
    seg.dist_coeffs = {dist_const_au};
    return seg;
}

// Star at 0.5 arcsec North of the asteroid at t_mid (same RA, Dec + ε)
static Star make_occulted_star() {
    constexpr double target_sep_arcsec = 0.5;
    Star star;
    star.ra  = RightAscension::from_deg(100.3);
    star.dec = Declination::from_deg(10.0 + target_sep_arcsec / 3600.0);
    // pm_ra_cosdec and pm_dec default to zero → has_proper_motion() = false
    return star;
}

static void verify_result(const std::vector<ClosestApproachResult>& results) {
    assert(results.size() == 1 && "FAIL: expected exactly 1 closest approach minimum");

    const auto& r            = results[0];
    const double sep_arcsec  = r.separation.to_arcsec();
    const double pa_deg      = r.position_angle_deg;

    std::cout << "t_ca           = JD " << r.t_ca.jd() << "\n";
    std::cout << "separation     = " << sep_arcsec << " arcsec\n";
    std::cout << "position_angle = " << pa_deg << " deg\n";
    std::cout << "relative_rate  = " << r.relative_rate_arcsec_min << " arcsec/min\n";

    assert(sep_arcsec < 1.0 &&
           "FAIL: minimum separation >= 1 arcsec");

    // Star is to the North: PA must be within 5° of 0° (i.e. < 5° or > 355°)
    const bool is_north_consistent = (pa_deg < 5.0 || pa_deg > 355.0);
    assert(is_north_consistent &&
           "FAIL: position angle not consistent with northern geometry");

    std::cout << "PASS: 1 minimum, sep=" << sep_arcsec
              << " arcsec < 1, PA=" << pa_deg << "° (North)\n";
}

int main() {
    const ChebyshevSegment segment = make_synthetic_segment();
    const Star             star    = make_occulted_star();

    const auto t1 = EpochTDB::from_jd(segment.t_start);
    const auto t2 = EpochTDB::from_jd(segment.t_end);

    const auto results = ClosestApproachFinder::find_in_segment(
        segment, star, t1, t2,
        Angle::from_arcsec(120.0),
        /*n_coarse=*/200);

    verify_result(results);
    return 0;
}
