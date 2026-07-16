/**
 * @file test_frame_types.cpp
 * @brief Frame type system: compile-time guards + numerical validation.
 *
 * The static_asserts are the real test. They do not check behaviour, they check
 * that the bugs we spent a day on CANNOT BE WRITTEN. If someone deletes Q(t)
 * from the chain again, this file stops compiling.
 */
#include <gtest/gtest.h>
#include "astdyn/math/Rotation.hpp"
#include "astdyn/coordinates/CelestialToTerrestrial.hpp"
#include "astdyn/coordinates/GeodeticTypes.hpp"
#include <type_traits>
#include <cmath>

using namespace astdyn;
using namespace astdyn::math;
using namespace astdyn::astrometry;
using astdyn::core::GCRF; using astdyn::core::CIRS;
using astdyn::core::TIRS; using astdyn::core::ITRF;

// GCC evaluates requires-expressions on non-dependent types eagerly, so the
// checks must go through a concept to be a proper SFINAE context.
template <class R, class V> concept CanApply   = requires(R r, V v) { r * v; };
template <class A, class B> concept CanCompose = requires(A a, B b) { a * b; };
template <class A, class B> concept CanAdd     = requires(A a, B b) { a + b; };
template <class A, class B> concept Assignable = requires(A a, B b) { a = b; };

// ---- bug #2: ERA applied to a GCRF vector (was: 16.6 km) ----
static_assert(!CanApply<Rotation<CIRS,TIRS>, Vec3<GCRF>>,
    "REGRESSION: applying ERA to a GCRF vector must be a type error");
static_assert(!CanApply<Rotation<CIRS,TIRS>, Direction<GCRF>>,
    "REGRESSION: same for directions");
static_assert(CanApply<Rotation<GCRF,ITRF>, Vec3<GCRF>>,
    "the full chain applied to a GCRF vector must compile");
static_assert(CanCompose<Rotation<CIRS,TIRS>, Rotation<GCRF,CIRS>>,
    "R(ERA) * Q(t) must compose");
static_assert(!CanCompose<Rotation<GCRF,CIRS>, Rotation<CIRS,TIRS>>,
    "Q * R is the reversed order and must not compile");
static_assert(std::is_same_v<
    decltype(std::declval<Rotation<TIRS,ITRF>>() * std::declval<Rotation<CIRS,TIRS>>()
           * std::declval<Rotation<GCRF,CIRS>>() * std::declval<Vec3<GCRF>>()), Vec3<ITRF>>,
    "W*R*Q on a GCRF vector must yield Vec3<ITRF>");
static_assert(!CanAdd<Vec3<GCRF>, Vec3<ITRF>>, "cross-frame addition must fail");
static_assert(CanAdd<Vec3<GCRF>, Vec3<GCRF>>, "same-frame addition must work");

// ---- bug #3: geocentric latitude used as geodetic (was: up to 21.4 km) ----
static_assert(!Assignable<coordinates::GeodeticLatitude, coordinates::GeocentricLatitude>,
    "REGRESSION: a geocentric latitude must not be usable as a geodetic one");
static_assert(!Assignable<coordinates::GeocentricLatitude, coordinates::GeodeticLatitude>,
    "and vice versa");

namespace {
struct Ref { double jd_tt; double Q[9]; };
const Ref kRefs[] = {
    {2451545.0, {9.99999999636946191e-01, 9.75665226388144902e-09, 2.69463804328460896e-05, -1.05112782367022817e-08, 9.99999999607867673e-01, 2.80047208916912528e-05, -2.69463801490472195e-05, -2.80047211647649307e-05, 9.99999999244814086e-01}},
    {2455197.5, {9.99999496686378331e-01, -1.40311849569640934e-08, -1.00330802329249061e-03, 1.50265280773886900e-09, 9.99999999922034588e-01, -1.24872248744164419e-05, 1.00330802338947792e-03, 1.24872170818024486e-05, 9.99999496608413141e-01}},
    {2458849.5, {9.99998173705995530e-01, 1.11552446733012567e-08, -1.91117363769368600e-03, 1.27594734942657051e-08, 9.99999999921711180e-01, 1.25130937306962971e-05, 1.91117363768364889e-03, -1.25130952636776094e-05, 9.99998173627706710e-01}},
    {2461236.5, {9.99996630952164844e-01, -4.96528091038783703e-09, -2.59578202473128520e-03, -7.05272207510154736e-08, 9.99999999577096843e-01, -2.90826669308392533e-05, 2.59578202377792411e-03, 2.90827520232350700e-05, 9.99996630529262576e-01}},
    {2462687.5, {9.99995508052246351e-01, -7.89248318956714545e-09, -2.99731135016042529e-03, 7.89945052877272236e-08, 9.99999999718634625e-01, 2.37218156384576585e-05, 2.99731134912986164e-03, -2.37219458524284866e-05, 9.99995507770882752e-01}},
    {2469807.5, {9.99988060822585845e-01, 2.46793752465418459e-08, -4.88653376985432545e-03, 2.36351365361098242e-07, 9.99999998573241067e-01, 5.34178071872025767e-05, 4.88653376420073761e-03, -5.34183243614538839e-05, 9.99988059395838014e-01}},
};
constexpr double R_EARTH_M = 6378137.0;
}  // namespace

TEST(FrameTypes, CipRotationMatchesErfa) {
    for (const auto& r : kRefs) {
        Eigen::Matrix3d Qref;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) Qref(i, j) = r.Q[3*i + j];
        const auto Q = coordinates::cip_rotation(time::EpochTT::from_jd(r.jd_tt));
        const Eigen::Matrix3d M = Q.matrix() * Qref.transpose();
        const double ang = std::acos(std::clamp((M.trace() - 1.0) / 2.0, -1.0, 1.0));
        EXPECT_LT(ang * R_EARTH_M, 1.0) << "JD " << r.jd_tt;
    }
}

TEST(FrameTypes, EraMatchesErfa) {
    // This test exists because it was missing. Q(t) was validated against ERFA
    // while the Earth Rotation Angle was assumed correct, and a spurious +0.5
    // day in the fractional-day term went unnoticed: exactly 180 degrees of
    // error, i.e. the shadow on the opposite side of the planet, with the
    // latitude left intact so that nothing looked obviously wrong.
    struct EraRef { double jd_ut1; double era_rad; };
    static const EraRef kEra[] = {
    {2451545.0000000, 4.89496121282375629e+00},
    {2455197.5000000, 1.75247638601582167e+00},
    {2458849.5000000, 1.74298312301005609e+00},
    {2461236.5000000, 5.10547392603933048e+00},
    {2461248.5478515, 5.61338307340346176e+00},
    {2462687.5000000, 4.93309526099962881e+00},
    {2469807.5000000, 1.74890769314321659e+00},
    };
    for (const auto& r : kEra) {
        const double era = coordinates::earth_rotation_angle(
            time::EpochUT1::from_jd(r.jd_ut1)).to_rad();
        double d = era - r.era_rad;
        while (d >  M_PI) d -= 2.0 * M_PI;
        while (d < -M_PI) d += 2.0 * M_PI;
        // 1 arcsec of ERA is ~31 m at the equator; require far better.
        EXPECT_LT(std::abs(d) * 6378137.0, 0.1)
            << "JD_UT1 " << r.jd_ut1 << ": " << d * 206264.806 << " arcsec";
    }
}

TEST(FrameTypes, GeodeticConversionIsExactOnTheEllipsoid) {
    const auto E = coordinates::Ellipsoid::wgs84();
    // A known pair: geocentric 45 deg -> geodetic 45.1924232 deg on WGS84.
    const auto gc = coordinates::GeocentricLatitude::from_angle(Angle::from_deg(45.0));
    EXPECT_NEAR(E.to_geodetic(gc).to_deg(), 45.1924232, 1e-6);
    // At the equator the two coincide.
    const auto eq = coordinates::GeocentricLatitude::from_angle(Angle::from_deg(0.0));
    EXPECT_NEAR(E.to_geodetic(eq).to_deg(), 0.0, 1e-12);
}

TEST(FrameTypes, ShadowAxisThroughGeocentreLandsOnEllipsoid) {
    const auto dir = Direction<GCRF>::from_xyz(std::cos(0.5)*std::cos(1.0),
                                               std::cos(0.5)*std::sin(1.0), std::sin(0.5));
    auto p = coordinates::shadow_point(physics::Distance::from_km(0.0),
                                       physics::Distance::from_km(0.0), dir,
                                       time::EpochTT::from_jd(2461236.5),
                                       time::EpochUT1::from_jd(2461236.5));
    ASSERT_TRUE(p.has_value());
    EXPECT_LT(std::abs(p->lat.to_deg()), 90.0);
}

TEST(FrameTypes, AxisMissingEarthReturnsNullopt) {
    const auto dir = Direction<GCRF>::from_xyz(std::cos(0.5)*std::cos(1.0),
                                               std::cos(0.5)*std::sin(1.0), std::sin(0.5));
    auto p = coordinates::shadow_point(physics::Distance::from_km(50000.0),
                                       physics::Distance::from_km(0.0), dir,
                                       time::EpochTT::from_jd(2461236.5),
                                       time::EpochUT1::from_jd(2461236.5));
    EXPECT_FALSE(p.has_value());
}
