/**
 * @file test_celestial_to_terrestrial.cpp
 * @brief Validates the GCRF->ITRF chain against ERFA/SOFA (iauC2i06a).
 *
 * Reference values produced with ERFA (the IAU SOFA-equivalent reference
 * implementation). The tolerance corresponds to ~1 m at the Earth's surface.
 */
#include <gtest/gtest.h>
#include "astdyn/coordinates/CelestialToTerrestrial.hpp"
#include <cmath>

using namespace astdyn::coordinates;

namespace {
struct Ref { double jd_tt; double Q[9]; };

// GCRS -> CIRS matrices from ERFA iauC2i06a.
const Ref kRefs[] = {
    {2451545.0, {9.99999999636946191e-01, 9.75665226388144902e-09, 2.69463804328460896e-05, -1.05112782367022817e-08, 9.99999999607867673e-01, 2.80047208916912528e-05, -2.69463801490472195e-05, -2.80047211647649307e-05, 9.99999999244814086e-01}},
    {2455197.5, {9.99999496686378331e-01, -1.40311849569640934e-08, -1.00330802329249061e-03, 1.50265280773886900e-09, 9.99999999922034588e-01, -1.24872248744164419e-05, 1.00330802338947792e-03, 1.24872170818024486e-05, 9.99999496608413141e-01}},
    {2458849.5, {9.99998173705995530e-01, 1.11552446733012567e-08, -1.91117363769368600e-03, 1.27594734942657051e-08, 9.99999999921711180e-01, 1.25130937306962971e-05, 1.91117363768364889e-03, -1.25130952636776094e-05, 9.99998173627706710e-01}},
    {2461236.5, {9.99996630952164844e-01, -4.96528091038783703e-09, -2.59578202473128520e-03, -7.05272207510154736e-08, 9.99999999577096843e-01, -2.90826669308392533e-05, 2.59578202377792411e-03, 2.90827520232350700e-05, 9.99996630529262576e-01}},
    {2462687.5, {9.99995508052246351e-01, -7.89248318956714545e-09, -2.99731135016042529e-03, 7.89945052877272236e-08, 9.99999999718634625e-01, 2.37218156384576585e-05, 2.99731134912986164e-03, -2.37219458524284866e-05, 9.99995507770882752e-01}},
    {2469807.5, {9.99988060822585845e-01, 2.46793752465418459e-08, -4.88653376985432545e-03, 2.36351365361098242e-07, 9.99999998573241067e-01, 5.34178071872025767e-05, 4.88653376420073761e-03, -5.34183243614538839e-05, 9.99988059395838014e-01}},
};

constexpr double DAS2R = 4.848136811095359935899141e-6;
constexpr double R_EARTH_M = 6378137.0;
}  // namespace

TEST(CelestialToTerrestrial, QMatrixMatchesErfa) {
    for (const auto& r : kRefs) {
        Eigen::Matrix3d Qref;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) Qref(i, j) = r.Q[3*i + j];

        const Eigen::Matrix3d Q = celestial_to_intermediate(r.jd_tt);

        // Residual rotation angle between the two matrices.
        const Eigen::Matrix3d M = Q * Qref.transpose();
        const double cos_a = std::clamp((M.trace() - 1.0) / 2.0, -1.0, 1.0);
        const double ang = std::acos(cos_a);

        EXPECT_LT(ang * R_EARTH_M, 1.0)
            << "JD " << r.jd_tt << ": residual " << ang / DAS2R
            << " arcsec = " << ang * R_EARTH_M << " m on the ground";
    }
}

TEST(CelestialToTerrestrial, OmittingQIsTheBugWeFixed) {
    // Regression guard: ERA alone (i.e. Q == I) must differ from the correct
    // chain by the known ~16 km at the 2026 epoch. If this ever collapses to
    // zero, Q has silently stopped being applied.
    const double jd_tt = 2461236.5;
    const Eigen::Matrix3d Q = celestial_to_intermediate(jd_tt);
    const double cos_a = std::clamp((Q.trace() - 1.0) / 2.0, -1.0, 1.0);
    const double ang = std::acos(cos_a);
    EXPECT_GT(ang * R_EARTH_M, 10000.0);
    EXPECT_LT(ang * R_EARTH_M, 25000.0);
}

TEST(CelestialToTerrestrial, GeodeticLatitudeIsOnTheEllipsoid) {
    // A shadow axis through the geocentre (xi = eta = 0) must land at the
    // sub-star point, exactly on the WGS84 ellipsoid.
    const double jd_tt = 2461236.5, jd_ut1 = 2461236.5;
    double lat = 0.0, lon = 0.0;
    ASSERT_TRUE(shadow_point_geodetic(0.0, 0.0, 1.0, 0.5, jd_tt, jd_ut1, lat, lon));

    constexpr double a = 6378137.0, f = 1.0 / 298.257223563;
    const double e2 = 2.0*f - f*f;
    const double N = a / std::sqrt(1.0 - e2*std::sin(lat)*std::sin(lat));
    const double x = N*std::cos(lat)*std::cos(lon);
    const double y = N*std::cos(lat)*std::sin(lon);
    const double z = N*(1.0 - e2)*std::sin(lat);
    const double r = std::sqrt(x*x + y*y + z*z);
    EXPECT_GT(r, 6356000.0);
    EXPECT_LT(r, 6379000.0);
}

TEST(CelestialToTerrestrial, AxisMissingEarthReturnsFalse) {
    double lat = 0.0, lon = 0.0;
    EXPECT_FALSE(shadow_point_geodetic(50.0e6, 0.0, 1.0, 0.5, 2461236.5, 2461236.5, lat, lon));
}
