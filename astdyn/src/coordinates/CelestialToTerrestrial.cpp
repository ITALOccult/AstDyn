#include "astdyn/coordinates/CelestialToTerrestrial.hpp"
#include <cmath>

namespace astdyn::coordinates {
namespace {

using math::Rotation;
using core::GCRF; using core::CIRS; using core::TIRS; using core::ITRF;

constexpr double DAS2R  = 4.848136811095359935899141e-6;
constexpr double D2PI   = 6.283185307179586476925287;
constexpr double TURNAS = 1296000.0;
constexpr double U2R    = DAS2R / 1e7;

// Elementary rotations, SOFA convention (they act on coordinates).
// Untagged on purpose: they are building blocks, never exposed.
inline Eigen::Matrix3d Rx(double a){ double s=std::sin(a),c=std::cos(a);
    Eigen::Matrix3d m; m << 1,0,0,  0,c,s,  0,-s,c; return m; }
inline Eigen::Matrix3d Ry(double a){ double s=std::sin(a),c=std::cos(a);
    Eigen::Matrix3d m; m << c,0,-s, 0,1,0,  s,0,c; return m; }
inline Eigen::Matrix3d Rz(double a){ double s=std::sin(a),c=std::cos(a);
    Eigen::Matrix3d m; m << c,s,0, -s,c,0,  0,0,1; return m; }

struct NutTerm { int nl,nlp,nf,nd,nom; double sp,spt,cp,ce,cet,se; };
constexpr NutTerm XLS[] = {
 { 0, 0, 0, 0,1, -172064161.0,-174666.0, 33386.0, 92052331.0, 9086.0, 15377.0},
 { 0, 0, 2,-2,2,  -13170906.0,  -1675.0,-13696.0,  5730336.0,-3015.0, -4587.0},
 { 0, 0, 2, 0,2,   -2276413.0,   -234.0,  2796.0,   978459.0, -485.0,  1374.0},
 { 0, 0, 0, 0,2,    2074554.0,    207.0,  -698.0,  -897492.0,  470.0,  -291.0},
 { 0, 1, 0, 0,0,    1475877.0,  -3633.0, 11817.0,    73871.0, -184.0, -1924.0},
 { 0, 1, 2,-2,2,    -516821.0,   1226.0,  -524.0,   224386.0, -677.0,  -174.0},
 { 1, 0, 0, 0,0,     711159.0,     73.0,  -872.0,    -6750.0,    0.0,   358.0},
 { 0, 0, 2, 0,1,    -387298.0,   -367.0,   380.0,   200728.0,   18.0,   318.0},
 { 1, 0, 2, 0,2,    -301461.0,    -36.0,   816.0,   129025.0,  -63.0,   367.0},
 { 0,-1, 2,-2,2,     215829.0,   -494.0,   111.0,   -95929.0,  299.0,   132.0},
 { 0, 0, 2,-2,1,     128227.0,    137.0,   181.0,   -68982.0,   -9.0,    39.0},
 {-1, 0, 2, 0,2,     123457.0,     11.0,    19.0,   -53311.0,   32.0,    -4.0},
 {-1, 0, 0, 2,0,     156994.0,     10.0,  -168.0,    -1235.0,    0.0,    82.0},
 { 1, 0, 0, 0,1,      63110.0,     63.0,    27.0,   -33228.0,    0.0,    -9.0},
 {-1, 0, 0, 0,1,     -57976.0,    -63.0,  -189.0,    31429.0,    0.0,   -75.0},
 {-1, 0, 2, 2,2,     -59641.0,    -11.0,   149.0,    25543.0,  -11.0,    66.0},
 { 1, 0, 2, 0,1,     -51613.0,    -42.0,   129.0,    26366.0,    0.0,    78.0},
 {-2, 0, 2, 0,1,      45893.0,     50.0,    31.0,   -24236.0,  -10.0,    20.0},
 { 0, 0, 0, 2,0,      63384.0,     11.0,  -150.0,    -1220.0,    0.0,    29.0},
 { 0, 0, 2, 2,2,     -38571.0,     -1.0,   158.0,    16452.0,  -11.0,    68.0},
};
constexpr double DPPLAN = -0.135e-3 * DAS2R;
constexpr double DEPLAN =  0.388e-3 * DAS2R;

void precession_fw06(double t, double& gamb, double& phib, double& psib, double& epsa) {
    gamb = (   -0.052928 + (10.556378 + (0.4932044 + (-0.00031238
             + (-0.000002788 + 0.0000000260*t)*t)*t)*t)*t) * DAS2R;
    phib = (84381.412819 + (-46.811016 + (0.0511268 + (0.00053289
             + (-0.000000440 - 0.0000000176*t)*t)*t)*t)*t) * DAS2R;
    psib = (   -0.041775 + (5038.481484 + (1.5584175 + (-0.00018522
             + (-0.000026452 - 0.0000000148*t)*t)*t)*t)*t) * DAS2R;
    epsa = (84381.406    + (-46.836769 + (-0.0001831 + (0.00200340
             + (-0.000000576 - 0.0000000434*t)*t)*t)*t)*t) * DAS2R;
}

void nutation_00b(double t, double& dpsi, double& deps) {
    const double el  = std::fmod(485868.249036 + 1717915923.2178*t, TURNAS) * DAS2R;
    const double elp = std::fmod(1287104.79305 +  129596581.0481*t, TURNAS) * DAS2R;
    const double f   = std::fmod(335779.526232 + 1739527262.8478*t, TURNAS) * DAS2R;
    const double d   = std::fmod(1072260.70369 + 1602961601.2090*t, TURNAS) * DAS2R;
    const double om  = std::fmod(450160.398036 -    6962890.5431*t, TURNAS) * DAS2R;
    double dp = 0.0, de = 0.0;
    for (int i = (int)(sizeof(XLS)/sizeof(XLS[0])) - 1; i >= 0; --i) {
        const NutTerm& x = XLS[i];
        const double arg = std::fmod(x.nl*el + x.nlp*elp + x.nf*f + x.nd*d + x.nom*om, D2PI);
        const double sarg = std::sin(arg), carg = std::cos(arg);
        dp += (x.sp + x.spt*t)*sarg + x.cp*carg;
        de += (x.ce + x.cet*t)*carg + x.se*sarg;
    }
    dpsi = dp*U2R + DPPLAN;
    deps = de*U2R + DEPLAN;
}

} // namespace

void cip_xy(time::EpochTT tt, astrometry::Angle& X, astrometry::Angle& Y) {
    const double t = (tt.jd() - 2451545.0) / 36525.0;
    double gamb, phib, psib, epsa, dpsi, deps;
    precession_fw06(t, gamb, phib, psib, epsa);
    nutation_00b(t, dpsi, deps);
    const double fj2 = -2.7774e-6 * t;
    dpsi += dpsi * (0.4697e-6 + fj2);
    deps += deps * fj2;
    const Eigen::Matrix3d rnpb =
        Rx(-(epsa + deps)) * Rz(-(psib + dpsi)) * Rx(phib) * Rz(gamb);
    X = astrometry::Angle::from_rad(rnpb(2,0));
    Y = astrometry::Angle::from_rad(rnpb(2,1));
}

math::Rotation<core::GCRF, core::CIRS> cip_rotation(time::EpochTT tt) {
    astrometry::Angle Xa, Ya;
    cip_xy(tt, Xa, Ya);
    const double X = Xa.to_rad(), Y = Ya.to_rad();
    const double s  = -X * Y / 2.0;
    const double r2 = X*X + Y*Y;
    const double e  = (r2 > 0.0) ? std::atan2(Y, X) : 0.0;
    const double d  = std::atan(std::sqrt(r2 / (1.0 - r2)));
    return math::Rotation<core::GCRF, core::CIRS>::from_matrix_unchecked(
        Rz(-(e + s)) * Ry(d) * Rz(e));
}

astrometry::Angle earth_rotation_angle(time::EpochUT1 ut1) {
    const double jd = ut1.jd();
    const double tu = jd - 2451545.0;
    const double f  = std::fmod(jd, 1.0);
    double era = D2PI * (f + 0.7790572732640 + 0.00273781191135448 * tu);
    era = std::fmod(era, D2PI);
    if (era < 0.0) era += D2PI;
    return astrometry::Angle::from_rad(era);
}

math::Rotation<core::CIRS, core::TIRS> era_rotation(time::EpochUT1 ut1) {
    return math::Rotation<core::CIRS, core::TIRS>::from_matrix_unchecked(
        Rz(earth_rotation_angle(ut1).to_rad()));
}

math::Rotation<core::TIRS, core::ITRF> polar_motion_rotation() {
    return math::Rotation<core::TIRS, core::ITRF>::identity();
}

math::Rotation<core::GCRF, core::ITRF> gcrf_to_itrf(time::EpochTT tt, time::EpochUT1 ut1) {
    // Reads exactly like the IAU formula: [ITRF] = W * R(ERA) * Q(t) * [GCRF].
    // Any other order is a compile error.
    return polar_motion_rotation() * era_rotation(ut1) * cip_rotation(tt);
}

} // namespace astdyn::coordinates

namespace astdyn::coordinates {

std::optional<SurfacePoint> shadow_point(
    physics::Distance xi, physics::Distance eta,
    const math::Direction<core::GCRF>& star_dir,
    time::EpochTT tt, time::EpochUT1 ut1,
    const Ellipsoid& E)
{
    using math::Direction; using math::Vec3;

    // Fundamental-plane basis, built in GCRF and tagged as such.
    const Direction<GCRF> k = star_dir;
    const Direction<GCRF> i = Direction<GCRF>::from_xyz(-k.y(), k.x(), 0.0);
    const Direction<GCRF> j = k.cross(i);

    // Rotate the basis into ITRF, where the ellipsoid is defined.
    // The compiler enforces that this is the full W*R*Q chain.
    const auto C = gcrf_to_itrf(tt, ut1);
    const Direction<ITRF> ke = C * k, ie = C * i, je = C * j;

    // Scaling z by 1/(1-f) maps the ellipsoid onto a sphere of radius a.
    const double u = 1.0 / (1.0 - E.f());
    auto sc = [u](const Eigen::Vector3d& v){ return Eigen::Vector3d(v.x(), v.y(), u*v.z()); };
    const Eigen::Vector3d Ks = sc(ke.eigen());
    const Eigen::Vector3d Q0 = xi.to_m()*sc(ie.eigen()) + eta.to_m()*sc(je.eigen());

    const double A = Ks.dot(Ks);
    const double B = 2.0 * Q0.dot(Ks);
    const double Cc = Q0.dot(Q0) - E.a_m()*E.a_m();
    const double disc = B*B - 4.0*A*Cc;
    if (disc < 0.0 || A <= 0.0) return std::nullopt;      // the axis misses the Earth

    const double s = (-B + std::sqrt(disc)) / (2.0*A);    // star-facing root
    const Vec3<ITRF> P = ie*xi.to_m() + je*eta.to_m() + ke*s;

    // Geocentric first -- that is what the geometry gives -- then an explicit,
    // ellipsoid-aware conversion. The type system forbids skipping this step.
    const auto phi_gc = GeocentricLatitude::from_angle(astrometry::Angle::from_rad(
        std::atan2(P.z_m(), std::hypot(P.x_m(), P.y_m()))));

    return SurfacePoint{ E.to_geodetic(phi_gc),
                         astrometry::Angle::from_rad(std::atan2(P.y_m(), P.x_m())) };
}

} // namespace astdyn::coordinates
