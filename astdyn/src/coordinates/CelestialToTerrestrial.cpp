#include "astdyn/coordinates/CelestialToTerrestrial.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>

namespace astdyn::coordinates {
namespace {

constexpr double DAS2R  = 4.848136811095359935899141e-6; // arcsec -> rad
constexpr double D2PI   = 6.283185307179586476925287;
constexpr double TURNAS = 1296000.0;
constexpr double U2R    = DAS2R / 1e7;  // 0.1 microarcsec -> rad

inline Eigen::Matrix3d Rx(double a){ double s=std::sin(a),c=std::cos(a);
    Eigen::Matrix3d m; m << 1,0,0,  0,c,s,  0,-s,c; return m; }
inline Eigen::Matrix3d Ry(double a){ double s=std::sin(a),c=std::cos(a);
    Eigen::Matrix3d m; m << c,0,-s, 0,1,0,  s,0,c; return m; }
inline Eigen::Matrix3d Rz(double a){ double s=std::sin(a),c=std::cos(a);
    Eigen::Matrix3d m; m << c,s,0, -s,c,0,  0,0,1; return m; }

// IAU 2000B luni-solar nutation, 20 largest terms.
// nl, nl', nF, nD, nOm, ps, ps_t, pc, ec, ec_t, es   [units of 0.1 uas]
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
constexpr double DPPLAN = -0.135e-3 * DAS2R;  // fixed offsets in lieu of planetary terms
constexpr double DEPLAN =  0.388e-3 * DAS2R;

} // namespace

void precession_fw06(double t, double& gamb, double& phib, double& psib, double& epsa)
{
    gamb = (   -0.052928 + (10.556378 + (0.4932044 + (-0.00031238
             + (-0.000002788 + 0.0000000260*t)*t)*t)*t)*t) * DAS2R;
    phib = (84381.412819 + (-46.811016 + (0.0511268 + (0.00053289
             + (-0.000000440 - 0.0000000176*t)*t)*t)*t)*t) * DAS2R;
    psib = (   -0.041775 + (5038.481484 + (1.5584175 + (-0.00018522
             + (-0.000026452 - 0.0000000148*t)*t)*t)*t)*t) * DAS2R;
    epsa = (84381.406    + (-46.836769 + (-0.0001831 + (0.00200340
             + (-0.000000576 - 0.0000000434*t)*t)*t)*t)*t) * DAS2R;
}

void nutation_00b(double t, double& dpsi, double& deps)
{
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

void cip_xy(double jd_tt, double& X, double& Y)
{
    const double t = (jd_tt - 2451545.0) / 36525.0;
    double gamb, phib, psib, epsa, dpsi, deps;
    precession_fw06(t, gamb, phib, psib, epsa);
    nutation_00b(t, dpsi, deps);
    const double fj2 = -2.7774e-6 * t;             // IAU 2006 adjustment of IAU 2000 nutation
    dpsi += dpsi * (0.4697e-6 + fj2);
    deps += deps * fj2;
    // Fukushima-Williams angles -> bias-precession-nutation matrix (cf. iauFw2m).
    const Eigen::Matrix3d rnpb =
        Rx(-(epsa + deps)) * Rz(-(psib + dpsi)) * Rx(phib) * Rz(gamb);
    X = rnpb(2,0);
    Y = rnpb(2,1);
}

Eigen::Matrix3d celestial_to_intermediate(double jd_tt)
{
    double X, Y;
    cip_xy(jd_tt, X, Y);
    const double s  = -X * Y / 2.0;   // CIO locator; residual < 1 m at the surface
    const double r2 = X*X + Y*Y;
    const double e  = (r2 > 0.0) ? std::atan2(Y, X) : 0.0;
    const double d  = std::atan(std::sqrt(r2 / (1.0 - r2)));
    return Rz(-(e + s)) * Ry(d) * Rz(e);
}

double earth_rotation_angle(double jd_ut1)
{
    const double tu = jd_ut1 - 2451545.0;
    const double f  = std::fmod(jd_ut1, 1.0) + 0.5;  // split improves precision
    double era = D2PI * (f + 0.7790572732640 + 0.00273781191135448 * tu);
    era = std::fmod(era, D2PI);
    return (era < 0.0) ? era + D2PI : era;
}

Eigen::Matrix3d gcrf_to_itrf(double jd_tt, double jd_ut1)
{
    return Rz(earth_rotation_angle(jd_ut1)) * celestial_to_intermediate(jd_tt);
}

bool shadow_point_geodetic(double xi_m, double eta_m,
                           double ra_rad, double dec_rad,
                           double jd_tt, double jd_ut1,
                           double& lat_rad, double& lon_rad)
{
    const double a  = constants::R_EARTH_EQUATORIAL * 1000.0;   // [m]
    const double fl = constants::EARTH_FLATTENING;
    const double e2 = 2.0*fl - fl*fl;

    // Fundamental-plane basis in GCRF: k toward the star, i East, j North.
    const double as = ra_rad, ds = dec_rad;
    const Eigen::Vector3d k(std::cos(ds)*std::cos(as), std::cos(ds)*std::sin(as), std::sin(ds));
    const Eigen::Vector3d i(-std::sin(as), std::cos(as), 0.0);
    const Eigen::Vector3d j = k.cross(i).normalized();

    // Rotate the basis into ITRF, where the ellipsoid is defined.
    const Eigen::Matrix3d C = gcrf_to_itrf(jd_tt, jd_ut1);
    const Eigen::Vector3d ke = C*k, ie = C*i, je = C*j;

    // Intersect P(s) = xi*ie + eta*je + s*ke with the ellipsoid.
    // Scaling z by 1/(1-f) turns the ellipsoid into a sphere of radius a.
    const double u = 1.0 / (1.0 - fl);
    auto sc = [u](const Eigen::Vector3d& v){ return Eigen::Vector3d(v.x(), v.y(), u*v.z()); };
    const Eigen::Vector3d Ks = sc(ke);
    const Eigen::Vector3d Q0 = xi_m*sc(ie) + eta_m*sc(je);

    const double A = Ks.dot(Ks);
    const double B = 2.0 * Q0.dot(Ks);
    const double Cc = Q0.dot(Q0) - a*a;
    const double disc = B*B - 4.0*A*Cc;
    if (disc < 0.0 || A <= 0.0) return false;      // the axis misses the Earth

    const double s = (-B + std::sqrt(disc)) / (2.0*A);   // star-facing root
    const Eigen::Vector3d P = xi_m*ie + eta_m*je + s*ke; // ITRF, unscaled

    // Geodetic latitude: exact for a point lying on the ellipsoid.
    lat_rad = std::atan2(P.z(), (1.0 - e2) * std::hypot(P.x(), P.y()));
    lon_rad = std::atan2(P.y(), P.x());
    return true;
}

} // namespace astdyn::coordinates
