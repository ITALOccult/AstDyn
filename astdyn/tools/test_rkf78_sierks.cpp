/**
 * Test Sierks usando il propagator RKF78 con stato JPL diretto
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <vector>

namespace astdyn {

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double k = 0.01720209895;
constexpr double k2 = k * k;
constexpr double AU_KM = 149597870.7;
constexpr double DAY_SEC = 86400.0;
constexpr double C_LIGHT = 173.1446326846693;
constexpr double c2 = C_LIGHT * C_LIGHT;

constexpr double GM_MERCURY = 4.9125474514508118e-11;
constexpr double GM_VENUS   = 7.2434524861627027e-10;
constexpr double GM_EMB     = 8.9970116036316091e-10;
constexpr double GM_MARS    = 9.5495351057792580e-11;
constexpr double GM_JUPITER = 2.8253458420837436e-07;
constexpr double GM_SATURN  = 8.4597151856806587e-08;
constexpr double GM_URANUS  = 1.2920249167819693e-08;
constexpr double GM_NEPTUNE = 1.5243589007842762e-08;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vec3 operator+(const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
    Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
    Vec3& operator+=(const Vec3& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    double norm2() const { return x*x + y*y + z*z; }
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
};

struct State {
    Vec3 r, v;
    State() = default;
    State(const Vec3& r_, const Vec3& v_) : r(r_), v(v_) {}
    State operator+(const State& s) const { return State(r+s.r, v+s.v); }
    State operator*(double k) const { return State(r*k, v*k); }
};

// Effemeridi pianeti ICRF
Vec3 getPlanetPosition(double jd, int planet) {
    double T = (jd - 2451545.0) / 36525.0;
    double a, e, I, L, omega_bar, Omega;
    
    switch(planet) {
        case 1: a=0.38709927; e=0.20563593; I=7.00497902; L=252.25032350+149472.67411175*T; omega_bar=77.45779628; Omega=48.33076593; break;
        case 2: a=0.72333566; e=0.00677672; I=3.39467605; L=181.97909950+58517.81538729*T; omega_bar=131.60246718; Omega=76.67984255; break;
        case 3: a=1.00000261; e=0.01671123; I=0; L=100.46457166+35999.37244981*T; omega_bar=102.93768193; Omega=0; break;
        case 4: a=1.52371034; e=0.09339410; I=1.84969142; L=-4.55343205+19140.30268499*T; omega_bar=-23.94362959; Omega=49.55953891; break;
        case 5: a=5.20288700; e=0.04838624; I=1.30439695; L=34.39644051+3034.74612775*T; omega_bar=14.72847983; Omega=100.47390909; break;
        case 6: a=9.53667594; e=0.05386179; I=2.48599187; L=49.95424423+1222.49362201*T; omega_bar=92.59887831; Omega=113.66242448; break;
        case 7: a=19.18916464; e=0.04725744; I=0.77263783; L=313.23810451+428.48202785*T; omega_bar=170.95427630; Omega=74.01692503; break;
        case 8: a=30.06992276; e=0.00859048; I=1.77004347; L=-55.12002969+218.45945325*T; omega_bar=44.96476227; Omega=131.78422574; break;
        default: return Vec3();
    }
    
    double i_rad = I * DEG2RAD;
    double Omega_rad = Omega * DEG2RAD;
    double omega = (omega_bar - Omega) * DEG2RAD;
    double M = (L - omega_bar) * DEG2RAD;
    M = std::fmod(M, 2*PI); if (M < 0) M += 2*PI;
    
    double E = M;
    for (int i = 0; i < 15; i++) {
        double dE = (M - E + e * std::sin(E)) / (1 - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < 1e-14) break;
    }
    
    double nu = 2.0 * std::atan2(std::sqrt(1+e) * std::sin(E/2), std::sqrt(1-e) * std::cos(E/2));
    double r = a * (1 - e * std::cos(E));
    double x_orb = r * std::cos(nu), y_orb = r * std::sin(nu);
    
    double cO = std::cos(Omega_rad), sO = std::sin(Omega_rad);
    double ci = std::cos(i_rad), si = std::sin(i_rad);
    double cw = std::cos(omega), sw = std::sin(omega);
    
    double x = (cO*cw - sO*sw*ci) * x_orb + (-cO*sw - sO*cw*ci) * y_orb;
    double y = (sO*cw + cO*sw*ci) * x_orb + (-sO*sw + cO*cw*ci) * y_orb;
    double z = (sw*si) * x_orb + (cw*si) * y_orb;
    
    double eps = 23.4392911 * DEG2RAD;
    return Vec3(x, std::cos(eps)*y - std::sin(eps)*z, std::sin(eps)*y + std::cos(eps)*z);
}

// Accelerazione
Vec3 acceleration(double t, const Vec3& r, const Vec3& v) {
    double r_norm = r.norm();
    double r3 = r_norm * r_norm * r_norm;
    Vec3 acc = r * (-k2 / r3);
    
    constexpr double GM[8] = {GM_MERCURY, GM_VENUS, GM_EMB, GM_MARS, GM_JUPITER, GM_SATURN, GM_URANUS, GM_NEPTUNE};
    for (int i = 0; i < 8; i++) {
        Vec3 r_planet = getPlanetPosition(t, i+1);
        Vec3 dr = r - r_planet;
        double dr_norm = dr.norm();
        double rp_norm = r_planet.norm();
        acc += dr * (-GM[i] / (dr_norm * dr_norm * dr_norm));
        acc += r_planet * (-GM[i] / (rp_norm * rp_norm * rp_norm));
    }
    return acc;
}

State deriv(double t, const State& y) {
    return State(y.v, acceleration(t, y.r, y.v));
}

// RKF78 - versione corretta
State step_rkf78(const State& y, double t, double& h, bool& accepted, double tol) {
    static constexpr double c[13] = {0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 1.0/2.0, 5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0, 1.0};
    
    State k[13];
    k[0] = deriv(t, y);
    
    State y1 = y + k[0] * (h * c[1]);
    k[1] = deriv(t + c[1]*h, y1);
    
    State y2 = y + (k[0] * (1.0/36.0) + k[1] * (3.0/36.0)) * h;
    k[2] = deriv(t + c[2]*h, y2);
    
    State y3 = y + (k[0] * (1.0/24.0) + k[2] * (3.0/24.0)) * h;
    k[3] = deriv(t + c[3]*h, y3);
    
    State y4 = y + (k[0] * (20.0/48.0) + k[2] * (-75.0/48.0) + k[3] * (75.0/48.0)) * h;
    k[4] = deriv(t + c[4]*h, y4);
    
    State y5 = y + (k[0] * (1.0/20.0) + k[3] * (5.0/20.0) + k[4] * (4.0/20.0)) * h;
    k[5] = deriv(t + c[5]*h, y5);
    
    State y6 = y + (k[0] * (-25.0/108.0) + k[3] * (125.0/108.0) + k[4] * (-260.0/108.0) + k[5] * (250.0/108.0)) * h;
    k[6] = deriv(t + c[6]*h, y6);
    
    State y7 = y + (k[0] * (31.0/300.0) + k[4] * (61.0/225.0) + k[5] * (-2.0/9.0) + k[6] * (13.0/900.0)) * h;
    k[7] = deriv(t + c[7]*h, y7);
    
    State y8 = y + (k[0] * 2.0 + k[3] * (-53.0/6.0) + k[4] * (704.0/45.0) + k[5] * (-107.0/9.0) + k[6] * (67.0/90.0) + k[7] * 3.0) * h;
    k[8] = deriv(t + c[8]*h, y8);
    
    State y9 = y + (k[0] * (-91.0/108.0) + k[3] * (23.0/108.0) + k[4] * (-976.0/135.0) + k[5] * (311.0/54.0) + k[6] * (-19.0/60.0) + k[7] * (17.0/6.0) + k[8] * (-1.0/12.0)) * h;
    k[9] = deriv(t + c[9]*h, y9);
    
    State y10 = y + (k[0] * (2383.0/4100.0) + k[3] * (-341.0/164.0) + k[4] * (4496.0/1025.0) + k[5] * (-301.0/82.0) + k[6] * (2133.0/4100.0) + k[7] * (45.0/82.0) + k[8] * (45.0/164.0) + k[9] * (18.0/41.0)) * h;
    k[10] = deriv(t + c[10]*h, y10);
    
    State y11 = y + (k[0] * (3.0/205.0) + k[5] * (-6.0/41.0) + k[6] * (-3.0/205.0) + k[7] * (-3.0/41.0) + k[8] * (3.0/41.0) + k[9] * (6.0/41.0)) * h;
    k[11] = deriv(t + c[11]*h, y11);
    
    State y12 = y + (k[0] * (-1777.0/4100.0) + k[3] * (-341.0/164.0) + k[4] * (4496.0/1025.0) + k[5] * (-289.0/82.0) + k[6] * (2193.0/4100.0) + k[7] * (51.0/82.0) + k[8] * (33.0/164.0) + k[9] * (12.0/41.0) + k[11] * 1.0) * h;
    k[12] = deriv(t + c[12]*h, y12);
    
    // Soluzione ordine 7
    static constexpr double b7[13] = {41.0/840.0, 0, 0, 0, 0, 34.0/105.0, 9.0/35.0, 9.0/35.0, 9.0/280.0, 9.0/280.0, 41.0/840.0, 0, 0};
    // Soluzione ordine 8
    static constexpr double b8[13] = {0, 0, 0, 0, 0, 34.0/105.0, 9.0/35.0, 9.0/35.0, 9.0/280.0, 9.0/280.0, 0, 41.0/840.0, 41.0/840.0};
    
    State y_new;
    State err;
    for (int i = 0; i < 13; i++) {
        y_new.r += k[i].r * (h * b8[i]);
        y_new.v += k[i].v * (h * b8[i]);
        err.r += k[i].r * (h * (b8[i] - b7[i]));
        err.v += k[i].v * (h * (b8[i] - b7[i]));
    }
    y_new.r = y_new.r + y.r;
    y_new.v = y_new.v + y.v;
    
    double error = std::max(err.r.norm(), err.v.norm() * 100);
    
    if (error < tol || std::abs(h) <= 1e-6) {
        accepted = true;
        if (error > 0) {
            double factor = 0.9 * std::pow(tol / error, 1.0/8.0);
            factor = std::max(0.1, std::min(5.0, factor));
            h *= factor;
        }
    } else {
        accepted = false;
        double factor = 0.9 * std::pow(tol / error, 1.0/8.0);
        factor = std::max(0.1, std::min(0.9, factor));
        h *= factor;
    }
    h = std::max(1e-6, std::min(1.0, std::abs(h))) * (h > 0 ? 1 : -1);
    
    return y_new;
}

State propagate(const State& y0, double t0, double t1, double tol = 1e-12) {
    double dt = t1 - t0;
    double h = dt > 0 ? std::min(0.1, dt/10) : std::max(-0.1, dt/10);
    double t = t0;
    State y = y0;
    int steps = 0;
    
    while ((dt > 0 && t < t1) || (dt < 0 && t > t1)) {
        if (dt > 0 && t + h > t1) h = t1 - t;
        if (dt < 0 && t + h < t1) h = t1 - t;
        
        bool accepted;
        State y_new = step_rkf78(y, t, h, accepted, tol);
        
        if (accepted) {
            t += h;
            y = y_new;
            steps++;
        }
    }
    std::cout << "  Steps: " << steps << "\n";
    return y;
}

} // namespace astdyn

int main() {
    using namespace astdyn;
    std::cout << std::fixed << std::setprecision(8);
    
    // Stato ICRF da JPL (Sierks, epoch 2461000.5)
    double jd0 = 2461000.5;
    State y0(
        Vec3(1.619720632845141E+08 / AU_KM, 4.290918361163563E+08 / AU_KM, 1.711012428418700E+08 / AU_KM),
        Vec3(-1.545063491937339E+01 * DAY_SEC / AU_KM, 4.181920917345073E+00 * DAY_SEC / AU_KM, 2.575783042670384E+00 * DAY_SEC / AU_KM)
    );
    
    double jd_target = 2461008.0913;
    
    std::cout << "=== TEST RKF78 CON STATO JPL ===\n\n";
    std::cout << "Stato iniziale:\n";
    std::cout << "  r = [" << y0.r.x << ", " << y0.r.y << ", " << y0.r.z << "]\n\n";
    
    State y1 = propagate(y0, jd0, jd_target);
    
    std::cout << "Stato finale:\n";
    std::cout << "  r = [" << y1.r.x << ", " << y1.r.y << ", " << y1.r.z << "]\n\n";
    
    Vec3 earth = getPlanetPosition(jd_target, 3);
    Vec3 geo = y1.r - earth;
    double delta = geo.norm();
    
    double lt = delta / C_LIGHT;
    geo = geo - y1.v * lt;
    delta = geo.norm();
    
    double ra = std::atan2(geo.y, geo.x);
    if (ra < 0) ra += 2*PI;
    double dec = std::asin(geo.z / delta);
    
    double ra_h = ra * 12.0 / PI;
    int ra_hh = (int)ra_h;
    double ra_m = (ra_h - ra_hh) * 60.0;
    int ra_mm = (int)ra_m;
    double ra_ss = (ra_m - ra_mm) * 60.0;
    
    double dec_d = std::abs(dec) * RAD2DEG;
    char sign = dec >= 0 ? '+' : '-';
    int dec_dd = (int)dec_d;
    double dec_m = (dec_d - dec_dd) * 60.0;
    int dec_mm = (int)dec_m;
    double dec_ss = (dec_m - dec_mm) * 60.0;
    
    std::cout << "POSIZIONE:\n";
    printf("  RA  = %02d %02d %06.3f\n", ra_hh, ra_mm, ra_ss);
    printf("  Dec = %c%02d %02d %05.2f\n", sign, dec_dd, dec_mm, dec_ss);
    std::cout << "\nJPL: RA = 04 53 11.25   Dec = +20 19 25.8\n";
    
    return 0;
}
