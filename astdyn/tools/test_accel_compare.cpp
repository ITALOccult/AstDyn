/**
 * Confronto accelerazioni tra i due metodi
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double k = 0.01720209895;
constexpr double k2 = k * k;
constexpr double AU_KM = 149597870.7;
constexpr double DAY_SEC = 86400.0;

constexpr double GM_MERCURY = 4.9125474514508118e-11;
constexpr double GM_VENUS   = 7.2434524861627027e-10;
constexpr double GM_EMB     = 8.9970116036316091e-10;
constexpr double GM_MARS    = 9.5495351057792580e-11;
constexpr double GM_JUPITER = 2.8253458420837436e-07;
constexpr double GM_SATURN  = 8.4597151856806587e-08;
constexpr double GM_URANUS  = 1.2920249167819693e-08;
constexpr double GM_NEPTUNE = 1.5243589007842762e-08;

using Vec3 = std::array<double, 3>;

double norm(const Vec3& v) { return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }

// Pianeta (semplificato - solo Terra, I=0)
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
        default: return {0, 0, 0};
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
    return {x, std::cos(eps)*y - std::sin(eps)*z, std::sin(eps)*y + std::cos(eps)*z};
}

// Metodo 1: come nel test che funziona (RK4)
Vec3 accel_v1(double t, const Vec3& r) {
    double r3 = std::pow(norm(r), 3);
    Vec3 acc = {-k2*r[0]/r3, -k2*r[1]/r3, -k2*r[2]/r3};
    
    double GM[] = {0, GM_MERCURY, GM_VENUS, GM_EMB, GM_MARS, GM_JUPITER, GM_SATURN, GM_URANUS, GM_NEPTUNE};
    for (int p = 1; p <= 8; p++) {
        Vec3 rp = getPlanetPosition(t, p);
        Vec3 d = {r[0]-rp[0], r[1]-rp[1], r[2]-rp[2]};
        double d3 = std::pow(norm(d), 3);
        double rp3 = std::pow(norm(rp), 3);
        acc[0] -= GM[p] * (d[0]/d3 + rp[0]/rp3);
        acc[1] -= GM[p] * (d[1]/d3 + rp[1]/rp3);
        acc[2] -= GM[p] * (d[2]/d3 + rp[2]/rp3);
    }
    return acc;
}

// Metodo 2: come nel propagator RKF78
Vec3 accel_v2(double t, const Vec3& r) {
    double r_norm = norm(r);
    double r3 = r_norm * r_norm * r_norm;
    Vec3 acc = {-k2*r[0]/r3, -k2*r[1]/r3, -k2*r[2]/r3};
    
    constexpr double GM[8] = {GM_MERCURY, GM_VENUS, GM_EMB, GM_MARS, GM_JUPITER, GM_SATURN, GM_URANUS, GM_NEPTUNE};
    for (int i = 0; i < 8; i++) {
        Vec3 r_planet = getPlanetPosition(t, i+1);
        Vec3 dr = {r[0]-r_planet[0], r[1]-r_planet[1], r[2]-r_planet[2]};
        double dr_norm = norm(dr);
        double rp_norm = norm(r_planet);
        double dr3 = dr_norm * dr_norm * dr_norm;
        double rp3 = rp_norm * rp_norm * rp_norm;
        acc[0] += -GM[i] * (dr[0]/dr3 + r_planet[0]/rp3);
        acc[1] += -GM[i] * (dr[1]/dr3 + r_planet[1]/rp3);
        acc[2] += -GM[i] * (dr[2]/dr3 + r_planet[2]/rp3);
    }
    return acc;
}

int main() {
    std::cout << std::scientific << std::setprecision(15);
    
    double jd = 2461000.5;
    Vec3 r = {1.619720632845141E+08 / AU_KM, 4.290918361163563E+08 / AU_KM, 1.711012428418700E+08 / AU_KM};
    
    Vec3 a1 = accel_v1(jd, r);
    Vec3 a2 = accel_v2(jd, r);
    
    std::cout << "Accelerazione v1 (RK4): [" << a1[0] << ", " << a1[1] << ", " << a1[2] << "]\n";
    std::cout << "Accelerazione v2 (RKF78): [" << a2[0] << ", " << a2[1] << ", " << a2[2] << "]\n";
    
    Vec3 diff = {a1[0]-a2[0], a1[1]-a2[1], a1[2]-a2[2]};
    std::cout << "Differenza: [" << diff[0] << ", " << diff[1] << ", " << diff[2] << "]\n";
    std::cout << "Norma diff: " << norm(diff) << "\n";
    
    return 0;
}
