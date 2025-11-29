/**
 * Confronto stato da elementi vs stato JPL
 */
#include <iostream>
#include <iomanip>
#include <cmath>

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double k = 0.01720209895;
constexpr double k2 = k * k;
constexpr double AU_KM = 149597870.7;
constexpr double DAY_SEC = 86400.0;

int main() {
    std::cout << std::fixed << std::setprecision(10);
    
    // Stato da JPL (ICRF)
    double jpl_x  = 1.619720632845141E+08 / AU_KM;
    double jpl_y  = 4.290918361163563E+08 / AU_KM;
    double jpl_z  = 1.711012428418700E+08 / AU_KM;
    double jpl_vx = -1.545063491937339E+01 * DAY_SEC / AU_KM;
    double jpl_vy =  4.181920917345073E+00 * DAY_SEC / AU_KM;
    double jpl_vz =  2.575783042670384E+00 * DAY_SEC / AU_KM;
    
    std::cout << "STATO JPL (ICRF, epoch 2461000.5):\n";
    std::cout << "  r = [" << jpl_x << ", " << jpl_y << ", " << jpl_z << "]\n";
    std::cout << "  v = [" << jpl_vx << ", " << jpl_vy << ", " << jpl_vz << "]\n\n";
    
    // Elementi orbitali (eclittici)
    double a = 3.1754733;
    double e = 0.0454207;
    double i = 2.9046 * DEG2RAD;
    double Omega = 104.16243 * DEG2RAD;
    double omega = 100.5141 * DEG2RAD;
    double M0 = 229.79088 * DEG2RAD;
    
    // Keplero
    double E = M0;
    for (int iter = 0; iter < 15; iter++) {
        double dE = (E - e * std::sin(E) - M0) / (1 - e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-14) break;
    }
    
    double sin_E = std::sin(E);
    double cos_E = std::cos(E);
    double sqrt_1_e2 = std::sqrt(1.0 - e * e);
    double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - e);
    double r = a * (1.0 - e * cos_E);
    
    double x_orb = r * std::cos(nu);
    double y_orb = r * std::sin(nu);
    double v_factor = std::sqrt(k2 * a) / r;
    double vx_orb = -v_factor * sin_E;
    double vy_orb = v_factor * sqrt_1_e2 * cos_E;
    
    double cO = std::cos(Omega), sO = std::sin(Omega);
    double cw = std::cos(omega), sw = std::sin(omega);
    double ci = std::cos(i), si = std::sin(i);
    
    double P11 = cO*cw - sO*sw*ci, P12 = -cO*sw - sO*cw*ci;
    double P21 = sO*cw + cO*sw*ci, P22 = -sO*sw + cO*cw*ci;
    double P31 = sw*si,            P32 = cw*si;
    
    double x_ecl = P11*x_orb + P12*y_orb;
    double y_ecl = P21*x_orb + P22*y_orb;
    double z_ecl = P31*x_orb + P32*y_orb;
    
    double vx_ecl = P11*vx_orb + P12*vy_orb;
    double vy_ecl = P21*vx_orb + P22*vy_orb;
    double vz_ecl = P31*vx_orb + P32*vy_orb;
    
    std::cout << "STATO DA ELEMENTI (Eclittico):\n";
    std::cout << "  r = [" << x_ecl << ", " << y_ecl << ", " << z_ecl << "]\n";
    std::cout << "  v = [" << vx_ecl << ", " << vy_ecl << ", " << vz_ecl << "]\n\n";
    
    // Rotazione → ICRF
    double eps = 23.4392911 * DEG2RAD;
    double ce = std::cos(eps), se = std::sin(eps);
    
    double x_icrf = x_ecl;
    double y_icrf = ce * y_ecl - se * z_ecl;
    double z_icrf = se * y_ecl + ce * z_ecl;
    
    double vx_icrf = vx_ecl;
    double vy_icrf = ce * vy_ecl - se * vz_ecl;
    double vz_icrf = se * vy_ecl + ce * vz_ecl;
    
    std::cout << "STATO DA ELEMENTI (ICRF dopo rotazione):\n";
    std::cout << "  r = [" << x_icrf << ", " << y_icrf << ", " << z_icrf << "]\n";
    std::cout << "  v = [" << vx_icrf << ", " << vy_icrf << ", " << vz_icrf << "]\n\n";
    
    // Differenze
    double dr = std::sqrt(std::pow(x_icrf-jpl_x,2) + std::pow(y_icrf-jpl_y,2) + std::pow(z_icrf-jpl_z,2));
    double dv = std::sqrt(std::pow(vx_icrf-jpl_vx,2) + std::pow(vy_icrf-jpl_vy,2) + std::pow(vz_icrf-jpl_vz,2));
    
    std::cout << "DIFFERENZE:\n";
    std::cout << "  Δr = " << dr << " AU = " << dr * AU_KM << " km\n";
    std::cout << "  Δv = " << dv << " AU/day = " << dv * AU_KM / DAY_SEC << " km/s\n";
    
    return 0;
}
