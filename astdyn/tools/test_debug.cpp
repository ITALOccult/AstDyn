/**
 * Test diagnostico per debug propagator Sierks
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double k = 0.01720209895;
constexpr double k2 = k * k;
constexpr double AU_KM = 149597870.7;
constexpr double C_LIGHT_AU = 173.1446326846693;

using Vec3 = std::array<double, 3>;

double norm(const Vec3& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Posizione Terra (semplificata)
Vec3 getEarth(double jd) {
    double T = (jd - 2451545.0) / 36525.0;
    double L = (100.46457166 + 35999.37244981*T) * DEG2RAD;
    double a = 1.00000261;
    double e = 0.01671123;
    double M = L - (102.93768193 + 0.32327364*T) * DEG2RAD;
    
    // Keplero
    double E = M;
    for (int i = 0; i < 15; i++) {
        double dE = (M - E + e * std::sin(E)) / (1 - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < 1e-14) break;
    }
    
    double nu = 2.0 * std::atan2(std::sqrt(1+e) * std::sin(E/2),
                                  std::sqrt(1-e) * std::cos(E/2));
    double r = a * (1 - e * std::cos(E));
    
    // Posizione eclittica
    double x_ecl = r * std::cos(L);
    double y_ecl = r * std::sin(L);
    double z_ecl = 0;
    
    // Rotazione eclittica → ICRF
    double eps = 23.4392911 * DEG2RAD;
    return {
        x_ecl,
        std::cos(eps) * y_ecl - std::sin(eps) * z_ecl,
        std::sin(eps) * y_ecl + std::cos(eps) * z_ecl
    };
}

// Elementi orbitali → stato ICRF
void elementsToICRF(double a, double e, double i_deg, double Omega_deg, 
                    double omega_deg, double M_deg, 
                    Vec3& pos, Vec3& vel) {
    double i = i_deg * DEG2RAD;
    double Omega = Omega_deg * DEG2RAD;
    double omega = omega_deg * DEG2RAD;
    double M0 = M_deg * DEG2RAD;
    
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
    
    // Matrice rotazione
    double cO = std::cos(Omega), sO = std::sin(Omega);
    double cw = std::cos(omega), sw = std::sin(omega);
    double ci = std::cos(i), si = std::sin(i);
    
    double P11 = cO*cw - sO*sw*ci, P12 = -cO*sw - sO*cw*ci;
    double P21 = sO*cw + cO*sw*ci, P22 = -sO*sw + cO*cw*ci;
    double P31 = sw*si,            P32 = cw*si;
    
    // Posizione/velocità eclittiche
    double x_ecl = P11*x_orb + P12*y_orb;
    double y_ecl = P21*x_orb + P22*y_orb;
    double z_ecl = P31*x_orb + P32*y_orb;
    
    double vx_ecl = P11*vx_orb + P12*vy_orb;
    double vy_ecl = P21*vx_orb + P22*vy_orb;
    double vz_ecl = P31*vx_orb + P32*vy_orb;
    
    // Rotazione → ICRF
    double eps = 23.4392911 * DEG2RAD;
    double ce = std::cos(eps), se = std::sin(eps);
    
    pos = {x_ecl, ce * y_ecl - se * z_ecl, se * y_ecl + ce * z_ecl};
    vel = {vx_ecl, ce * vy_ecl - se * vz_ecl, se * vy_ecl + ce * vz_ecl};
}

int main() {
    std::cout << std::fixed << std::setprecision(10);
    
    // Elementi Sierks (da MPC/JPL)
    double a = 3.1754733;
    double e = 0.0454207;
    double i = 2.9046;
    double Omega = 104.16243;
    double omega = 100.5141;
    double M = 229.79088;
    double epoch = 2461000.5;
    double jd_target = 2461008.0913;
    
    std::cout << "=== TEST DIAGNOSTICO SIERKS ===\n\n";
    
    // Calcola stato all'epoca
    Vec3 r0, v0;
    elementsToICRF(a, e, i, Omega, omega, M, r0, v0);
    
    std::cout << "Stato ICRF all'epoca " << epoch << ":\n";
    std::cout << "  r = [" << r0[0] << ", " << r0[1] << ", " << r0[2] << "]\n";
    std::cout << "  v = [" << v0[0] << ", " << v0[1] << ", " << v0[2] << "]\n\n";
    
    // Propaga con moto kepleriano semplice (per diagnostica)
    double dt = jd_target - epoch;
    double n = k / std::pow(a, 1.5);  // moto medio
    double M_new = M * DEG2RAD + n * dt;
    
    // Nuova posizione (trascurando perturbazioni)
    Vec3 r1, v1;
    double M_new_deg = std::fmod(M_new * RAD2DEG, 360.0);
    if (M_new_deg < 0) M_new_deg += 360.0;
    elementsToICRF(a, e, i, Omega, omega, M_new_deg, r1, v1);
    
    std::cout << "Stato ICRF a " << jd_target << " (Kepleriano puro):\n";
    std::cout << "  r = [" << r1[0] << ", " << r1[1] << ", " << r1[2] << "]\n\n";
    
    // Terra
    Vec3 earth = getEarth(jd_target);
    std::cout << "Terra ICRF: [" << earth[0] << ", " << earth[1] << ", " << earth[2] << "]\n\n";
    
    // Vettore geocentrico
    Vec3 geo = {r1[0] - earth[0], r1[1] - earth[1], r1[2] - earth[2]};
    double delta = norm(geo);
    
    // RA/Dec
    double ra_rad = std::atan2(geo[1], geo[0]);
    if (ra_rad < 0) ra_rad += 2*PI;
    double dec_rad = std::asin(geo[2] / delta);
    
    double ra_h = ra_rad * 12.0 / PI;
    int ra_hh = (int)ra_h;
    double ra_m = (ra_h - ra_hh) * 60.0;
    int ra_mm = (int)ra_m;
    double ra_ss = (ra_m - ra_mm) * 60.0;
    
    double dec_d = std::abs(dec_rad) * RAD2DEG;
    char sign = dec_rad >= 0 ? '+' : '-';
    int dec_dd = (int)dec_d;
    double dec_m = (dec_d - dec_dd) * 60.0;
    int dec_mm = (int)dec_m;
    double dec_ss = (dec_m - dec_mm) * 60.0;
    
    std::cout << "POSIZIONE (Kepleriano puro, no perturbazioni):\n";
    printf("  RA  = %02d %02d %06.3f\n", ra_hh, ra_mm, ra_ss);
    printf("  Dec = %c%02d %02d %05.2f\n", sign, dec_dd, dec_mm, dec_ss);
    std::cout << "  Δ   = " << delta << " AU\n\n";
    
    std::cout << "JPL HORIZONS:\n";
    std::cout << "  RA  = 04 53 11.25\n";
    std::cout << "  Dec = +20 19 25.8\n";
    std::cout << "  Δ   = 2.290657 AU\n";
    
    return 0;
}
