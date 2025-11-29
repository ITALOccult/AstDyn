/**
 * @file test_jpl_full.cpp
 * @brief Test propagazione COMPLETA da elementi JPL (6+ anni)
 * 
 * - Stato ICRF diretto da JPL
 * - Perturbazioni planetarie complete (8 pianeti)
 * - Integrazione RKF78
 * - Confronto con JPL Horizons a multiple epoche
 * 
 * Asteroide: (11234) 1999 JS82
 * Epoca elementi: 2019-Jan-26 (JD 2458509.5)
 * Target: 2025-Oct-22 e 2025-Dec-21
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <vector>
#include <fstream>
#include <sstream>

// ============================================================================
// COSTANTI
// ============================================================================

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double AU_KM = 149597870.7;
constexpr double C_LIGHT = 299792.458;  // km/s
constexpr double DAY_SEC = 86400.0;

// GM in AU³/day²
constexpr double GM_SUN     = 2.9591220828559093e-04;
constexpr double GM_MERCURY = 4.9125474514508118e-11;
constexpr double GM_VENUS   = 7.2434524861627027e-10;
constexpr double GM_EARTH   = 8.8876925870231834e-10;  // Terra + Luna
constexpr double GM_MARS    = 9.5495351057792580e-11;
constexpr double GM_JUPITER = 2.8253458420837619e-07;
constexpr double GM_SATURN  = 8.4597151856806587e-08;
constexpr double GM_URANUS  = 1.2920249167819693e-08;
constexpr double GM_NEPTUNE = 1.5243589007842762e-08;

// Stato ICRF iniziale da JPL per (11234) - epoca 2019-Jan-26
constexpr double JPL_EPOCH_JD = 2458509.5;
constexpr double JPL_X  =  2.015534527930346;
constexpr double JPL_Y  =  1.560170291279843;
constexpr double JPL_Z  =  0.07755625121716653;
constexpr double JPL_VX = -0.006439826187731527;
constexpr double JPL_VY =  0.007976810840048847;
constexpr double JPL_VZ =  0.004075596542667446;

using Vec3 = std::array<double, 3>;
using State = std::array<double, 6>;

// ============================================================================
// OPERAZIONI VETTORIALI
// ============================================================================

Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}
Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
Vec3 operator*(double s, const Vec3& v) {
    return {s*v[0], s*v[1], s*v[2]};
}
double dot(const Vec3& a, const Vec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

// ============================================================================
// EFFEMERIDI PLANETARIE (VSOP87 semplificato + termini principali)
// ============================================================================

// Elementi medi J2000.0 per i pianeti (eclittica)
// a, e, i, Omega, omega_bar, L0 (tutti in AU o gradi)
// Rate: da/dt, de/dt, etc. (per secolo)
struct PlanetElements {
    double a0, a1;      // a + da/dt * T
    double e0, e1;
    double i0, i1;
    double Om0, Om1;    // Omega
    double w0, w1;      // omega_bar (longitude of perihelion)
    double L0, L1;      // Mean longitude
};

// Dati da Standish (1992) - JPL
static const PlanetElements PLANETS[] = {
    // Mercury
    {0.38709927, 0.00000037, 0.20563593, 0.00001906, 7.00497902, -0.00594749,
     48.33076593, -0.12534081, 77.45779628, 0.16047689, 252.25032350, 149472.67411175},
    // Venus
    {0.72333566, 0.00000390, 0.00677672, -0.00004107, 3.39467605, -0.00078890,
     76.67984255, -0.27769418, 131.60246718, 0.00268329, 181.97909950, 58517.81538729},
    // Earth-Moon Barycenter
    {1.00000261, 0.00000562, 0.01671123, -0.00004392, -0.00001531, -0.01294668,
     0.0, 0.0, 102.93768193, 0.32327364, 100.46457166, 35999.37244981},
    // Mars
    {1.52371034, 0.00001847, 0.09339410, 0.00007882, 1.84969142, -0.00813131,
     49.55953891, -0.29257343, -23.94362959, 0.44441088, -4.55343205, 19140.30268499},
    // Jupiter
    {5.20288700, -0.00011607, 0.04838624, -0.00013253, 1.30439695, -0.00183714,
     100.47390909, 0.20469106, 14.72847983, 0.21252668, 34.39644051, 3034.74612775},
    // Saturn
    {9.53667594, -0.00125060, 0.05386179, -0.00050991, 2.48599187, 0.00193609,
     113.66242448, -0.28867794, 92.59887831, -0.41897216, 49.95424423, 1222.49362201},
    // Uranus
    {19.18916464, -0.00196176, 0.04725744, -0.00004397, 0.77263783, -0.00242939,
     74.01692503, 0.04240589, 170.95427630, 0.40805281, 313.23810451, 428.48202785},
    // Neptune
    {30.06992276, 0.00026291, 0.00859048, 0.00005105, 1.77004347, 0.00035372,
     131.78422574, -0.00508664, 44.96476227, -0.32241464, -55.12002969, 218.45945325}
};

double solveKepler(double M, double e) {
    double E = M;
    for (int i = 0; i < 30; ++i) {
        double dE = (M - E + e * std::sin(E)) / (1.0 - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < 1e-14) break;
    }
    return E;
}

// Posizione pianeta in ECLITTICA J2000, poi converti a ICRF
Vec3 getPlanetPositionICRF(int planet, double jd) {
    const double J2000 = 2451545.0;
    double T = (jd - J2000) / 36525.0;  // Secoli da J2000
    
    const auto& p = PLANETS[planet];
    
    // Elementi osculanti
    double a = p.a0 + p.a1 * T;
    double e = p.e0 + p.e1 * T;
    double i = (p.i0 + p.i1 * T) * DEG2RAD;
    double Om = (p.Om0 + p.Om1 * T) * DEG2RAD;
    double w_bar = (p.w0 + p.w1 * T) * DEG2RAD;
    double L = (p.L0 + p.L1 * T) * DEG2RAD;
    
    double omega = w_bar - Om;  // Argomento del perielio
    double M = L - w_bar;       // Anomalia media
    
    // Normalizza M
    while (M < 0) M += 2*PI;
    while (M > 2*PI) M -= 2*PI;
    
    // Risolvi Keplero
    double E = solveKepler(M, e);
    
    // Posizione nel piano orbitale
    double cosE = std::cos(E);
    double sinE = std::sin(E);
    double x_orb = a * (cosE - e);
    double y_orb = a * std::sqrt(1.0 - e*e) * sinE;
    
    // Rotazione: piano orbitale -> eclittica
    double cosO = std::cos(Om);
    double sinO = std::sin(Om);
    double cosw = std::cos(omega);
    double sinw = std::sin(omega);
    double cosi = std::cos(i);
    double sini = std::sin(i);
    
    double P1 = cosO * cosw - sinO * sinw * cosi;
    double P2 = -cosO * sinw - sinO * cosw * cosi;
    double Q1 = sinO * cosw + cosO * sinw * cosi;
    double Q2 = -sinO * sinw + cosO * cosw * cosi;
    double W1 = sinw * sini;
    double W2 = cosw * sini;
    
    double x_ecl = P1 * x_orb + P2 * y_orb;
    double y_ecl = Q1 * x_orb + Q2 * y_orb;
    double z_ecl = W1 * x_orb + W2 * y_orb;
    
    // Rotazione: eclittica -> equatoriale (ICRF)
    constexpr double eps = 23.4392911 * DEG2RAD;  // Obliquità J2000
    double c = std::cos(eps);
    double s = std::sin(eps);
    
    return {x_ecl,
            c * y_ecl - s * z_ecl,
            s * y_ecl + c * z_ecl};
}

// ============================================================================
// ACCELERAZIONE CON PERTURBAZIONI
// ============================================================================

State computeAcceleration(const State& state, double jd) {
    Vec3 r = {state[0], state[1], state[2]};
    Vec3 v = {state[3], state[4], state[5]};
    
    double r_mag = norm(r);
    double r3 = r_mag * r_mag * r_mag;
    
    // Accelerazione centrale (Sole)
    Vec3 acc = {-GM_SUN * r[0] / r3,
                -GM_SUN * r[1] / r3,
                -GM_SUN * r[2] / r3};
    
    // Perturbazioni planetarie
    static const double GM[] = {GM_MERCURY, GM_VENUS, GM_EARTH, GM_MARS,
                                 GM_JUPITER, GM_SATURN, GM_URANUS, GM_NEPTUNE};
    
    for (int i = 0; i < 8; ++i) {
        Vec3 rp = getPlanetPositionICRF(i, jd);
        Vec3 d = r - rp;  // Asteroide - pianeta
        
        double d_mag = norm(d);
        double rp_mag = norm(rp);
        
        double d3 = d_mag * d_mag * d_mag;
        double rp3 = rp_mag * rp_mag * rp_mag;
        
        // Termine diretto + indiretto
        acc[0] -= GM[i] * (d[0]/d3 + rp[0]/rp3);
        acc[1] -= GM[i] * (d[1]/d3 + rp[1]/rp3);
        acc[2] -= GM[i] * (d[2]/d3 + rp[2]/rp3);
    }
    
    return {v[0], v[1], v[2], acc[0], acc[1], acc[2]};
}

// ============================================================================
// INTEGRATORE RKF78
// ============================================================================

class RKF78Integrator {
private:
    // Coefficienti RKF78 (Fehlberg 1968)
    static constexpr double c[13] = {
        0.0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 1.0/2.0,
        5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0, 1.0
    };
    
    // Matrice dei coefficienti a[i][j] per RKF78
    static constexpr double a[13][12] = {
        {0},
        {2.0/27.0},
        {1.0/36.0, 1.0/12.0},
        {1.0/24.0, 0, 1.0/8.0},
        {5.0/12.0, 0, -25.0/16.0, 25.0/16.0},
        {1.0/20.0, 0, 0, 1.0/4.0, 1.0/5.0},
        {-25.0/108.0, 0, 0, 125.0/108.0, -65.0/27.0, 125.0/54.0},
        {31.0/300.0, 0, 0, 0, 61.0/225.0, -2.0/9.0, 13.0/900.0},
        {2.0, 0, 0, -53.0/6.0, 704.0/45.0, -107.0/9.0, 67.0/90.0, 3.0},
        {-91.0/108.0, 0, 0, 23.0/108.0, -976.0/135.0, 311.0/54.0, -19.0/60.0, 17.0/6.0, -1.0/12.0},
        {2383.0/4100.0, 0, 0, -341.0/164.0, 4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0},
        {3.0/205.0, 0, 0, 0, 0, -6.0/41.0, -3.0/205.0, -3.0/41.0, 3.0/41.0, 6.0/41.0, 0},
        {-1777.0/4100.0, 0, 0, -341.0/164.0, 4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0, 0, 1.0}
    };
    
    // Pesi ordine 8
    static constexpr double b8[13] = {
        41.0/840.0, 0, 0, 0, 0, 34.0/105.0, 9.0/35.0, 9.0/35.0,
        9.0/280.0, 9.0/280.0, 41.0/840.0, 0, 0
    };
    
    // Pesi ordine 7 (per stima errore)
    static constexpr double b7[13] = {
        0, 0, 0, 0, 0, 34.0/105.0, 9.0/35.0, 9.0/35.0,
        9.0/280.0, 9.0/280.0, 0, 41.0/840.0, 41.0/840.0
    };
    
    double tol_;
    double h_min_, h_max_;
    
public:
    RKF78Integrator(double tol = 1e-12, double h_min = 0.001, double h_max = 10.0)
        : tol_(tol), h_min_(h_min), h_max_(h_max) {}
    
    State step(const State& y, double t, double& h, double& err) {
        std::array<State, 13> k;
        
        // Calcola tutti gli stage
        for (int stage = 0; stage < 13; ++stage) {
            State y_stage = y;
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < stage; ++j) {
                    y_stage[i] += h * a[stage][j] * k[j][i];
                }
            }
            k[stage] = computeAcceleration(y_stage, t + c[stage] * h);
        }
        
        // Soluzione ordine 8
        State y8 = y;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 13; ++j) {
                y8[i] += h * b8[j] * k[j][i];
            }
        }
        
        // Soluzione ordine 7
        State y7 = y;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 13; ++j) {
                y7[i] += h * b7[j] * k[j][i];
            }
        }
        
        // Stima errore
        err = 0.0;
        for (int i = 0; i < 6; ++i) {
            double e = std::abs(y8[i] - y7[i]);
            double scale = std::abs(y8[i]) + std::abs(h * k[0][i]) + 1e-10;
            err = std::max(err, e / scale);
        }
        
        return y8;
    }
    
    State integrate(const State& y0, double t0, double t_end, int& total_steps) {
        State y = y0;
        double t = t0;
        double h = (t_end - t0) / 100.0;  // Passo iniziale
        
        if (h > h_max_) h = h_max_;
        if (std::abs(h) < h_min_) h = (h > 0) ? h_min_ : -h_min_;
        
        total_steps = 0;
        
        while ((h > 0 && t < t_end) || (h < 0 && t > t_end)) {
            // Non superare t_end
            if ((h > 0 && t + h > t_end) || (h < 0 && t + h < t_end)) {
                h = t_end - t;
            }
            
            double err;
            State y_new = step(y, t, h, err);
            
            if (err < tol_ || std::abs(h) <= h_min_) {
                // Accetta il passo
                y = y_new;
                t += h;
                total_steps++;
                
                // Adatta passo
                if (err > 0) {
                    double factor = 0.9 * std::pow(tol_ / err, 1.0/8.0);
                    factor = std::max(0.1, std::min(4.0, factor));
                    h *= factor;
                }
            } else {
                // Rifiuta e riduci
                h *= 0.5;
            }
            
            // Limita h
            if (std::abs(h) > h_max_) h = (h > 0) ? h_max_ : -h_max_;
            if (std::abs(h) < h_min_) h = (h > 0) ? h_min_ : -h_min_;
        }
        
        return y;
    }
};

// ============================================================================
// CONVERSIONE A RA/DEC
// ============================================================================

void stateToRADec(const State& s, double& ra_h, double& ra_m, double& ra_s,
                  int& dec_d, int& dec_m, double& dec_s) {
    // Stato già in ICRF (equatoriale)
    double ra = std::atan2(s[1], s[0]);
    if (ra < 0) ra += 2*PI;
    double dec = std::atan2(s[2], std::sqrt(s[0]*s[0] + s[1]*s[1]));
    
    // RA in ore
    double ra_deg = ra * RAD2DEG;
    double ra_hours = ra_deg / 15.0;
    ra_h = std::floor(ra_hours);
    ra_m = std::floor((ra_hours - ra_h) * 60.0);
    ra_s = ((ra_hours - ra_h) * 60.0 - ra_m) * 60.0;
    
    // Dec in gradi
    double dec_deg = dec * RAD2DEG;
    int sign = (dec_deg >= 0) ? 1 : -1;
    dec_deg = std::abs(dec_deg);
    dec_d = sign * static_cast<int>(std::floor(dec_deg));
    dec_m = static_cast<int>(std::floor((dec_deg - std::abs(dec_d)) * 60.0));
    dec_s = ((dec_deg - std::abs(dec_d)) * 60.0 - dec_m) * 60.0;
}

double computeDeltaArcsec(double ra1_h, double ra1_m, double ra1_s,
                          int dec1_d, int dec1_m, double dec1_s,
                          double ra2_h, double ra2_m, double ra2_s,
                          int dec2_d, int dec2_m, double dec2_s) {
    double ra1 = (ra1_h + ra1_m/60.0 + ra1_s/3600.0) * 15.0;  // gradi
    double ra2 = (ra2_h + ra2_m/60.0 + ra2_s/3600.0) * 15.0;
    
    double dec1 = std::abs(dec1_d) + dec1_m/60.0 + dec1_s/3600.0;
    if (dec1_d < 0) dec1 = -dec1;
    double dec2 = std::abs(dec2_d) + dec2_m/60.0 + dec2_s/3600.0;
    if (dec2_d < 0) dec2 = -dec2;
    
    double dra = (ra1 - ra2) * 3600.0;   // arcsec
    double ddec = (dec1 - dec2) * 3600.0;
    
    // Correzione coseno declinazione per RA
    double cos_dec = std::cos((dec1 + dec2) / 2.0 * DEG2RAD);
    dra *= cos_dec;
    
    return std::sqrt(dra*dra + ddec*ddec);
}

// ============================================================================
// MAIN
// ============================================================================

int main() {
    std::cout << std::fixed;
    
    std::cout << "================================================================\n";
    std::cout << "  Test Propagazione COMPLETA (6.9 anni)\n";
    std::cout << "  Asteroide (11234) 1999 JS82\n";
    std::cout << "  Integratore: RKF78 con perturbazioni planetarie\n";
    std::cout << "================================================================\n\n";
    
    // Stato iniziale ICRF da JPL
    State s0 = {JPL_X, JPL_Y, JPL_Z, JPL_VX, JPL_VY, JPL_VZ};
    
    std::cout << "1. STATO INIZIALE (JD " << std::setprecision(1) << JPL_EPOCH_JD 
              << " = 2019-Jan-26)\n";
    std::cout << std::setprecision(12);
    std::cout << "   r = [" << s0[0] << ", " << s0[1] << ", " << s0[2] << "] AU\n";
    std::cout << "   v = [" << s0[3] << ", " << s0[4] << ", " << s0[5] << "] AU/day\n";
    std::cout << "   |r| = " << norm({s0[0], s0[1], s0[2]}) << " AU\n\n";
    
    // Target dates
    struct Target {
        double jd;
        const char* date;
        double jpl_ra_h, jpl_ra_m, jpl_ra_s;
        int jpl_dec_d, jpl_dec_m;
        double jpl_dec_s;
        double jpl_delta;
    };
    
    std::vector<Target> targets = {
        {2458540.5, "2019-Feb-26", 2, 59, 36.94, 4, 33, 53.2, 2.5484},   // +31 days
        {2460970.5, "2025-Oct-22", 15, 18, 51.75, -6, 25, 31.9, 2.8112}, // +2461 days
        {2461030.5, "2025-Dec-21", 16, 5, 27.35, -10, 38, 8.1, 2.8046}   // +2521 days
    };
    
    RKF78Integrator integrator(1e-12, 0.01, 5.0);
    
    std::cout << "2. PROPAGAZIONE E CONFRONTO\n\n";
    
    for (const auto& t : targets) {
        double dt = t.jd - JPL_EPOCH_JD;
        
        std::cout << "================================================================\n";
        std::cout << "  " << t.date << " (JD " << std::setprecision(1) << t.jd 
                  << ", Δt = " << std::setprecision(0) << dt << " giorni)\n";
        std::cout << "================================================================\n";
        
        int steps;
        State s_prop = integrator.integrate(s0, JPL_EPOCH_JD, t.jd, steps);
        
        std::cout << std::setprecision(9);
        std::cout << "   r = [" << s_prop[0] << ", " << s_prop[1] << ", " << s_prop[2] << "] AU\n";
        std::cout << "   |r| = " << norm({s_prop[0], s_prop[1], s_prop[2]}) << " AU\n";
        std::cout << "   Passi RKF78: " << steps << "\n";
        
        double ra_h, ra_m, ra_s;
        int dec_d, dec_m;
        double dec_s;
        stateToRADec(s_prop, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s);
        
        std::cout << "\n   AstDyn:\n";
        std::cout << std::setprecision(2);
        std::cout << "   RA  = " << std::setw(2) << (int)ra_h << " " 
                  << std::setw(2) << (int)ra_m << " " 
                  << std::setprecision(2) << std::setw(5) << ra_s << "\n";
        std::cout << "   Dec = " << std::showpos << std::setw(3) << dec_d << std::noshowpos
                  << " " << std::setw(2) << dec_m << " " 
                  << std::setprecision(1) << std::setw(4) << dec_s << "\n";
        
        std::cout << "\n   JPL Horizons:\n";
        std::cout << std::setprecision(2);
        std::cout << "   RA  = " << std::setw(2) << (int)t.jpl_ra_h << " " 
                  << std::setw(2) << (int)t.jpl_ra_m << " " 
                  << std::setprecision(2) << std::setw(5) << t.jpl_ra_s << "\n";
        std::cout << "   Dec = " << std::showpos << std::setw(3) << t.jpl_dec_d << std::noshowpos
                  << " " << std::setw(2) << t.jpl_dec_m << " " 
                  << std::setprecision(1) << std::setw(4) << t.jpl_dec_s << "\n";
        
        double total_err = computeDeltaArcsec(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s,
                                               t.jpl_ra_h, t.jpl_ra_m, t.jpl_ra_s,
                                               t.jpl_dec_d, t.jpl_dec_m, t.jpl_dec_s);
        
        std::cout << std::setprecision(1);
        std::cout << "\n   ERRORE TOTALE: " << total_err << " arcsec";
        if (total_err < 10) {
            std::cout << " ✓✓✓ ECCELLENTE\n";
        } else if (total_err < 60) {
            std::cout << " ✓✓ BUONO (< 1')\n";
        } else if (total_err < 300) {
            std::cout << " ✓ ACCETTABILE (< 5')\n";
        } else {
            std::cout << " ✗ DA MIGLIORARE\n";
        }
        std::cout << std::endl;
    }
    
    // Test round-trip
    std::cout << "================================================================\n";
    std::cout << "  TEST ROUND-TRIP (precisione integratore)\n";
    std::cout << "================================================================\n";
    
    double jd_mid = 2460000.0;  // ~2023
    int steps1, steps2;
    State s_mid = integrator.integrate(s0, JPL_EPOCH_JD, jd_mid, steps1);
    State s_back = integrator.integrate(s_mid, jd_mid, JPL_EPOCH_JD, steps2);
    
    double dr = norm({s_back[0]-s0[0], s_back[1]-s0[1], s_back[2]-s0[2]});
    double dv = norm({s_back[3]-s0[3], s_back[4]-s0[4], s_back[5]-s0[5]});
    
    std::cout << std::setprecision(3);
    std::cout << "   Andata: " << steps1 << " passi, Ritorno: " << steps2 << " passi\n";
    std::cout << std::setprecision(6);
    std::cout << "   Errore posizione: " << dr << " AU = " << dr * AU_KM * 1000.0 << " m\n";
    std::cout << "   Errore velocità:  " << dv << " AU/day\n";
    
    std::cout << "\n================================================================\n";
    std::cout << "  RIEPILOGO\n";
    std::cout << "================================================================\n";
    std::cout << "  Asteroide: (11234) 1999 JS82\n";
    std::cout << "  Epoca elementi: 2019-Jan-26 (JPL)\n";
    std::cout << "  Intervallo massimo: 6.9 anni\n";
    std::cout << "  Perturbazioni: 8 pianeti (Mercurio-Nettuno)\n";
    std::cout << "  Integratore: RKF78 (ordine 7/8 adattivo)\n";
    std::cout << "  Tolleranza: 1e-12\n";
    std::cout << "================================================================\n";
    
    return 0;
}
