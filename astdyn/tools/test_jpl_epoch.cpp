/**
 * @file test_jpl_epoch.cpp
 * @brief Test propagazione da elementi JPL con confronto Horizons
 * 
 * Elementi JPL per (11234) 1999 JS82:
 * EPOCH = 2458509.5 (2019-Jan-26.00 TDB)
 * EC = 0.0494518850848228
 * QR = 2.548364809395927 AU (perihelion)
 * TP = 2458547.7031539734 (time of perihelion)
 * OM = 112.6647729938423° (ascending node)
 * W  = 292.4561825062523° (argument of perihelion)
 * IN = 12.77459773563021° (inclination)
 * A  = 2.6809424682551 AU
 * MA = 351.4222763330189° (mean anomaly at epoch)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <vector>

// Costanti
constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double AU_KM = 149597870.7;
constexpr double GM_SUN = 2.9591220828559093e-04;  // AU³/day²

// Elementi orbitali JPL per (11234)
constexpr double JPL_EPOCH_JD = 2458509.5;  // 2019-Jan-26.00 TDB
constexpr double JPL_A = 2.6809424682551;    // AU
constexpr double JPL_E = 0.0494518850848228;
constexpr double JPL_I = 12.77459773563021;  // deg
constexpr double JPL_OM = 112.6647729938423; // Omega (node) deg
constexpr double JPL_W = 292.4561825062523;  // omega (peri) deg
constexpr double JPL_MA = 351.4222763330189; // Mean anomaly deg

// GM pianeti (AU³/day²)
constexpr double GM_MERCURY = 4.9125474514508118e-11;
constexpr double GM_VENUS   = 7.2434524861627027e-10;
constexpr double GM_EARTH   = 8.8876925870231834e-10;
constexpr double GM_MARS    = 9.5495351057792580e-11;
constexpr double GM_JUPITER = 2.8253458420837619e-07;
constexpr double GM_SATURN  = 8.4597151856806587e-08;
constexpr double GM_URANUS  = 1.2920249167819693e-08;
constexpr double GM_NEPTUNE = 1.5243589007842762e-08;

using Vec3 = std::array<double, 3>;
using State = std::array<double, 6>;

// Operazioni vettoriali
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

// Risoluzione equazione di Keplero
double solveKepler(double M, double e, double tol = 1e-14) {
    double E = M;
    for (int i = 0; i < 50; ++i) {
        double dE = (M - E + e * std::sin(E)) / (1.0 - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < tol) break;
    }
    return E;
}

// Elementi kepleriani -> stato cartesiano (eclittica J2000)
State keplerToState(double a, double e, double i, double Omega, double omega, double M) {
    // Risolvi Keplero
    double E = solveKepler(M, e);
    
    // Posizione nel piano orbitale
    double cosE = std::cos(E);
    double sinE = std::sin(E);
    double sqrtae = std::sqrt(1.0 - e*e);
    
    double x_orb = a * (cosE - e);
    double y_orb = a * sqrtae * sinE;
    
    // Velocità nel piano orbitale
    double n = std::sqrt(GM_SUN / (a*a*a));  // mean motion
    double r = a * (1.0 - e * cosE);
    double vx_orb = -a * n * sinE / (1.0 - e * cosE);
    double vy_orb = a * n * sqrtae * cosE / (1.0 - e * cosE);
    
    // Rotazione al frame eclittico
    double cosO = std::cos(Omega);
    double sinO = std::sin(Omega);
    double cosw = std::cos(omega);
    double sinw = std::sin(omega);
    double cosi = std::cos(i);
    double sini = std::sin(i);
    
    // Matrice di rotazione
    double P1 = cosO * cosw - sinO * sinw * cosi;
    double P2 = -cosO * sinw - sinO * cosw * cosi;
    double Q1 = sinO * cosw + cosO * sinw * cosi;
    double Q2 = -sinO * sinw + cosO * cosw * cosi;
    double W1 = sinw * sini;
    double W2 = cosw * sini;
    
    State s;
    s[0] = P1 * x_orb + P2 * y_orb;
    s[1] = Q1 * x_orb + Q2 * y_orb;
    s[2] = W1 * x_orb + W2 * y_orb;
    s[3] = P1 * vx_orb + P2 * vy_orb;
    s[4] = Q1 * vx_orb + Q2 * vy_orb;
    s[5] = W1 * vx_orb + W2 * vy_orb;
    
    return s;
}

// Posizione pianeta (semplificata ma accurata per test)
Vec3 getPlanetPosition(int planet, double jd) {
    // Elementi medi J2000 per i pianeti
    struct PlanetElements {
        double a, e, i, Om, w, L0, n;  // n = mean motion (deg/day)
    };
    
    static const PlanetElements elements[] = {
        {0.38709927, 0.20563593, 7.00497902, 48.33076593, 77.45779628, 252.25032350, 4.09233445},   // Mercury
        {0.72333566, 0.00677672, 3.39467605, 76.67984255, 131.60246718, 181.97909950, 1.60213034},  // Venus
        {1.00000261, 0.01671123, -0.00001531, 0.0, 102.93768193, 100.46457166, 0.98560028},         // Earth
        {1.52371034, 0.09339410, 1.84969142, 49.55953891, -23.94362959, -4.55343205, 0.52402068},   // Mars
        {5.20288700, 0.04838624, 1.30439695, 100.47390909, 14.72847983, 34.39644051, 0.08308529},   // Jupiter
        {9.53667594, 0.05386179, 2.48599187, 113.66242448, 92.59887831, 49.95424423, 0.03349791},   // Saturn
        {19.18916464, 0.04725744, 0.77263783, 74.01692503, 170.95427630, 313.23810451, 0.01176904}, // Uranus
        {30.06992276, 0.00859048, 1.77004347, 131.78422574, 44.96476227, -55.12002969, 0.00606020}  // Neptune
    };
    
    const double J2000 = 2451545.0;
    double dt = jd - J2000;
    
    const auto& p = elements[planet];
    double L = (p.L0 + p.n * dt) * DEG2RAD;
    double M = L - p.w * DEG2RAD;
    
    // Normalizza M
    while (M < 0) M += 2*PI;
    while (M > 2*PI) M -= 2*PI;
    
    double E = solveKepler(M, p.e);
    
    double cosE = std::cos(E);
    double sinE = std::sin(E);
    double x_orb = p.a * (cosE - p.e);
    double y_orb = p.a * std::sqrt(1 - p.e*p.e) * sinE;
    
    double Om = p.Om * DEG2RAD;
    double w = p.w * DEG2RAD;
    double inc = p.i * DEG2RAD;
    
    double cosO = std::cos(Om);
    double sinO = std::sin(Om);
    double cosw = std::cos(w);
    double sinw = std::sin(w);
    double cosi = std::cos(inc);
    double sini = std::sin(inc);
    
    double P1 = cosO * cosw - sinO * sinw * cosi;
    double P2 = -cosO * sinw - sinO * cosw * cosi;
    double Q1 = sinO * cosw + cosO * sinw * cosi;
    double Q2 = -sinO * sinw + cosO * cosw * cosi;
    double W1 = sinw * sini;
    double W2 = cosw * sini;
    
    return {P1 * x_orb + P2 * y_orb,
            Q1 * x_orb + Q2 * y_orb,
            W1 * x_orb + W2 * y_orb};
}

// Accelerazione gravitazionale
State computeAcceleration(const State& s, double jd) {
    Vec3 r = {s[0], s[1], s[2]};
    double r_mag = norm(r);
    double r3 = r_mag * r_mag * r_mag;
    
    // Accelerazione centrale (Sole)
    Vec3 acc = {-GM_SUN * r[0] / r3,
                -GM_SUN * r[1] / r3,
                -GM_SUN * r[2] / r3};
    
    // Perturbazioni planetarie
    double gm[] = {GM_MERCURY, GM_VENUS, GM_EARTH, GM_MARS, 
                   GM_JUPITER, GM_SATURN, GM_URANUS, GM_NEPTUNE};
    
    for (int i = 0; i < 8; ++i) {
        Vec3 rp = getPlanetPosition(i, jd);
        Vec3 d = r - rp;  // asteroide - pianeta
        double d_mag = norm(d);
        double rp_mag = norm(rp);
        
        double d3 = d_mag * d_mag * d_mag;
        double rp3 = rp_mag * rp_mag * rp_mag;
        
        // Termine diretto + indiretto
        acc[0] -= gm[i] * (d[0]/d3 + rp[0]/rp3);
        acc[1] -= gm[i] * (d[1]/d3 + rp[1]/rp3);
        acc[2] -= gm[i] * (d[2]/d3 + rp[2]/rp3);
    }
    
    return {s[3], s[4], s[5], acc[0], acc[1], acc[2]};
}

// Integratore RKF78
class RKF78 {
public:
    static constexpr int nStages = 13;
    
    // Coefficienti RKF78 (Fehlberg)
    static constexpr double c[13] = {
        0.0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 1.0/2.0, 
        5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0, 1.0
    };
    
    // Pesi ordine 8
    static constexpr double b8[13] = {
        41.0/840.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0, 
        9.0/35.0, 9.0/280.0, 9.0/280.0, 41.0/840.0, 0.0, 0.0
    };
    
    // Pesi ordine 7 (per stima errore)
    static constexpr double b7[13] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0,
        9.0/35.0, 9.0/280.0, 9.0/280.0, 0.0, 41.0/840.0, 41.0/840.0
    };
    
    State step(const State& y, double t, double h) {
        std::array<State, 13> k;
        
        // Stage 1
        k[0] = computeAcceleration(y, t);
        
        // Stages 2-13 (semplificato per test)
        for (int stage = 1; stage < 13; ++stage) {
            State y_temp = y;
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < stage; ++j) {
                    y_temp[i] += h * 0.1 * k[j][i]; // Coefficienti semplificati
                }
            }
            k[stage] = computeAcceleration(y_temp, t + c[stage] * h);
        }
        
        // Combina con pesi ordine 8
        State y_new = y;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 13; ++j) {
                y_new[i] += h * b8[j] * k[j][i];
            }
        }
        
        return y_new;
    }
};

// Propagazione con RK4 (più semplice e affidabile per test)
State propagateRK4(const State& y0, double jd0, double jd_target, int steps) {
    double h = (jd_target - jd0) / steps;
    State y = y0;
    double jd = jd0;
    
    for (int i = 0; i < steps; ++i) {
        State k1 = computeAcceleration(y, jd);
        
        State y2;
        for (int j = 0; j < 6; ++j) y2[j] = y[j] + 0.5*h*k1[j];
        State k2 = computeAcceleration(y2, jd + 0.5*h);
        
        State y3;
        for (int j = 0; j < 6; ++j) y3[j] = y[j] + 0.5*h*k2[j];
        State k3 = computeAcceleration(y3, jd + 0.5*h);
        
        State y4;
        for (int j = 0; j < 6; ++j) y4[j] = y[j] + h*k3[j];
        State k4 = computeAcceleration(y4, jd + h);
        
        for (int j = 0; j < 6; ++j) {
            y[j] += h * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]) / 6.0;
        }
        jd += h;
    }
    
    return y;
}

// Calcola RA/Dec eliocentrica
void computeRADec(const State& s, double& ra_h, double& ra_m, double& ra_s,
                  double& dec_d, double& dec_m, double& dec_s) {
    // Converti eclittica -> equatoriale (obliquità J2000)
    constexpr double eps = 23.4392911 * DEG2RAD;
    double x_eq = s[0];
    double y_eq = s[1] * std::cos(eps) - s[2] * std::sin(eps);
    double z_eq = s[1] * std::sin(eps) + s[2] * std::cos(eps);
    
    double ra = std::atan2(y_eq, x_eq);
    if (ra < 0) ra += 2*PI;
    double dec = std::atan2(z_eq, std::sqrt(x_eq*x_eq + y_eq*y_eq));
    
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
    dec_d = sign * std::floor(dec_deg);
    dec_m = std::floor((dec_deg - std::abs(dec_d)) * 60.0);
    dec_s = ((dec_deg - std::abs(dec_d)) * 60.0 - dec_m) * 60.0;
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "================================================================\n";
    std::cout << "  Test Propagazione da Elementi JPL\n";
    std::cout << "  Asteroide (11234) 1999 JS82\n";
    std::cout << "================================================================\n\n";
    
    // Elementi JPL
    std::cout << "1. ELEMENTI ORBITALI JPL\n";
    std::cout << "   Epoca: JD " << std::setprecision(1) << JPL_EPOCH_JD 
              << " (2019-Jan-26.00 TDB)\n";
    std::cout << std::setprecision(10);
    std::cout << "   a  = " << JPL_A << " AU\n";
    std::cout << "   e  = " << JPL_E << "\n";
    std::cout << "   i  = " << JPL_I << "°\n";
    std::cout << "   Ω  = " << JPL_OM << "°\n";
    std::cout << "   ω  = " << JPL_W << "°\n";
    std::cout << "   M  = " << JPL_MA << "°\n\n";
    
    // Converti a radianti e calcola stato iniziale
    State s0 = keplerToState(JPL_A, JPL_E, 
                              JPL_I * DEG2RAD, 
                              JPL_OM * DEG2RAD, 
                              JPL_W * DEG2RAD, 
                              JPL_MA * DEG2RAD);
    
    std::cout << "2. STATO ALL'EPOCA JPL\n";
    std::cout << std::setprecision(9);
    std::cout << "   r = [" << s0[0] << ", " << s0[1] << ", " << s0[2] << "] AU\n";
    std::cout << "   v = [" << s0[3] << ", " << s0[4] << ", " << s0[5] << "] AU/day\n";
    std::cout << "   |r| = " << norm({s0[0], s0[1], s0[2]}) << " AU\n\n";
    
    // Epoche target
    double jd_oct22 = 2460970.5;  // 2025-Oct-22
    double jd_dec21 = 2461030.5;  // 2025-Dec-21
    
    double dt_oct = jd_oct22 - JPL_EPOCH_JD;
    double dt_dec = jd_dec21 - JPL_EPOCH_JD;
    
    std::cout << "3. PROPAGAZIONE\n";
    std::cout << "   Epoca JPL: JD " << std::setprecision(1) << JPL_EPOCH_JD << "\n";
    std::cout << "   Target 1: JD " << jd_oct22 << " (2025-Oct-22) Δt = " 
              << std::setprecision(0) << dt_oct << " giorni\n";
    std::cout << "   Target 2: JD " << jd_dec21 << " (2025-Dec-21) Δt = " 
              << dt_dec << " giorni\n\n";
    
    // Propagazione a 2025-Oct-22
    int steps_per_day = 10;
    int steps_oct = static_cast<int>(std::abs(dt_oct)) * steps_per_day;
    
    std::cout << "================================================================\n";
    std::cout << "  PROPAGAZIONE A 2025-Oct-22 (Δt = " << dt_oct << " giorni)\n";
    std::cout << "================================================================\n";
    
    State s_oct = propagateRK4(s0, JPL_EPOCH_JD, jd_oct22, steps_oct);
    
    std::cout << std::setprecision(9);
    std::cout << "   r = [" << s_oct[0] << ", " << s_oct[1] << ", " << s_oct[2] << "] AU\n";
    std::cout << "   |r| = " << norm({s_oct[0], s_oct[1], s_oct[2]}) << " AU\n";
    
    double ra_h, ra_m, ra_s, dec_d, dec_m, dec_s;
    computeRADec(s_oct, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s);
    
    std::cout << "\n   POSIZIONE ELIOCENTRICA:\n";
    std::cout << std::setprecision(2);
    std::cout << "   RA  = " << std::setw(2) << (int)ra_h << " " 
              << std::setw(2) << (int)ra_m << " " 
              << std::setprecision(3) << std::setw(6) << ra_s << "\n";
    std::cout << "   Dec = " << std::showpos << std::setw(3) << (int)dec_d << std::noshowpos
              << " " << std::setw(2) << (int)dec_m << " " 
              << std::setprecision(2) << std::setw(5) << dec_s << "\n";
    
    std::cout << "\n   JPL HORIZONS (2025-Oct-22, eliocentrico):\n";
    std::cout << "   RA  = 15 18 51.75\n";
    std::cout << "   Dec = -06 25 31.9\n";
    std::cout << "   Δ   = 2.8112 AU\n";
    
    // Calcola differenze
    double jpl_ra_oct = 15.0 + 18.0/60.0 + 51.75/3600.0;
    double astdyn_ra_oct = ra_h + ra_m/60.0 + ra_s/3600.0;
    double delta_ra_oct = (astdyn_ra_oct - jpl_ra_oct) * 15.0 * 3600.0; // arcsec
    
    double jpl_dec_oct = -(6.0 + 25.0/60.0 + 31.9/3600.0);
    double astdyn_dec_oct = dec_d + (dec_d >= 0 ? 1 : -1) * (dec_m/60.0 + dec_s/3600.0);
    double delta_dec_oct = (astdyn_dec_oct - jpl_dec_oct) * 3600.0; // arcsec
    
    std::cout << std::setprecision(1);
    std::cout << "\n   DIFFERENZA:\n";
    std::cout << "   ΔRA  = " << delta_ra_oct << " arcsec\n";
    std::cout << "   ΔDec = " << delta_dec_oct << " arcsec\n";
    
    // Propagazione a 2025-Dec-21
    int steps_dec = static_cast<int>(std::abs(dt_dec)) * steps_per_day;
    
    std::cout << "\n================================================================\n";
    std::cout << "  PROPAGAZIONE A 2025-Dec-21 (Δt = " << dt_dec << " giorni)\n";
    std::cout << "================================================================\n";
    
    State s_dec = propagateRK4(s0, JPL_EPOCH_JD, jd_dec21, steps_dec);
    
    std::cout << std::setprecision(9);
    std::cout << "   r = [" << s_dec[0] << ", " << s_dec[1] << ", " << s_dec[2] << "] AU\n";
    std::cout << "   |r| = " << norm({s_dec[0], s_dec[1], s_dec[2]}) << " AU\n";
    
    computeRADec(s_dec, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s);
    
    std::cout << "\n   POSIZIONE ELIOCENTRICA:\n";
    std::cout << std::setprecision(2);
    std::cout << "   RA  = " << std::setw(2) << (int)ra_h << " " 
              << std::setw(2) << (int)ra_m << " " 
              << std::setprecision(3) << std::setw(6) << ra_s << "\n";
    std::cout << "   Dec = " << std::showpos << std::setw(3) << (int)dec_d << std::noshowpos
              << " " << std::setw(2) << " " << std::setw(2) << (int)dec_m << " " 
              << std::setprecision(2) << std::setw(5) << dec_s << "\n";
    
    std::cout << "\n   JPL HORIZONS (2025-Dec-21, eliocentrico):\n";
    std::cout << "   RA  = 16 05 27.35\n";
    std::cout << "   Dec = -10 38 08.1\n";
    std::cout << "   Δ   = 2.8046 AU\n";
    
    double jpl_ra_dec = 16.0 + 5.0/60.0 + 27.35/3600.0;
    double astdyn_ra_dec = ra_h + ra_m/60.0 + ra_s/3600.0;
    double delta_ra_dec = (astdyn_ra_dec - jpl_ra_dec) * 15.0 * 3600.0;
    
    double jpl_dec_dec = -(10.0 + 38.0/60.0 + 8.1/3600.0);
    double astdyn_dec_dec = dec_d + (dec_d >= 0 ? 1 : -1) * (dec_m/60.0 + dec_s/3600.0);
    double delta_dec_dec = (astdyn_dec_dec - jpl_dec_dec) * 3600.0;
    
    std::cout << std::setprecision(1);
    std::cout << "\n   DIFFERENZA:\n";
    std::cout << "   ΔRA  = " << delta_ra_dec << " arcsec\n";
    std::cout << "   ΔDec = " << delta_dec_dec << " arcsec\n";
    
    std::cout << "\n================================================================\n";
    std::cout << "  RIEPILOGO CONFRONTO\n";
    std::cout << "================================================================\n";
    std::cout << "  Propagazione da elementi JPL (epoca 2019-Jan-26)\n";
    std::cout << "  Intervallo: ~" << std::setprecision(1) << dt_dec << " giorni (~" 
              << dt_dec/365.25 << " anni)\n\n";
    
    std::cout << "  2025-Oct-22: ΔRA = " << delta_ra_oct << "\", ΔDec = " << delta_dec_oct << "\"\n";
    std::cout << "  2025-Dec-21: ΔRA = " << delta_ra_dec << "\", ΔDec = " << delta_dec_dec << "\"\n";
    std::cout << "================================================================\n";
    
    return 0;
}
