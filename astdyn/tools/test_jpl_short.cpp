/**
 * @file test_jpl_short.cpp
 * @brief Test propagazione BREVE da elementi JPL (1 mese)
 * 
 * Verifica precisione su breve termine partendo dall'epoca degli elementi
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>

// Costanti
constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double GM_SUN = 2.9591220828559093e-04;  // AU³/day²

// Elementi orbitali JPL per (11234) - ECLITTICA J2000
constexpr double JPL_EPOCH_JD = 2458509.5;  // 2019-Jan-26.00 TDB
constexpr double JPL_A = 2.6809424682551;
constexpr double JPL_E = 0.0494518850848228;
constexpr double JPL_I = 12.77459773563021;  // deg - eclittica
constexpr double JPL_OM = 112.6647729938423; // Omega
constexpr double JPL_W = 292.4561825062523;  // omega
constexpr double JPL_MA = 351.4222763330189; // Mean anomaly

// Stato ICRF iniziale da JPL (equatoriale!)
constexpr double JPL_X = 2.015534527930346;
constexpr double JPL_Y = 1.560170291279843;
constexpr double JPL_Z = 0.07755625121716653;
constexpr double JPL_VX = -0.006439826187731527;
constexpr double JPL_VY = 0.007976810840048847;
constexpr double JPL_VZ = 0.004075596542667446;

using Vec3 = std::array<double, 3>;
using State = std::array<double, 6>;

double norm(const Vec3& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Risoluzione equazione di Keplero
double solveKepler(double M, double e) {
    double E = M;
    for (int i = 0; i < 50; ++i) {
        double dE = (M - E + e * std::sin(E)) / (1.0 - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < 1e-14) break;
    }
    return E;
}

// Elementi kepleriani -> stato cartesiano
// frame: 0 = eclittica, 1 = equatoriale (ICRF)
State keplerToState(double a, double e, double i, double Omega, double omega, double M, int frame = 0) {
    double E = solveKepler(M, e);
    
    double cosE = std::cos(E);
    double sinE = std::sin(E);
    double sqrtae = std::sqrt(1.0 - e*e);
    
    double x_orb = a * (cosE - e);
    double y_orb = a * sqrtae * sinE;
    
    double n = std::sqrt(GM_SUN / (a*a*a));
    double vx_orb = -a * n * sinE / (1.0 - e * cosE);
    double vy_orb = a * n * sqrtae * cosE / (1.0 - e * cosE);
    
    double cosO = std::cos(Omega);
    double sinO = std::sin(Omega);
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
    
    State s;
    s[0] = P1 * x_orb + P2 * y_orb;
    s[1] = Q1 * x_orb + Q2 * y_orb;
    s[2] = W1 * x_orb + W2 * y_orb;
    s[3] = P1 * vx_orb + P2 * vy_orb;
    s[4] = Q1 * vx_orb + Q2 * vy_orb;
    s[5] = W1 * vx_orb + W2 * vy_orb;
    
    // Se richiesto frame equatoriale, ruota da eclittica
    if (frame == 1) {
        constexpr double eps = 23.4392911 * DEG2RAD;
        double c = std::cos(eps);
        double ss = std::sin(eps);
        
        double y_new = c * s[1] + ss * s[2];
        double z_new = -ss * s[1] + c * s[2];
        s[1] = y_new;
        s[2] = z_new;
        
        double vy_new = c * s[4] + ss * s[5];
        double vz_new = -ss * s[4] + c * s[5];
        s[4] = vy_new;
        s[5] = vz_new;
    }
    
    return s;
}

// Accelerazione (solo Sole per test breve)
State computeAcceleration(const State& s, double jd) {
    double r = std::sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
    double r3 = r * r * r;
    
    return {s[3], s[4], s[5],
            -GM_SUN * s[0] / r3,
            -GM_SUN * s[1] / r3,
            -GM_SUN * s[2] / r3};
}

// RK4
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

// Calcola RA/Dec da stato equatoriale
void computeRADec(const State& s, double& ra_h, double& ra_m, double& ra_s,
                  double& dec_d, double& dec_m, double& dec_s) {
    double ra = std::atan2(s[1], s[0]);
    if (ra < 0) ra += 2*PI;
    double dec = std::atan2(s[2], std::sqrt(s[0]*s[0] + s[1]*s[1]));
    
    double ra_deg = ra * RAD2DEG;
    double ra_hours = ra_deg / 15.0;
    ra_h = std::floor(ra_hours);
    ra_m = std::floor((ra_hours - ra_h) * 60.0);
    ra_s = ((ra_hours - ra_h) * 60.0 - ra_m) * 60.0;
    
    double dec_deg = dec * RAD2DEG;
    int sign = (dec_deg >= 0) ? 1 : -1;
    dec_deg = std::abs(dec_deg);
    dec_d = sign * std::floor(dec_deg);
    dec_m = std::floor((dec_deg - std::abs(dec_d)) * 60.0);
    dec_s = ((dec_deg - std::abs(dec_d)) * 60.0 - dec_m) * 60.0;
}

int main() {
    std::cout << std::fixed;
    
    std::cout << "================================================================\n";
    std::cout << "  Test Propagazione BREVE da Elementi JPL\n";
    std::cout << "  Asteroide (11234) 1999 JS82\n";
    std::cout << "================================================================\n\n";
    
    // METODO 1: Stato ICRF diretto da JPL
    std::cout << "1. STATO ICRF DA JPL (EQUATORIALE)\n";
    State s_icrf = {JPL_X, JPL_Y, JPL_Z, JPL_VX, JPL_VY, JPL_VZ};
    std::cout << std::setprecision(12);
    std::cout << "   r = [" << s_icrf[0] << ", " << s_icrf[1] << ", " << s_icrf[2] << "] AU\n";
    std::cout << "   v = [" << s_icrf[3] << ", " << s_icrf[4] << ", " << s_icrf[5] << "] AU/day\n";
    std::cout << "   |r| = " << norm({s_icrf[0], s_icrf[1], s_icrf[2]}) << " AU\n\n";
    
    // METODO 2: Elementi kepleriani -> stato eclittico -> equatoriale
    std::cout << "2. ELEMENTI KEPLERIANI -> STATO EQUATORIALE\n";
    State s_kep = keplerToState(JPL_A, JPL_E, 
                                 JPL_I * DEG2RAD, 
                                 JPL_OM * DEG2RAD, 
                                 JPL_W * DEG2RAD, 
                                 JPL_MA * DEG2RAD, 
                                 1);  // frame equatoriale
    std::cout << "   r = [" << s_kep[0] << ", " << s_kep[1] << ", " << s_kep[2] << "] AU\n";
    std::cout << "   v = [" << s_kep[3] << ", " << s_kep[4] << ", " << s_kep[5] << "] AU/day\n";
    std::cout << "   |r| = " << norm({s_kep[0], s_kep[1], s_kep[2]}) << " AU\n\n";
    
    // Differenza
    std::cout << "3. CONFRONTO STATI ALL'EPOCA\n";
    double dr = norm({s_icrf[0]-s_kep[0], s_icrf[1]-s_kep[1], s_icrf[2]-s_kep[2]});
    double dv = norm({s_icrf[3]-s_kep[3], s_icrf[4]-s_kep[4], s_icrf[5]-s_kep[5]});
    std::cout << std::setprecision(6);
    std::cout << "   Δr = " << dr << " AU = " << dr * 149597870.7 << " km\n";
    std::cout << "   Δv = " << dv << " AU/day = " << dv * 149597870.7 / 86400.0 << " km/s\n\n";
    
    // Propaga 1 mese (31 giorni) - JD 2458509.5 -> 2458540.5 (2019-Feb-26)
    double jd_target = 2458540.5;  // 2019-Feb-26
    double dt = jd_target - JPL_EPOCH_JD;
    
    std::cout << "================================================================\n";
    std::cout << "  PROPAGAZIONE +31 GIORNI (ICRF diretto)\n";
    std::cout << "================================================================\n";
    
    State s_prop = propagateRK4(s_icrf, JPL_EPOCH_JD, jd_target, 310);
    
    std::cout << std::setprecision(9);
    std::cout << "   Epoca: JD " << std::setprecision(1) << jd_target << " (2019-Feb-26)\n";
    std::cout << std::setprecision(9);
    std::cout << "   r = [" << s_prop[0] << ", " << s_prop[1] << ", " << s_prop[2] << "] AU\n";
    std::cout << "   |r| = " << norm({s_prop[0], s_prop[1], s_prop[2]}) << " AU\n";
    
    double ra_h, ra_m, ra_s, dec_d, dec_m, dec_s;
    computeRADec(s_prop, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s);
    
    std::cout << "\n   POSIZIONE ASTROMETRICA (ICRF):\n";
    std::cout << std::setprecision(2);
    std::cout << "   RA  = " << std::setw(2) << (int)ra_h << " " 
              << std::setw(2) << (int)ra_m << " " 
              << std::setprecision(3) << std::setw(6) << ra_s << "\n";
    std::cout << "   Dec = " << std::showpos << std::setw(3) << (int)dec_d << std::noshowpos
              << " " << std::setw(2) << (int)dec_m << " " 
              << std::setprecision(2) << std::setw(5) << dec_s << "\n";
    
    std::cout << "\n   JPL HORIZONS (2019-Feb-26, eliocentrico):\n";
    std::cout << "   RA  =  2 59 36.94\n";
    std::cout << "   Dec = +04 33 53.2\n";
    std::cout << "   Δ   = 2.5484 AU\n";
    
    // Calcola differenze
    double jpl_ra = 2.0 + 59.0/60.0 + 36.94/3600.0;
    double astdyn_ra = ra_h + ra_m/60.0 + ra_s/3600.0;
    double delta_ra = (astdyn_ra - jpl_ra) * 15.0 * 3600.0; // arcsec
    
    double jpl_dec = 4.0 + 33.0/60.0 + 53.2/3600.0;
    double astdyn_dec = std::abs(dec_d) + dec_m/60.0 + dec_s/3600.0;
    if (dec_d < 0) astdyn_dec = -astdyn_dec;
    double delta_dec = (astdyn_dec - jpl_dec) * 3600.0; // arcsec
    
    std::cout << std::setprecision(2);
    std::cout << "\n   DIFFERENZA:\n";
    std::cout << "   ΔRA  = " << delta_ra << " arcsec\n";
    std::cout << "   ΔDec = " << delta_dec << " arcsec\n";
    
    double total_err = std::sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
    std::cout << "   Totale = " << total_err << " arcsec = " 
              << total_err/60.0 << " arcmin\n";
    
    std::cout << "\n================================================================\n";
    std::cout << "  ANALISI\n";
    std::cout << "================================================================\n";
    std::cout << "  Intervallo propagazione: " << dt << " giorni\n";
    std::cout << "  Modello: Solo gravità solare (kepleriano)\n";
    std::cout << "  Frame: ICRF (equatoriale J2000)\n";
    
    if (total_err < 60) {
        std::cout << "  RISULTATO: ✓ Errore < 1 arcmin (accettabile per 1 mese)\n";
    } else if (total_err < 300) {
        std::cout << "  RISULTATO: ◯ Errore < 5 arcmin (perturbazioni mancanti)\n";
    } else {
        std::cout << "  RISULTATO: ✗ Errore > 5 arcmin (problema sistematico)\n";
    }
    std::cout << "================================================================\n";
    
    return 0;
}
