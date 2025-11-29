/**
 * @file test_propagation_jpl.cpp
 * @brief Test propagazione ±1 mese con confronto JPL Horizons
 * 
 * Propaga elementi AstDyS di (11234) a:
 * - Epoca - 30 giorni (1 mese indietro)
 * - Epoca + 30 giorni (1 mese avanti)
 * 
 * Confronta con effemeridi JPL Horizons
 * 
 * @author AstDyS Team
 * @date 29 Novembre 2025
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <array>
#include <vector>
#include <string>

//=============================================================================
// COSTANTI
//=============================================================================

namespace constants {
    constexpr double PI = 3.14159265358979323846;
    constexpr double TWO_PI = 2.0 * PI;
    constexpr double DEG2RAD = PI / 180.0;
    constexpr double RAD2DEG = 180.0 / PI;
    constexpr double ARCSEC2RAD = PI / (180.0 * 3600.0);
    constexpr double RAD2ARCSEC = 180.0 * 3600.0 / PI;
    
    constexpr double k = 0.01720209895;
    constexpr double k2 = k * k;
    constexpr double c_light = 173.1446326846693;
    constexpr double c2 = c_light * c_light;
    constexpr double AU_km = 149597870.7;
    
    constexpr double OBLIQUITY_J2000 = 23.439291111 * DEG2RAD;
    constexpr double JD_J2000 = 2451545.0;
    
    // GM pianeti [AU³/day²]
    constexpr double GM_planets[8] = {
        4.9125474514508118e-11,  // Mercury
        7.2434524861627027e-10,  // Venus
        8.9970116036316091e-10,  // EMB
        9.5495351057792580e-11,  // Mars
        2.8253458420837436e-07,  // Jupiter
        8.4597151856806587e-08,  // Saturn
        1.2920249167819693e-08,  // Uranus
        1.5243589007842762e-08   // Neptune
    };
}

using namespace constants;

//=============================================================================
// STRUTTURE
//=============================================================================

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

struct KeplerElements {
    std::string name;
    double a, e, i, Omega, omega, M;  // i, Omega, omega, M in radianti
    double epoch;  // MJD
};

//=============================================================================
// PARSER AstDyS
//=============================================================================

KeplerElements readAstDySElements(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open: " + filename);
    }
    
    std::string line;
    double a, h, k_elem, p, q, lambda;
    double epoch = 0;
    
    while (std::getline(file, line)) {
        if (line.find("EQU") != std::string::npos) {
            std::istringstream iss(line);
            std::string equ;
            iss >> equ >> a >> h >> k_elem >> p >> q >> lambda;
        }
        if (line.find("MJD") != std::string::npos) {
            std::istringstream iss(line);
            std::string mjd_str;
            iss >> mjd_str >> epoch;
            break;
        }
    }
    
    KeplerElements kep;
    kep.name = "(11234)";
    kep.a = a;
    kep.epoch = epoch;
    
    // Converti equinoziali → kepleriani
    kep.e = std::sqrt(h*h + k_elem*k_elem);
    double LP = std::atan2(h, k_elem);
    
    double tan_i_half = std::sqrt(p*p + q*q);
    kep.i = 2.0 * std::atan(tan_i_half);
    kep.Omega = std::atan2(p, q);
    kep.omega = LP - kep.Omega;
    kep.M = lambda * DEG2RAD - LP;
    
    // Normalizza
    while (kep.omega < 0) kep.omega += TWO_PI;
    while (kep.omega >= TWO_PI) kep.omega -= TWO_PI;
    while (kep.M < 0) kep.M += TWO_PI;
    while (kep.M >= TWO_PI) kep.M -= TWO_PI;
    
    return kep;
}

//=============================================================================
// EFFEMERIDI PLANETARIE
//=============================================================================

Vec3 getPlanetPosition(double jd, int planet) {
    double T = (jd - JD_J2000) / 36525.0;
    
    double a, e, I, L, omega_bar, Omega;
    
    switch(planet) {
        case 1: a=0.38709927+0.00000037*T; e=0.20563593+0.00001906*T; I=7.00497902-0.00594749*T;
                L=252.25032350+149472.67411175*T; omega_bar=77.45779628+0.16047689*T; Omega=48.33076593-0.12534081*T; break;
        case 2: a=0.72333566+0.00000390*T; e=0.00677672-0.00004107*T; I=3.39467605-0.00078890*T;
                L=181.97909950+58517.81538729*T; omega_bar=131.60246718+0.00268329*T; Omega=76.67984255-0.27769418*T; break;
        case 3: a=1.00000261+0.00000562*T; e=0.01671123-0.00004392*T; I=-0.00001531-0.01294668*T;
                L=100.46457166+35999.37244981*T; omega_bar=102.93768193+0.32327364*T; Omega=0.0; break;
        case 4: a=1.52371034+0.00001847*T; e=0.09339410+0.00007882*T; I=1.84969142-0.00813131*T;
                L=-4.55343205+19140.30268499*T; omega_bar=-23.94362959+0.44441088*T; Omega=49.55953891-0.29257343*T; break;
        case 5: a=5.20288700-0.00011607*T; e=0.04838624-0.00013253*T; I=1.30439695-0.00183714*T;
                L=34.39644051+3034.74612775*T; omega_bar=14.72847983+0.21252668*T; Omega=100.47390909+0.20469106*T; break;
        case 6: a=9.53667594-0.00125060*T; e=0.05386179-0.00050991*T; I=2.48599187+0.00193609*T;
                L=49.95424423+1222.49362201*T; omega_bar=92.59887831-0.41897216*T; Omega=113.66242448-0.28867794*T; break;
        case 7: a=19.18916464-0.00196176*T; e=0.04725744-0.00004397*T; I=0.77263783-0.00242939*T;
                L=313.23810451+428.48202785*T; omega_bar=170.95427630+0.40805281*T; Omega=74.01692503+0.04240589*T; break;
        case 8: a=30.06992276+0.00026291*T; e=0.00859048+0.00005105*T; I=1.77004347+0.00035372*T;
                L=-55.12002969+218.45945325*T; omega_bar=44.96476227-0.32241464*T; Omega=131.78422574-0.00508664*T; break;
        default: return Vec3();
    }
    
    double i_rad = I * DEG2RAD;
    double O_rad = Omega * DEG2RAD;
    double w = (omega_bar - Omega) * DEG2RAD;
    double M = std::fmod((L - omega_bar) * DEG2RAD, TWO_PI);
    if (M < 0) M += TWO_PI;
    
    double E = M;
    for (int iter = 0; iter < 15; iter++) {
        double dE = (M - E + e * std::sin(E)) / (1 - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < 1e-14) break;
    }
    
    double nu = 2.0 * std::atan2(std::sqrt(1+e) * std::sin(E/2), std::sqrt(1-e) * std::cos(E/2));
    double r = a * (1 - e * std::cos(E));
    
    double x_orb = r * std::cos(nu), y_orb = r * std::sin(nu);
    
    double cO = std::cos(O_rad), sO = std::sin(O_rad);
    double ci = std::cos(i_rad), si = std::sin(i_rad);
    double cw = std::cos(w), sw = std::sin(w);
    
    return Vec3((cO*cw-sO*sw*ci)*x_orb+(-cO*sw-sO*cw*ci)*y_orb,
                (sO*cw+cO*sw*ci)*x_orb+(-sO*sw+cO*cw*ci)*y_orb,
                sw*si*x_orb+cw*si*y_orb);
}

//=============================================================================
// PROPAGATORE CON RK4
//=============================================================================

class Propagator {
public:
    int steps_count = 0;
    
    State elementsToState(const KeplerElements& kep) {
        double a = kep.a, e = kep.e;
        double inc = kep.i, Omega = kep.Omega, omega = kep.omega, M0 = kep.M;
        
        double E = M0;
        for (int iter = 0; iter < 15; iter++) {
            double dE = (E - e * std::sin(E) - M0) / (1 - e * std::cos(E));
            E -= dE;
            if (std::abs(dE) < 1e-14) break;
        }
        
        double sin_E = std::sin(E), cos_E = std::cos(E);
        double sqrt_1_e2 = std::sqrt(1.0 - e * e);
        double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - e);
        double r = a * (1.0 - e * cos_E);
        
        double x_orb = r * std::cos(nu), y_orb = r * std::sin(nu);
        double v_factor = std::sqrt(k2 * a) / r;
        double vx_orb = -v_factor * sin_E, vy_orb = v_factor * sqrt_1_e2 * cos_E;
        
        double cO = std::cos(Omega), sO = std::sin(Omega);
        double cw = std::cos(omega), sw = std::sin(omega);
        double ci = std::cos(inc), si = std::sin(inc);
        
        double P11 = cO*cw - sO*sw*ci, P12 = -cO*sw - sO*cw*ci;
        double P21 = sO*cw + cO*sw*ci, P22 = -sO*sw + cO*cw*ci;
        double P31 = sw*si,            P32 = cw*si;
        
        return State(Vec3(P11*x_orb+P12*y_orb, P21*x_orb+P22*y_orb, P31*x_orb+P32*y_orb),
                     Vec3(P11*vx_orb+P12*vy_orb, P21*vx_orb+P22*vy_orb, P31*vx_orb+P32*vy_orb));
    }
    
    State propagate(const State& y0, double jd0, double jd1) {
        double dt = jd1 - jd0;
        double h = dt > 0 ? std::min(0.1, dt/10) : std::max(-0.1, dt/10);
        double t = jd0;
        State y = y0;
        steps_count = 0;
        
        while ((dt > 0 && t < jd1 - 1e-10) || (dt < 0 && t > jd1 + 1e-10)) {
            if (dt > 0 && t + h > jd1) h = jd1 - t;
            if (dt < 0 && t + h < jd1) h = jd1 - t;
            
            State k1 = deriv(t, y);
            State k2 = deriv(t + h/2, y + k1*(h/2));
            State k3 = deriv(t + h/2, y + k2*(h/2));
            State k4 = deriv(t + h, y + k3*h);
            
            y = y + (k1 + k2*2 + k3*2 + k4) * (h/6);
            t += h;
            steps_count++;
        }
        return y;
    }
    
    void computeApparentPosition(const State& state, double jd,
                                  double& ra, double& dec, double& delta) {
        Vec3 earth = getPlanetPosition(jd, 3);
        Vec3 geo = state.r - earth;
        delta = geo.norm();
        
        // Light-time correction
        for (int iter = 0; iter < 3; iter++) {
            double lt = delta / c_light;
            Vec3 r_retard = state.r - state.v * lt;
            geo = r_retard - earth;
            delta = geo.norm();
        }
        
        // Eclittica → Equatoriale
        double cos_eps = std::cos(OBLIQUITY_J2000), sin_eps = std::sin(OBLIQUITY_J2000);
        double x_eq = geo.x;
        double y_eq = geo.y * cos_eps - geo.z * sin_eps;
        double z_eq = geo.y * sin_eps + geo.z * cos_eps;
        
        ra = std::atan2(y_eq, x_eq);
        if (ra < 0) ra += TWO_PI;
        dec = std::asin(z_eq / delta);
    }
    
private:
    State deriv(double t, const State& y) {
        Vec3 acc;
        double r_norm = y.r.norm();
        double r3 = r_norm * r_norm * r_norm;
        
        // Kepleriano
        acc = y.r * (-k2 / r3);
        
        // Perturbazioni planetarie
        for (int i = 0; i < 8; i++) {
            Vec3 rp = getPlanetPosition(t, i+1);
            Vec3 dr = rp - y.r;
            double dr3 = std::pow(dr.norm(), 3);
            double rp3 = std::pow(rp.norm(), 3);
            acc += dr * (GM_planets[i] / dr3) - rp * (GM_planets[i] / rp3);
        }
        
        // Correzione relativistica (Schwarzschild)
        double v2 = y.v.norm2();
        double rv = y.r.dot(y.v);
        double factor = k2 / (c2 * r3);
        acc += y.r * (factor * (4*k2/r_norm - v2)) + y.v * (factor * 4 * rv);
        
        return State(y.v, acc);
    }
};

//=============================================================================
// UTILITY
//=============================================================================

std::string formatRA(double ra_rad) {
    double h = ra_rad * 12.0 / PI;
    if (h < 0) h += 24;
    int hh = (int)h;
    double m = (h - hh) * 60;
    int mm = (int)m;
    double ss = (m - mm) * 60;
    char buf[32];
    std::snprintf(buf, 32, "%02d %02d %06.3f", hh, mm, ss);
    return buf;
}

std::string formatDec(double dec_rad) {
    double d = dec_rad * RAD2DEG;
    char sign = d >= 0 ? '+' : '-';
    d = std::abs(d);
    int dd = (int)d;
    double m = (d - dd) * 60;
    int mm = (int)m;
    double ss = (m - mm) * 60;
    char buf[32];
    std::snprintf(buf, 32, "%c%02d %02d %05.2f", sign, dd, mm, ss);
    return buf;
}

void mjdToDate(double mjd, int& year, int& month, int& day) {
    int jd = (int)(mjd + 2400001);
    int l = jd + 68569;
    int n = 4 * l / 146097;
    l = l - (146097 * n + 3) / 4;
    int i = 4000 * (l + 1) / 1461001;
    l = l - 1461 * i / 4 + 31;
    int j = 80 * l / 2447;
    day = l - 2447 * j / 80;
    l = j / 11;
    month = j + 2 - 12 * l;
    year = 100 * (n - 49) + i + l;
}

//=============================================================================
// MAIN
//=============================================================================

int main() {
    std::cout << std::fixed;
    std::cout << "================================================================\n";
    std::cout << "  Test Propagazione ±1 Mese con Confronto JPL\n";
    std::cout << "  Asteroide (11234)\n";
    std::cout << "================================================================\n\n";
    
    // Leggi elementi
    std::cout << "1. ELEMENTI ORBITALI (AstDyS)\n";
    KeplerElements kep;
    try {
        kep = readAstDySElements("../data/11234.eq1");
    } catch (...) {
        // Fallback: elementi hardcoded
        kep.name = "(11234)";
        kep.a = 2.6808535916678031;
        // h = e*sin(LP) = 0.032872036471001
        // k = e*cos(LP) = 0.036254405825130
        double h = 0.032872036471001;
        double k_elem = 0.036254405825130;
        kep.e = std::sqrt(h*h + k_elem*k_elem);  // = 0.04894
        double LP = std::atan2(h, k_elem);
        
        // p = tan(i/2)*sin(LN) = 0.103391596538937
        // q = tan(i/2)*cos(LN) = -0.042907901689093
        double p = 0.103391596538937;
        double q = -0.042907901689093;
        double tan_i_half = std::sqrt(p*p + q*q);
        kep.i = 2.0 * std::atan(tan_i_half);
        kep.Omega = std::atan2(p, q);
        kep.omega = LP - kep.Omega;
        kep.M = 235.8395861037268 * DEG2RAD - LP;
        
        while (kep.omega < 0) kep.omega += TWO_PI;
        while (kep.M < 0) kep.M += TWO_PI;
        
        kep.epoch = 61000.0;  // MJD
    }
    
    int y, m, d;
    mjdToDate(kep.epoch, y, m, d);
    
    std::cout << "   Epoca: MJD " << std::setprecision(1) << kep.epoch 
              << " (" << y << "-" << std::setfill('0') << std::setw(2) << m 
              << "-" << std::setw(2) << d << ")\n" << std::setfill(' ');
    std::cout << "   a = " << std::setprecision(7) << kep.a << " AU\n";
    std::cout << "   e = " << std::setprecision(7) << kep.e << "\n";
    std::cout << "   i = " << std::setprecision(4) << kep.i * RAD2DEG << "°\n";
    std::cout << "   Ω = " << kep.Omega * RAD2DEG << "°\n";
    std::cout << "   ω = " << kep.omega * RAD2DEG << "°\n";
    std::cout << "   M = " << kep.M * RAD2DEG << "°\n\n";
    
    // Propagatore
    Propagator prop;
    State y0 = prop.elementsToState(kep);
    double jd0 = kep.epoch + 2400000.5;
    
    std::cout << "2. STATO ALL'EPOCA\n";
    std::cout << "   r = [" << std::setprecision(9) << y0.r.x << ", " 
              << y0.r.y << ", " << y0.r.z << "] AU\n";
    std::cout << "   v = [" << y0.v.x << ", " << y0.v.y << ", " << y0.v.z << "] AU/day\n";
    std::cout << "   |r| = " << y0.r.norm() << " AU\n\n";
    
    // Epoche target
    double dt_month = 30.0;  // 1 mese ≈ 30 giorni
    
    double jd_back = jd0 - dt_month;   // 1 mese indietro
    double jd_fwd = jd0 + dt_month;    // 1 mese avanti
    
    int y1, m1, d1, y2, m2, d2;
    mjdToDate(jd_back - 2400000.5, y1, m1, d1);
    mjdToDate(jd_fwd - 2400000.5, y2, m2, d2);
    
    std::cout << "3. PROPAGAZIONE\n";
    std::cout << "   Epoca iniziale: JD " << std::setprecision(1) << jd0 << "\n";
    std::cout << "   1 mese indietro: JD " << jd_back << " (" << y1 << "-" 
              << std::setfill('0') << std::setw(2) << m1 << "-" << std::setw(2) << d1 << ")\n";
    std::cout << "   1 mese avanti:   JD " << jd_fwd << " (" << y2 << "-" 
              << std::setw(2) << m2 << "-" << std::setw(2) << d2 << ")\n\n" << std::setfill(' ');
    
    // Propaga indietro
    std::cout << "================================================================\n";
    std::cout << "  PROPAGAZIONE 1 MESE INDIETRO (Δt = -30 giorni)\n";
    std::cout << "================================================================\n";
    
    State y_back = prop.propagate(y0, jd0, jd_back);
    std::cout << "   Passi RK4: " << prop.steps_count << "\n";
    std::cout << "   r = [" << std::setprecision(9) << y_back.r.x << ", " 
              << y_back.r.y << ", " << y_back.r.z << "] AU\n";
    
    double ra_back, dec_back, delta_back;
    prop.computeApparentPosition(y_back, jd_back, ra_back, dec_back, delta_back);
    
    std::cout << "\n   POSIZIONE ASTROMETRICA:\n";
    std::cout << "   RA  = " << formatRA(ra_back) << "\n";
    std::cout << "   Dec = " << formatDec(dec_back) << "\n";
    std::cout << "   Δ   = " << std::setprecision(6) << delta_back << " AU\n";
    
    // Propaga avanti
    std::cout << "\n================================================================\n";
    std::cout << "  PROPAGAZIONE 1 MESE AVANTI (Δt = +30 giorni)\n";
    std::cout << "================================================================\n";
    
    State y_fwd = prop.propagate(y0, jd0, jd_fwd);
    std::cout << "   Passi RK4: " << prop.steps_count << "\n";
    std::cout << "   r = [" << std::setprecision(9) << y_fwd.r.x << ", " 
              << y_fwd.r.y << ", " << y_fwd.r.z << "] AU\n";
    
    double ra_fwd, dec_fwd, delta_fwd;
    prop.computeApparentPosition(y_fwd, jd_fwd, ra_fwd, dec_fwd, delta_fwd);
    
    std::cout << "\n   POSIZIONE ASTROMETRICA:\n";
    std::cout << "   RA  = " << formatRA(ra_fwd) << "\n";
    std::cout << "   Dec = " << formatDec(dec_fwd) << "\n";
    std::cout << "   Δ   = " << std::setprecision(6) << delta_fwd << " AU\n";
    
    // Test round-trip
    std::cout << "\n================================================================\n";
    std::cout << "  TEST ROUND-TRIP\n";
    std::cout << "================================================================\n";
    
    // Avanti e poi indietro
    State y_rt = prop.propagate(y_fwd, jd_fwd, jd0);
    Vec3 dr = y_rt.r - y0.r;
    Vec3 dv = y_rt.v - y0.v;
    
    std::cout << "   Propagazione: epoca → +30d → epoca\n";
    std::cout << "   Errore posizione: " << std::scientific << std::setprecision(3) 
              << dr.norm() << " AU = " << dr.norm() * AU_km * 1000 << " m\n";
    std::cout << "   Errore velocità:  " << dv.norm() << " AU/day\n";
    
    // Indietro e poi avanti
    State y_rt2 = prop.propagate(y_back, jd_back, jd0);
    Vec3 dr2 = y_rt2.r - y0.r;
    
    std::cout << "\n   Propagazione: epoca → -30d → epoca\n";
    std::cout << "   Errore posizione: " << dr2.norm() << " AU = " 
              << dr2.norm() * AU_km * 1000 << " m\n";
    
    // Riepilogo per confronto JPL
    std::cout << std::fixed << "\n================================================================\n";
    std::cout << "  RIEPILOGO PER CONFRONTO JPL HORIZONS\n";
    std::cout << "================================================================\n";
    std::cout << "  Oggetto: (11234)\n\n";
    
    std::cout << "  EPOCA -30 GIORNI (JD " << std::setprecision(1) << jd_back << "):\n";
    std::cout << "    Data: " << y1 << "-" << std::setfill('0') << std::setw(2) << m1 
              << "-" << std::setw(2) << d1 << " 00:00 TDB\n" << std::setfill(' ');
    std::cout << "    AstDyn RA:  " << formatRA(ra_back) << "\n";
    std::cout << "    AstDyn Dec: " << formatDec(dec_back) << "\n";
    std::cout << "    AstDyn Δ:   " << std::setprecision(6) << delta_back << " AU\n\n";
    
    std::cout << "  EPOCA +30 GIORNI (JD " << std::setprecision(1) << jd_fwd << "):\n";
    std::cout << "    Data: " << y2 << "-" << std::setfill('0') << std::setw(2) << m2 
              << "-" << std::setw(2) << d2 << " 00:00 TDB\n" << std::setfill(' ');
    std::cout << "    AstDyn RA:  " << formatRA(ra_fwd) << "\n";
    std::cout << "    AstDyn Dec: " << formatDec(dec_fwd) << "\n";
    std::cout << "    AstDyn Δ:   " << std::setprecision(6) << delta_fwd << " AU\n\n";
    
    std::cout << "  Per verificare, usa JPL Horizons:\n";
    std::cout << "  https://ssd.jpl.nasa.gov/horizons/\n";
    std::cout << "  Target: 11234\n";
    std::cout << "  Observer: @sun (eliocentrico) o 500 (geocentrico)\n";
    std::cout << "================================================================\n";
    
    return 0;
}
