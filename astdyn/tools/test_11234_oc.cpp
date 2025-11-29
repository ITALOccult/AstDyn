/**
 * @file test_11234_oc.cpp
 * @brief Test O-C con dati reali AstDyS per (11234)
 * 
 * Dati scaricati da:
 * - Elementi: https://newton.spacedys.com/~astdys2/epoch/numbered/11/11234.eq1
 * - Osservazioni: https://newton.spacedys.com/~astdys2/mpcobs/numbered/11/11234.rwo
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
#include <algorithm>
#include <map>

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
    double a;       // Semiasse maggiore [AU]
    double e;       // Eccentricità
    double i;       // Inclinazione [rad]
    double Omega;   // Longitudine nodo [rad]
    double omega;   // Argomento perielio [rad]
    double M;       // Anomalia media [rad]
    double epoch;   // Epoca [MJD]
};

struct Observation {
    double mjd;
    double ra_obs;      // [rad]
    double dec_obs;     // [rad]
    double ra_sigma;    // [arcsec]
    double dec_sigma;   // [arcsec]
    std::string obs_code;
    double ra_residual;  // [arcsec]
    double dec_residual; // [arcsec]
    double chi;
    bool is_outlier;
};

struct FitStats {
    int n_obs;
    int n_outliers;
    double rms_ra;
    double rms_dec;
    double rms_total;
    double chi2_reduced;
    double time_span_years;
};

//=============================================================================
// PARSER AstDyS
//=============================================================================

/**
 * @brief Legge elementi equinoziali da file .eq1 AstDyS e converte in Kepler
 */
KeplerElements readAstDySElements(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open: " + filename);
    }
    
    std::string line;
    double a, h, k_elem, p, q, lambda;
    double epoch = 0;
    
    while (std::getline(file, line)) {
        // Cerca riga EQU
        if (line.find("EQU") != std::string::npos) {
            std::istringstream iss(line);
            std::string equ;
            iss >> equ >> a >> h >> k_elem >> p >> q >> lambda;
        }
        // Cerca epoca MJD
        if (line.find("MJD") != std::string::npos) {
            std::istringstream iss(line);
            std::string mjd_str;
            iss >> mjd_str >> epoch;
            break;
        }
    }
    
    // Converti equinoziali → kepleriani
    KeplerElements kep;
    kep.a = a;
    kep.epoch = epoch;
    
    // e*sin(LP), e*cos(LP) → e, omega+Omega
    kep.e = std::sqrt(h*h + k_elem*k_elem);
    double LP = std::atan2(h, k_elem);  // Long. perielio [rad]
    
    // tan(i/2)*sin(LN), tan(i/2)*cos(LN) → i, Omega
    double tan_i_half = std::sqrt(p*p + q*q);
    kep.i = 2.0 * std::atan(tan_i_half);
    kep.Omega = std::atan2(p, q);  // Long. nodo [rad]
    
    // omega = LP - Omega
    kep.omega = LP - kep.Omega;
    
    // lambda = M + LP → M = lambda - LP
    kep.M = lambda * DEG2RAD - LP;
    
    // Normalizza angoli
    while (kep.omega < 0) kep.omega += TWO_PI;
    while (kep.omega >= TWO_PI) kep.omega -= TWO_PI;
    while (kep.M < 0) kep.M += TWO_PI;
    while (kep.M >= TWO_PI) kep.M -= TWO_PI;
    
    return kep;
}

/**
 * @brief Legge osservazioni da file .rwo AstDyS
 */
std::vector<Observation> readAstDySRWO(const std::string& filename, int max_obs = 0) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open: " + filename);
    }
    
    std::vector<Observation> obs;
    std::string line;
    bool in_header = true;
    
    while (std::getline(file, line)) {
        // Skip header
        if (line.find("END_OF_HEADER") != std::string::npos) {
            in_header = false;
            continue;
        }
        if (in_header || line.empty() || line[0] == '!') continue;
        
        // Parse RWO format (very complex, simplified parsing)
        // Format: designation flags date ra dec rms_ra rms_dec ... obs_code ...
        
        // Estrai campi base
        if (line.length() < 100) continue;
        
        try {
            Observation o;
            
            // Designation (primi 10 caratteri circa)
            std::string des = line.substr(0, 10);
            
            // Flags e data
            // Cerchiamo il pattern della data: YYYY MM DD.dddddd
            size_t pos = 14;  // Dopo "11234     O A"
            
            // Cerca l'anno (4 cifre che iniziano con 19 o 20)
            size_t year_pos = line.find("19", pos);
            if (year_pos == std::string::npos) year_pos = line.find("20", pos);
            if (year_pos == std::string::npos) continue;
            
            // Parse data
            int year, month;
            double day;
            std::istringstream date_ss(line.substr(year_pos, 20));
            date_ss >> year >> month >> day;
            
            // Converti in MJD
            int a = (14 - month) / 12;
            int y = year + 4800 - a;
            int m = month + 12*a - 3;
            int jdn = (int)day + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045;
            o.mjd = jdn - 2400000.5 + (day - (int)day);
            
            // Cerca RA: HH MM SS.sss
            size_t ra_pos = year_pos + 20;
            while (ra_pos < line.length() && !std::isdigit(line[ra_pos])) ra_pos++;
            
            int ra_h, ra_m;
            double ra_s;
            std::istringstream ra_ss(line.substr(ra_pos, 15));
            ra_ss >> ra_h >> ra_m >> ra_s;
            o.ra_obs = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0 * DEG2RAD;
            
            // Cerca Dec: ±DD MM SS.ss
            size_t dec_pos = ra_pos + 15;
            while (dec_pos < line.length() && line[dec_pos] != '+' && line[dec_pos] != '-') dec_pos++;
            if (dec_pos >= line.length()) continue;
            
            char dec_sign = line[dec_pos];
            int dec_d, dec_m;
            double dec_s;
            std::istringstream dec_ss(line.substr(dec_pos + 1, 12));
            dec_ss >> dec_d >> dec_m >> dec_s;
            o.dec_obs = (dec_d + dec_m/60.0 + dec_s/3600.0) * DEG2RAD;
            if (dec_sign == '-') o.dec_obs = -o.dec_obs;
            
            // Default errori
            o.ra_sigma = 0.5;
            o.dec_sigma = 0.5;
            
            // Cerca codice osservatorio (3 caratteri dopo la magnitudine)
            size_t code_pos = line.length() - 20;
            while (code_pos < line.length() && !std::isalnum(line[code_pos])) code_pos++;
            if (code_pos + 3 <= line.length()) {
                o.obs_code = line.substr(code_pos, 3);
            }
            
            o.ra_residual = 0;
            o.dec_residual = 0;
            o.chi = 0;
            o.is_outlier = false;
            
            obs.push_back(o);
            
            if (max_obs > 0 && (int)obs.size() >= max_obs) break;
            
        } catch (...) {
            continue;  // Skip malformed lines
        }
    }
    
    // Ordina per epoca
    std::sort(obs.begin(), obs.end(), 
              [](const Observation& a, const Observation& b) { return a.mjd < b.mjd; });
    
    return obs;
}

//=============================================================================
// EFFEMERIDI
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
// PROPAGATORE
//=============================================================================

class Propagator {
public:
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
        
        while ((dt > 0 && t < jd1 - 1e-10) || (dt < 0 && t > jd1 + 1e-10)) {
            if (dt > 0 && t + h > jd1) h = jd1 - t;
            if (dt < 0 && t + h < jd1) h = jd1 - t;
            
            State k1 = deriv(t, y);
            State k2 = deriv(t + h/2, y + k1*(h/2));
            State k3 = deriv(t + h/2, y + k2*(h/2));
            State k4 = deriv(t + h, y + k3*h);
            
            y = y + (k1 + k2*2 + k3*2 + k4) * (h/6);
            t += h;
        }
        return y;
    }
    
    void computeApparentPosition(const State& state, double jd,
                                  double& ra, double& dec, double& delta) {
        Vec3 earth = getPlanetPosition(jd, 3);
        Vec3 geo = state.r - earth;
        delta = geo.norm();
        
        for (int iter = 0; iter < 3; iter++) {
            double lt = delta / c_light;
            Vec3 r_retard = state.r - state.v * lt;
            geo = r_retard - earth;
            delta = geo.norm();
        }
        
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
        
        acc = y.r * (-k2 / r3);
        
        for (int i = 0; i < 8; i++) {
            Vec3 rp = getPlanetPosition(t, i+1);
            Vec3 dr = rp - y.r;
            double dr3 = std::pow(dr.norm(), 3);
            double rp3 = std::pow(rp.norm(), 3);
            acc += dr * (GM_planets[i] / dr3) - rp * (GM_planets[i] / rp3);
        }
        
        double v2 = y.v.norm2();
        double rv = y.r.dot(y.v);
        double factor = k2 / (c2 * r3);
        acc += y.r * (factor * (4*k2/r_norm - v2)) + y.v * (factor * 4 * rv);
        
        return State(y.v, acc);
    }
};

//=============================================================================
// MAIN
//=============================================================================

int main(int argc, char* argv[]) {
    std::cout << std::fixed;
    std::cout << "================================================================\n";
    std::cout << "  Test O-C con dati reali AstDyS\n";
    std::cout << "  Asteroide (11234)\n";
    std::cout << "================================================================\n\n";
    
    std::string data_dir = "../data/";
    if (argc > 1) data_dir = argv[1];
    
    // Leggi elementi
    std::cout << "1. LETTURA ELEMENTI ORBITALI\n";
    KeplerElements kep;
    try {
        kep = readAstDySElements(data_dir + "11234.eq1");
        std::cout << "   File: 11234.eq1\n";
        std::cout << "   Epoca: MJD " << kep.epoch << "\n";
        std::cout << "   a = " << kep.a << " AU\n";
        std::cout << "   e = " << kep.e << "\n";
        std::cout << "   i = " << kep.i * RAD2DEG << "°\n";
        std::cout << "   Ω = " << kep.Omega * RAD2DEG << "°\n";
        std::cout << "   ω = " << kep.omega * RAD2DEG << "°\n";
        std::cout << "   M = " << kep.M * RAD2DEG << "°\n\n";
    } catch (const std::exception& e) {
        std::cerr << "   ERRORE: " << e.what() << "\n";
        return 1;
    }
    
    // Leggi osservazioni
    std::cout << "2. LETTURA OSSERVAZIONI RWO\n";
    std::vector<Observation> obs;
    try {
        obs = readAstDySRWO(data_dir + "11234.rwo", 500);  // Limita a 500
        std::cout << "   File: 11234.rwo\n";
        std::cout << "   Osservazioni lette: " << obs.size() << "\n";
        if (!obs.empty()) {
            std::cout << "   Prima: MJD " << obs.front().mjd << "\n";
            std::cout << "   Ultima: MJD " << obs.back().mjd << "\n";
            std::cout << "   Arco: " << (obs.back().mjd - obs.front().mjd) / 365.25 << " anni\n\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "   ERRORE: " << e.what() << "\n";
        return 1;
    }
    
    // Propagatore
    Propagator prop;
    State y0 = prop.elementsToState(kep);
    double jd0 = kep.epoch + 2400000.5;
    
    std::cout << "3. STATO INIZIALE\n";
    std::cout << "   r = [" << y0.r.x << ", " << y0.r.y << ", " << y0.r.z << "] AU\n";
    std::cout << "   v = [" << y0.v.x << ", " << y0.v.y << ", " << y0.v.z << "] AU/day\n\n";
    
    // Calcola residui
    std::cout << "4. CALCOLO RESIDUI O-C\n";
    
    double sum_ra2 = 0, sum_dec2 = 0;
    int n_computed = 0;
    int n_outliers = 0;
    double max_chi = 0;
    
    // Statistiche per osservatorio
    std::map<std::string, std::vector<double>> obs_stats;
    
    for (size_t i = 0; i < obs.size(); i++) {
        double jd_obs = obs[i].mjd + 2400000.5;
        
        // Propaga
        State y = prop.propagate(y0, jd0, jd_obs);
        
        // Posizione apparente
        double ra_calc, dec_calc, delta;
        prop.computeApparentPosition(y, jd_obs, ra_calc, dec_calc, delta);
        
        // Residui
        double cos_dec = std::cos(obs[i].dec_obs);
        obs[i].ra_residual = (obs[i].ra_obs - ra_calc) * cos_dec * RAD2ARCSEC;
        obs[i].dec_residual = (obs[i].dec_obs - dec_calc) * RAD2ARCSEC;
        
        // Chi
        double chi_ra = obs[i].ra_residual / obs[i].ra_sigma;
        double chi_dec = obs[i].dec_residual / obs[i].dec_sigma;
        obs[i].chi = std::sqrt(chi_ra*chi_ra + chi_dec*chi_dec);
        
        obs[i].is_outlier = (obs[i].chi > 5.0);
        
        if (!obs[i].is_outlier) {
            sum_ra2 += obs[i].ra_residual * obs[i].ra_residual;
            sum_dec2 += obs[i].dec_residual * obs[i].dec_residual;
            n_computed++;
            
            obs_stats[obs[i].obs_code].push_back(obs[i].chi);
        } else {
            n_outliers++;
        }
        
        if (obs[i].chi > max_chi) max_chi = obs[i].chi;
        
        // Progress
        if ((i + 1) % 100 == 0) {
            std::cout << "   " << (i+1) << "/" << obs.size() << " completate\r" << std::flush;
        }
    }
    std::cout << "   " << obs.size() << "/" << obs.size() << " completate\n\n";
    
    // Statistiche
    FitStats stats;
    stats.n_obs = obs.size();
    stats.n_outliers = n_outliers;
    stats.rms_ra = std::sqrt(sum_ra2 / n_computed);
    stats.rms_dec = std::sqrt(sum_dec2 / n_computed);
    stats.rms_total = std::sqrt((sum_ra2 + sum_dec2) / (2 * n_computed));
    stats.chi2_reduced = (sum_ra2/0.25 + sum_dec2/0.25) / (2 * n_computed - 6);
    stats.time_span_years = (obs.back().mjd - obs.front().mjd) / 365.25;
    
    std::cout << "================================================================\n";
    std::cout << "  RISULTATI FIT O-C\n";
    std::cout << "================================================================\n";
    std::cout << "  Osservazioni totali:  " << stats.n_obs << "\n";
    std::cout << "  Outlier (>5σ):        " << stats.n_outliers << "\n";
    std::cout << "  Usate nel fit:        " << n_computed << "\n";
    std::cout << "  Arco temporale:       " << stats.time_span_years << " anni\n\n";
    
    std::cout << "  RMS residui:\n";
    std::cout << "    RA*cos(δ): " << std::setprecision(3) << stats.rms_ra << " arcsec\n";
    std::cout << "    Dec:       " << stats.rms_dec << " arcsec\n";
    std::cout << "    Totale:    " << stats.rms_total << " arcsec\n\n";
    
    std::cout << "  χ² ridotto:  " << std::setprecision(2) << stats.chi2_reduced << "\n";
    std::cout << "  Chi max:     " << max_chi << "\n\n";
    
    // Top 10 osservatori
    std::cout << "  Osservatori principali:\n";
    std::vector<std::pair<std::string, int>> obs_count;
    for (const auto& p : obs_stats) {
        obs_count.push_back({p.first, (int)p.second.size()});
    }
    std::sort(obs_count.begin(), obs_count.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    for (int i = 0; i < std::min(10, (int)obs_count.size()); i++) {
        std::cout << "    " << obs_count[i].first << ": " << obs_count[i].second << " obs\n";
    }
    
    // Sample residui
    std::cout << "\n  Primi 10 residui:\n";
    std::cout << "  " << std::string(70, '-') << "\n";
    std::cout << "  MJD          RA res\"   Dec res\"   Chi    Cod\n";
    std::cout << "  " << std::string(70, '-') << "\n";
    
    for (int i = 0; i < std::min(10, (int)obs.size()); i++) {
        const auto& o = obs[i];
        std::cout << "  " << std::setprecision(2) << o.mjd
                  << "  " << std::setw(8) << std::setprecision(3) << o.ra_residual
                  << "  " << std::setw(8) << o.dec_residual
                  << "  " << std::setw(6) << std::setprecision(2) << o.chi
                  << "  " << o.obs_code
                  << (o.is_outlier ? " *" : "")
                  << "\n";
    }
    
    std::cout << "\n================================================================\n";
    
    return 0;
}
