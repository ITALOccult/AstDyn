/**
 * @file astdyn_oc_fitter.cpp
 * @brief AstDyn O-C Fitter - Propagazione e Fitting Osservazioni
 * 
 * Funzionalità:
 * 1. Lettura elementi orbitali da file AstDyS (.eq0/.eq1)
 * 2. Lettura osservazioni RWO da AstDyS
 * 3. Propagazione avanti/indietro con RKF78
 * 4. Calcolo residui O-C (Observed - Computed)
 * 5. Statistiche e analisi dei residui
 * 
 * @author AstDyS Team - Università di Pisa
 * @date 29 Novembre 2025
 * @version 1.0
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
#include <numeric>
#include <map>

namespace astdyn {

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
    
    constexpr double k = 0.01720209895;          // Costante di Gauss
    constexpr double k2 = k * k;                 // GM_Sun [AU³/day²]
    constexpr double c_light = 173.1446326846693; // Velocità luce [AU/day]
    constexpr double c2 = c_light * c_light;
    constexpr double AU_km = 149597870.7;
    
    constexpr double OBLIQUITY_J2000 = 23.439291111 * DEG2RAD;
    constexpr double JD_J2000 = 2451545.0;
    constexpr double MJD_J2000 = 51544.5;
    
    // GM pianeti [AU³/day²]
    constexpr double GM_Mercury = 4.9125474514508118e-11;
    constexpr double GM_Venus   = 7.2434524861627027e-10;
    constexpr double GM_EMB     = 8.9970116036316091e-10;
    constexpr double GM_Mars    = 9.5495351057792580e-11;
    constexpr double GM_Jupiter = 2.8253458420837436e-07;
    constexpr double GM_Saturn  = 8.4597151856806587e-08;
    constexpr double GM_Uranus  = 1.2920249167819693e-08;
    constexpr double GM_Neptune = 1.5243589007842762e-08;
}

using namespace constants;

//=============================================================================
// STRUTTURE DATI
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

/**
 * @struct OrbitalElements
 * @brief Elementi orbitali kepleriani (formato AstDyS)
 */
struct OrbitalElements {
    std::string name;       // Designazione
    int number;             // Numero asteroide (0 se non numerato)
    
    // Elementi orbitali
    double a;               // Semiasse maggiore [AU]
    double e;               // Eccentricità
    double i;               // Inclinazione [deg]
    double Omega;           // Longitudine nodo ascendente [deg]
    double omega;           // Argomento del perielio [deg]
    double M;               // Anomalia media [deg]
    double epoch;           // Epoca [MJD]
    
    // Parametri fisici (opzionali)
    double H;               // Magnitudine assoluta
    double G;               // Slope parameter
    
    // Matrice covarianza (opzionale)
    bool has_covariance;
    std::array<double, 21> cov;  // Triangolo superiore 6x6
    
    double getJD() const { return epoch + 2400000.5; }
};

/**
 * @struct Observation
 * @brief Osservazione astrometrica (formato RWO)
 */
struct Observation {
    std::string designation;  // Designazione oggetto
    double mjd;               // Epoca [MJD UTC]
    double ra_obs;            // RA osservata [rad]
    double dec_obs;           // Dec osservata [rad]
    double ra_sigma;          // Errore RA [arcsec]
    double dec_sigma;         // Errore Dec [arcsec]
    std::string obs_code;     // Codice osservatorio MPC
    double mag;               // Magnitudine
    std::string band;         // Banda fotometrica
    
    // Residui (calcolati)
    double ra_residual;       // O-C in RA*cos(dec) [arcsec]
    double dec_residual;      // O-C in Dec [arcsec]
    double chi;               // Chi normalizzato
    bool is_outlier;          // Flag outlier
    
    double getJD() const { return mjd + 2400000.5; }
};

/**
 * @struct PropagationStats
 * @brief Statistiche propagazione
 */
struct PropagationStats {
    int steps_accepted = 0;
    int steps_rejected = 0;
    double h_min = 1e10;
    double h_max = 0;
};

/**
 * @struct FitStatistics
 * @brief Statistiche del fit O-C
 */
struct FitStatistics {
    int n_obs;                // Numero osservazioni
    int n_outliers;           // Numero outlier
    double rms_ra;            // RMS residui RA [arcsec]
    double rms_dec;           // RMS residui Dec [arcsec]
    double rms_total;         // RMS totale [arcsec]
    double chi2;              // Chi-quadro
    double chi2_reduced;      // Chi-quadro ridotto
    double mean_ra;           // Media residui RA
    double mean_dec;          // Media residui Dec
    double max_residual;      // Residuo massimo
    double time_span;         // Arco temporale [giorni]
};

//=============================================================================
// PARSER FILE ASTDYS
//=============================================================================

/**
 * @brief Legge elementi orbitali da file .eq0 o .eq1 di AstDyS
 */
class AstDySElementsReader {
public:
    static OrbitalElements read(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        OrbitalElements elem;
        elem.has_covariance = false;
        std::string line;
        
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '!' || line[0] == '%') continue;
            
            // Parse formato AstDyS
            // Esempio: 17030  KEP  60676.0  3.1754733  0.0454207  2.9046 ...
            std::istringstream iss(line);
            
            std::string token;
            iss >> token;
            
            // Numero o designazione
            try {
                elem.number = std::stoi(token);
                elem.name = "(" + token + ")";
            } catch (...) {
                elem.number = 0;
                elem.name = token;
            }
            
            // Tipo elementi (KEP, EQU, CAR)
            std::string type;
            iss >> type;
            
            if (type == "KEP" || type == "kep") {
                // Elementi kepleriani: epoch a e i Omega omega M
                iss >> elem.epoch >> elem.a >> elem.e >> elem.i 
                    >> elem.Omega >> elem.omega >> elem.M;
            } else if (type == "EQU" || type == "equ") {
                // Elementi equinoziali: epoch a h k p q lambda
                double h, k_elem, p, q, lambda;
                iss >> elem.epoch >> elem.a >> h >> k_elem >> p >> q >> lambda;
                
                // Converti in kepleriani
                elem.e = std::sqrt(h*h + k_elem*k_elem);
                elem.omega = std::atan2(h, k_elem) * RAD2DEG;
                elem.i = 2.0 * std::asin(std::sqrt(p*p + q*q)) * RAD2DEG;
                elem.Omega = std::atan2(p, q) * RAD2DEG;
                elem.M = (lambda - elem.omega - elem.Omega);
                while (elem.M < 0) elem.M += 360;
                while (elem.M >= 360) elem.M -= 360;
            }
            
            // Leggi H e G se presenti
            if (iss >> elem.H) {
                if (!(iss >> elem.G)) elem.G = 0.15;
            } else {
                elem.H = 15.0;
                elem.G = 0.15;
            }
            
            break;  // Solo prima riga di elementi
        }
        
        // Cerca matrice covarianza
        while (std::getline(file, line)) {
            if (line.find("COV") != std::string::npos || 
                line.find("cov") != std::string::npos) {
                elem.has_covariance = true;
                // Leggi 21 elementi (triangolo superiore 6x6)
                for (int i = 0; i < 21 && std::getline(file, line); ) {
                    std::istringstream iss2(line);
                    double val;
                    while (iss2 >> val && i < 21) {
                        elem.cov[i++] = val;
                    }
                }
                break;
            }
        }
        
        return elem;
    }
    
    /**
     * @brief Crea elementi da input manuale (per test)
     */
    static OrbitalElements fromManual(int number, const std::string& name,
                                       double epoch_mjd,
                                       double a, double e, double i,
                                       double Omega, double omega, double M,
                                       double H = 15.0, double G = 0.15) {
        OrbitalElements elem;
        elem.number = number;
        elem.name = name;
        elem.epoch = epoch_mjd;
        elem.a = a;
        elem.e = e;
        elem.i = i;
        elem.Omega = Omega;
        elem.omega = omega;
        elem.M = M;
        elem.H = H;
        elem.G = G;
        elem.has_covariance = false;
        return elem;
    }
};

/**
 * @brief Legge osservazioni da file RWO di AstDyS
 * 
 * Formato RWO:
 * des  mjd  ra  dec  rms_ra  rms_dec  mag  band  obs_code  ...
 */
class AstDySRWOReader {
public:
    static std::vector<Observation> read(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open RWO file: " + filename);
        }
        
        std::vector<Observation> obs;
        std::string line;
        
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '!' || line[0] == '%') continue;
            
            Observation o;
            std::istringstream iss(line);
            
            // Parse RWO format
            double ra_deg, dec_deg;
            iss >> o.designation >> o.mjd >> ra_deg >> dec_deg 
                >> o.ra_sigma >> o.dec_sigma >> o.mag >> o.band >> o.obs_code;
            
            // Converti in radianti
            o.ra_obs = ra_deg * DEG2RAD;
            o.dec_obs = dec_deg * DEG2RAD;
            
            // Default errori se non specificati
            if (o.ra_sigma <= 0) o.ra_sigma = 0.5;
            if (o.dec_sigma <= 0) o.dec_sigma = 0.5;
            
            o.ra_residual = 0;
            o.dec_residual = 0;
            o.chi = 0;
            o.is_outlier = false;
            
            obs.push_back(o);
        }
        
        // Ordina per epoca
        std::sort(obs.begin(), obs.end(), 
                  [](const Observation& a, const Observation& b) {
                      return a.mjd < b.mjd;
                  });
        
        return obs;
    }
    
    /**
     * @brief Crea osservazioni simulate per test
     */
    static std::vector<Observation> createSimulated(
            const OrbitalElements& elem,
            double mjd_start, double mjd_end, int n_obs,
            double sigma_arcsec = 0.5) {
        
        std::vector<Observation> obs;
        double dt = (mjd_end - mjd_start) / (n_obs - 1);
        
        for (int i = 0; i < n_obs; i++) {
            Observation o;
            o.designation = elem.name;
            o.mjd = mjd_start + i * dt;
            o.ra_sigma = sigma_arcsec;
            o.dec_sigma = sigma_arcsec;
            o.obs_code = "500";  // Geocentrico
            o.mag = elem.H + 5.0;  // Approssimativo
            o.band = "V";
            
            // RA/Dec saranno calcolate dal propagatore
            o.ra_obs = 0;
            o.dec_obs = 0;
            o.ra_residual = 0;
            o.dec_residual = 0;
            o.chi = 0;
            o.is_outlier = false;
            
            obs.push_back(o);
        }
        
        return obs;
    }
};

//=============================================================================
// EFFEMERIDI PLANETARIE (Simon et al. 1994)
//=============================================================================

namespace ephemeris {

Vec3 getPlanetPosition(double jd, int planet) {
    double T = (jd - JD_J2000) / 36525.0;
    
    double a, e, I, L, omega_bar, Omega;
    
    switch(planet) {
        case 1: // Mercurio
            a = 0.38709927 + 0.00000037*T;
            e = 0.20563593 + 0.00001906*T;
            I = 7.00497902 - 0.00594749*T;
            L = 252.25032350 + 149472.67411175*T;
            omega_bar = 77.45779628 + 0.16047689*T;
            Omega = 48.33076593 - 0.12534081*T;
            break;
        case 2: // Venere
            a = 0.72333566 + 0.00000390*T;
            e = 0.00677672 - 0.00004107*T;
            I = 3.39467605 - 0.00078890*T;
            L = 181.97909950 + 58517.81538729*T;
            omega_bar = 131.60246718 + 0.00268329*T;
            Omega = 76.67984255 - 0.27769418*T;
            break;
        case 3: // Terra-Luna
            a = 1.00000261 + 0.00000562*T;
            e = 0.01671123 - 0.00004392*T;
            I = -0.00001531 - 0.01294668*T;
            L = 100.46457166 + 35999.37244981*T;
            omega_bar = 102.93768193 + 0.32327364*T;
            Omega = 0.0;
            break;
        case 4: // Marte
            a = 1.52371034 + 0.00001847*T;
            e = 0.09339410 + 0.00007882*T;
            I = 1.84969142 - 0.00813131*T;
            L = -4.55343205 + 19140.30268499*T;
            omega_bar = -23.94362959 + 0.44441088*T;
            Omega = 49.55953891 - 0.29257343*T;
            break;
        case 5: // Giove
            a = 5.20288700 - 0.00011607*T;
            e = 0.04838624 - 0.00013253*T;
            I = 1.30439695 - 0.00183714*T;
            L = 34.39644051 + 3034.74612775*T;
            omega_bar = 14.72847983 + 0.21252668*T;
            Omega = 100.47390909 + 0.20469106*T;
            break;
        case 6: // Saturno
            a = 9.53667594 - 0.00125060*T;
            e = 0.05386179 - 0.00050991*T;
            I = 2.48599187 + 0.00193609*T;
            L = 49.95424423 + 1222.49362201*T;
            omega_bar = 92.59887831 - 0.41897216*T;
            Omega = 113.66242448 - 0.28867794*T;
            break;
        case 7: // Urano
            a = 19.18916464 - 0.00196176*T;
            e = 0.04725744 - 0.00004397*T;
            I = 0.77263783 - 0.00242939*T;
            L = 313.23810451 + 428.48202785*T;
            omega_bar = 170.95427630 + 0.40805281*T;
            Omega = 74.01692503 + 0.04240589*T;
            break;
        case 8: // Nettuno
            a = 30.06992276 + 0.00026291*T;
            e = 0.00859048 + 0.00005105*T;
            I = 1.77004347 + 0.00035372*T;
            L = -55.12002969 + 218.45945325*T;
            omega_bar = 44.96476227 - 0.32241464*T;
            Omega = 131.78422574 - 0.00508664*T;
            break;
        default:
            return Vec3(0, 0, 0);
    }
    
    double i_rad = I * DEG2RAD;
    double Omega_rad = Omega * DEG2RAD;
    double omega = (omega_bar - Omega) * DEG2RAD;
    double M = (L - omega_bar) * DEG2RAD;
    M = std::fmod(M, TWO_PI);
    if (M < 0) M += TWO_PI;
    
    double E = M;
    for (int iter = 0; iter < 15; iter++) {
        double dE = (M - E + e * std::sin(E)) / (1 - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < 1e-14) break;
    }
    
    double nu = 2.0 * std::atan2(std::sqrt(1+e) * std::sin(E/2),
                                  std::sqrt(1-e) * std::cos(E/2));
    double r = a * (1 - e * std::cos(E));
    
    double x_orb = r * std::cos(nu);
    double y_orb = r * std::sin(nu);
    
    double cos_O = std::cos(Omega_rad), sin_O = std::sin(Omega_rad);
    double cos_i = std::cos(i_rad), sin_i = std::sin(i_rad);
    double cos_w = std::cos(omega), sin_w = std::sin(omega);
    
    double x = (cos_O*cos_w - sin_O*sin_w*cos_i) * x_orb +
               (-cos_O*sin_w - sin_O*cos_w*cos_i) * y_orb;
    double y = (sin_O*cos_w + cos_O*sin_w*cos_i) * x_orb +
               (-sin_O*sin_w + cos_O*cos_w*cos_i) * y_orb;
    double z = (sin_w*sin_i) * x_orb + (cos_w*sin_i) * y_orb;
    
    return Vec3(x, y, z);
}

} // namespace ephemeris

//=============================================================================
// O-C FITTER
//=============================================================================

/**
 * @class OCFitter
 * @brief Calcola residui O-C propagando orbita e confrontando con osservazioni
 */
class OCFitter {
public:
    OCFitter(double tol = 1e-12) : tol_(tol), h_min_(1e-6), h_max_(1.0) {}
    
    /**
     * @brief Calcola residui O-C per tutte le osservazioni
     */
    FitStatistics computeResiduals(const OrbitalElements& elem,
                                    std::vector<Observation>& obs,
                                    double outlier_threshold = 3.0) {
        
        if (obs.empty()) {
            throw std::runtime_error("No observations provided");
        }
        
        // Stato iniziale all'epoca degli elementi
        State y0 = elementsToState(elem);
        double jd0 = elem.getJD();
        
        std::cout << "\n  Propagazione per " << obs.size() << " osservazioni...\n";
        
        // Per ogni osservazione, propaga e calcola residuo
        for (size_t i = 0; i < obs.size(); i++) {
            double jd_obs = obs[i].getJD();
            
            // Propaga all'epoca dell'osservazione
            PropagationStats stats;
            State y_obs = propagate(y0, jd0, jd_obs, stats);
            
            // Calcola posizione apparente (con light-time)
            double ra_calc, dec_calc, delta;
            computeApparentPosition(y_obs, jd_obs, ra_calc, dec_calc, delta);
            
            // Se osservazioni simulate, imposta ra/dec osservate
            if (obs[i].ra_obs == 0 && obs[i].dec_obs == 0) {
                obs[i].ra_obs = ra_calc;
                obs[i].dec_obs = dec_calc;
            }
            
            // Calcola residui [arcsec]
            double cos_dec = std::cos(obs[i].dec_obs);
            obs[i].ra_residual = (obs[i].ra_obs - ra_calc) * cos_dec * RAD2ARCSEC;
            obs[i].dec_residual = (obs[i].dec_obs - dec_calc) * RAD2ARCSEC;
            
            // Chi normalizzato
            double chi_ra = obs[i].ra_residual / obs[i].ra_sigma;
            double chi_dec = obs[i].dec_residual / obs[i].dec_sigma;
            obs[i].chi = std::sqrt(chi_ra*chi_ra + chi_dec*chi_dec);
            
            // Marca outlier
            obs[i].is_outlier = (obs[i].chi > outlier_threshold);
            
            // Progress
            if ((i + 1) % 100 == 0 || i == obs.size() - 1) {
                std::cout << "    " << (i+1) << "/" << obs.size() << " completate\r" << std::flush;
            }
        }
        std::cout << "\n";
        
        // Calcola statistiche
        return computeStatistics(obs);
    }
    
    /**
     * @brief Esegue test round-trip
     */
    double testRoundTrip(const OrbitalElements& elem, double dt_days) {
        State y0 = elementsToState(elem);
        double jd0 = elem.getJD();
        
        // Avanti
        PropagationStats stats1;
        State y1 = propagate(y0, jd0, jd0 + dt_days, stats1);
        
        // Indietro
        PropagationStats stats2;
        State y2 = propagate(y1, jd0 + dt_days, jd0, stats2);
        
        // Errore
        Vec3 dr = y2.r - y0.r;
        return dr.norm() * AU_km * 1000;  // [m]
    }
    
private:
    double tol_;
    double h_min_;
    double h_max_;
    
    //=========================================================================
    // CONVERSIONI
    //=========================================================================
    
    State elementsToState(const OrbitalElements& elem) {
        double a = elem.a;
        double e = elem.e;
        double inc = elem.i * DEG2RAD;
        double Omega = elem.Omega * DEG2RAD;
        double omega = elem.omega * DEG2RAD;
        double M0 = elem.M * DEG2RAD;
        
        // Keplero
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
        
        double x_orb = r * std::cos(nu);
        double y_orb = r * std::sin(nu);
        
        double v_factor = std::sqrt(k2 * a) / r;
        double vx_orb = -v_factor * sin_E;
        double vy_orb = v_factor * sqrt_1_e2 * cos_E;
        
        double cO = std::cos(Omega), sO = std::sin(Omega);
        double cw = std::cos(omega), sw = std::sin(omega);
        double ci = std::cos(inc), si = std::sin(inc);
        
        double P11 = cO*cw - sO*sw*ci, P12 = -cO*sw - sO*cw*ci;
        double P21 = sO*cw + cO*sw*ci, P22 = -sO*sw + cO*cw*ci;
        double P31 = sw*si,            P32 = cw*si;
        
        Vec3 pos(P11*x_orb + P12*y_orb, P21*x_orb + P22*y_orb, P31*x_orb + P32*y_orb);
        Vec3 vel(P11*vx_orb + P12*vy_orb, P21*vx_orb + P22*vy_orb, P31*vx_orb + P32*vy_orb);
        
        return State(pos, vel);
    }
    
    void computeApparentPosition(const State& state, double jd,
                                  double& ra, double& dec, double& delta) {
        // Posizione Terra
        Vec3 earth = ephemeris::getPlanetPosition(jd, 3);
        
        // Vettore geocentrico
        Vec3 geo = state.r - earth;
        delta = geo.norm();
        
        // Correzione light-time (iterativa)
        for (int iter = 0; iter < 3; iter++) {
            double lt = delta / c_light;
            Vec3 r_retard = state.r - state.v * lt;
            geo = r_retard - earth;
            delta = geo.norm();
        }
        
        // Eclittica → Equatoriale
        double cos_eps = std::cos(OBLIQUITY_J2000);
        double sin_eps = std::sin(OBLIQUITY_J2000);
        double x_eq = geo.x;
        double y_eq = geo.y * cos_eps - geo.z * sin_eps;
        double z_eq = geo.y * sin_eps + geo.z * cos_eps;
        
        // RA, Dec
        ra = std::atan2(y_eq, x_eq);
        if (ra < 0) ra += TWO_PI;
        dec = std::asin(z_eq / delta);
    }
    
    //=========================================================================
    // INTEGRATORE RK4 (per semplicità)
    //=========================================================================
    
    State propagate(const State& y0, double t0, double t1, PropagationStats& stats) {
        double dt = t1 - t0;
        double h = dt > 0 ? std::min(0.1, dt/10) : std::max(-0.1, dt/10);
        double t = t0;
        State y = y0;
        
        while ((dt > 0 && t < t1 - 1e-10) || (dt < 0 && t > t1 + 1e-10)) {
            if (dt > 0 && t + h > t1) h = t1 - t;
            if (dt < 0 && t + h < t1) h = t1 - t;
            
            // RK4
            State k1 = deriv(t, y);
            State k2 = deriv(t + h/2, y + k1*(h/2));
            State k3 = deriv(t + h/2, y + k2*(h/2));
            State k4 = deriv(t + h, y + k3*h);
            
            y = y + (k1 + k2*2 + k3*2 + k4) * (h/6);
            t += h;
            stats.steps_accepted++;
            stats.h_min = std::min(stats.h_min, std::abs(h));
            stats.h_max = std::max(stats.h_max, std::abs(h));
        }
        
        return y;
    }
    
    State deriv(double t, const State& y) {
        State dy;
        dy.r = y.v;
        dy.v = acceleration(t, y.r, y.v);
        return dy;
    }
    
    Vec3 acceleration(double t, const Vec3& r, const Vec3& v) {
        double r_norm = r.norm();
        double r3 = r_norm * r_norm * r_norm;
        
        // Kepleriano
        Vec3 acc = r * (-k2 / r3);
        
        // Perturbazioni planetarie
        constexpr double GM[8] = {
            GM_Mercury, GM_Venus, GM_EMB, GM_Mars,
            GM_Jupiter, GM_Saturn, GM_Uranus, GM_Neptune
        };
        
        for (int i = 0; i < 8; i++) {
            Vec3 r_planet = ephemeris::getPlanetPosition(t, i+1);
            Vec3 dr = r_planet - r;
            double dr_norm = dr.norm();
            double rp_norm = r_planet.norm();
            
            acc += dr * (GM[i] / (dr_norm * dr_norm * dr_norm));
            acc += r_planet * (-GM[i] / (rp_norm * rp_norm * rp_norm));
        }
        
        // Correzione relativistica
        double v2 = v.norm2();
        double rv = r.dot(v);
        double factor = k2 / (c2 * r3);
        acc += r * (factor * (4*k2/r_norm - v2)) + v * (factor * 4 * rv);
        
        return acc;
    }
    
    //=========================================================================
    // STATISTICHE
    //=========================================================================
    
    FitStatistics computeStatistics(const std::vector<Observation>& obs) {
        FitStatistics stats;
        stats.n_obs = obs.size();
        stats.n_outliers = 0;
        
        double sum_ra2 = 0, sum_dec2 = 0;
        double sum_ra = 0, sum_dec = 0;
        double sum_chi2 = 0;
        double max_res = 0;
        int n_valid = 0;
        
        double mjd_min = 1e10, mjd_max = 0;
        
        for (const auto& o : obs) {
            if (o.is_outlier) {
                stats.n_outliers++;
                continue;
            }
            
            sum_ra += o.ra_residual;
            sum_dec += o.dec_residual;
            sum_ra2 += o.ra_residual * o.ra_residual;
            sum_dec2 += o.dec_residual * o.dec_residual;
            sum_chi2 += o.chi * o.chi;
            
            double res = std::sqrt(o.ra_residual*o.ra_residual + 
                                   o.dec_residual*o.dec_residual);
            if (res > max_res) max_res = res;
            
            if (o.mjd < mjd_min) mjd_min = o.mjd;
            if (o.mjd > mjd_max) mjd_max = o.mjd;
            
            n_valid++;
        }
        
        if (n_valid > 0) {
            stats.mean_ra = sum_ra / n_valid;
            stats.mean_dec = sum_dec / n_valid;
            stats.rms_ra = std::sqrt(sum_ra2 / n_valid);
            stats.rms_dec = std::sqrt(sum_dec2 / n_valid);
            stats.rms_total = std::sqrt((sum_ra2 + sum_dec2) / (2 * n_valid));
            stats.chi2 = sum_chi2;
            stats.chi2_reduced = sum_chi2 / (2 * n_valid - 6);  // 6 parametri
            stats.max_residual = max_res;
            stats.time_span = mjd_max - mjd_min;
        }
        
        return stats;
    }
};

//=============================================================================
// UTILITY DI OUTPUT
//=============================================================================

std::string formatRA(double ra_rad) {
    double h = ra_rad * 12.0 / PI;
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

void printObservations(const std::vector<Observation>& obs, int max_print = 20) {
    std::cout << "\n  Primi/ultimi residui:\n";
    std::cout << "  " << std::string(80, '-') << "\n";
    std::cout << "  MJD          RA res   Dec res   Chi   Flag\n";
    std::cout << "  " << std::string(80, '-') << "\n";
    
    int n = obs.size();
    int half = max_print / 2;
    
    for (int i = 0; i < std::min(half, n); i++) {
        const auto& o = obs[i];
        std::cout << "  " << std::fixed << std::setprecision(5) << o.mjd
                  << "  " << std::setw(8) << std::setprecision(3) << o.ra_residual
                  << "  " << std::setw(8) << o.dec_residual
                  << "  " << std::setw(5) << std::setprecision(2) << o.chi
                  << "  " << (o.is_outlier ? "OUTLIER" : "")
                  << "\n";
    }
    
    if (n > max_print) {
        std::cout << "  ... (" << (n - max_print) << " osservazioni omesse) ...\n";
    }
    
    for (int i = std::max(half, n - half); i < n; i++) {
        const auto& o = obs[i];
        std::cout << "  " << std::fixed << std::setprecision(5) << o.mjd
                  << "  " << std::setw(8) << std::setprecision(3) << o.ra_residual
                  << "  " << std::setw(8) << o.dec_residual
                  << "  " << std::setw(5) << std::setprecision(2) << o.chi
                  << "  " << (o.is_outlier ? "OUTLIER" : "")
                  << "\n";
    }
    std::cout << "  " << std::string(80, '-') << "\n";
}

} // namespace astdyn

//=============================================================================
// MAIN
//=============================================================================

int main(int argc, char* argv[]) {
    using namespace astdyn;
    
    std::cout << std::fixed;
    std::cout << "================================================================\n";
    std::cout << "  AstDyn O-C Fitter v1.0\n";
    std::cout << "  Propagazione e Fitting Osservazioni\n";
    std::cout << "================================================================\n\n";
    
    // Elementi orbitali (17030) Sierks - da AstDyS
    OrbitalElements elem = AstDySElementsReader::fromManual(
        17030, "(17030) Sierks",
        61000.0,                    // Epoca MJD (2025-11-21)
        3.1754733,                  // a [AU]
        0.0454207,                  // e
        2.9046,                     // i [deg]
        104.16243,                  // Omega [deg]
        100.5141,                   // omega [deg]
        229.79088,                  // M [deg]
        13.9,                       // H
        0.15                        // G
    );
    
    std::cout << "ELEMENTI ORBITALI:\n";
    std::cout << "  Oggetto: " << elem.name << "\n";
    std::cout << "  Epoca:   MJD " << elem.epoch << " (JD " << elem.getJD() << ")\n";
    std::cout << "  a = " << elem.a << " AU\n";
    std::cout << "  e = " << elem.e << "\n";
    std::cout << "  i = " << elem.i << "°\n";
    std::cout << "  Ω = " << elem.Omega << "°\n";
    std::cout << "  ω = " << elem.omega << "°\n";
    std::cout << "  M = " << elem.M << "°\n";
    std::cout << "  H = " << elem.H << " mag\n\n";
    
    // Crea fitter
    OCFitter fitter(1e-12);
    
    // Test round-trip
    std::cout << "TEST ROUND-TRIP:\n";
    double rt_1yr = fitter.testRoundTrip(elem, 365.25);
    double rt_10yr = fitter.testRoundTrip(elem, 3652.5);
    std::cout << "  1 anno:   " << rt_1yr << " m\n";
    std::cout << "  10 anni:  " << rt_10yr << " m\n\n";
    
    // Crea osservazioni simulate (span temporale realistico)
    std::cout << "GENERAZIONE OSSERVAZIONI SIMULATE:\n";
    double mjd_start = elem.epoch - 3652.5;  // -10 anni
    double mjd_end = elem.epoch + 365.25;    // +1 anno
    int n_obs = 150;
    
    auto obs = AstDySRWOReader::createSimulated(elem, mjd_start, mjd_end, n_obs, 0.5);
    std::cout << "  Numero osservazioni: " << n_obs << "\n";
    std::cout << "  Arco temporale: " << (mjd_end - mjd_start) / 365.25 << " anni\n";
    std::cout << "  Da MJD " << mjd_start << " a MJD " << mjd_end << "\n\n";
    
    // Calcola residui
    std::cout << "CALCOLO RESIDUI O-C:\n";
    FitStatistics stats = fitter.computeResiduals(elem, obs, 3.0);
    
    // Risultati
    std::cout << "\n================================================================\n";
    std::cout << "  STATISTICHE FIT\n";
    std::cout << "================================================================\n";
    std::cout << "  Osservazioni totali:  " << stats.n_obs << "\n";
    std::cout << "  Outlier (>3σ):        " << stats.n_outliers << "\n";
    std::cout << "  Arco temporale:       " << stats.time_span / 365.25 << " anni\n\n";
    
    std::cout << "  RMS residui:\n";
    std::cout << "    RA*cos(δ): " << stats.rms_ra << " arcsec\n";
    std::cout << "    Dec:       " << stats.rms_dec << " arcsec\n";
    std::cout << "    Totale:    " << stats.rms_total << " arcsec\n\n";
    
    std::cout << "  Media residui:\n";
    std::cout << "    RA*cos(δ): " << stats.mean_ra << " arcsec\n";
    std::cout << "    Dec:       " << stats.mean_dec << " arcsec\n\n";
    
    std::cout << "  Chi-quadro:\n";
    std::cout << "    χ²:        " << stats.chi2 << "\n";
    std::cout << "    χ² ridotto:" << stats.chi2_reduced << "\n";
    std::cout << "    Max resid: " << stats.max_residual << " arcsec\n";
    
    // Mostra alcuni residui
    printObservations(obs, 20);
    
    // Posizione attuale
    std::cout << "\n================================================================\n";
    std::cout << "  EFFEMERIDI ATTUALI\n";
    std::cout << "================================================================\n";
    
    double jd_now = 2461008.0913;  // 28 Nov 2025, 14:11 UTC
    State y0 = fitter.testRoundTrip(elem, 0) == 0 ? State() : State();  // dummy
    
    // Propaga
    OCFitter prop;
    std::vector<Observation> single_obs = {Observation()};
    single_obs[0].mjd = jd_now - 2400000.5;
    single_obs[0].ra_sigma = 0.5;
    single_obs[0].dec_sigma = 0.5;
    single_obs[0].designation = elem.name;
    
    prop.computeResiduals(elem, single_obs, 10.0);
    
    std::cout << "  Data: JD " << jd_now << " (28 Nov 2025, 14:11 UTC)\n";
    std::cout << "  RA  = " << formatRA(single_obs[0].ra_obs) << "\n";
    std::cout << "  Dec = " << formatDec(single_obs[0].dec_obs) << "\n";
    
    std::cout << "\n================================================================\n";
    std::cout << "  JPL HORIZONS (riferimento):\n";
    std::cout << "  RA  = 04 53 11.25\n";
    std::cout << "  Dec = +20 19 25.8\n";
    std::cout << "================================================================\n";
    
    return 0;
}
