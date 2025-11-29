/**
 * @file test_astdys_full.cpp
 * @brief Test propagazione da elementi AstDyS (equinoziali)
 * 
 * Elementi AstDyS per (11234):
 * Format: OEF2.0 Equinoctial, ECLM J2000
 * Epoca: MJD 61000.0 TDT (JD 2461000.5)
 * 
 * Elementi equinoziali:
 *   a = 2.6808535916678031 AU
 *   h = e*sin(LP) = 0.032872036471001    (LP = longitude of perihelion)
 *   k = e*cos(LP) = 0.036254405825130
 *   p = tan(i/2)*sin(LN) = 0.103391596538937  (LN = longitude of node)
 *   q = tan(i/2)*cos(LN) = -0.042907901689093
 *   λ = mean longitude = 235.8395861037268°
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <vector>

// ============================================================================
// COSTANTI
// ============================================================================

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double AU_KM = 149597870.7;

// GM in AU³/day²
constexpr double GM_SUN     = 2.9591220828559093e-04;
constexpr double GM_MERCURY = 4.9125474514508118e-11;
constexpr double GM_VENUS   = 7.2434524861627027e-10;
constexpr double GM_EARTH   = 8.8876925870231834e-10;
constexpr double GM_MARS    = 9.5495351057792580e-11;
constexpr double GM_JUPITER = 2.8253458420837619e-07;
constexpr double GM_SATURN  = 8.4597151856806587e-08;
constexpr double GM_URANUS  = 1.2920249167819693e-08;
constexpr double GM_NEPTUNE = 1.5243589007842762e-08;

// Elementi equinoziali AstDyS per (11234)
constexpr double ASTDYS_EPOCH_MJD = 61000.0;
constexpr double ASTDYS_EPOCH_JD = 2461000.5;  // MJD + 2400000.5

constexpr double ASTDYS_A      = 2.6808535916678031;    // AU
constexpr double ASTDYS_H      = 0.032872036471001;     // e*sin(LP)
constexpr double ASTDYS_K      = 0.036254405825130;     // e*cos(LP)
constexpr double ASTDYS_P      = 0.103391596538937;     // tan(i/2)*sin(LN)
constexpr double ASTDYS_Q      = -0.042907901689093;    // tan(i/2)*cos(LN)
constexpr double ASTDYS_LAMBDA = 235.8395861037268;     // mean longitude (deg)

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
// CONVERSIONE ELEMENTI EQUINOZIALI -> STATO
// ============================================================================

/**
 * Converte elementi equinoziali a elementi kepleriani classici
 * 
 * Input (equinoziali AstDyS):
 *   a = semiasse maggiore
 *   h = e*sin(ϖ)  dove ϖ = Ω + ω (longitude of perihelion)
 *   k = e*cos(ϖ)
 *   p = tan(i/2)*sin(Ω)
 *   q = tan(i/2)*cos(Ω)
 *   λ = mean longitude = Ω + ω + M
 * 
 * Output (kepleriani):
 *   a, e, i, Ω, ω, M
 */
void equinoctialToKeplerian(double a, double h, double k, double p, double q, double lambda,
                             double& e, double& i, double& Omega, double& omega, double& M) {
    // Eccentricità
    e = std::sqrt(h*h + k*k);
    
    // Inclinazione
    double tan_i_2 = std::sqrt(p*p + q*q);
    i = 2.0 * std::atan(tan_i_2);
    
    // Longitudine del nodo ascendente
    Omega = std::atan2(p, q);
    if (Omega < 0) Omega += 2*PI;
    
    // Longitudine del perielio
    double LP = std::atan2(h, k);  // ϖ = Ω + ω
    if (LP < 0) LP += 2*PI;
    
    // Argomento del perielio
    omega = LP - Omega;
    while (omega < 0) omega += 2*PI;
    while (omega > 2*PI) omega -= 2*PI;
    
    // Anomalia media
    M = lambda * DEG2RAD - LP;
    while (M < 0) M += 2*PI;
    while (M > 2*PI) M -= 2*PI;
}

double solveKepler(double M, double e) {
    double E = M;
    for (int i = 0; i < 30; ++i) {
        double dE = (M - E + e * std::sin(E)) / (1.0 - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < 1e-14) break;
    }
    return E;
}

/**
 * Elementi kepleriani (eclittica) -> Stato cartesiano ICRF (equatoriale)
 */
State keplerToStateICRF(double a, double e, double i, double Omega, double omega, double M) {
    // Risolvi Keplero
    double E = solveKepler(M, e);
    
    // Posizione nel piano orbitale
    double cosE = std::cos(E);
    double sinE = std::sin(E);
    double sqrtae = std::sqrt(1.0 - e*e);
    
    double x_orb = a * (cosE - e);
    double y_orb = a * sqrtae * sinE;
    
    // Velocità nel piano orbitale
    double n = std::sqrt(GM_SUN / (a*a*a));
    double r = a * (1.0 - e * cosE);
    double vx_orb = -a * n * sinE / (1.0 - e * cosE);
    double vy_orb = a * n * sqrtae * cosE / (1.0 - e * cosE);
    
    // Rotazione: piano orbitale -> eclittica
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
    
    double x_ecl = P1 * x_orb + P2 * y_orb;
    double y_ecl = Q1 * x_orb + Q2 * y_orb;
    double z_ecl = W1 * x_orb + W2 * y_orb;
    
    double vx_ecl = P1 * vx_orb + P2 * vy_orb;
    double vy_ecl = Q1 * vx_orb + Q2 * vy_orb;
    double vz_ecl = W1 * vx_orb + W2 * vy_orb;
    
    // Rotazione: eclittica -> equatoriale (ICRF)
    constexpr double eps = 23.4392911 * DEG2RAD;
    double c = std::cos(eps);
    double s = std::sin(eps);
    
    State state;
    state[0] = x_ecl;
    state[1] = c * y_ecl - s * z_ecl;
    state[2] = s * y_ecl + c * z_ecl;
    state[3] = vx_ecl;
    state[4] = c * vy_ecl - s * vz_ecl;
    state[5] = s * vy_ecl + c * vz_ecl;
    
    return state;
}

// ============================================================================
// EFFEMERIDI PLANETARIE
// ============================================================================

struct PlanetElements {
    double a0, a1;
    double e0, e1;
    double i0, i1;
    double Om0, Om1;
    double w0, w1;
    double L0, L1;
};

static const PlanetElements PLANETS[] = {
    {0.38709927, 0.00000037, 0.20563593, 0.00001906, 7.00497902, -0.00594749,
     48.33076593, -0.12534081, 77.45779628, 0.16047689, 252.25032350, 149472.67411175},
    {0.72333566, 0.00000390, 0.00677672, -0.00004107, 3.39467605, -0.00078890,
     76.67984255, -0.27769418, 131.60246718, 0.00268329, 181.97909950, 58517.81538729},
    {1.00000261, 0.00000562, 0.01671123, -0.00004392, -0.00001531, -0.01294668,
     0.0, 0.0, 102.93768193, 0.32327364, 100.46457166, 35999.37244981},
    {1.52371034, 0.00001847, 0.09339410, 0.00007882, 1.84969142, -0.00813131,
     49.55953891, -0.29257343, -23.94362959, 0.44441088, -4.55343205, 19140.30268499},
    {5.20288700, -0.00011607, 0.04838624, -0.00013253, 1.30439695, -0.00183714,
     100.47390909, 0.20469106, 14.72847983, 0.21252668, 34.39644051, 3034.74612775},
    {9.53667594, -0.00125060, 0.05386179, -0.00050991, 2.48599187, 0.00193609,
     113.66242448, -0.28867794, 92.59887831, -0.41897216, 49.95424423, 1222.49362201},
    {19.18916464, -0.00196176, 0.04725744, -0.00004397, 0.77263783, -0.00242939,
     74.01692503, 0.04240589, 170.95427630, 0.40805281, 313.23810451, 428.48202785},
    {30.06992276, 0.00026291, 0.00859048, 0.00005105, 1.77004347, 0.00035372,
     131.78422574, -0.00508664, 44.96476227, -0.32241464, -55.12002969, 218.45945325}
};

Vec3 getPlanetPositionICRF(int planet, double jd) {
    const double J2000 = 2451545.0;
    double T = (jd - J2000) / 36525.0;
    
    const auto& p = PLANETS[planet];
    
    double a = p.a0 + p.a1 * T;
    double e = p.e0 + p.e1 * T;
    double inc = (p.i0 + p.i1 * T) * DEG2RAD;
    double Om = (p.Om0 + p.Om1 * T) * DEG2RAD;
    double w_bar = (p.w0 + p.w1 * T) * DEG2RAD;
    double L = (p.L0 + p.L1 * T) * DEG2RAD;
    
    double omega = w_bar - Om;
    double M = L - w_bar;
    while (M < 0) M += 2*PI;
    while (M > 2*PI) M -= 2*PI;
    
    double E = solveKepler(M, e);
    
    double cosE = std::cos(E);
    double sinE = std::sin(E);
    double x_orb = a * (cosE - e);
    double y_orb = a * std::sqrt(1.0 - e*e) * sinE;
    
    double cosO = std::cos(Om);
    double sinO = std::sin(Om);
    double cosw = std::cos(omega);
    double sinw = std::sin(omega);
    double cosi = std::cos(inc);
    double sini = std::sin(inc);
    
    double P1 = cosO * cosw - sinO * sinw * cosi;
    double P2 = -cosO * sinw - sinO * cosw * cosi;
    double Q1 = sinO * cosw + cosO * sinw * cosi;
    double Q2 = -sinO * sinw + cosO * cosw * cosi;
    double W1 = sinw * sini;
    double W2 = cosw * sini;
    
    double x_ecl = P1 * x_orb + P2 * y_orb;
    double y_ecl = Q1 * x_orb + Q2 * y_orb;
    double z_ecl = W1 * x_orb + W2 * y_orb;
    
    constexpr double eps = 23.4392911 * DEG2RAD;
    double c = std::cos(eps);
    double s = std::sin(eps);
    
    return {x_ecl, c * y_ecl - s * z_ecl, s * y_ecl + c * z_ecl};
}

// ============================================================================
// ACCELERAZIONE
// ============================================================================

State computeAcceleration(const State& state, double jd) {
    Vec3 r = {state[0], state[1], state[2]};
    Vec3 v = {state[3], state[4], state[5]};
    
    double r_mag = norm(r);
    double r3 = r_mag * r_mag * r_mag;
    
    Vec3 acc = {-GM_SUN * r[0] / r3,
                -GM_SUN * r[1] / r3,
                -GM_SUN * r[2] / r3};
    
    static const double GM[] = {GM_MERCURY, GM_VENUS, GM_EARTH, GM_MARS,
                                 GM_JUPITER, GM_SATURN, GM_URANUS, GM_NEPTUNE};
    
    for (int i = 0; i < 8; ++i) {
        Vec3 rp = getPlanetPositionICRF(i, jd);
        Vec3 d = r - rp;
        
        double d_mag = norm(d);
        double rp_mag = norm(rp);
        
        double d3 = d_mag * d_mag * d_mag;
        double rp3 = rp_mag * rp_mag * rp_mag;
        
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
    static constexpr double c[13] = {
        0.0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 1.0/2.0,
        5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0, 1.0
    };
    
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
    
    static constexpr double b8[13] = {
        41.0/840.0, 0, 0, 0, 0, 34.0/105.0, 9.0/35.0, 9.0/35.0,
        9.0/280.0, 9.0/280.0, 41.0/840.0, 0, 0
    };
    
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
        
        for (int stage = 0; stage < 13; ++stage) {
            State y_stage = y;
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < stage; ++j) {
                    y_stage[i] += h * a[stage][j] * k[j][i];
                }
            }
            k[stage] = computeAcceleration(y_stage, t + c[stage] * h);
        }
        
        State y8 = y;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 13; ++j) {
                y8[i] += h * b8[j] * k[j][i];
            }
        }
        
        State y7 = y;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 13; ++j) {
                y7[i] += h * b7[j] * k[j][i];
            }
        }
        
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
        double h = (t_end - t0) / 100.0;
        
        if (h > h_max_) h = h_max_;
        if (std::abs(h) < h_min_) h = (h > 0) ? h_min_ : -h_min_;
        
        total_steps = 0;
        
        while ((h > 0 && t < t_end) || (h < 0 && t > t_end)) {
            if ((h > 0 && t + h > t_end) || (h < 0 && t + h < t_end)) {
                h = t_end - t;
            }
            
            double err;
            State y_new = step(y, t, h, err);
            
            if (err < tol_ || std::abs(h) <= h_min_) {
                y = y_new;
                t += h;
                total_steps++;
                
                if (err > 0) {
                    double factor = 0.9 * std::pow(tol_ / err, 1.0/8.0);
                    factor = std::max(0.1, std::min(4.0, factor));
                    h *= factor;
                }
            } else {
                h *= 0.5;
            }
            
            if (std::abs(h) > h_max_) h = (h > 0) ? h_max_ : -h_max_;
            if (std::abs(h) < h_min_) h = (h > 0) ? h_min_ : -h_min_;
        }
        
        return y;
    }
};

// ============================================================================
// CONVERSIONE RA/DEC
// ============================================================================

void stateToRADec(const State& s, double& ra_h, double& ra_m, double& ra_s,
                  int& dec_d, int& dec_m, double& dec_s) {
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
    dec_d = sign * static_cast<int>(std::floor(dec_deg));
    dec_m = static_cast<int>(std::floor((dec_deg - std::abs(dec_d)) * 60.0));
    dec_s = ((dec_deg - std::abs(dec_d)) * 60.0 - dec_m) * 60.0;
}

double computeDeltaArcsec(double ra1_h, double ra1_m, double ra1_s,
                          int dec1_d, int dec1_m, double dec1_s,
                          double ra2_h, double ra2_m, double ra2_s,
                          int dec2_d, int dec2_m, double dec2_s) {
    double ra1 = (ra1_h + ra1_m/60.0 + ra1_s/3600.0) * 15.0;
    double ra2 = (ra2_h + ra2_m/60.0 + ra2_s/3600.0) * 15.0;
    
    double dec1 = std::abs(dec1_d) + dec1_m/60.0 + dec1_s/3600.0;
    if (dec1_d < 0) dec1 = -dec1;
    double dec2 = std::abs(dec2_d) + dec2_m/60.0 + dec2_s/3600.0;
    if (dec2_d < 0) dec2 = -dec2;
    
    double dra = (ra1 - ra2) * 3600.0;
    double ddec = (dec1 - dec2) * 3600.0;
    
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
    std::cout << "  Test Propagazione da Elementi AstDyS\n";
    std::cout << "  Asteroide (11234) 1999 JS82\n";
    std::cout << "  Formato: OEF2.0 Equinoctial, ECLM J2000\n";
    std::cout << "================================================================\n\n";
    
    // 1. Mostra elementi equinoziali
    std::cout << "1. ELEMENTI EQUINOZIALI AstDyS\n";
    std::cout << "   Epoca: MJD " << std::setprecision(1) << ASTDYS_EPOCH_MJD 
              << " (JD " << ASTDYS_EPOCH_JD << " = 2025-Nov-21)\n";
    std::cout << std::setprecision(12);
    std::cout << "   a = " << ASTDYS_A << " AU\n";
    std::cout << "   h = e*sin(ϖ) = " << ASTDYS_H << "\n";
    std::cout << "   k = e*cos(ϖ) = " << ASTDYS_K << "\n";
    std::cout << "   p = tan(i/2)*sin(Ω) = " << ASTDYS_P << "\n";
    std::cout << "   q = tan(i/2)*cos(Ω) = " << ASTDYS_Q << "\n";
    std::cout << "   λ = " << ASTDYS_LAMBDA << "°\n\n";
    
    // 2. Converti a kepleriani
    double e, i, Omega, omega, M;
    equinoctialToKeplerian(ASTDYS_A, ASTDYS_H, ASTDYS_K, ASTDYS_P, ASTDYS_Q, ASTDYS_LAMBDA,
                           e, i, Omega, omega, M);
    
    std::cout << "2. ELEMENTI KEPLERIANI (conversione)\n";
    std::cout << std::setprecision(10);
    std::cout << "   a = " << ASTDYS_A << " AU\n";
    std::cout << "   e = " << e << "\n";
    std::cout << "   i = " << i * RAD2DEG << "°\n";
    std::cout << "   Ω = " << Omega * RAD2DEG << "°\n";
    std::cout << "   ω = " << omega * RAD2DEG << "°\n";
    std::cout << "   M = " << M * RAD2DEG << "°\n\n";
    
    // 3. Converti a stato ICRF
    State s0 = keplerToStateICRF(ASTDYS_A, e, i, Omega, omega, M);
    
    std::cout << "3. STATO ICRF ALL'EPOCA\n";
    std::cout << std::setprecision(12);
    std::cout << "   r = [" << s0[0] << ", " << s0[1] << ", " << s0[2] << "] AU\n";
    std::cout << "   v = [" << s0[3] << ", " << s0[4] << ", " << s0[5] << "] AU/day\n";
    std::cout << "   |r| = " << norm({s0[0], s0[1], s0[2]}) << " AU\n\n";
    
    // 4. Propagazione a epoche di test
    struct Target {
        double jd;
        const char* date;
        double jpl_ra_h, jpl_ra_m, jpl_ra_s;
        int jpl_dec_d, jpl_dec_m;
        double jpl_dec_s;
    };
    
    std::vector<Target> targets = {
        {2460970.5, "2025-Oct-22", 15, 18, 51.75, -6, 25, 31.9},
        {2461030.5, "2025-Dec-21", 16, 5, 27.35, -10, 38, 8.1}
    };
    
    RKF78Integrator integrator(1e-12, 0.01, 5.0);
    
    std::cout << "4. PROPAGAZIONE E CONFRONTO CON JPL\n\n";
    
    for (const auto& t : targets) {
        double dt = t.jd - ASTDYS_EPOCH_JD;
        
        std::cout << "================================================================\n";
        std::cout << "  " << t.date << " (Δt = " << std::setprecision(0) << dt << " giorni)\n";
        std::cout << "================================================================\n";
        
        int steps;
        State s_prop = integrator.integrate(s0, ASTDYS_EPOCH_JD, t.jd, steps);
        
        std::cout << std::setprecision(9);
        std::cout << "   r = [" << s_prop[0] << ", " << s_prop[1] << ", " << s_prop[2] << "] AU\n";
        std::cout << "   |r| = " << norm({s_prop[0], s_prop[1], s_prop[2]}) << " AU\n";
        std::cout << "   Passi RKF78: " << steps << "\n";
        
        double ra_h, ra_m, ra_s;
        int dec_d, dec_m;
        double dec_s;
        stateToRADec(s_prop, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s);
        
        std::cout << "\n   AstDyn (da AstDyS):\n";
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
            std::cout << " ◯ DA VERIFICARE\n";
        }
        std::cout << std::endl;
    }
    
    // Test round-trip
    std::cout << "================================================================\n";
    std::cout << "  TEST ROUND-TRIP\n";
    std::cout << "================================================================\n";
    
    double jd_far = 2461030.5;
    int steps1, steps2;
    State s_far = integrator.integrate(s0, ASTDYS_EPOCH_JD, jd_far, steps1);
    State s_back = integrator.integrate(s_far, jd_far, ASTDYS_EPOCH_JD, steps2);
    
    double dr = norm({s_back[0]-s0[0], s_back[1]-s0[1], s_back[2]-s0[2]});
    
    std::cout << std::setprecision(6);
    std::cout << "   Errore round-trip (" << steps1+steps2 << " passi): " 
              << dr * AU_KM * 1000.0 << " m\n";
    
    std::cout << "\n================================================================\n";
    std::cout << "  RIEPILOGO\n";
    std::cout << "================================================================\n";
    std::cout << "  Sorgente elementi: AstDyS (OEF2.0 Equinoctial)\n";
    std::cout << "  Epoca elementi: MJD 61000.0 (2025-Nov-21)\n";
    std::cout << "  Frame: ECLM J2000 -> ICRF\n";
    std::cout << "  Perturbazioni: 8 pianeti\n";
    std::cout << "  Integratore: RKF78 adattivo\n";
    std::cout << "================================================================\n";
    
    return 0;
}
