/**
 * Confronta propagazione propagator vs test_sierks_jpl
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double k = 0.01720209895;
constexpr double k2 = k * k;
constexpr double GM_SUN = k2;
constexpr double AU_KM = 149597870.7;
constexpr double DAY_SEC = 86400.0;

constexpr double GM_MERCURY = 4.9125474514508118e-11;
constexpr double GM_VENUS   = 7.2434524861627027e-10;
constexpr double GM_EARTH   = 8.9970116036316091e-10;
constexpr double GM_MARS    = 9.5495351057792580e-11;
constexpr double GM_JUPITER = 2.8253458420837436e-07;
constexpr double GM_SATURN  = 8.4597151856806587e-08;
constexpr double GM_URANUS  = 1.2920249167819693e-08;
constexpr double GM_NEPTUNE = 1.5243589007842762e-08;

using Vec3 = std::array<double, 3>;
using State = std::array<double, 6>;

Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
Vec3 operator*(double s, const Vec3& v) { return {s*v[0], s*v[1], s*v[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }

// Effemeridi pianeti ICRF (versione test)
Vec3 getPlanetICRF_test(double jd, int planet) {
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
    double x_orb = r * std::cos(nu);
    double y_orb = r * std::sin(nu);
    
    double cO = std::cos(Omega_rad), sO = std::sin(Omega_rad);
    double ci = std::cos(i_rad), si = std::sin(i_rad);
    double cw = std::cos(omega), sw = std::sin(omega);
    
    double x = (cO*cw - sO*sw*ci) * x_orb + (-cO*sw - sO*cw*ci) * y_orb;
    double y = (sO*cw + cO*sw*ci) * x_orb + (-sO*sw + cO*cw*ci) * y_orb;
    double z = (sw*si) * x_orb + (cw*si) * y_orb;
    
    double eps = 23.4392911 * DEG2RAD;
    return {x, std::cos(eps)*y - std::sin(eps)*z, std::sin(eps)*y + std::cos(eps)*z};
}

// Accelerazione (Sole + 8 pianeti) - identica al test
State deriv(double t, const State& s) {
    Vec3 r = {s[0], s[1], s[2]};
    Vec3 v = {s[3], s[4], s[5]};
    double r3 = std::pow(norm(r), 3);
    
    Vec3 acc = {-GM_SUN*r[0]/r3, -GM_SUN*r[1]/r3, -GM_SUN*r[2]/r3};
    
    double GM[] = {0, GM_MERCURY, GM_VENUS, GM_EARTH, GM_MARS, GM_JUPITER, GM_SATURN, GM_URANUS, GM_NEPTUNE};
    for (int p = 1; p <= 8; p++) {
        Vec3 rp = getPlanetICRF_test(t, p);
        Vec3 d = r - rp;
        double d3 = std::pow(norm(d), 3);
        double rp3 = std::pow(norm(rp), 3);
        acc[0] -= GM[p] * (d[0]/d3 + rp[0]/rp3);
        acc[1] -= GM[p] * (d[1]/d3 + rp[1]/rp3);
        acc[2] -= GM[p] * (d[2]/d3 + rp[2]/rp3);
    }
    
    return {v[0], v[1], v[2], acc[0], acc[1], acc[2]};
}

// RK4 
State rk4_step(double t, const State& y, double h) {
    State k1 = deriv(t, y);
    State y1; for(int i=0; i<6; i++) y1[i] = y[i] + 0.5*h*k1[i];
    State k2 = deriv(t + 0.5*h, y1);
    State y2; for(int i=0; i<6; i++) y2[i] = y[i] + 0.5*h*k2[i];
    State k3 = deriv(t + 0.5*h, y2);
    State y3; for(int i=0; i<6; i++) y3[i] = y[i] + h*k3[i];
    State k4 = deriv(t + h, y3);
    
    State yn;
    for(int i=0; i<6; i++) yn[i] = y[i] + h*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
    return yn;
}

int main() {
    std::cout << std::fixed << std::setprecision(10);
    
    // Stato iniziale da elementi Sierks convertiti a ICRF
    // (stessi valori che usa il propagator)
    double jd0 = 2461000.5;
    State y0 = {
        1.0827162342,  // X
        2.8683018724,  // Y
        1.1437412744,  // Z
        -0.0089234884, // VX
        0.0024152608,  // VY
        0.0014876383   // VZ
    };
    
    double jd_target = 2461008.0913;
    
    // Propaga con RK4
    double dt = 0.01;
    double t = jd0;
    State y = y0;
    
    while (t < jd_target) {
        if (t + dt > jd_target) dt = jd_target - t;
        y = rk4_step(t, y, dt);
        t += dt;
    }
    
    std::cout << "RISULTATO PROPAGAZIONE (RK4, 8 pianeti):\n";
    std::cout << "  r = [" << y[0] << ", " << y[1] << ", " << y[2] << "]\n\n";
    
    std::cout << "RISULTATO ASTDYN_PROPAGATOR:\n";
    std::cout << "  r = [1.022805, 2.883920, 1.153479]\n\n";
    
    double dx = y[0] - 1.022805;
    double dy = y[1] - 2.883920;
    double dz = y[2] - 1.153479;
    double diff = std::sqrt(dx*dx + dy*dy + dz*dz);
    std::cout << "DIFFERENZA: " << diff << " AU = " << diff * AU_KM << " km\n";
    
    return 0;
}
