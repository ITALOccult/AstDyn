#include <iostream>
#include <iomanip>
#include <cmath>

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;

// Versione propagator
namespace v1 {
    void getEarth(double jd, double& x, double& y, double& z) {
        double T = (jd - 2451545.0) / 36525.0;
        double a = 1.00000261 + 0.00000562*T;
        double e = 0.01671123 - 0.00004392*T;
        double I = -0.00001531 - 0.01294668*T;
        double L = 100.46457166 + 35999.37244981*T;
        double omega_bar = 102.93768193 + 0.32327364*T;
        double Omega = 0.0;
        
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
        
        double x_ecl = (cO*cw - sO*sw*ci) * x_orb + (-cO*sw - sO*cw*ci) * y_orb;
        double y_ecl = (sO*cw + cO*sw*ci) * x_orb + (-sO*sw + cO*cw*ci) * y_orb;
        double z_ecl = (sw*si) * x_orb + (cw*si) * y_orb;
        
        double eps = 23.4392911 * DEG2RAD;
        x = x_ecl;
        y = std::cos(eps)*y_ecl - std::sin(eps)*z_ecl;
        z = std::sin(eps)*y_ecl + std::cos(eps)*z_ecl;
    }
}

// Versione test_sierks_jpl
namespace v2 {
    void getEarth(double jd, double& x, double& y, double& z) {
        double T = (jd - 2451545.0) / 36525.0;
        double a = 1.00000261;  // Senza variazione temporale!
        double e = 0.01671123;
        double I = 0;  // Inclinazione 0
        double L = 100.46457166 + 35999.37244981*T;
        double omega_bar = 102.93768193;  // Senza variazione temporale!
        double Omega = 0;
        
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
        
        double cO = std::cos(Omega * DEG2RAD), sO = std::sin(Omega * DEG2RAD);
        double ci = std::cos(I * DEG2RAD), si = std::sin(I * DEG2RAD);
        double cw = std::cos(omega), sw = std::sin(omega);
        
        double x_ecl = (cO*cw - sO*sw*ci) * x_orb + (-cO*sw - sO*cw*ci) * y_orb;
        double y_ecl = (sO*cw + cO*sw*ci) * x_orb + (-sO*sw + cO*cw*ci) * y_orb;
        double z_ecl = (sw*si) * x_orb + (cw*si) * y_orb;
        
        double eps = 23.4392911 * DEG2RAD;
        x = x_ecl;
        y = std::cos(eps)*y_ecl - std::sin(eps)*z_ecl;
        z = std::sin(eps)*y_ecl + std::cos(eps)*z_ecl;
    }
}

int main() {
    double jd = 2461008.0913;
    double x1, y1, z1, x2, y2, z2;
    
    v1::getEarth(jd, x1, y1, z1);
    v2::getEarth(jd, x2, y2, z2);
    
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Propagator: [" << x1 << ", " << y1 << ", " << z1 << "]\n";
    std::cout << "Test JPL:   [" << x2 << ", " << y2 << ", " << z2 << "]\n";
    
    double dx = x1-x2, dy = y1-y2, dz = z1-z2;
    double diff = std::sqrt(dx*dx + dy*dy + dz*dz);
    std::cout << "Differenza: " << diff << " AU = " << diff * 149597870.7 << " km\n";
    
    return 0;
}
