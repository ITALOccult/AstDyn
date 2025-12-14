/**
 * @file orbital_conversions.cpp
 * @brief Implementazione conversioni orbitali
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 */

#include "orbital_conversions.h"
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace ioccultcalc {

// ============================================================================
// Conversione Equinoziale → Kepleriano
// ============================================================================

KeplerianElements OrbitalConversions::equinoctialToKeplerian(
    const EquinoctialElements& eq) {
    
    KeplerianElements kep;
    kep.name = eq.name;
    kep.a = eq.a;
    kep.epoch_jd = eq.getEpochJD();
    
    // 1. Eccentricità: e = sqrt(h² + k²)
    kep.e = std::sqrt(eq.h*eq.h + eq.k*eq.k);
    
    // 2. Inclinazione: i = 2·atan(sqrt(p² + q²))
    double tan_i_2 = std::sqrt(eq.p*eq.p + eq.q*eq.q);
    kep.i = 2.0 * std::atan(tan_i_2);
    
    // 3. Longitudine nodo ascendente: Ω = atan2(p, q)
    double Omega_rad = std::atan2(eq.p, eq.q);
    if (Omega_rad < 0.0) Omega_rad += TWO_PI;
    kep.Omega = Omega_rad;
    
    // 4. Longitudine del perielio: ϖ = atan2(h, k)
    double LP_rad = std::atan2(eq.h, eq.k);
    if (LP_rad < 0.0) LP_rad += TWO_PI;
    
    // 5. Argomento del perielio: ω = ϖ - Ω
    double omega_rad = LP_rad - Omega_rad;
    kep.omega = normalizeAngle(omega_rad);
    
    // 6. Anomalia media: M = λ - ϖ
#ifdef ASTDYN_EQ1_LAMBDA_RADIANS
    double lambda_rad = eq.lambda;
#else
    double lambda_rad = eq.lambda * DEG_TO_RAD;
#endif
    double M_rad = lambda_rad - LP_rad;
    kep.M = normalizeAngle(M_rad);
    
    return kep;
}

// ============================================================================
// Conversione Kepleriano → Cartesiano (Eclittico)
// ============================================================================

CartesianState OrbitalConversions::keplerianToCartesian(
    const KeplerianElements& kep) {
    
    // Validazione input
    if (!kep.isValid()) {
        throw std::runtime_error("Invalid Keplerian elements");
    }
    
    // 1. Risolvi equazione di Keplero per anomalia eccentrica E
    double E = solveKeplerEquation(kep.M, kep.e);
    
    // 2. Calcola anomalia vera ν
    double sin_E = std::sin(E);
    double cos_E = std::cos(E);
    double sqrt_1_e2 = std::sqrt(1.0 - kep.e * kep.e);
    
    // ν = atan2(sqrt(1-e²)·sin(E), cos(E) - e)
    double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - kep.e);
    
    // 3. Raggio vettore: r = a(1 - e·cos(E))
    double r = kep.a * (1.0 - kep.e * cos_E);
    
    // 4. Posizione nel piano orbitale
    double x_orb = r * std::cos(nu);
    double y_orb = r * std::sin(nu);
    // z_orb = 0 (piano orbitale)
    
    // 5. Velocità nel piano orbitale
    // v = sqrt(μ/a) / r
    double v_factor = std::sqrt(GM_SUN * kep.a) / r;
    double vx_orb = -v_factor * sin_E;
    double vy_orb = v_factor * sqrt_1_e2 * cos_E;
    // vz_orb = 0
    
    // 6. Matrice di Gauss (piano orbitale → eclittico J2000)
    // Rotazioni: Rz(Ω) · Rx(i) · Rz(ω)
    double cO = std::cos(kep.Omega);
    double sO = std::sin(kep.Omega);
    double cw = std::cos(kep.omega);
    double sw = std::sin(kep.omega);
    double ci = std::cos(kep.i);
    double si = std::sin(kep.i);
    
    // Elementi matrice di trasformazione
    double P11 = cO*cw - sO*sw*ci;
    double P12 = -cO*sw - sO*cw*ci;
    double P21 = sO*cw + cO*sw*ci;
    double P22 = -sO*sw + cO*cw*ci;
    double P31 = sw*si;
    double P32 = cw*si;
    
    // 7. Applica rotazione
    CartesianState state;
    state.position(0) = P11*x_orb + P12*y_orb;
    state.position(1) = P21*x_orb + P22*y_orb;
    state.position(2) = P31*x_orb + P32*y_orb;
    
    state.velocity(0) = P11*vx_orb + P12*vy_orb;
    state.velocity(1) = P21*vx_orb + P22*vy_orb;
    state.velocity(2) = P31*vx_orb + P32*vy_orb;
    
    state.epoch_jd = kep.epoch_jd;
    
    return state;
}

// ============================================================================
// Conversione Eclittico → ICRF
// ============================================================================

CartesianState OrbitalConversions::eclipticToICRF(
    const CartesianState& ecliptic) {
    
    // Rotazione attorno asse X con obliquità ε = 23.439291°
    // ICRF: x' = x
    //       y' = cos(ε)·y - sin(ε)·z
    //       z' = sin(ε)·y + cos(ε)·z
    
    double c_eps = std::cos(OBLIQUITY_J2000);
    double s_eps = std::sin(OBLIQUITY_J2000);
    
    CartesianState icrf;
    icrf.epoch_jd = ecliptic.epoch_jd;
    
    // Posizione
    icrf.position(0) = ecliptic.position(0);
    icrf.position(1) = c_eps * ecliptic.position(1) - 
                       s_eps * ecliptic.position(2);
    icrf.position(2) = s_eps * ecliptic.position(1) + 
                       c_eps * ecliptic.position(2);
    
    // Velocità
    icrf.velocity(0) = ecliptic.velocity(0);
    icrf.velocity(1) = c_eps * ecliptic.velocity(1) - 
                       s_eps * ecliptic.velocity(2);
    icrf.velocity(2) = s_eps * ecliptic.velocity(1) + 
                       c_eps * ecliptic.velocity(2);
    
    return icrf;
}

CartesianState OrbitalConversions::icrfToEcliptic(
    const CartesianState& icrf) {
    
    // Rotazione inversa (trasposta)
    double c_eps = std::cos(OBLIQUITY_J2000);
    double s_eps = std::sin(OBLIQUITY_J2000);
    
    CartesianState ecliptic;
    ecliptic.epoch_jd = icrf.epoch_jd;
    
    // Posizione
    ecliptic.position(0) = icrf.position(0);
    ecliptic.position(1) = c_eps * icrf.position(1) + 
                           s_eps * icrf.position(2);
    ecliptic.position(2) = -s_eps * icrf.position(1) + 
                            c_eps * icrf.position(2);
    
    // Velocità
    ecliptic.velocity(0) = icrf.velocity(0);
    ecliptic.velocity(1) = c_eps * icrf.velocity(1) + 
                           s_eps * icrf.velocity(2);
    ecliptic.velocity(2) = -s_eps * icrf.velocity(1) + 
                            c_eps * icrf.velocity(2);
    
    return ecliptic;
}

// ============================================================================
// Risoluzione Equazione di Keplero
// ============================================================================

double OrbitalConversions::solveKeplerEquation(
    double M, double e, double tol, int max_iter) {
    
    // Validazione
    if (e < 0.0 || e >= 1.0) {
        throw std::runtime_error("Invalid eccentricity for Kepler equation");
    }
    
    // Prima approssimazione
    double E = M;
    
    // Metodo di Newton-Raphson
    for (int iter = 0; iter < max_iter; iter++) {
        double f = E - e * std::sin(E) - M;
        double fp = 1.0 - e * std::cos(E);
        
        // Protezione divisione per zero
        if (std::abs(fp) < 1e-15) {
            throw std::runtime_error("Derivative near zero in Kepler solver");
        }
        
        double dE = f / fp;
        E -= dE;
        
        // Check convergenza
        if (std::abs(dE) < tol) {
            return E;
        }
    }
    
    // Non convergente
    std::ostringstream oss;
    oss << "Kepler equation did not converge after " << max_iter 
        << " iterations (M=" << M << ", e=" << e << ")";
    throw std::runtime_error(oss.str());
}

// ============================================================================
// Normalizzazione Angoli
// ============================================================================

double OrbitalConversions::normalizeAngle(double angle) {
    // Normalizza in [0, 2π)
    angle = std::fmod(angle, TWO_PI);
    if (angle < 0.0) {
        angle += TWO_PI;
    }
    return angle;
}

double OrbitalConversions::normalizeAngleSigned(double angle) {
    // Normalizza in [-π, π)
    angle = std::fmod(angle + PI, TWO_PI);
    if (angle < 0.0) {
        angle += TWO_PI;
    }
    return angle - PI;
}

// ============================================================================
// Validazione
// ============================================================================

bool OrbitalConversions::validateICRF(const CartesianState& state) {
    // 1. Verifica valori finiti
    if (!std::isfinite(state.position.norm()) ||
        !std::isfinite(state.velocity.norm())) {
        return false;
    }
    
    // 2. Verifica range posizione [0.1, 100 AU]
    double r = state.position.norm();
    if (r < 0.1 || r > 100.0) {
        return false;
    }
    
    // 3. Verifica velocità ragionevole [< 100 AU/day]
    double v = state.velocity.norm();
    if (v > 100.0) {
        return false;
    }
    
    // 4. Verifica epoca valida
    if (state.epoch_jd < 2400000.0 || state.epoch_jd > 2500000.0) {
        return false;
    }
    
    return true;
}

// ============================================================================
// Conversione Cartesiano → Kepleriano
// ============================================================================

KeplerianElements OrbitalConversions::cartesianToKeplerian(
    const CartesianState& state, double mu) {
    
    KeplerianElements kep;
    kep.epoch_jd = state.epoch_jd;
    
    const Eigen::Vector3d& r = state.position;
    const Eigen::Vector3d& v = state.velocity;
    
    double r_mag = r.norm();
    double v_mag = v.norm();
    
    // Momento angolare specifico
    Eigen::Vector3d h = r.cross(v);
    double h_mag = h.norm();
    
    // Vettore di Laplace-Runge-Lenz (eccentricità)
    Eigen::Vector3d e_vec = v.cross(h) / mu - r / r_mag;
    kep.e = e_vec.norm();
    
    // Energia specifica
    double xi = v_mag*v_mag / 2.0 - mu / r_mag;
    
    // Semiasse maggiore
    kep.a = -mu / (2.0 * xi);
    
    // Inclinazione
    kep.i = std::acos(h(2) / h_mag);
    
    // Vettore nodo ascendente
    Eigen::Vector3d n(-h(1), h(0), 0.0);
    double n_mag = n.norm();
    
    // Longitudine nodo ascendente
    if (n_mag > 1e-10) {
        kep.Omega = std::atan2(n(1), n(0));
        if (kep.Omega < 0.0) kep.Omega += TWO_PI;
    } else {
        kep.Omega = 0.0;  // Orbita equatoriale
    }
    
    // Argomento del perielio
    if (n_mag > 1e-10 && kep.e > 1e-10) {
        double cos_omega = n.dot(e_vec) / (n_mag * kep.e);
        cos_omega = std::max(-1.0, std::min(1.0, cos_omega));
        kep.omega = std::acos(cos_omega);
        if (e_vec(2) < 0.0) {
            kep.omega = TWO_PI - kep.omega;
        }
    } else {
        kep.omega = 0.0;
    }
    
    // Anomalia vera
    double nu;
    if (kep.e > 1e-10) {
        double cos_nu = e_vec.dot(r) / (kep.e * r_mag);
        cos_nu = std::max(-1.0, std::min(1.0, cos_nu));
        nu = std::acos(cos_nu);
        if (r.dot(v) < 0.0) {
            nu = TWO_PI - nu;
        }
    } else {
        // Orbita circolare
        nu = std::atan2(r(1), r(0));
        if (nu < 0.0) nu += TWO_PI;
    }
    
    // Anomalia eccentrica
    double E = 2.0 * std::atan(std::tan(nu / 2.0) / 
                               std::sqrt((1.0 + kep.e) / (1.0 - kep.e)));
    if (E < 0.0) E += TWO_PI;
    
    // Anomalia media
    kep.M = E - kep.e * std::sin(E);
    if (kep.M < 0.0) kep.M += TWO_PI;
    
    return kep;
}

} // namespace ioccultcalc
