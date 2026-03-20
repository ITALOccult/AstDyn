#pragma once
#include "constants.h"
#include "linalg.h"
#include <cmath>
#include <stdexcept>

namespace orbfit {

// ─────────────────────────────────────────────────────────────
//  Orbital elements (Keplerian, heliocentric J2000 ecliptic)
// ─────────────────────────────────────────────────────────────
struct OrbitalElements {
    double a;     // semi-major axis  [AU]
    double e;     // eccentricity     [dimensionless]
    double i;     // inclination      [rad]
    double Omega; // longitude of ascending node [rad]
    double omega; // argument of perihelion       [rad]
    double M0;    // mean anomaly at epoch        [rad]
    double epoch; // reference epoch  [JD]

    // Derived
    double n() const {  // mean motion [rad/day]
        return std::sqrt(GM_SUN / (a*a*a));
    }
    double T() const {  // orbital period [day]
        return 2.0*M_PI / n();
    }
    double q() const { return a*(1.0-e); }  // perihelion distance [AU]
    double Q() const { return a*(1.0+e); }  // aphelion distance   [AU]
};

// ─────────────────────────────────────────────────────────────
//  State vector: position + velocity (heliocentric, AU, AU/day)
// ─────────────────────────────────────────────────────────────
struct StateVector {
    Vec3   pos;   // [AU]
    Vec3   vel;   // [AU/day]
    double epoch; // [JD]
};

// ─────────────────────────────────────────────────────────────
//  Solve Kepler's equation  M = E - e sin E  for E
//  Uses Halley's method — converges in ~3 iterations
// ─────────────────────────────────────────────────────────────
inline double solve_kepler(double M, double e, double tol=1e-12) {
    // Normalise M to [0, 2pi)
    M = std::fmod(M, 2.0*M_PI);
    if(M < 0) M += 2.0*M_PI;

    // Initial guess
    double E = (e < 0.8) ? M : M_PI;

    for(int it=0; it<50; it++){
        double sinE = std::sin(E);
        double cosE = std::cos(E);
        double f    = E - e*sinE - M;
        double fp   = 1.0 - e*cosE;
        double fpp  = e*sinE;
        double dE   = -f*fp / (fp*fp - 0.5*f*fpp);   // Halley
        E += dE;
        if(std::abs(dE) < tol) break;
    }
    return E;
}

// ─────────────────────────────────────────────────────────────
//  Keplerian elements → state vector
// ─────────────────────────────────────────────────────────────
inline StateVector elements_to_state(const OrbitalElements& el, double t) {
    double dt = t - el.epoch;
    double M  = el.M0 + el.n() * dt;
    double E  = solve_kepler(M, el.e);

    double sinE = std::sin(E);
    double cosE = std::cos(E);
    double sqrt1e2 = std::sqrt(1.0 - el.e*el.e);

    // Position in orbital plane
    double r    = el.a * (1.0 - el.e*cosE);
    double x_op = el.a * (cosE - el.e);
    double y_op = el.a * sqrt1e2 * sinE;

    // Velocity in orbital plane
    double fac   = el.a * el.n() / (1.0 - el.e*cosE);
    double vx_op = -fac * sinE;
    double vy_op =  fac * sqrt1e2 * cosE;

    // Rotation matrices: Rz(-Omega) Rx(-i) Rz(-omega)
    double cO=std::cos(el.Omega), sO=std::sin(el.Omega);
    double ci=std::cos(el.i),     si=std::sin(el.i);
    double co=std::cos(el.omega), so=std::sin(el.omega);

    // Direction cosines (Peri-focal to ecliptic)
    double Px = cO*co - sO*so*ci;
    double Py = sO*co + cO*so*ci;
    double Pz = so*si;
    double Qx = -cO*so - sO*co*ci;
    double Qy = -sO*so + cO*co*ci;
    double Qz =  co*si;

    StateVector sv;
    sv.epoch = t;
    sv.pos = Vec3(Px*x_op + Qx*y_op,
                  Py*x_op + Qy*y_op,
                  Pz*x_op + Qz*y_op);
    sv.vel = Vec3(Px*vx_op + Qx*vy_op,
                  Py*vx_op + Qy*vy_op,
                  Pz*vx_op + Qz*vy_op);
    return sv;
}

// ─────────────────────────────────────────────────────────────
//  State vector → Keplerian elements  (two-body, heliocentric)
// ─────────────────────────────────────────────────────────────
inline OrbitalElements state_to_elements(const StateVector& sv) {
    const Vec3& r = sv.pos;
    const Vec3& v = sv.vel;
    double rv = r.norm();
    double vv = v.norm();

    // Angular momentum
    Vec3   h    = r.cross(v);
    double habs = h.norm();

    // Eccentricity vector
    double mu = GM_SUN;
    Vec3   ecc = (r*(vv*vv - mu/rv) - v*(r.dot(v))) / mu;
    double e   = ecc.norm();

    // Node vector
    Vec3 k(0,0,1);
    Vec3 N = k.cross(h);
    double Nabs = N.norm();

    // Semi-major axis
    double xi = 0.5*vv*vv - mu/rv;   // specific energy
    double a  = -0.5*mu/xi;

    // Inclination
    double i = std::acos(std::max(-1.0, std::min(1.0, h.z/habs)));

    // Longitude of ascending node
    double Omega = 0.0;
    if(Nabs > 1e-15) {
        Omega = std::acos(std::max(-1.0, std::min(1.0, N.x/Nabs)));
        if(N.y < 0) Omega = 2.0*M_PI - Omega;
    }

    // Argument of perihelion
    double omega = 0.0;
    if(Nabs > 1e-15 && e > 1e-10) {
        omega = std::acos(std::max(-1.0, std::min(1.0, N.dot(ecc)/(Nabs*e))));
        if(ecc.z < 0) omega = 2.0*M_PI - omega;
    }

    // True anomaly
    double f = 0.0;
    if(e > 1e-10) {
        f = std::acos(std::max(-1.0, std::min(1.0, ecc.dot(r)/(e*rv))));
        if(r.dot(v) < 0) f = 2.0*M_PI - f;
    }

    // Mean anomaly from true anomaly
    double cosf = std::cos(f), sinf = std::sin(f);
    double tanE2 = std::sin(f)*std::sqrt(1.0-e*e) / (1.0+e*cosf);
    double E  = 2.0*std::atan(std::sqrt((1.0-e)/(1.0+e))*std::tan(f/2.0));
    if(E < 0) E += 2.0*M_PI;
    double M0 = E - e*std::sin(E);

    OrbitalElements el;
    el.a     = a;
    el.e     = e;
    el.i     = i;
    el.Omega = Omega;
    el.omega = omega;
    el.M0    = M0;
    el.epoch = sv.epoch;
    return el;
}

// ─────────────────────────────────────────────────────────────
//  Partial derivatives  ∂(r,v)/∂(a,e,i,Ω,ω,M0)  — 6×6 Jacobian
//  Computed by finite differences (step: 1e-6 relative)
// ─────────────────────────────────────────────────────────────
inline Matrix elements_jacobian(const OrbitalElements& el, double t) {
    const double EPS = 1e-6;
    Matrix J(6,6);

    // Pack elements into array
    auto pack = [](const OrbitalElements& e) -> std::array<double,6> {
        return {e.a, e.e, e.i, e.Omega, e.omega, e.M0};
    };
    auto unpack = [&](std::array<double,6> p, OrbitalElements base) -> OrbitalElements {
        base.a=p[0]; base.e=p[1]; base.i=p[2];
        base.Omega=p[3]; base.omega=p[4]; base.M0=p[5];
        return base;
    };

    auto sv0 = elements_to_state(el, t);
    std::array<double,6> p0 = pack(el);

    for(int j=0;j<6;j++){
        auto pp = p0; pp[j] += EPS * std::max(std::abs(p0[j]), 1e-10);
        auto pm = p0; pm[j] -= EPS * std::max(std::abs(p0[j]), 1e-10);
        double h2 = pp[j] - pm[j];

        auto svp = elements_to_state(unpack(pp,el), t);
        auto svm = elements_to_state(unpack(pm,el), t);

        J(0,j) = (svp.pos.x - svm.pos.x)/h2;
        J(1,j) = (svp.pos.y - svm.pos.y)/h2;
        J(2,j) = (svp.pos.z - svm.pos.z)/h2;
        J(3,j) = (svp.vel.x - svm.vel.x)/h2;
        J(4,j) = (svp.vel.y - svm.vel.y)/h2;
        J(5,j) = (svp.vel.z - svm.vel.z)/h2;
    }
    return J;
}

} // namespace orbfit
