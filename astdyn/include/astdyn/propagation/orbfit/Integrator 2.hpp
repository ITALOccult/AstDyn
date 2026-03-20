#pragma once
#include "constants.h"
#include "linalg.h"
#include "orbital_elements.h"
#include <vector>
#include <cmath>
#include <functional>

namespace orbfit {

// ─────────────────────────────────────────────────────────────
//  Planetary ephemeris — simple analytical approximation
//  (VSOP87 truncated; suitable for NEO/asteroid work < 100yr)
//  Returns heliocentric ecliptic J2000 position [AU] at JD t
// ─────────────────────────────────────────────────────────────

struct PlanetState {
    Vec3 pos; // AU, heliocentric ecliptic J2000
    Vec3 vel; // AU/day
    double gm;
};

// Keplerian mean elements at J2000 + secular rates (deg, deg/century)
// Source: Standish 1992, simplified
struct PlanetKep {
    double a0, adot;   // AU, AU/cy
    double e0, edot;   // -, -/cy
    double i0, idot;   // deg, deg/cy
    double L0, Ldot;   // mean longitude deg, deg/cy
    double lp0, lpdot; // long. perihelion deg, deg/cy
    double Om0, Omdot; // long. node deg, deg/cy
    double gm;
};

static const PlanetKep PLANET_KEP[9] = {
    // Mercury
    {0.38709927, 3.7e-7, 0.20563593, 1.906e-5,
     7.00497902, -0.00594749, 252.25032350, 149472.67411175,
     77.45779628,  0.16047689, 48.33076593, -0.12534081, GM_MERCURY},
    // Venus
    {0.72333566, 3.9e-6, 0.00677672, -4.107e-5,
     3.39467605, -0.0007889,  181.97909950, 58517.81538729,
     131.60246718, 0.00268329, 76.67984255, -0.27769418, GM_VENUS},
    // Earth-Moon barycentre
    {1.00000261, 1.23e-6, 0.01671123, -4.392e-5,
    -0.00001531, -0.01294668, 100.46457166, 35999.37244981,
    102.93768193,  0.32327364, 0.0, 0.0, GM_EARTH+GM_MOON},
    // Mars
    {1.52371034, 1.847e-4, 0.09339410, 7.882e-5,
     1.84969142, -0.00813131, -4.55343205, 19140.30268499,
    -23.94362959,  0.44441088, 49.55953891, -0.29257343, GM_MARS},
    // Jupiter
    {5.20288700, -1.1607e-3, 0.04838624, -1.3253e-4,
     1.30439695, -0.00183714,  34.39644051, 3034.74612775,
     14.72847983,  0.21252668, 100.47390909,  0.20469106, GM_JUPITER},
    // Saturn
    {9.53667594, -1.25060e-2, 0.05386179, -5.0991e-4,
     2.48599187,  0.00193609,  49.95424423, 1222.49362201,
     92.59887831, -0.41897216, 113.66242448, -0.28867794, GM_SATURN},
    // Uranus
    {19.18916464, -1.96176e-2, 0.04725744, -4.397e-5,
     0.77263783, -0.00242939, 313.23810451,  428.48202785,
    170.95427630,  0.40805281, 74.01692503,  0.04240589, GM_URANUS},
    // Neptune
    {30.06992276,  2.62159e-3, 0.00859048,  5.105e-5,
     1.77004347,  3.343e-5,  -55.12002969,  218.45945325,
     44.96476227, -0.32241464, 131.78422574, -0.00508664, GM_NEPTUNE},
    // Pluto (for completeness, GM tiny)
    {39.48211675, -3.12250e-3, 0.24882730, 5.1701e-4,
     17.14001206,  1.4004e-4, 238.92903833,  145.20780515,
    224.06891629, -0.04062942, 110.30393684, -0.01183482, 7.35e-9}
};

inline PlanetState planet_state(int pid, double t_jd) {
    // centuries since J2000
    double T = (t_jd - JD_J2000) / 36525.0;
    const auto& p = PLANET_KEP[pid];

    double a     = p.a0  + p.adot *T;
    double e     = p.e0  + p.edot *T;
    double i     = (p.i0  + p.idot *T)*DEG2RAD;
    double L     = (p.L0  + p.Ldot *T)*DEG2RAD;
    double lp    = (p.lp0 + p.lpdot*T)*DEG2RAD;
    double Omega = (p.Om0 + p.Omdot*T)*DEG2RAD;
    double omega = lp - Omega;          // arg. perihelion
    double M     = L  - lp;             // mean anomaly
    // Normalise
    M = std::fmod(M, 2.0*M_PI);
    if(M<0) M+=2.0*M_PI;

    OrbitalElements el;
    el.a=a; el.e=e; el.i=i; el.Omega=Omega; el.omega=omega;
    el.M0=M; el.epoch=t_jd;

    auto sv = elements_to_state(el, t_jd);
    return {sv.pos, sv.vel, p.gm};
}

// ─────────────────────────────────────────────────────────────
//  Gravitational acceleration  ẍ = -GM_sun r/r³ + Σ planet perturbations
//  + relativistic correction (1PN, Soffel 2003)
// ─────────────────────────────────────────────────────────────

// State packed as [x,y,z, vx,vy,vz]
using State6 = std::array<double,6>;

inline State6 eom(const State6& s, double t,
                  bool perturbers=true,
                  bool relativistic=false)
{
    Vec3 r{s[0],s[1],s[2]};
    Vec3 v{s[3],s[4],s[5]};
    double rv = r.norm();
    double rv3 = rv*rv*rv;

    // Central body
    Vec3 acc = r * (-GM_SUN / rv3);

    // N-body perturbations (8 planets)
    if(perturbers){
        for(int pid=0; pid<8; pid++){
            auto ps = planet_state(pid, t);
            Vec3 dp  = r - ps.pos;        // target - planet
            double dp3 = std::pow(dp.norm(), 3.0);
            double rp  = ps.pos.norm();
            double rp3 = rp*rp*rp;
            // Direct + indirect terms
            acc -= dp * (ps.gm / dp3);
            acc -= ps.pos * (ps.gm / rp3);
        }
    }

    // Post-Newtonian (1PN) correction
    if(relativistic){
        double vv2 = v.norm2();
        double rdotv = r.dot(v);
        double c2  = CLIGHT*CLIGHT;
        Vec3 pn = r*(-GM_SUN/rv3) *
                  (4.0*GM_SUN/rv - vv2 + 4.0*rdotv*rdotv/rv/rv*(3.0/2.0))
                + v*(4.0*rdotv/rv);
        acc += pn * (1.0/c2);
    }

    return {v.x, v.y, v.z, acc.x, acc.y, acc.z};
}

// ─────────────────────────────────────────────────────────────
//  RK4 integrator (fixed step)
// ─────────────────────────────────────────────────────────────
inline State6 rk4_step(const State6& y, double t, double h,
                        bool perturbers=true, bool relativ=false)
{
    auto f = [&](const State6& s, double tt){ return eom(s,tt,perturbers,relativ); };

    auto add = [](const State6& a, const State6& b, double s) -> State6 {
        State6 r;
        for(int i=0;i<6;i++) r[i]=a[i]+b[i]*s;
        return r;
    };

    auto k1 = f(y,          t        );
    auto k2 = f(add(y,k1,h/2), t+h/2);
    auto k3 = f(add(y,k2,h/2), t+h/2);
    auto k4 = f(add(y,k3,h),   t+h  );

    State6 yn;
    for(int i=0;i<6;i++)
        yn[i] = y[i] + h/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    return yn;
}

// ─────────────────────────────────────────────────────────────
//  Dormand-Prince RK45 (adaptive) — used as main integrator
// ─────────────────────────────────────────────────────────────

struct IntegratorOptions {
    double rtol       = 1e-10;
    double atol       = 1e-12;
    double h_init     = 0.1;    // initial step [day]
    double h_max      = 5.0;    // max step [day]
    double h_min      = 1e-6;   // min step [day]
    bool   perturbers = true;
    bool   relativistic = false;
};

inline StateVector propagate(const StateVector& sv0, double t_end,
                              const IntegratorOptions& opt = {})
{
    // Dormand-Prince tableau
    static const double a21=1.0/5,
        a31=3.0/40, a32=9.0/40,
        a41=44.0/45, a42=-56.0/15, a43=32.0/9,
        a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729,
        a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247,
            a64=49.0/176, a65=-5103.0/18656;

    static const double b1=35.0/384, b3=500.0/1113, b4=125.0/192,
                        b5=-2187.0/6784, b6=11.0/84;
    static const double e1=71.0/57600, e3=-71.0/16695, e4=71.0/1920,
                        e5=-17253.0/339200, e6=22.0/525, e7=-1.0/40;

    State6 y{sv0.pos.x, sv0.pos.y, sv0.pos.z,
             sv0.vel.x, sv0.vel.y, sv0.vel.z};
    double t = sv0.epoch;
    double h = (t_end > t) ? opt.h_init : -opt.h_init;
    double sign = (t_end > t) ? 1.0 : -1.0;

    auto F = [&](const State6& s, double tt){
        return eom(s, tt, opt.perturbers, opt.relativistic);
    };
    auto ax = [](const State6& a, const State6& b, double s)->State6{
        State6 r; for(int i=0;i<6;i++) r[i]=a[i]+b[i]*s; return r;
    };
    auto axby = [](const State6& a, double sa,
                   const State6& b, double sb)->State6{
        State6 r; for(int i=0;i<6;i++) r[i]=a[i]*sa+b[i]*sb; return r;
    };
    auto addN = [](std::initializer_list<std::pair<const State6*,double>> terms)->State6{
        State6 r{}; r.fill(0.0);
        for(auto& [v,c]: terms) for(int i=0;i<6;i++) r[i]+=(*v)[i]*c;
        return r;
    };

    int maxstep = 1000000;
    for(int step=0; step<maxstep; step++){
        // Check if we've reached t_end
        if(sign*(t_end - t) <= 0.0) break;
        // Clip step to not overshoot
        if(std::abs(h) > std::abs(t_end - t))
            h = t_end - t;

        auto k1 = F(y, t);
        auto k2 = F(ax(y,k1, h*a21), t + h*0.2);
        auto k3 = F(addN({{&y,1},{&k1,h*a31},{&k2,h*a32}}), t+h*0.3);
        auto k4 = F(addN({{&y,1},{&k1,h*a41},{&k2,h*a42},{&k3,h*a43}}), t+h*0.8);
        auto k5 = F(addN({{&y,1},{&k1,h*a51},{&k2,h*a52},{&k3,h*a53},{&k4,h*a54}}), t+h);
        auto k6 = F(addN({{&y,1},{&k1,h*a61},{&k2,h*a62},{&k3,h*a63},{&k4,h*a64},{&k5,h*a65}}), t+h);

        // 5th order solution
        State6 y5;
        for(int i=0;i<6;i++)
            y5[i] = y[i] + h*(b1*k1[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);

        auto k7 = F(y5, t+h);

        // Error estimate (difference of 5th and 4th order)
        double err = 0.0;
        for(int i=0;i<6;i++){
            double ei = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
            double sc = opt.atol + opt.rtol*std::max(std::abs(y[i]),std::abs(y5[i]));
            err += (ei/sc)*(ei/sc);
        }
        err = std::sqrt(err/6.0);

        if(err <= 1.0){
            // Accept step
            t += h;
            y  = y5;
        }

        // New step size (PI controller)
        double factor = 0.9 * std::pow(1.0/(err+1e-10), 0.2);
        factor = std::max(0.2, std::min(factor, 10.0));
        h *= factor;
        h  = std::copysign(std::min(std::abs(h), opt.h_max), sign);
        h  = std::copysign(std::max(std::abs(h), opt.h_min), sign);
    }

    StateVector sv_out;
    sv_out.epoch = t_end;
    sv_out.pos = Vec3(y[0],y[1],y[2]);
    sv_out.vel = Vec3(y[3],y[4],y[5]);
    return sv_out;
}

// ─────────────────────────────────────────────────────────────
//  Variational equations: propagate 6×6 state transition matrix
//  STM: Φ satisfies  dΦ/dt = A(t) Φ,  Φ(t0)=I
//  where A = ∂f/∂y (Jacobian of EOM)
// ─────────────────────────────────────────────────────────────
struct StateSTM {
    State6            y;    // 6 state variables
    std::array<double,36> Phi; // 6x6 STM, row-major
};

inline StateSTM propagate_stm(const StateVector& sv0, double t_end,
                               const IntegratorOptions& opt = {})
{
    // EOM for augmented state (y, Phi)  -> 42 equations
    // A matrix computed by finite differences
    auto eom_stm = [&](const StateSTM& s, double t) -> StateSTM {
        StateSTM ds;
        // state derivative
        auto dy = eom(s.y, t, opt.perturbers, opt.relativistic);
        ds.y = dy;

        // A = ∂f/∂y  (6×6, finite diff)
        const double eps = 1e-7;
        std::array<std::array<double,6>,6> A{};
        for(int j=0;j<6;j++){
            State6 yp=s.y, ym=s.y;
            yp[j]+=eps; ym[j]-=eps;
            auto fp=eom(yp,t,opt.perturbers,opt.relativistic);
            auto fm=eom(ym,t,opt.perturbers,opt.relativistic);
            for(int i=0;i<6;i++) A[i][j]=(fp[i]-fm[i])/(2*eps);
        }

        // dPhi/dt = A * Phi
        for(int i=0;i<6;i++)
            for(int j=0;j<6;j++){
                double v=0;
                for(int k=0;k<6;k++) v+=A[i][k]*s.Phi[k*6+j];
                ds.Phi[i*6+j]=v;
            }
        return ds;
    };

    // RK4 for augmented system
    StateSTM s0;
    s0.y = {sv0.pos.x, sv0.pos.y, sv0.pos.z,
            sv0.vel.x, sv0.vel.y, sv0.vel.z};
    // Identity STM
    s0.Phi.fill(0.0);
    for(int i=0;i<6;i++) s0.Phi[i*6+i]=1.0;

    double t = sv0.epoch;
    double h = (t_end>t) ? opt.h_init : -opt.h_init;
    double sign = (t_end>t)?1.0:-1.0;

    auto add_stm = [](const StateSTM& a, const StateSTM& b, double sc)->StateSTM{
        StateSTM r;
        for(int i=0;i<6;i++)  r.y[i]=a.y[i]+b.y[i]*sc;
        for(int i=0;i<36;i++) r.Phi[i]=a.Phi[i]+b.Phi[i]*sc;
        return r;
    };

    auto s = s0;
    int maxstep = 100000;
    for(int step=0; step<maxstep && sign*(t_end-t)>0; step++){
        if(std::abs(h) > std::abs(t_end-t)) h=t_end-t;
        auto k1=eom_stm(s,t);
        auto k2=eom_stm(add_stm(s,k1,h/2),t+h/2);
        auto k3=eom_stm(add_stm(s,k2,h/2),t+h/2);
        auto k4=eom_stm(add_stm(s,k3,h),  t+h  );
        for(int i=0;i<6;i++)
            s.y[i]+= h/6*(k1.y[i]+2*k2.y[i]+2*k3.y[i]+k4.y[i]);
        for(int i=0;i<36;i++)
            s.Phi[i]+=h/6*(k1.Phi[i]+2*k2.Phi[i]+2*k3.Phi[i]+k4.Phi[i]);
        t+=h;
    }
    return s;
}

} // namespace orbfit
