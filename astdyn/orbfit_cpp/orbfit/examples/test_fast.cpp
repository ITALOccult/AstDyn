/**
 * OrbFit C++ — Unit Tests (fast, two-body only)
 * Build: g++ -O2 -std=c++17 -I../include -o test_fast test_fast.cpp
 */

#include "orbfit.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <random>

using namespace orbfit;

bool nearly_equal(double a, double b, double tol=1e-6){
    return std::abs(a-b) < tol*(1+std::abs(a));
}

void test_kepler(){
    std::cout << "[TEST] Kepler equation solver...\n";
    for(double e : {0.0, 0.1, 0.5, 0.9, 0.99}){
        for(double M : {0.0, 0.5, 1.0, 2.0, 3.14159}){
            double E = solve_kepler(M, e);
            double M2= E - e*std::sin(E);
            // Normalise
            M = std::fmod(M,2*M_PI); if(M<0) M+=2*M_PI;
            assert(std::abs(M2-M)<1e-11 && "Kepler solver failed");
        }
    }
    std::cout << "  PASS\n";
}

void test_roundtrip(){
    std::cout << "[TEST] Elements ↔ State roundtrip...\n";
    OrbitalElements el;
    el.a=2.5; el.e=0.15; el.i=10*DEG2RAD;
    el.Omega=120*DEG2RAD; el.omega=250*DEG2RAD; el.M0=30*DEG2RAD;
    el.epoch=JD_J2000;

    auto sv  = elements_to_state(el, JD_J2000);
    auto el2 = state_to_elements(sv);

    assert(nearly_equal(el.a, el2.a, 1e-8));
    assert(nearly_equal(el.e, el2.e, 1e-8));
    assert(nearly_equal(el.i, el2.i, 1e-7));
    assert(nearly_equal(el.Omega, el2.Omega, 1e-6));
    assert(nearly_equal(el.omega, el2.omega, 1e-6));
    std::cout << "  PASS  (a err=" << std::abs(el.a-el2.a) << ")\n";
}

void test_propagation(){
    std::cout << "[TEST] Two-body energy conservation over 1 year...\n";
    OrbitalElements el;
    el.a=1.5; el.e=0.2; el.i=5*DEG2RAD;
    el.Omega=30*DEG2RAD; el.omega=100*DEG2RAD; el.M0=0.0;
    el.epoch=JD_J2000;

    auto sv0 = elements_to_state(el, JD_J2000);

    IntegratorOptions opt;
    opt.perturbers=false;
    opt.rtol=1e-12; opt.atol=1e-14;
    auto sv1 = propagate(sv0, JD_J2000 + 365.25, opt);

    double E0 = 0.5*sv0.vel.norm2() - GM_SUN/sv0.pos.norm();
    double E1 = 0.5*sv1.vel.norm2() - GM_SUN/sv1.pos.norm();
    double dE = std::abs(E1-E0)/std::abs(E0);
    std::cout << "  ΔE/E = " << std::scientific << dE;
    assert(dE < 1e-9 && "Energy not conserved");
    std::cout << "  PASS\n";
}

void test_dc(){
    std::cout << "[TEST] Differential corrections recovery...\n";

    // True orbit
    OrbitalElements true_el;
    true_el.a=2.0; true_el.e=0.1; true_el.i=15*DEG2RAD;
    true_el.Omega=80*DEG2RAD; true_el.omega=130*DEG2RAD; true_el.M0=60*DEG2RAD;
    true_el.epoch=JD_J2000;

    // Generate clean observations (no noise, 2-body)
    IntegratorOptions iopt; iopt.perturbers=false;
    std::vector<Observation> obs;
    double t0 = JD_J2000 + 200;

    for(int k=0; k<12; k++){
        double t = t0 + k*2.0;
        auto earth = planet_state(2, t);
        auto sv = propagate(elements_to_state(true_el, true_el.epoch), t, iopt);

        Vec3 rho_ecl = sv.pos - earth.pos;
        double rho = rho_ecl.norm();
        const double eps=23.4392911*DEG2RAD;
        double ce=std::cos(eps),se=std::sin(eps);
        double xeq=rho_ecl.x;
        double yeq=ce*rho_ecl.y-se*rho_ecl.z;
        double zeq=se*rho_ecl.y+ce*rho_ecl.z;
        double ra=std::atan2(yeq,xeq); if(ra<0) ra+=2*M_PI;
        double dec=std::asin(zeq/rho);

        Observation o;
        o.t_obs=t; o.ra=ra; o.dec=dec;
        o.sigma_ra=o.sigma_dec=1*ARCSEC2RAD;
        o.obs_pos=earth.pos;
        obs.push_back(o);
    }

    // IOD on obs 0,5,11
    auto iod = gauss_iod(obs[0], obs[5], obs[11]);
    std::cout << "  IOD: a=" << std::setprecision(4) << iod.elements.a
              << "  true=" << true_el.a << "\n";

    // DC
    DCOptions dco; dco.max_iter=25; dco.verbose=false; dco.iopt=iopt;
    auto fit = differential_corrections(iod.elements, obs, dco);

    double da = std::abs(fit.elements.a - true_el.a);
    double de = std::abs(fit.elements.e - true_el.e);
    std::cout << "  Fit: a=" << fit.elements.a << " (err=" << da << ")"
              << "  e=" << fit.elements.e << " (err=" << de << ")\n";
    std::cout << "  RMS=" << fit.rms << " arcsec  iter=" << fit.iterations
              << "  converged=" << fit.converged << "\n";

    assert(da < 1e-4 && "Semi-major axis not recovered");
    assert(de < 1e-4 && "Eccentricity not recovered");
    std::cout << "  PASS\n";
}

void test_moid(){
    std::cout << "[TEST] MOID estimation...\n";
    // Earth-crossing orbit should have MOID near 0
    OrbitalElements neo;
    neo.a=1.0; neo.e=0.05; neo.i=0*DEG2RAD;
    neo.Omega=0; neo.omega=0; neo.M0=0; neo.epoch=JD_J2000;
    double m = moid_estimate(neo, JD_J2000, 90);
    std::cout << "  Earth-like MOID = " << std::setprecision(5) << m << " AU\n";
    assert(m < 0.2 && "MOID too large for Earth-like orbit");
    std::cout << "  PASS\n";
}

int main(){
    std::cout << "╔═══════════════════════════════╗\n";
    std::cout << "║  OrbFit C++ — Unit Tests      ║\n";
    std::cout << "╚═══════════════════════════════╝\n\n";
    test_kepler();
    test_roundtrip();
    test_propagation();
    test_dc();
    test_moid();
    std::cout << "\n✓  All tests passed.\n";
    return 0;
}
