/**
 * OrbFit C++ — Example: fit a synthetic asteroid
 *
 * Generates synthetic astrometric observations for a known orbit,
 * runs the full IOD + DC pipeline, and checks recovery accuracy.
 *
 * Build:
 *   g++ -O2 -std=c++17 -I../include -o demo demo.cpp
 */

#include "orbfit.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

using namespace orbfit;

// ─────────────────────────────────────────────────────────────
//  Generate synthetic observations for a known orbit
// ─────────────────────────────────────────────────────────────
std::vector<Observation> generate_obs(const OrbitalElements& true_el,
                                       double t_start,
                                       int n_obs,
                                       double dt_days,
                                       double sigma_arcsec = 0.5,
                                       unsigned seed = 42)
{
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, 1.0);
    double sigma_rad = sigma_arcsec * ARCSEC2RAD;

    // Observer at geocentre (simplified: use Earth heliocentric position)
    std::vector<Observation> obs;

    for(int k=0; k<n_obs; k++){
        double t = t_start + k*dt_days;
        // Get Earth heliocentric pos
        auto earth = planet_state(2, t);

        // True target state (no perturbations for ground truth)
        IntegratorOptions iopt_true;
        iopt_true.perturbers=false;
        auto sv = propagate(elements_to_state(true_el, true_el.epoch), t, iopt_true);

        // Topocentric vector (ecliptic)
        Vec3 rho_ecl = sv.pos - earth.pos;
        double rho   = rho_ecl.norm();

        // Convert to equatorial
        const double eps = 23.4392911 * DEG2RAD;
        double ce=std::cos(eps), se=std::sin(eps);
        double xeq = rho_ecl.x;
        double yeq = ce*rho_ecl.y - se*rho_ecl.z;
        double zeq = se*rho_ecl.y + ce*rho_ecl.z;

        double ra  = std::atan2(yeq, xeq);
        if(ra<0) ra+=2.0*M_PI;
        double dec = std::asin(zeq/rho);

        // Add Gaussian noise
        double cos_dec = std::cos(dec);
        double ra_noisy  = ra  + noise(rng)*sigma_rad/cos_dec;
        double dec_noisy = dec + noise(rng)*sigma_rad;

        Observation o;
        o.t_obs    = t;
        o.ra       = ra_noisy;
        o.dec      = dec_noisy;
        o.sigma_ra = sigma_rad;
        o.sigma_dec= sigma_rad;
        o.obs_pos  = earth.pos;
        o.station_code = "500";  // geocentre
        obs.push_back(o);
    }
    return obs;
}

// ─────────────────────────────────────────────────────────────
//  Print comparison: true vs fitted elements
// ─────────────────────────────────────────────────────────────
void compare_elements(const OrbitalElements& truth,
                       const FitResult& fit)
{
    auto s = fit.sigma();
    std::cout << "\n=== Element Recovery ===\n";
    std::cout << std::setw(8) << "Param"
              << std::setw(16) << "True"
              << std::setw(16) << "Fitted"
              << std::setw(16) << "Residual"
              << std::setw(14) << "Pull (σ)\n";
    std::cout << std::string(70,'-') << "\n";

    auto row = [&](const char* name, double tr, double fi, double sig, double fac=1.0){
        double res=std::abs(fi-tr)*fac;
        double pull=(sig>0)?res/sig:0;
        std::cout << std::setw(8) << name
                  << std::setw(16) << std::setprecision(8) << tr*fac
                  << std::setw(16) << fi*fac
                  << std::setw(16) << res
                  << std::setw(14) << std::setprecision(4) << pull << "\n";
    };
    std::cout << std::fixed;
    row("a [AU]",   truth.a,     fit.elements.a,     s[0]);
    row("e",        truth.e,     fit.elements.e,     s[1]);
    row("i [deg]",  truth.i,     fit.elements.i,     s[2], RAD2DEG);
    row("Ω [deg]",  truth.Omega, fit.elements.Omega, s[3], RAD2DEG);
    row("ω [deg]",  truth.omega, fit.elements.omega, s[4], RAD2DEG);
    row("M0 [deg]", truth.M0,    fit.elements.M0,    s[5], RAD2DEG);
}

// ─────────────────────────────────────────────────────────────
//  Main
// ─────────────────────────────────────────────────────────────
int main(){
    std::cout << "╔══════════════════════════════════════════════╗\n";
    std::cout << "║  OrbFit C++ — Orbital Fitting Pipeline       ║\n";
    std::cout << "╚══════════════════════════════════════════════╝\n\n";

    // ── True orbital elements (MBA-like asteroid) ──────────────
    OrbitalElements true_el;
    true_el.a     = 2.5634;              // AU
    true_el.e     = 0.1423;
    true_el.i     = 8.745 * DEG2RAD;    // rad
    true_el.Omega = 112.34 * DEG2RAD;
    true_el.omega = 273.11 * DEG2RAD;
    true_el.M0    = 47.22  * DEG2RAD;
    true_el.epoch = JD_J2000;           // J2000.0

    std::cout << "True orbit:\n";
    std::cout << "  a=" << true_el.a << " AU,  e=" << true_el.e
              << ",  i=" << true_el.i*RAD2DEG << "°\n";
    std::cout << "  Ω=" << true_el.Omega*RAD2DEG << "°"
              << ",  ω=" << true_el.omega*RAD2DEG << "°"
              << ",  M0=" << true_el.M0*RAD2DEG << "°\n\n";

    // ── Generate observations (30-day arc, 1 obs/day) ──────────
    int    n_obs   = 20;
    double dt      = 1.5;        // [day]
    double t_start = JD_J2000 + 100.0;
    double sigma   = 1.0;        // [arcsec]

    std::cout << "Generating " << n_obs << " synthetic observations "
              << "(σ=" << sigma << "\", arc=" << (n_obs-1)*dt << " days)...\n";
    auto obs = generate_obs(true_el, t_start, n_obs, dt, sigma);

    // ── Step 1: Initial Orbit Determination ───────────────────
    std::cout << "\n[1] Initial Orbit Determination (Gauss method)...\n";
    auto [i1,i2,i3] = select_iod_obs(obs);
    GaussSolution iod = gauss_iod(obs[i1], obs[i2], obs[i3]);
    std::cout << "  IOD RMS = " << std::setprecision(4) << iod.rms_init << " arcsec\n";
    std::cout << "  a_iod = " << iod.elements.a << " AU,  e_iod = " << iod.elements.e << "\n";

    // ── Step 2: Differential Corrections ──────────────────────
    std::cout << "\n[2] Differential Corrections (Gauss-Newton/LM)...\n";
    DCOptions dco;
    dco.max_iter  = 30;
    dco.verbose   = true;
    dco.iopt.perturbers  = false;  // two-body for speed in demo
    dco.iopt.relativistic= false;

    FitResult fit = differential_corrections(iod.elements, obs, dco);

    // ── Results ───────────────────────────────────────────────
    std::cout << "\n";
    fit.print();
    compare_elements(true_el, fit);

    // ── Step 3: LOV sampling (3σ ellipsoid) ───────────────────
    std::cout << "\n[3] Sampling Line of Variations (3σ, 11 VOs)...\n";
    ConfidenceEllipsoid ell;
    ell.covariance  = fit.covariance;
    ell.sigma_level = 3.0;
    auto lov_pts = ell.lov_sample(fit.elements, 11);
    std::cout << "  Generated " << lov_pts.size() << " Virtual Objects\n";
    std::cout << "  a range: ["
              << std::setprecision(6)
              << lov_pts.front().a << ", " << lov_pts.back().a << "] AU\n";

    // ── Step 4: Close approach check (simple demo) ────────────
    std::cout << "\n[4] MOID estimate...\n";
    double moid = moid_estimate(fit.elements);
    std::cout << "  MOID ≈ " << std::setprecision(5) << moid << " AU"
              << "  (" << std::setprecision(1) << moid/EARTH_RADIUS << " Earth radii)\n";

    // ── Test: NEO-like orbit ───────────────────────────────────
    std::cout << "\n══════════════════════════════════\n";
    std::cout << "TEST 2: NEO-like orbit (a=1.2 AU, e=0.45)\n";
    std::cout << "══════════════════════════════════\n";

    OrbitalElements neo_el;
    neo_el.a     = 1.2;
    neo_el.e     = 0.45;
    neo_el.i     = 23.0 * DEG2RAD;
    neo_el.Omega = 60.0 * DEG2RAD;
    neo_el.omega = 200.0 * DEG2RAD;
    neo_el.M0    = 10.0 * DEG2RAD;
    neo_el.epoch = JD_J2000;

    auto obs_neo = generate_obs(neo_el, t_start, 15, 2.0, 0.5, 123);

    auto [j1,j2,j3] = select_iod_obs(obs_neo);
    auto iod_neo = gauss_iod(obs_neo[j1], obs_neo[j2], obs_neo[j3]);

    DCOptions dco_neo; dco_neo.max_iter=30; dco_neo.verbose=false;
    dco_neo.iopt.perturbers=false;
    auto fit_neo = differential_corrections(iod_neo.elements, obs_neo, dco_neo);

    std::cout << "NEO fit: a=" << std::setprecision(5) << fit_neo.elements.a
              << " AU (true=" << neo_el.a << ")"
              << "  e=" << fit_neo.elements.e
              << " (true=" << neo_el.e << ")\n";
    std::cout << "RMS=" << std::setprecision(4) << fit_neo.rms
              << " arcsec  converged=" << fit_neo.converged << "\n";

    double moid_neo = moid_estimate(fit_neo.elements);
    std::cout << "MOID ≈ " << std::setprecision(5) << moid_neo
              << " AU  ("  << moid_neo/EARTH_RADIUS << " R⊕)\n";
    if(moid_neo < 0.05)
        std::cout << "⚠  MOID < 0.05 AU — potential Potentially Hazardous Asteroid!\n";

    std::cout << "\nDone.\n";
    return 0;
}
