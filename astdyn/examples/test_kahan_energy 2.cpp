/**
 * @file test_kahan_energy.cpp
 * @brief Two-part verification of Kahan compensated summation.
 *
 * Part 1 — Summation accuracy over 365 steps (5x improvement):
 *   Models 1 year of daily integration steps: one large "background" term
 *   + 365 small "perturbation" increments + cancelling neg term.
 *   Each small increment is below ULP(large), so naive summation loses all
 *   365 of them; Kahan accumulates the compensation until it exceeds ULP/2
 *   and recovers the sum. Expected improvement: >> 5x.
 *
 * Part 2 — Energy conservation (1 year, 2.5 AU, two-body, no ephemeris):
 *   |ΔE/E| < 1e-8 after 365-day RKF78 integration.
 */

#include "astdyn/math/kahan_sum.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/PropagatorSettings.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/core/Constants.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::constants;

// ============================================================================
// Part 1: Summation accuracy — 365-term sequence
// ============================================================================
//
// Sequence: [large, δ, δ, ..., δ (N=365 times), -large]
// Each δ = {1/365, 1/365, 1/365}  →  exact_sum = {1, 1, 1}
//
// Naive: δ << ULP(large) ≈ 0.125, so each δ rounds to 0 → sum = 0. Error = 1.
// Kahan: compensation grows by δ each step; after ~23 steps exceeds ULP/2,
//        triggers a carry into acc, and the full sum is recovered.

static constexpr int    N_STEPS   = 365;
static constexpr double LARGE_VAL = 1e15;

static Eigen::Vector3d exact_sum() {
    return Eigen::Vector3d::Ones();  // N * (1/N) = 1 per component
}

static double naive_summation_error() {
    const Eigen::Vector3d large = Eigen::Vector3d::Constant(LARGE_VAL);
    const Eigen::Vector3d delta = Eigen::Vector3d::Constant(1.0 / N_STEPS);
    const Eigen::Vector3d neg   = Eigen::Vector3d::Constant(-LARGE_VAL);

    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    acc += large;
    for (int k = 0; k < N_STEPS; ++k) acc += delta;  // each lost: ULP too large
    acc += neg;
    return (acc - exact_sum()).norm() / exact_sum().norm();
}

static double kahan_summation_error() {
    const Eigen::Vector3d large = Eigen::Vector3d::Constant(LARGE_VAL);
    const Eigen::Vector3d delta = Eigen::Vector3d::Constant(1.0 / N_STEPS);
    const Eigen::Vector3d neg   = Eigen::Vector3d::Constant(-LARGE_VAL);

    Eigen::Vector3d acc  = Eigen::Vector3d::Zero();
    Eigen::Vector3d comp = Eigen::Vector3d::Zero();
    math::kahan_add(acc, comp, large);
    for (int k = 0; k < N_STEPS; ++k) math::kahan_add(acc, comp, delta);
    math::kahan_add(acc, comp, neg);
    return (acc - exact_sum()).norm() / exact_sum().norm();
}

static void verify_kahan_vs_naive() {
    const double err_naive   = naive_summation_error();
    const double err_kahan   = kahan_summation_error();
    const double improvement = (err_kahan > 0.0) ? err_naive / err_kahan : 1e20;

    std::cout << "naive_error  = " << err_naive  << " (relative, 365-step sequence)\n";
    std::cout << "kahan_error  = " << err_kahan  << " (relative, 365-step sequence)\n";
    std::cout << "improvement  = " << improvement << "x\n";

    assert(err_naive > 0.0 &&
           "naive summation unexpectedly exact — check LARGE_VAL and N_STEPS");
    assert(err_kahan * 5.0 <= err_naive &&
           "FAIL: Kahan not 5x better than naive over 365 steps");
    std::cout << "PASS: Kahan " << improvement << "x better than naive\n";
}

// ============================================================================
// Part 2: Energy conservation over 1 year (two-body, no ephemeris)
// ============================================================================

// Two-body: circular orbit at 2.5 AU — position [AU], velocity [AU/day]
static Eigen::VectorXd make_circular_orbit_25au() {
    const double r_au  = 2.5;
    const double v_aud = std::sqrt(GMS / r_au);
    Eigen::VectorXd y0(6);
    y0 << r_au, 0.0, 0.0,
          0.0,  v_aud, 0.0;
    return y0;
}

// Specific mechanical energy E = ½v² - GMS/r  [AU²/day²]
static double specific_energy_au_day2(const Eigen::VectorXd& y_au) {
    const double r_au  = y_au.head<3>().norm();
    const double v2_au = y_au.tail<3>().squaredNorm();
    return 0.5 * v2_au - GMS / r_au;
}

static std::shared_ptr<Propagator> make_two_body_propagator() {
    PropagatorSettings settings;
    settings.include_planets    = false;
    settings.include_moon       = false;
    settings.include_asteroids  = false;
    settings.include_relativity = false;
    settings.include_earth_j2   = false;
    settings.include_sun_j2     = false;
    return std::make_shared<Propagator>(
        std::make_shared<RKF78Integrator>(0.1),
        std::make_shared<ephemeris::PlanetaryEphemeris>(),
        settings);
}

static void verify_energy_1year(Propagator& propagator) {
    const Eigen::VectorXd y0_au  = make_circular_orbit_25au();
    const double          e0     = specific_energy_au_day2(y0_au);
    const double          t0_mjd = MJD2000;
    const double          tf_mjd = MJD2000 + DAYS_PER_YEAR;

    const Eigen::VectorXd yf_au = propagator.integrate_raw_au(y0_au, t0_mjd, tf_mjd);
    const double          ef    = specific_energy_au_day2(yf_au);
    const double          de_e  = std::abs((ef - e0) / e0);

    std::cout << "E_0     = " << e0  << " AU²/day²\n";
    std::cout << "E_final = " << ef  << " AU²/day²\n";
    std::cout << "|ΔE/E|  = " << de_e << "  (1 year, 2.5 AU, two-body)\n";

    assert(de_e < 1e-8 &&
           "FAIL: |ΔE/E| >= 1e-8 after 1-year two-body integration");
    std::cout << "PASS: |ΔE/E| = " << de_e << " < 1e-8\n";
}

// ============================================================================
// main
// ============================================================================

int main() {
    verify_kahan_vs_naive();

    auto propagator = make_two_body_propagator();
    verify_energy_1year(*propagator);

    return 0;
}
