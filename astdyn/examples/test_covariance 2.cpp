/**
 * @file test_covariance.cpp
 * @brief Minimal verification of CovariancePropagator.
 *
 * Property tested:
 *   Without process noise (Q = 0), trace(P_f) >= trace(P_0).
 *
 * Physical basis: for typical heliocentric orbits, positional uncertainty
 * grows along the orbit track (along-track spreading), so the overall
 * trace of the covariance increases over time.
 *
 * Usage: this example requires no external ephemeris file because planets
 * are disabled in PropagatorSettings (pure two-body test).
 */

#include "astdyn/propagation/CovariancePropagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/PropagatorSettings.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/core/Types.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::constants;

// Circular orbit at 3 AU (typical MBA), in ECLIPJ2000
static physics::CartesianStateTyped<core::ECLIPJ2000> make_circular_mba_state() {
    const double t0_mjd = constants::MJD2000;
    const double r_au = 3.0;
    const double r_m  = r_au * constants::AU * 1000.0;

    // Circular velocity: v = sqrt(GM_sun_SI / r)
    const double v_ms = std::sqrt(constants::GM_SUN * 1e9 / r_m);

    return physics::CartesianStateTyped<core::ECLIPJ2000>::from_si(
        time::EpochTDB::from_mjd(t0_mjd),
        r_m, 0.0, 0.0,          // position [m]
        0.0, v_ms, 0.0,          // velocity [m/s]
        constants::GM_SUN * 1e9  // GM [m³/s²]
    );
}

// Isotropic diagonal covariance: 150 km in position, tiny velocity uncertainty
static physics::CovarianceMatrixAU<core::ECLIPJ2000> make_initial_covariance() {
    const double sigma_pos_au  = 150.0 / constants::AU_TO_KM;   // 150 km in AU
    const double sigma_vel_aud = 1e-8;                           // [AU/day]

    Matrix6d p0 = Matrix6d::Zero();
    p0.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * (sigma_pos_au  * sigma_pos_au);
    p0.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * (sigma_vel_aud * sigma_vel_aud);
    return physics::CovarianceMatrixAU<core::ECLIPJ2000>(p0);
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
        settings
    );
}

static void print_covariance_summary(
    const physics::CovarianceMatrixAU<core::ECLIPJ2000>& p0,
    const CovariancePropagator::Result& result)
{
    const double trace_p0 = p0.matrix().trace();
    const double trace_pf = result.covariance.matrix().trace();

    std::cout << "trace(P_0) = " << trace_p0 << " AU²(+AU²/day²)\n";
    std::cout << "trace(P_f) = " << trace_pf << " AU²(+AU²/day²)\n";
    std::cout << "sigma_pos  = " << result.covariance.sigma_pos_km()
              << " km   (initial: " << p0.sigma_pos_km() << " km)\n";
    std::cout << "sigma_vel  = " << result.covariance.sigma_vel_ms()
              << " m/s  (initial: " << p0.sigma_vel_ms() << " m/s)\n";

    // Without Q, uncertainty should not decrease (along-track spreading dominates)
    assert(trace_pf >= trace_p0 && "FAIL: covariance trace decreased without process noise");
    std::cout << "PASS: trace(P_f) >= trace(P_0) with Q = 0\n";
}

int main() {
    const auto initial_state      = make_circular_mba_state();
    const auto initial_covariance = make_initial_covariance();
    const auto process_noise      = physics::CovarianceMatrixAU<core::ECLIPJ2000>{}; // Q = 0

    const auto propagator   = make_two_body_propagator();
    auto cov_propagator     = CovariancePropagator::make(propagator);
    const double dt_days    = 30.0;
    const auto target_epoch = time::EpochTDB::from_mjd(
        initial_state.epoch.mjd() + dt_days);

    const auto result = cov_propagator.propagate_with_covariance(
        initial_state, initial_covariance, target_epoch, process_noise);

    print_covariance_summary(initial_covariance, result);
    return 0;
}
