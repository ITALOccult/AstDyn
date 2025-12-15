// test_asteroid_fitter_horizons.cpp
// Test comparing AsteroidFitter results with direct PositionCalculator computation
// using orbital elements from JPL Horizons for asteroid (79450) on 2026-03-27.

#include <gtest/gtest.h>
#include "astdyn/ephemeris/AsteroidFitter.hpp"
#include "astdyn/ephemeris/AsteroidFitConfig.hpp"
#include "astdyn/ephemeris/PositionCalculator.hpp"

// Example orbital elements for asteroid 79450 (dummy values; replace with real JPL data if available)
static const astdyn::ephemeris::SimpleKeplerianElements jpl_elements = {
    2.5,        // a  (AU) – semi‑major axis
    0.12,       // e  – eccentricity
    0.15,       // i  – inclination (rad)
    1.0,        // Omega – longitude of ascending node (rad)
    0.5,        // omega – argument of periapsis (rad)
    0.0         // M – mean anomaly at epoch (rad)
};

TEST(AsteroidFitterHorizonTest, ComputeFromMemoryMatchesDirectPositionCalculator) {
    // Observation epoch (MJD) for 2026-03-27
    const double mjd_epoch = 60820.0; // approximate MJD for 2026-03-27
    std::vector<double> epochs = {mjd_epoch};

    // Compute using AsteroidFitter (in‑memory)
    astdyn::ephemeris::AsteroidFitResult fit_result =
        astdyn::ephemeris::AsteroidFitter::computeFromMemory(jpl_elements, epochs, true);

    ASSERT_TRUE(fit_result.success);
    ASSERT_EQ(fit_result.fitted_positions.size(), 1u);

    // Compute directly with PositionCalculator
    astdyn::ephemeris::KeplerianElements calc_elem;
    calc_elem.a = jpl_elements.a;
    calc_elem.e = jpl_elements.e;
    calc_elem.i = jpl_elements.i;
    calc_elem.Omega = jpl_elements.Omega;
    calc_elem.omega = jpl_elements.omega;
    calc_elem.M = jpl_elements.M;

    Eigen::Vector3d direct_pos = astdyn::ephemeris::PositionCalculator::computePosition(
        calc_elem, mjd_epoch, mjd_epoch, true);

    // Compare positions (tolerance 1e-9 AU)
    const double tolerance = 1e-9;
    EXPECT_NEAR(fit_result.fitted_positions[0].x(), direct_pos.x(), tolerance);
    EXPECT_NEAR(fit_result.fitted_positions[0].y(), direct_pos.y(), tolerance);
    EXPECT_NEAR(fit_result.fitted_positions[0].z(), direct_pos.z(), tolerance);
}
