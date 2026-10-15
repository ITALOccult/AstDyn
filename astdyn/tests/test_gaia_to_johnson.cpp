/**
 * @file test_gaia_to_johnson.cpp
 * @brief Johnson V and Cousins R from Gaia G and BP-RP (Riello et al. 2021).
 */

#include <gtest/gtest.h>
#include "astdyn/catalog/GaiaToJohnson.hpp"

using astdyn::catalog::cousins_r_from_gaia;
using astdyn::catalog::johnson_v_from_gaia;

TEST(GaiaToJohnson, ZeroColourGivesTheConstantTerms) {
    EXPECT_NEAR(*johnson_v_from_gaia(10.0, 0.0), 10.02704, 1e-9);
    EXPECT_NEAR(*cousins_r_from_gaia(10.0, 0.0), 10.02275, 1e-9);
}

TEST(GaiaToJohnson, SolarColourStarIsFainterInVThanInG) {
    // BP-RP = 0.82 (G2V): G - V = -0.152, V - R about 0.35 as for the Sun.
    const double v = *johnson_v_from_gaia(12.0, 0.82);
    const double r = *cousins_r_from_gaia(12.0, 0.82);
    EXPECT_NEAR(v - 12.0, 0.1525, 1e-3);
    EXPECT_NEAR(v - r, 0.35, 0.03);
}

TEST(GaiaToJohnson, OutsideTheCalibratedColourRangeNothingIsReturned) {
    EXPECT_FALSE(johnson_v_from_gaia(12.0, 5.2).has_value());
    EXPECT_FALSE(cousins_r_from_gaia(12.0, 3.0).has_value());
    EXPECT_FALSE(johnson_v_from_gaia(12.0, -0.7).has_value());
}
