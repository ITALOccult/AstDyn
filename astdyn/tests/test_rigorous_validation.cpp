/**
 * @file test_rigorous_validation.cpp
 * @brief Rigorous validation tests for AstDyn strong types and transformations.
 */

#include <gtest/gtest.h>
#include <astdyn/core/Types.hpp>
#include <astdyn/core/Constants.hpp>
#include "src/types/vectors.hpp"
#include "src/core/frame_tags.hpp"
#include "src/core/units.hpp"
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/observations/ObservatoryDatabase.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <type_traits>

using namespace astdyn;
using namespace astdyn::core;
using namespace astdyn::types;
using namespace astdyn::propagation;

// SFINAE helper to check if two Vector3 can be added
template <typename V1, typename V2, typename = void>
struct is_addable : std::false_type {};

template <typename V1, typename V2>
struct is_addable<V1, V2, std::void_t<decltype(std::declval<V1>() + std::declval<V2>())>> : std::true_type {};

/**
 * @brief Test 1: Unit Type Safety (Static)
 * Verifies that vectors of different frames cannot be added.
 */
TEST(RigorousValidation, UnitTypeSafety) {
    using GCRF_Meter = Vector3<GCRF, Meter>;
    using ITRF_Meter = Vector3<ITRF, Meter>;

    // Same frame, same unit -> Addable
    static_assert(is_addable<GCRF_Meter, GCRF_Meter>::value, "Identical vectors should be addable");
    
    // Different frame -> NOT addable (compilation error if attempted)
    static_assert(!is_addable<GCRF_Meter, ITRF_Meter>::value, "Different frames should NOT be addable");
    
    SUCCEED();
}

/**
 * @brief Test 2: Transformation Test
 * Verify geodetic to GCRF conversion at a specific Instant.
 */
TEST(RigorousValidation, TransformationTest) {
    auto& db = observations::ObservatoryDatabase::getInstance();
    db.loadDefaultObservatories();
    auto mk_opt = db.getObservatory("568"); // Mauna Kea
    ASSERT_TRUE(mk_opt.has_value());
    const auto& mk = *mk_opt;

    // Epoch: J2000.0 (2000-01-01 12:00:00 UTC)
    time::EpochUTC t = time::EpochUTC::from_mjd(51544.5);
    
    // Convert to GCRF
    auto pos_gcrf = mk.getPositionGCRF(t);
    
    // Verify distance remains Earth-radius scale (~6378 km)
    double dist = std::sqrt(pos_gcrf.x * pos_gcrf.x + pos_gcrf.y * pos_gcrf.y + pos_gcrf.z * pos_gcrf.z);
    EXPECT_NEAR(dist, 6378137.0 + 4200.0, 30000.0); // Approx radius + alt, accounting for oblateness
    
    // Check non-zero components
    EXPECT_NE(pos_gcrf.x, 0.0);
    EXPECT_NE(pos_gcrf.y, 0.0);
}

/**
 * @brief Test 3: Propagation Test
 * Propagate a circular orbit for one period.
 */
TEST(RigorousValidation, PropagationTest) {
    // Earth-like circular orbit in meters
    double a = 7000000.0; // 7000 km
    double mu = constants::GM_EARTH * 1e9; // km^3/s^2 to m^3/s^2
    
    auto initial = physics::KeplerianStateTyped<core::GCRF>::from_traditional(
        time::EpochTDB::from_mjd(51544.5),
        a / (constants::AU * 1000.0), 0.0, 0.0, 0.0, 0.0, 0.0,
        physics::GravitationalParameter::from_si(mu)
    );
    
    double a_au = initial.a.to_au();
    double mu_au_d2 = initial.gm.to_au3_d2();
    double n_rad_day = std::sqrt(mu_au_d2 / (a_au * a_au * a_au));
    double period_days = constants::TWO_PI / n_rad_day;
    
    time::EpochTDB t_start = time::EpochTDB::from_mjd(51544.5);
    time::EpochTDB t_end = time::EpochTDB::from_mjd(51544.5 + period_days);
    
    auto final = TwoBodyPropagator::propagate(initial, t_end);
    
    // Period should be exactly recovered
    EXPECT_NEAR(final.a.to_au(), initial.a.to_au(), 1e-7);
    EXPECT_LT(final.e, constants::EPSILON * 1e6); // Tight tolerance
    
    // Cartesian recovery
    auto cart_init = propagation::keplerian_to_cartesian(initial);
    auto cart_final = propagation::keplerian_to_cartesian(final);
    
    double dist_err = std::sqrt(
        std::pow(cart_init.position.x_si() - cart_final.position.x_si(), 2) +
        std::pow(cart_init.position.y_si() - cart_final.position.y_si(), 2) +
        std::pow(cart_init.position.z_si() - cart_final.position.z_si(), 2)
    );
    
    EXPECT_LT(dist_err, 1e-2); // < 1 centimeter
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
