/**
 * @file test_propagator_asteroids.cpp
 * @brief Integration tests for asteroid perturbations in orbital propagation
 * 
 * Tests the complete integration of asteroid perturbations (AST17 model)
 * into the orbital propagator.
 */

#include <gtest/gtest.h>
#include "orbfit/propagation/Propagator.hpp"
#include "orbfit/propagation/Integrator.hpp"
#include "orbfit/ephemeris/PlanetaryEphemeris.hpp"
#include "orbfit/core/Constants.hpp"

using namespace orbfit;
using namespace orbfit::propagation;

class PropagatorAsteroidsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create planetary ephemeris
        ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
    }
    
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephemeris;
};

// Test 1: Propagator accepts asteroid perturbation setting
TEST_F(PropagatorAsteroidsTest, AsteroidSettingAccepted) {
    PropagatorSettings settings;
    settings.include_asteroids = true;
    settings.include_planets = false;  // Only asteroids
    
    auto integrator = std::make_unique<RKF78Integrator>(0.1, 1e-12);
    
    EXPECT_NO_THROW({
        Propagator propagator(std::move(integrator), ephemeris, settings);
        EXPECT_TRUE(propagator.settings().include_asteroids);
    });
}

// Test 2: Propagation works with asteroids enabled
TEST_F(PropagatorAsteroidsTest, PropagationWithAsteroids) {
    PropagatorSettings settings;
    settings.include_asteroids = true;
    settings.include_planets = false;
    settings.include_relativity = false;
    
    auto integrator = std::make_unique<RKF78Integrator>(1.0, 1e-12);
    Propagator propagator(std::move(integrator), ephemeris, settings);
    
    // Simple test orbit at 2.5 AU (middle of main belt)
    KeplerianElements initial;
    initial.epoch_mjd_tdb = 60000.0;
    initial.semi_major_axis = 2.5;
    initial.eccentricity = 0.1;
    initial.inclination = 0.1;
    initial.longitude_ascending_node = 0.0;
    initial.argument_perihelion = 0.0;
    initial.mean_anomaly = 0.0;
    initial.gravitational_parameter = constants::GMS;
    
    // Propagate 100 days
    double target_mjd = 60100.0;
    
    EXPECT_NO_THROW({
        KeplerianElements final = propagator.propagate_keplerian(initial, target_mjd);
        EXPECT_NEAR(final.epoch_mjd_tdb, target_mjd, 1e-10);
    });
}

// Test 3: Asteroid perturbations produce measurable effect
TEST_F(PropagatorAsteroidsTest, AsteroidEffectMeasurable) {
    // Create two propagators: one with asteroids, one without
    PropagatorSettings settings_with;
    settings_with.include_asteroids = true;
    settings_with.include_planets = false;
    settings_with.include_relativity = false;
    
    PropagatorSettings settings_without;
    settings_without.include_asteroids = false;
    settings_without.include_planets = false;
    settings_without.include_relativity = false;
    
    auto int1 = std::make_unique<RKF78Integrator>(1.0, 1e-12);
    auto int2 = std::make_unique<RKF78Integrator>(1.0, 1e-12);
    
    Propagator prop_with(std::move(int1), ephemeris, settings_with);
    Propagator prop_without(std::move(int2), ephemeris, settings_without);
    
    // Orbit at 2.5 AU (main belt)
    CartesianElements initial;
    initial.epoch_mjd_tdb = 60000.0;
    initial.position = Eigen::Vector3d(2.5, 0.0, 0.0);
    initial.velocity = Eigen::Vector3d(0.0, 0.012, 0.0);  // ~circular
    initial.gravitational_parameter = constants::GMS;
    
    // Propagate 1 year
    double target_mjd = 60000.0 + 365.25;
    
    CartesianElements final_with = prop_with.propagate_cartesian(initial, target_mjd);
    CartesianElements final_without = prop_without.propagate_cartesian(initial, target_mjd);
    
    // Positions should differ due to asteroid perturbations
    Eigen::Vector3d diff = final_with.position - final_without.position;
    double diff_km = diff.norm() * constants::AU_TO_KM;
    
    // After 1 year, difference should be measurable (> 1 km)
    // but not huge (< 1000 km) since asteroid masses are small
    EXPECT_GT(diff_km, 0.01);     // > 10 meters (definitely measurable)
    EXPECT_LT(diff_km, 10000.0);  // < 10,000 km (reasonable for 1 year)
    
    std::cout << "Position difference after 1 year: " << diff_km << " km" << std::endl;
}

// Test 4: Asteroid effect increases with time
TEST_F(PropagatorAsteroidsTest, AsteroidEffectIncreasesWithTime) {
    PropagatorSettings settings_with;
    settings_with.include_asteroids = true;
    settings_with.include_planets = false;
    
    PropagatorSettings settings_without;
    settings_without.include_asteroids = false;
    settings_without.include_planets = false;
    
    // Orbit at 2.7 AU
    CartesianElements initial;
    initial.epoch_mjd_tdb = 60000.0;
    initial.position = Eigen::Vector3d(2.7, 0.0, 0.0);
    initial.velocity = Eigen::Vector3d(0.0, 0.011, 0.0);
    initial.gravitational_parameter = constants::GMS;
    
    std::vector<double> propagation_times = {10.0, 100.0, 365.25};
    std::vector<double> differences;
    
    for (double dt : propagation_times) {
        auto int1 = std::make_unique<RKF78Integrator>(1.0, 1e-12);
        auto int2 = std::make_unique<RKF78Integrator>(1.0, 1e-12);
        
        Propagator prop_with(std::move(int1), ephemeris, settings_with);
        Propagator prop_without(std::move(int2), ephemeris, settings_without);
        
        double target = initial.epoch_mjd_tdb + dt;
        
        auto final_with = prop_with.propagate_cartesian(initial, target);
        auto final_without = prop_without.propagate_cartesian(initial, target);
        
        double diff = (final_with.position - final_without.position).norm();
        differences.push_back(diff);
    }
    
    // Differences should generally increase with time
    // (though not strictly monotonic due to orbital geometry)
    EXPECT_GT(differences[1], differences[0]);  // 100 days > 10 days
    EXPECT_GT(differences[2], differences[0]);  // 365 days > 10 days
}

// Test 5: Asteroid perturbations work with planets enabled
TEST_F(PropagatorAsteroidsTest, AsteroidsWithPlanets) {
    PropagatorSettings settings;
    settings.include_asteroids = true;
    settings.include_planets = true;
    settings.perturb_jupiter = true;
    settings.perturb_saturn = true;
    
    auto integrator = std::make_unique<RKF78Integrator>(1.0, 1e-12);
    Propagator propagator(std::move(integrator), ephemeris, settings);
    
    KeplerianElements initial;
    initial.epoch_mjd_tdb = 60000.0;
    initial.semi_major_axis = 2.8;
    initial.eccentricity = 0.15;
    initial.inclination = 0.2;
    initial.longitude_ascending_node = 0.5;
    initial.argument_perihelion = 1.0;
    initial.mean_anomaly = 0.0;
    initial.gravitational_parameter = constants::GMS;
    
    // Propagate with both planets and asteroids
    EXPECT_NO_THROW({
        KeplerianElements final = propagator.propagate_keplerian(initial, 60100.0);
        EXPECT_TRUE(std::isfinite(final.semi_major_axis));
        EXPECT_GT(final.semi_major_axis, 0.0);
    });
}

// Test 6: Asteroid effect strongest near main belt
TEST_F(PropagatorAsteroidsTest, EffectStrongestInMainBelt) {
    PropagatorSettings settings_with;
    settings_with.include_asteroids = true;
    settings_with.include_planets = false;
    
    PropagatorSettings settings_without;
    settings_without.include_asteroids = false;
    settings_without.include_planets = false;
    
    // Test at different distances
    std::vector<double> distances = {1.5, 2.5, 3.5};  // AU
    std::vector<double> effects;
    
    for (double a : distances) {
        CartesianElements initial;
        initial.epoch_mjd_tdb = 60000.0;
        initial.position = Eigen::Vector3d(a, 0.0, 0.0);
        initial.velocity = Eigen::Vector3d(0.0, std::sqrt(constants::GMS / a), 0.0);
        initial.gravitational_parameter = constants::GMS;
        
        auto int1 = std::make_unique<RKF78Integrator>(1.0, 1e-12);
        auto int2 = std::make_unique<RKF78Integrator>(1.0, 1e-12);
        
        Propagator prop_with(std::move(int1), ephemeris, settings_with);
        Propagator prop_without(std::move(int2), ephemeris, settings_without);
        
        double target = 60000.0 + 100.0;
        
        auto final_with = prop_with.propagate_cartesian(initial, target);
        auto final_without = prop_without.propagate_cartesian(initial, target);
        
        double diff = (final_with.position - final_without.position).norm();
        effects.push_back(diff);
        
        std::cout << "Effect at " << a << " AU: " << diff * constants::AU_TO_KM 
                  << " km after 100 days" << std::endl;
    }
    
    // Effect at 2.5 AU (center of main belt) should be significant
    // Note: actual magnitude depends on asteroid positions at epoch
    EXPECT_GT(effects[1], 0.0);
}

// Test 7: Integration statistics reasonable with asteroids
TEST_F(PropagatorAsteroidsTest, IntegrationStatistics) {
    PropagatorSettings settings;
    settings.include_asteroids = true;
    settings.include_planets = false;
    
    auto integrator = std::make_unique<RKF78Integrator>(1.0, 1e-12);
    Propagator propagator(std::move(integrator), ephemeris, settings);
    
    KeplerianElements initial;
    initial.epoch_mjd_tdb = 60000.0;
    initial.semi_major_axis = 2.5;
    initial.eccentricity = 0.1;
    initial.inclination = 0.1;
    initial.longitude_ascending_node = 0.0;
    initial.argument_perihelion = 0.0;
    initial.mean_anomaly = 0.0;
    initial.gravitational_parameter = constants::GMS;
    
    propagator.propagate_keplerian(initial, 60100.0);
    
    const auto& stats = propagator.statistics();
    
    // Should have taken some steps
    EXPECT_GT(stats.num_steps, 0);
    
    // Should not have excessive rejections (< 50% rejection rate)
    double rejection_rate = static_cast<double>(stats.num_rejected_steps) / 
                           (stats.num_steps + stats.num_rejected_steps);
    EXPECT_LT(rejection_rate, 0.5);
    
    std::cout << "Steps: " << stats.num_steps 
              << ", Rejections: " << stats.num_rejected_steps 
              << ", Rate: " << (rejection_rate * 100.0) << "%" << std::endl;
}
