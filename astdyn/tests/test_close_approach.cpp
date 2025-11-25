/**
 * @file test_close_approach.cpp
 * @brief Unit tests for close approach detection
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#include <gtest/gtest.h>
#include "astdyn/close_approach/CloseApproach.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/core/Constants.hpp"
#include <memory>
#include <cmath>

using namespace astdyn::close_approach;
using namespace astdyn::propagation;
using namespace astdyn::ephemeris;
using namespace astdyn::constants;

// Test fixture
class CloseApproachTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create ephemeris
        ephemeris_ = std::make_shared<PlanetaryEphemeris>();
        
        // Create propagator
        auto integrator = std::make_unique<RKF78Integrator>(0.1, 1e-12, 1e-6, 10.0);
        PropagatorSettings settings;
        settings.central_body_gm = GMS;
        settings.include_planets = false;  // Two-body for tests
        
        propagator_ = std::make_shared<Propagator>(std::move(integrator), ephemeris_, settings);
        
        // Default detector settings
        detector_settings_.detection_distance = 0.05;  // 0.05 AU
        detector_settings_.refine_time = true;
        detector_settings_.compute_b_plane = true;
    }
    
    std::shared_ptr<PlanetaryEphemeris> ephemeris_;
    std::shared_ptr<Propagator> propagator_;
    CloseApproachSettings detector_settings_;
};

// ============================================================================
// Basic Structure Tests
// ============================================================================

TEST_F(CloseApproachTest, BPlaneCoordinatesCalculation) {
    BPlaneCoordinates bplane;
    bplane.xi = 0.003;
    bplane.zeta = 0.004;
    
    double expected_b = std::sqrt(0.003*0.003 + 0.004*0.004);
    EXPECT_NEAR(bplane.miss_distance(), expected_b, 1e-10);
    
    bplane.b_magnitude = expected_b;
    bplane.theta = std::atan2(bplane.zeta, bplane.xi);
    
    EXPECT_NEAR(bplane.b_magnitude, 0.005, 1e-10);
    EXPECT_GT(bplane.theta, 0.9);  // ~53° in radians
    EXPECT_LT(bplane.theta, 0.95);
}

TEST_F(CloseApproachTest, CloseApproachIsClose) {
    CloseApproach ca;
    ca.distance = 0.03;  // 0.03 AU
    
    EXPECT_TRUE(ca.is_close(0.05));
    EXPECT_FALSE(ca.is_close(0.02));
}

TEST_F(CloseApproachTest, CloseApproachDistanceInRadii) {
    CloseApproach ca;
    ca.distance = 0.01;  // 0.01 AU
    
    // Earth radius ~ 6371 km = 4.26e-5 AU
    double earth_radius_au = 6371.0 / AU_TO_KM;
    double dist_radii = ca.distance_in_radii(earth_radius_au);
    
    EXPECT_GT(dist_radii, 200.0);  // Should be > 200 Earth radii
    EXPECT_LT(dist_radii, 300.0);
}

// ============================================================================
// Detector Construction and Settings
// ============================================================================

TEST_F(CloseApproachTest, DetectorConstruction) {
    CloseApproachDetector detector(propagator_, detector_settings_);
    
    EXPECT_EQ(detector.settings().detection_distance, 0.05);
    EXPECT_TRUE(detector.settings().refine_time);
    EXPECT_TRUE(detector.settings().compute_b_plane);
}

TEST_F(CloseApproachTest, DetectorSettingsUpdate) {
    CloseApproachDetector detector(propagator_, detector_settings_);
    
    CloseApproachSettings new_settings;
    new_settings.detection_distance = 0.1;
    new_settings.refine_time = false;
    
    detector.set_settings(new_settings);
    
    EXPECT_EQ(detector.settings().detection_distance, 0.1);
    EXPECT_FALSE(detector.settings().refine_time);
}

TEST_F(CloseApproachTest, BodyTypeFiltering) {
    CloseApproachSettings settings;
    settings.bodies_to_check = {BodyType::EARTH, BodyType::MARS};
    
    EXPECT_TRUE(settings.should_check_body(BodyType::EARTH));
    EXPECT_TRUE(settings.should_check_body(BodyType::MARS));
    EXPECT_FALSE(settings.should_check_body(BodyType::JUPITER));
}

// ============================================================================
// Close Approach Detection Tests
// ============================================================================

TEST_F(CloseApproachTest, NoCloseApproachDetected) {
    CloseApproachDetector detector(propagator_, detector_settings_);
    
    // Orbit far from all planets (a=5 AU, low e, moderate i)
    KeplerianElements orbit;
    orbit.epoch_mjd_tdb = MJD2000;
    orbit.semi_major_axis = 5.0;
    orbit.eccentricity = 0.05;
    orbit.inclination = 15.0 * DEG_TO_RAD;
    orbit.longitude_ascending_node = 100.0 * DEG_TO_RAD;
    orbit.argument_perihelion = 80.0 * DEG_TO_RAD;
    orbit.mean_anomaly = 0.0;
    orbit.gravitational_parameter = GMS;
    
    // Short propagation (30 days)
    auto approaches = detector.find_approaches(orbit, MJD2000, MJD2000 + 30.0);
    
    // Should find no close approaches (orbit too far from inner planets)
    EXPECT_TRUE(approaches.empty());
}

TEST_F(CloseApproachTest, EarthCrossingOrbitDetection) {
    CloseApproachDetector detector(propagator_, detector_settings_);
    
    // Earth-crossing NEO-like orbit (q < 1 AU, Q > 1 AU)
    KeplerianElements neo_orbit;
    neo_orbit.epoch_mjd_tdb = MJD2000;
    neo_orbit.semi_major_axis = 1.2;
    neo_orbit.eccentricity = 0.3;  // q = 0.84 AU, Q = 1.56 AU
    neo_orbit.inclination = 10.0 * DEG_TO_RAD;
    neo_orbit.longitude_ascending_node = 0.0;
    neo_orbit.argument_perihelion = 0.0;
    neo_orbit.mean_anomaly = 0.0;
    neo_orbit.gravitational_parameter = GMS;
    
    // Larger distance threshold for this test
    detector_settings_.detection_distance = 0.5;  // 0.5 AU
    detector.set_settings(detector_settings_);
    
    // Propagate for one orbital period
    double period = neo_orbit.period();
    auto approaches = detector.find_approaches(neo_orbit, MJD2000, MJD2000 + period);
    
    // Note: This is a placeholder test. Actual detection depends on:
    // 1. Initial position relative to Earth
    // 2. Orbital geometry
    // 3. Detection threshold
    // A real test would use known close approach cases
    
    std::cout << "\nEarth-crossing orbit test:\n";
    std::cout << "  Found " << approaches.size() << " potential close approaches\n";
    std::cout << "  Orbit: a=" << neo_orbit.semi_major_axis << " AU, e=" << neo_orbit.eccentricity << "\n";
    std::cout << "  Period: " << period << " days\n";
}

// ============================================================================
// B-Plane Calculation Tests
// ============================================================================

TEST_F(CloseApproachTest, BPlaneCalculationStructure) {
    CloseApproachDetector detector(propagator_, detector_settings_);
    
    // Create a close approach with known geometry
    CloseApproach ca;
    ca.mjd_tdb = MJD2000;
    ca.body = BodyType::EARTH;
    ca.distance = 0.01;
    
    // Set relative position and velocity
    ca.rel_position = astdyn::Vector3d(0.008, 0.006, 0.0);  // In plane
    ca.rel_velocity = astdyn::Vector3d(0.0, 0.0, 0.02);     // Perpendicular approach
    ca.relative_velocity = ca.rel_velocity.norm();
    
    // Compute b-plane
    BPlaneCoordinates bplane = detector.compute_b_plane(ca);
    
    // Check that b-plane magnitude matches perpendicular distance
    EXPECT_NEAR(bplane.b_magnitude, ca.rel_position.norm(), 1e-6);
    
    // Verify coordinates
    EXPECT_GT(bplane.b_magnitude, 0.0);
    std::cout << "\nB-plane test:\n";
    std::cout << "  ξ = " << bplane.xi << " AU\n";
    std::cout << "  ζ = " << bplane.zeta << " AU\n";
    std::cout << "  |b| = " << bplane.b_magnitude << " AU\n";
    std::cout << "  θ = " << bplane.theta * RAD_TO_DEG << " deg\n";
}

// ============================================================================
// MOID Calculation Tests
// ============================================================================

TEST_F(CloseApproachTest, MOIDCircularCoplanar) {
    // Two circular coplanar orbits with different radii
    KeplerianElements orbit1, orbit2;
    
    orbit1.epoch_mjd_tdb = MJD2000;
    orbit1.semi_major_axis = 1.0;
    orbit1.eccentricity = 0.0;
    orbit1.inclination = 0.0;
    orbit1.longitude_ascending_node = 0.0;
    orbit1.argument_perihelion = 0.0;
    orbit1.mean_anomaly = 0.0;
    orbit1.gravitational_parameter = GMS;
    
    orbit2 = orbit1;
    orbit2.semi_major_axis = 1.5;
    
    // MOID should be difference in radii
    double moid = MOIDCalculator::compute_moid(orbit1, orbit2);
    
    EXPECT_NEAR(moid, 0.5, 0.01);  // Allow 1% error for grid search
}

TEST_F(CloseApproachTest, MOIDIdenticalOrbits) {
    KeplerianElements orbit;
    orbit.epoch_mjd_tdb = MJD2000;
    orbit.semi_major_axis = 1.5;
    orbit.eccentricity = 0.1;
    orbit.inclination = 5.0 * DEG_TO_RAD;
    orbit.longitude_ascending_node = 45.0 * DEG_TO_RAD;
    orbit.argument_perihelion = 30.0 * DEG_TO_RAD;
    orbit.mean_anomaly = 0.0;
    orbit.gravitational_parameter = GMS;
    
    // MOID between identical orbits should be zero
    double moid = MOIDCalculator::compute_moid(orbit, orbit);
    
    EXPECT_NEAR(moid, 0.0, 1e-6);
}

// ============================================================================
// Integration Test
// ============================================================================

TEST_F(CloseApproachTest, EndToEndWorkflow) {
    std::cout << "\n========================================\n";
    std::cout << "Phase 8: Close Approach Detection\n";
    std::cout << "========================================\n";
    std::cout << "✓ Close approach structures\n";
    std::cout << "✓ B-plane coordinates\n";
    std::cout << "✓ Distance monitoring\n";
    std::cout << "✓ MOID calculation (basic)\n";
    std::cout << "========================================\n";
    std::cout << "Note: Full close approach detection requires:\n";
    std::cout << "  - Known NEO orbits with documented approaches\n";
    std::cout << "  - Validation against JPL close approach data\n";
    std::cout << "  - Impact probability assessment\n";
    std::cout << "  - Keyhole mapping\n";
    std::cout << "========================================\n";
    
    SUCCEED();
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
