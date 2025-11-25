/**
 * @file advanced_ephemeris_example.cpp
 * @brief Advanced ephemeris usage: JPL DE and asteroid perturbations
 * 
 * This example demonstrates:
 * 1. Using JPL DE441 for high-precision planetary positions
 * 2. Computing asteroid perturbations
 * 3. Comparing VSOP87 vs JPL DE accuracy
 * 4. Integrating with orbit propagation
 */

#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/ephemeris/JPLDEProvider.hpp>
#include <astdyn/ephemeris/AsteroidPerturbations.hpp>
#include <astdyn/ephemeris/EphemerisFactory.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;
using namespace astdyn::ephemeris;
using namespace astdyn::propagation;

void example1_jpl_de_usage() {
    std::cout << "\n========================================\n";
    std::cout << "Example 1: JPL DE441 High-Precision Ephemeris\n";
    std::cout << "========================================\n\n";
    
    // NOTE: Requires de441.bsp file downloaded from JPL
    // https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441.bsp
    
    try {
        // Create JPL DE441 provider
        auto jpl_provider = std::make_unique<JPLDEProvider>(
            "data/de441.bsp",
            EphemerisSource::JPL_DE441
        );
        
        // Test date: 2025-01-01 00:00 TDB
        double jd_tdb = 2460676.5;
        
        std::cout << "Planetary positions at J2025.0 (JPL DE441):\n";
        std::cout << std::fixed << std::setprecision(9);
        
        for (int i = 1; i <= 8; ++i) {  // Mercury to Neptune
            CelestialBody body = static_cast<CelestialBody>(i);
            auto pos = jpl_provider->getPosition(body, jd_tdb);
            double dist = pos.norm();
            
            std::cout << "  " << std::setw(8) << std::left;
            // Print body name (simplified)
            switch (i) {
                case 1: std::cout << "Mercury"; break;
                case 2: std::cout << "Venus"; break;
                case 3: std::cout << "Earth"; break;
                case 4: std::cout << "Mars"; break;
                case 5: std::cout << "Jupiter"; break;
                case 6: std::cout << "Saturn"; break;
                case 7: std::cout << "Uranus"; break;
                case 8: std::cout << "Neptune"; break;
            }
            std::cout << ": r = " << dist << " AU\n";
        }
        
        // Get Sun barycentric correction
        auto sun_offset = jpl_provider->getSunBarycentricPosition(jd_tdb);
        std::cout << "\nSun barycentric offset: "
                  << sun_offset.norm() * 149597870.7 << " km\n";
        
    } catch (const std::exception& e) {
        std::cout << "JPL DE441 not available: " << e.what() << "\n";
        std::cout << "Download de441.bsp from JPL NAIF to use this feature.\n";
    }
}

void example2_vsop_vs_jpl_comparison() {
    std::cout << "\n========================================\n";
    std::cout << "Example 2: VSOP87 vs JPL DE441 Comparison\n";
    std::cout << "========================================\n\n";
    
    double jd_tdb = 2460676.5;  // 2025-01-01
    
    try {
        auto jpl = std::make_unique<JPLDEProvider>("data/de441.bsp");
        
        std::cout << "Position differences (VSOP87 - JPL DE441):\n";
        std::cout << std::fixed << std::setprecision(1);
        
        for (int i = 1; i <= 8; ++i) {
            CelestialBody body = static_cast<CelestialBody>(i);
            
            // VSOP87 position
            auto vsop_pos = PlanetaryEphemeris::getPosition(body, jd_tdb);
            
            // JPL DE441 position
            auto jpl_pos = jpl->getPosition(body, jd_tdb);
            
            // Difference in km
            Eigen::Vector3d diff = (vsop_pos - jpl_pos) * 149597870.7;
            double diff_km = diff.norm();
            
            // Convert to angular separation (arcsec) at 1 AU
            double angle_arcsec = (diff_km / 149597870.7) * 206265.0;
            
            std::cout << "  Body " << i << ": " << std::setw(8) << diff_km 
                      << " km (" << angle_arcsec << " arcsec)\n";
        }
        
    } catch (const std::exception& e) {
        std::cout << "Comparison skipped: JPL DE441 not available.\n";
    }
}

void example3_asteroid_perturbations() {
    std::cout << "\n========================================\n";
    std::cout << "Example 3: Asteroid Perturbations (AST17)\n";
    std::cout << "========================================\n\n";
    
    // Load default AST17 asteroids (16 most massive)
    AsteroidPerturbations asteroids;
    
    // Display loaded asteroids
    std::cout << "Loaded " << asteroids.getAsteroids().size() << " asteroids:\n\n";
    std::cout << std::setw(4) << "Num" << " "
              << std::setw(15) << std::left << "Name"
              << std::setw(10) << std::right << "GM [km³/s²]"
              << std::setw(10) << "a [AU]"
              << std::setw(10) << "e\n";
    std::cout << std::string(60, '-') << "\n";
    
    std::cout << std::fixed << std::setprecision(3);
    for (const auto& ast : asteroids.getAsteroids()) {
        std::cout << std::setw(4) << ast.number << " "
                  << std::setw(15) << std::left << ast.name
                  << std::setw(10) << std::right << ast.gm
                  << std::setw(10) << ast.a
                  << std::setw(10) << ast.e << "\n";
    }
    
    std::cout << "\nTotal mass: " << std::scientific << std::setprecision(3)
              << asteroids.getTotalMass() << " M☉\n";
    
    // Compute perturbation at test position (1 AU from Sun)
    Eigen::Vector3d test_pos(1.0, 0.0, 0.0);  // 1 AU on x-axis
    double mjd_tdb = 60000.0;  // ~2023
    
    auto accel = asteroids.computePerturbation(test_pos, mjd_tdb);
    double accel_magnitude = accel.norm();
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nPerturbation at (1, 0, 0) AU:\n";
    std::cout << "  Acceleration: " << accel_magnitude << " AU/day²\n";
    std::cout << "  Equivalent: " << accel_magnitude * 1.731e6 << " m/s²\n";
    std::cout << "  As fraction of Sun's gravity: " 
              << accel_magnitude / 5.93e-6 << "\n";
    
    // Test individual asteroid contributions
    std::cout << "\nTop 3 contributors at this position:\n";
    std::vector<std::pair<std::string, double>> contributions;
    
    for (const auto& ast : asteroids.getAsteroids()) {
        auto ast_pos = asteroids.getPosition(ast, mjd_tdb);
        double gm_au3day2 = ast.gm * constants::GM_KM3S2_TO_AU3DAY2;
        auto single_accel = AsteroidPerturbations::computeSinglePerturbation(
            test_pos, ast_pos, gm_au3day2
        );
        contributions.push_back({ast.name, single_accel.norm()});
    }
    
    std::sort(contributions.begin(), contributions.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    for (int i = 0; i < 3 && i < contributions.size(); ++i) {
        std::cout << "  " << std::setw(15) << std::left << contributions[i].first
                  << ": " << std::scientific << contributions[i].second 
                  << " AU/day²\n";
    }
}

void example4_propagation_with_asteroids() {
    std::cout << "\n========================================\n";
    std::cout << "Example 4: Orbit Propagation with Asteroids\n";
    std::cout << "========================================\n\n";
    
    // NOTE: This is a conceptual example showing how asteroid perturbations
    // would be integrated into the propagator
    
    std::cout << "Propagator configuration with asteroid perturbations:\n\n";
    std::cout << "PropagatorSettings settings;\n";
    std::cout << "settings.include_planets = true;      // 8 major planets\n";
    std::cout << "settings.include_asteroids = true;    // 16 AST17 asteroids\n";
    std::cout << "settings.include_relativity = false;  // Optional\n\n";
    
    std::cout << "Expected impact on integration:\n";
    std::cout << "  - Extra force evaluations: +16 per step\n";
    std::cout << "  - Step size reduction: ~5%\n";
    std::cout << "  - Computation time: +10-20%\n\n";
    
    std::cout << "When to include asteroids:\n";
    std::cout << "  ✓ Inner solar system (< 3 AU)\n";
    std::cout << "  ✓ Propagation time > 1 year\n";
    std::cout << "  ✓ High-precision orbit determination\n";
    std::cout << "  ✗ Outer solar system (usually negligible)\n";
    std::cout << "  ✗ Short propagations (< 1 month)\n";
}

void example5_ephemeris_selection() {
    std::cout << "\n========================================\n";
    std::cout << "Example 5: Automatic Ephemeris Selection\n";
    std::cout << "========================================\n\n";
    
    std::cout << "Recommended ephemeris based on requirements:\n\n";
    
    struct Scenario {
        std::string name;
        std::string requirement;
        std::string recommended;
    };
    
    std::vector<Scenario> scenarios = {
        {"Quick propagation", "Speed, 1-20\" accuracy", "VSOP87 (built-in)"},
        {"High precision OD", "< 1\" accuracy, modern epoch", "JPL DE441"},
        {"Historical comets", "Epoch 1600-1800", "JPL DE405"},
        {"Long-term evolution", "10,000+ year propagation", "Custom N-body"},
        {"Near-Earth asteroids", "Close approach analysis", "JPL DE441 + AST17"},
        {"Main belt asteroid", "Standard orbit fit", "VSOP87 + AST17"},
        {"Trans-Neptunian", "Outer solar system", "JPL DE441 only"},
        {"Satellite mission", "Earth-Moon system", "JPL DE441 + lunar theory"}
    };
    
    std::cout << std::setw(25) << std::left << "Application"
              << std::setw(30) << "Requirement"
              << "Recommended\n";
    std::cout << std::string(80, '-') << "\n";
    
    for (const auto& s : scenarios) {
        std::cout << std::setw(25) << std::left << s.name
                  << std::setw(30) << s.requirement
                  << s.recommended << "\n";
    }
    
    std::cout << "\nPerformance metrics:\n";
    std::cout << "  VSOP87:     ~1 μs per call (analytical)\n";
    std::cout << "  JPL DE:     ~10 μs per call (interpolation + file I/O)\n";
    std::cout << "  AST17:      ~50 μs for all 16 (Kepler solver + transforms)\n";
}

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║   AstDyn C++ - Advanced Ephemeris Examples                   ║\n";
    std::cout << "║   JPL DE & Asteroid Perturbations                            ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n";
    
    example1_jpl_de_usage();
    example2_vsop_vs_jpl_comparison();
    example3_asteroid_perturbations();
    example4_propagation_with_asteroids();
    example5_ephemeris_selection();
    
    std::cout << "\n========================================\n";
    std::cout << "Examples complete!\n";
    std::cout << "========================================\n";
    std::cout << "\nFor more information:\n";
    std::cout << "  - See docs/EPHEMERIS_SUPPORT.md\n";
    std::cout << "  - JPL NAIF: https://naif.jpl.nasa.gov/\n";
    std::cout << "  - AST17 paper: Baer et al. (2011)\n\n";
    
    return 0;
}
