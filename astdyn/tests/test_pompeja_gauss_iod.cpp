/**
 * @file test_pompeja_gauss_iod.cpp
 * @brief Test Gauss Initial Orbit Determination using asteroid (203) Pompeja
 * 
 * Uses real observations from AstDyS to test the Gauss IOD implementation.
 * Compares computed orbit with reference orbit from OrbFit Fortran.
 * 
 * Data sources:
 * - Observations: https://newton.spacedys.com/~astdys2/mpcobs/numbered/0/203.rwo
 * - Orbital elements: https://newton.spacedys.com/~astdys2/epoch/numbered/0/203.eq1
 */

#include <gtest/gtest.h>
#include <astdyn/orbit_determination/GaussIOD.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/observations/Observation.hpp>
#include <astdyn/io/AstDynConfig.hpp>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace astdyn;
using namespace astdyn::orbit_determination;
using namespace astdyn::ephemeris;
using namespace astdyn::constants;
using namespace astdyn::observations;
using namespace astdyn::config;

/**
 * @brief Helper to download RWO file if not present
 */
static std::string get_pompeja_rwo_file() {
    std::string test_data_dir = std::string(std::getenv("HOME")) + 
        "/VisualStudio Code/GitHub/ITALOccultLibrary/astdyn/tests/data";
    std::string rwo_file = test_data_dir + "/pompeja_203.rwo";
    
    // Check if file exists
    std::ifstream check(rwo_file);
    if (!check.good()) {
        std::cout << "Downloading Pompeja observations from AstDyS...\n";
        std::string cmd = "mkdir -p '" + test_data_dir + "' && "
                         "curl -s 'https://newton.spacedys.com/~astdys2/mpcobs/numbered/0/203.rwo' -o '" 
                         + rwo_file + "'";
        int ret = std::system(cmd.c_str());
        if (ret != 0) {
            throw std::runtime_error("Failed to download RWO file");
        }
        std::cout << "✓ Downloaded to: " << rwo_file << "\n";
    }
    
    return rwo_file;
}

/**
 * @brief Test fixture for Pompeja Gauss IOD
 */
class PompejaGaussIODTest : public ::testing::Test {
protected:
    std::vector<OpticalObservation> observations;
    std::string rwo_file;
    
    void SetUp() override {
        // Download RWO file if needed
        rwo_file = get_pompeja_rwo_file();
        
        // Parse RWO file using RWOFileHandler
        auto rwo_observations = RWOFileHandler::read(rwo_file);
        
        // Extract optical observations
        for (const auto& rwo_obs : rwo_observations) {
            observations.push_back(rwo_obs.observation);
        }
        
        std::cout << "\n✓ Loaded " << observations.size() << " observations from RWO file\n";
    }
};

TEST_F(PompejaGaussIODTest, BasicGaussIOD) {
    std::cout << "\n========================================\n";
    std::cout << "  Gauss IOD Test: Pompeja (203)\n";
    std::cout << "========================================\n";
    std::cout << std::fixed << std::setprecision(6);
    
    // Print observations
    std::cout << "\nTotal observations: " << observations.size() << "\n";
    std::cout << "First observation: MJD " << observations.front().mjd_utc << "\n";
    std::cout << "Last observation:  MJD " << observations.back().mjd_utc << "\n";
    std::cout << "Time span: " << (observations.back().mjd_utc - observations.front().mjd_utc) 
              << " days ≈ " << (observations.back().mjd_utc - observations.front().mjd_utc) / 365.25 
              << " years\n";
    
    // Configure Gauss IOD
    GaussIODSettings settings;
    settings.max_iterations = 50;
    settings.tolerance = 1e-8;
    settings.use_light_time = true;
    settings.verbose = true;
    
    // Create Gauss IOD object
    GaussIOD gauss_iod(settings);
    
    // Run Gauss IOD using the 3 observations directly
    std::cout << "\nRunning Gauss IOD...\n";
    auto result = gauss_iod.compute_from_three(
        observations[0], observations[1], observations[2]);
    
    // Check result
    std::cout << "\n========================================\n";
    if (result.success) {
        std::cout << "✓ Gauss IOD SUCCEEDED\n";
        std::cout << "========================================\n\n";
        
        std::cout << "Convergence:\n";
        std::cout << "  Iterations: " << result.iterations << "\n";
        std::cout << "  Selected observations:\n";
        std::cout << "    Obs 1: index " << result.obs_index_1 
                  << ", MJD " << observations[result.obs_index_1].mjd_utc << "\n";
        std::cout << "    Obs 2: index " << result.obs_index_2 
                  << ", MJD " << observations[result.obs_index_2].mjd_utc << "\n";
        std::cout << "    Obs 3: index " << result.obs_index_3 
                  << ", MJD " << observations[result.obs_index_3].mjd_utc << "\n\n";
        
        std::cout << "Heliocentric Cartesian State (J2000 Ecliptic):\n";
        std::cout << std::setprecision(9);
        std::cout << "  Position [AU]:\n";
        std::cout << "    x = " << result.state.position(0) << "\n";
        std::cout << "    y = " << result.state.position(1) << "\n";
        std::cout << "    z = " << result.state.position(2) << "\n";
        std::cout << "    |r| = " << result.state.position.norm() << " AU\n\n";
        
        std::cout << "  Velocity [AU/day]:\n";
        std::cout << "    vx = " << result.state.velocity(0) << "\n";
        std::cout << "    vy = " << result.state.velocity(1) << "\n";
        std::cout << "    vz = " << result.state.velocity(2) << "\n";
        std::cout << "    |v| = " << result.state.velocity.norm() << " AU/day\n\n";
        
        std::cout << "Slant ranges [AU]:\n";
        std::cout << "    ρ₁ = " << result.slant_range_1 << "\n";
        std::cout << "    ρ₂ = " << result.slant_range_2 << "\n";
        std::cout << "    ρ₃ = " << result.slant_range_3 << "\n\n";
        
        // Compute orbital elements
        Eigen::Vector3d r = result.state.position;
        Eigen::Vector3d v = result.state.velocity;
        double r_norm = r.norm();
        double v_norm = v.norm();
        
        // Specific orbital energy
        double energy = 0.5 * v_norm * v_norm - GMS / r_norm;
        double a = -GMS / (2.0 * energy);
        
        // Angular momentum
        Eigen::Vector3d h = r.cross(v);
        double h_norm = h.norm();
        
        // Eccentricity vector
        Eigen::Vector3d e_vec = ((v_norm * v_norm - GMS / r_norm) * r - r.dot(v) * v) / GMS;
        double e = e_vec.norm();
        
        // Inclination
        double i = std::acos(h(2) / h_norm);
        
        // Node vector
        Eigen::Vector3d n(-h(1), h(0), 0.0);
        double n_norm = n.norm();
        
        // Longitude of ascending node
        double Omega = 0.0;
        if (n_norm > 1e-10) {
            Omega = std::acos(n(0) / n_norm);
            if (n(1) < 0.0) {
                Omega = TWO_PI - Omega;
            }
        }
        
        // Argument of perihelion
        double omega = 0.0;
        if (n_norm > 1e-10 && e > 1e-10) {
            omega = std::acos(n.dot(e_vec) / (n_norm * e));
            if (e_vec(2) < 0.0) {
                omega = TWO_PI - omega;
            }
        }
        
        // True anomaly
        double nu = 0.0;
        if (e > 1e-10) {
            nu = std::acos(e_vec.dot(r) / (e * r_norm));
            if (r.dot(v) < 0.0) {
                nu = TWO_PI - nu;
            }
        }
        
        // Mean anomaly
        double E = 2.0 * std::atan(std::sqrt((1.0 - e) / (1.0 + e)) * std::tan(nu / 2.0));
        double M = E - e * std::sin(E);
        if (M < 0.0) M += TWO_PI;
        
        // Period
        double P = TWO_PI * std::sqrt(a * a * a / GMS);
        
        std::cout << "Keplerian Elements:\n";
        std::cout << std::setprecision(6);
        std::cout << "  a = " << a << " AU\n";
        std::cout << "  e = " << e << "\n";
        std::cout << "  i = " << i * RAD_TO_DEG << "°\n";
        std::cout << "  Ω = " << Omega * RAD_TO_DEG << "°\n";
        std::cout << "  ω = " << omega * RAD_TO_DEG << "°\n";
        std::cout << "  M = " << M * RAD_TO_DEG << "°\n";
        std::cout << "  P = " << P << " days ≈ " << P/365.25 << " years\n";
        std::cout << "========================================\n\n";
        
        // Assertions for reasonable orbit
        EXPECT_GT(a, 2.0);  // Should be in main belt
        EXPECT_LT(a, 4.0);  // Not too far out
        EXPECT_GT(e, 0.0);  // Positive eccentricity
        EXPECT_LT(e, 0.3);  // Not too eccentric for main belt
        EXPECT_LT(i * RAD_TO_DEG, 30.0);  // Reasonable inclination
        EXPECT_GT(r_norm, 1.5);  // Not too close to Sun
        EXPECT_LT(r_norm, 5.0);  // Not too far
        
    } else {
        std::cout << "✗ Gauss IOD FAILED\n";
        std::cout << "========================================\n";
        std::cout << "Error: " << result.error_message << "\n";
        std::cout << "========================================\n\n";
        FAIL() << "Gauss IOD failed: " << result.error_message;
    }
}

TEST_F(PompejaGaussIODTest, CompareWithReferenceOrbit) {
    std::cout << "\n========================================\n";
    std::cout << "  Compare with Reference Orbit\n";
    std::cout << "========================================\n";
    
    // Reference orbit from AstDyS (equinoctial elements at MJD 61000.0)
    // Converting to approximate values at epoch of first observation (MJD 10708)
    // Note: This is a rough comparison since epochs differ significantly
    
    std::cout << "\nReference orbit (AstDyS @ MJD 61000.0):\n";
    std::cout << "  a = 2.7385 AU\n";
    std::cout << "  e = 0.0611\n";
    std::cout << "  i = 1.56°\n\n";
    
    // Run Gauss IOD
    GaussIODSettings settings;
    settings.max_iterations = 50;
    settings.tolerance = 1e-8;
    settings.use_light_time = true;
    settings.verbose = false;
    
    GaussIOD gauss_iod(settings);
    auto result = gauss_iod.compute_from_three(
        observations[0], observations[1], observations[2]);
    
    ASSERT_TRUE(result.success) << "Gauss IOD failed: " << result.error_message;
    
    // Compute orbital elements from result
    Eigen::Vector3d r = result.state.position;
    Eigen::Vector3d v = result.state.velocity;
    double r_norm = r.norm();
    double v_norm = v.norm();
    
    double energy = 0.5 * v_norm * v_norm - GMS / r_norm;
    double a = -GMS / (2.0 * energy);
    
    Eigen::Vector3d h = r.cross(v);
    double h_norm = h.norm();
    
    Eigen::Vector3d e_vec = ((v_norm * v_norm - GMS / r_norm) * r - r.dot(v) * v) / GMS;
    double e = e_vec.norm();
    
    double i = std::acos(h(2) / h_norm);
    
    std::cout << "Gauss IOD result:\n";
    std::cout << std::setprecision(6);
    std::cout << "  a = " << a << " AU\n";
    std::cout << "  e = " << e << "\n";
    std::cout << "  i = " << i * RAD_TO_DEG << "°\n\n";
    
    // Note: We expect some difference due to:
    // 1. Different epochs (1879 vs 2026)
    // 2. Only 3 observations used
    // 3. Perturbations not accounted for in Gauss method
    // 4. Reference orbit is from differential correction with many observations
    
    std::cout << "Differences:\n";
    std::cout << "  Δa = " << (a - 2.7385) << " AU (" 
              << 100.0 * (a - 2.7385) / 2.7385 << "%)\n";
    std::cout << "  Δe = " << (e - 0.0611) << " (" 
              << 100.0 * (e - 0.0611) / 0.0611 << "%)\n";
    std::cout << "  Δi = " << (i * RAD_TO_DEG - 1.56) << "°\n\n";
    
    std::cout << "Note: Significant differences expected due to:\n";
    std::cout << "  - Different epochs (150 years apart)\n";
    std::cout << "  - Only 3 observations vs full orbit fit\n";
    std::cout << "  - Two-body approximation in Gauss method\n";
    std::cout << "========================================\n\n";
    
    // Very loose tolerances since we're comparing different epochs
    EXPECT_GT(a, 2.0);
    EXPECT_LT(a, 4.0);
    EXPECT_GT(e, 0.0);
    EXPECT_LT(e, 0.5);
    EXPECT_LT(i * RAD_TO_DEG, 30.0);
}

TEST_F(PompejaGaussIODTest, RobustnessCheck) {
    std::cout << "\n========================================\n";
    std::cout << "  Robustness Check\n";
    std::cout << "========================================\n\n";
    
    // Test with different settings
    std::vector<int> max_iters = {10, 30, 50, 100};
    std::vector<double> tolerances = {1e-6, 1e-8, 1e-10};
    
    int successful_runs = 0;
    int total_runs = 0;
    
    for (auto max_iter : max_iters) {
        for (auto tol : tolerances) {
            total_runs++;
            
            GaussIODSettings settings;
            settings.max_iterations = max_iter;
            settings.tolerance = tol;
            settings.use_light_time = true;
            settings.verbose = false;
            
            GaussIOD gauss_iod(settings);
            auto result = gauss_iod.compute_from_three(
                observations[0], observations[1], observations[2]);
            
            if (result.success) {
                successful_runs++;
                
                // Check sanity of result
                double r_norm = result.state.position.norm();
                if (r_norm > 1.5 && r_norm < 5.0) {
                    std::cout << "✓ max_iter=" << max_iter 
                              << ", tol=" << tol 
                              << " → SUCCESS (|r|=" << std::setprecision(3) << r_norm << " AU)\n";
                } else {
                    std::cout << "⚠ max_iter=" << max_iter 
                              << ", tol=" << tol 
                              << " → SUSPICIOUS (|r|=" << std::setprecision(3) << r_norm << " AU)\n";
                }
            } else {
                std::cout << "✗ max_iter=" << max_iter 
                          << ", tol=" << tol 
                          << " → FAILED\n";
            }
        }
    }
    
    double success_rate = 100.0 * successful_runs / total_runs;
    std::cout << "\n========================================\n";
    std::cout << "Success rate: " << successful_runs << "/" << total_runs 
              << " (" << std::setprecision(1) << success_rate << "%)\n";
    std::cout << "========================================\n\n";
    
    EXPECT_GT(success_rate, 75.0) << "Success rate too low";
}
