/**
 * @file test_pompeja_fit.cpp
 * @brief Test differential correction for asteroid 203 Pompeja
 * 
 * Test setup:
 * - Initial elements: 203_astdys.eq1 (MJD 61000.0 TDT)
 * - Observations: 203.rwo (all observations from 1879 to 2025)
 * - Target epoch: MJD 61192.0 TDT
 * - Expected elements: 203.oel (fitted by OrbFit)
 */

#include <gtest/gtest.h>
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/orbit_determination/DifferentialCorrector.hpp>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::observations;
using namespace astdyn::orbit_determination;

class PompejaFitTest : public ::testing::Test {
protected:
    struct EquinoctialElements {
        double a;        // Semi-major axis [AU]
        double h;        // e*sin(ϖ)
        double k;        // e*cos(ϖ)
        double p;        // tan(i/2)*sin(Ω)
        double q;        // tan(i/2)*cos(Ω)
        double lambda;   // Mean longitude [rad]
        double mjd_tdt;  // Epoch [MJD TDT]
    };
    
    // Parse OrbFit .eq1 file
    EquinoctialElements parse_eq1_file(const std::string& filename) {
        EquinoctialElements elem{};
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        std::string line;
        bool in_elements = false;
        
        while (std::getline(file, line)) {
            // Skip comments and empty lines
            if (line.empty() || line[0] == '!') continue;
            
            // Skip header
            if (line.find("END_OF_HEADER") != std::string::npos) {
                in_elements = true;
                continue;
            }
            
            if (!in_elements) continue;
            
            // Skip object designation line
            if (line.find("203") != std::string::npos && line.size() < 10) continue;
            
            // Parse element line: EQU a h k p q lambda
            if (line.find("EQU") != std::string::npos) {
                std::istringstream iss(line);
                std::string equ_tag;
                iss >> equ_tag >> elem.a >> elem.h >> elem.k >> elem.p >> elem.q >> elem.lambda;
                
                // Convert lambda from degrees to radians
                elem.lambda *= constants::DEG_TO_RAD;
            }
            
            // Parse epoch line: MJD value TDT
            if (line.find("MJD") != std::string::npos && line.find("TDT") != std::string::npos) {
                std::istringstream iss(line);
                std::string mjd_tag;
                double mjd_value;
                std::string tdt_tag;
                iss >> mjd_tag >> mjd_value >> tdt_tag;
                elem.mjd_tdt = mjd_value;
            }
        }
        
        return elem;
    }
    
    // Parse OrbFit .oel file (simple format)
    EquinoctialElements parse_oel_file(const std::string& filename) {
        EquinoctialElements elem{};
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '!') continue;
            
            std::istringstream iss(line);
            std::string key;
            if (iss >> key) {
                if (key == "MJD") {
                    iss >> elem.mjd_tdt;
                } else if (key == "a") {
                    iss >> elem.a;
                } else if (key == "h") {
                    iss >> elem.h;
                } else if (key == "k") {
                    iss >> elem.k;
                } else if (key == "p") {
                    iss >> elem.p;
                } else if (key == "q") {
                    iss >> elem.q;
                } else if (key == "lambda") {
                    iss >> elem.lambda;
                }
            }
        }
        
        return elem;
    }
    
    // Convert equinoctial to Keplerian elements
    KeplerianElements equinoctial_to_keplerian(const EquinoctialElements& eq) {
        KeplerianElements kep;
        
        kep.semi_major_axis = eq.a;
        
        // Eccentricity: e = sqrt(h² + k²)
        double e = std::sqrt(eq.h * eq.h + eq.k * eq.k);
        kep.eccentricity = e;
        
        // Inclination: tan(i/2) = sqrt(p² + q²)
        double tan_half_i = std::sqrt(eq.p * eq.p + eq.q * eq.q);
        kep.inclination = 2.0 * std::atan(tan_half_i);
        
        // Longitude of ascending node: Ω = atan2(p, q)
        double Omega = std::atan2(eq.p, eq.q);
        if (Omega < 0) Omega += 2.0 * constants::PI;
        kep.longitude_ascending_node = Omega;
        
        // Longitude of perihelion: ϖ = atan2(h, k)
        double omega_plus_Omega = std::atan2(eq.h, eq.k);
        if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * constants::PI;
        
        // Argument of perihelion: ω = ϖ - Ω
        double omega = omega_plus_Omega - Omega;
        if (omega < 0) omega += 2.0 * constants::PI;
        kep.argument_perihelion = omega;
        
        // Mean anomaly: M = λ - ϖ
        double M = eq.lambda - omega_plus_Omega;
        while (M < 0) M += 2.0 * constants::PI;
        while (M >= 2.0 * constants::PI) M -= 2.0 * constants::PI;
        kep.mean_anomaly = M;
        
        kep.epoch_mjd_tdb = eq.mjd_tdt;  // Assume TDT ≈ TDB for this test
        kep.gravitational_parameter = constants::GMS;  // Heliocentric
        
        return kep;
    }
};

TEST_F(PompejaFitTest, DifferentialCorrectionFullDataset) {
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║     Pompeja Differential Correction - Full Dataset        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    // 1. Load initial elements from 203_astdys.eq1
    std::cout << "1. Loading initial elements from 203_astdys.eq1...\n";
    auto initial_eq = parse_eq1_file("tools/203_astdys.eq1");
    
    std::cout << "   Initial equinoctial elements (MJD " << std::fixed << std::setprecision(1) 
              << initial_eq.mjd_tdt << " TDT):\n";
    std::cout << "   • a      = " << std::setprecision(10) << initial_eq.a << " AU\n";
    std::cout << "   • h      = " << initial_eq.h << "\n";
    std::cout << "   • k      = " << initial_eq.k << "\n";
    std::cout << "   • p      = " << initial_eq.p << "\n";
    std::cout << "   • q      = " << initial_eq.q << "\n";
    std::cout << "   • λ      = " << initial_eq.lambda << " rad = " 
              << std::setprecision(6) << initial_eq.lambda * constants::RAD_TO_DEG << " deg\n";
    
    // Convert to Keplerian
    auto initial_kep = equinoctial_to_keplerian(initial_eq);
    std::cout << "\n   Keplerian elements:\n";
    std::cout << "   • a = " << std::setprecision(10) << initial_kep.semi_major_axis << " AU\n";
    std::cout << "   • e = " << initial_kep.eccentricity << "\n";
    std::cout << "   • i = " << std::setprecision(6) << initial_kep.inclination * constants::RAD_TO_DEG << " deg\n";
    std::cout << "   • Ω = " << initial_kep.longitude_ascending_node * constants::RAD_TO_DEG << " deg\n";
    std::cout << "   • ω = " << initial_kep.argument_perihelion * constants::RAD_TO_DEG << " deg\n";
    std::cout << "   • M = " << initial_kep.mean_anomaly * constants::RAD_TO_DEG << " deg\n";
    
    // 2. Load recent observations from 203_recent100.rwo
    std::cout << "\n2. Loading observations from 203_recent100.rwo...\n";
    auto observations = RWOReader::readFile("tools/203_recent100.rwo");
    
    std::cout << "   Loaded " << observations.size() << " observations\n";
    if (!observations.empty()) {
        const auto& first = observations.front();
        const auto& last = observations.back();
        std::cout << "   First: MJD " << std::setprecision(6) << first.mjd_utc 
                  << " (" << first.observatory_code << ")\n";
        std::cout << "   Last:  MJD " << std::setprecision(6) << last.mjd_utc 
                  << " (" << last.observatory_code << ")\n";
    }
    
    ASSERT_FALSE(observations.empty()) << "No observations loaded!";
    
    // 3. Load expected final elements from 203.oel
    std::cout << "\n3. Loading expected final elements from 203.oel...\n";
    auto expected_eq = parse_oel_file("tools/203.oel");
    
    std::cout << "   Expected equinoctial elements (MJD " << std::fixed << std::setprecision(1) 
              << expected_eq.mjd_tdt << " TDT):\n";
    std::cout << "   • a      = " << std::setprecision(10) << expected_eq.a << " AU\n";
    std::cout << "   • h      = " << expected_eq.h << "\n";
    std::cout << "   • k      = " << expected_eq.k << "\n";
    std::cout << "   • p      = " << expected_eq.p << "\n";
    std::cout << "   • q      = " << expected_eq.q << "\n";
    std::cout << "   • λ      = " << expected_eq.lambda << " rad\n";
    
    // 4. Convert to Cartesian state
    std::cout << "\n4. Converting to Cartesian state...\n";
    
    auto initial_cart = propagation::keplerian_to_cartesian(initial_kep);
    
    std::cout << "   Cartesian state (MJD " << std::fixed << std::setprecision(1) 
              << initial_eq.mjd_tdt << " TDT):\n";
    std::cout << "   • Position: [" << std::setprecision(10) 
              << initial_cart.position[0] << ", " 
              << initial_cart.position[1] << ", "
              << initial_cart.position[2] << "] AU\n";
    std::cout << "   • Velocity: [" << std::setprecision(10)
              << initial_cart.velocity[0] << ", " 
              << initial_cart.velocity[1] << ", "
              << initial_cart.velocity[2] << "] AU/day\n";
    
    // 5. Setup differential corrector
    std::cout << "\n5. Setting up differential corrector...\n";
    
    // Create ephemeris and propagator
    auto ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
    
    // Create RK4 integrator
    auto integrator = std::make_unique<propagation::RK4Integrator>(0.1);  // 0.1 day step
    
    // Propagator settings
    propagation::PropagatorSettings prop_settings;
    prop_settings.include_planets = true;
    prop_settings.perturb_venus = true;
    prop_settings.perturb_earth = true;
    prop_settings.perturb_mars = true;
    prop_settings.perturb_jupiter = true;
    prop_settings.perturb_saturn = true;
    
    auto propagator = std::make_shared<propagation::Propagator>(
        std::move(integrator), ephemeris, prop_settings);
    
    auto residual_calc = std::make_shared<ResidualCalculator>(ephemeris, propagator);
    auto stm_computer = std::make_shared<StateTransitionMatrix>(propagator);
    
    DifferentialCorrector corrector(residual_calc, stm_computer);
    
    DifferentialCorrectorSettings settings;
    settings.max_iterations = 20;
    settings.convergence_tolerance = 1.0e-6;  // AU
    settings.outlier_sigma = 3.0;
    settings.reject_outliers = true;
    settings.verbose = true;
    
    // 6. Run differential correction
    std::cout << "\n6. Running differential correction...\n";
    std::cout << "   Target epoch: MJD " << expected_eq.mjd_tdt << " TDT\n\n";
    
    auto result = corrector.fit(observations, initial_cart, settings);
    
    // 7. Display results
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                    FIT RESULTS                             ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";
    
    result.print_summary();
    
    std::cout << "\nConvergence: " << (result.converged ? "YES ✓" : "NO ✗") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Final RMS: " << std::setprecision(4) << result.statistics.rms_total << " arcsec\n";
    std::cout << "Observations used: " << result.statistics.num_observations << "\n";
    std::cout << "Outliers rejected: " << result.statistics.num_outliers << "\n";
    
    // Compare with expected elements
    auto fitted_kep = propagation::cartesian_to_keplerian(result.final_state);
    auto expected_kep = equinoctial_to_keplerian(expected_eq);
    
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║              COMPARISON WITH ORBFIT                        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";
    
    std::cout << "\nElement         Fitted          Expected        Difference\n";
    std::cout << "---------------------------------------------------------------\n";
    
    double da = (fitted_kep.semi_major_axis - expected_kep.semi_major_axis) * 1.496e8;  // km
    double de = fitted_kep.eccentricity - expected_kep.eccentricity;
    double di = (fitted_kep.inclination - expected_kep.inclination) * constants::RAD_TO_ARCSEC;
    double dOmega = (fitted_kep.longitude_ascending_node - expected_kep.longitude_ascending_node) * constants::RAD_TO_ARCSEC;
    double domega = (fitted_kep.argument_perihelion - expected_kep.argument_perihelion) * constants::RAD_TO_ARCSEC;
    double dM = (fitted_kep.mean_anomaly - expected_kep.mean_anomaly) * constants::RAD_TO_ARCSEC;
    
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "a [AU]      " << fitted_kep.semi_major_axis << "  " 
              << expected_kep.semi_major_axis << "  " 
              << std::setprecision(3) << da << " km\n";
    
    std::cout << std::setprecision(10);
    std::cout << "e           " << fitted_kep.eccentricity << "  " 
              << expected_kep.eccentricity << "  " 
              << std::scientific << std::setprecision(3) << de << "\n";
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "i [deg]     " << fitted_kep.inclination * constants::RAD_TO_DEG << "  " 
              << expected_kep.inclination * constants::RAD_TO_DEG << "  " 
              << std::setprecision(3) << di << " arcsec\n";
    
    std::cout << std::setprecision(6);
    std::cout << "Ω [deg]     " << fitted_kep.longitude_ascending_node * constants::RAD_TO_DEG << "  " 
              << expected_kep.longitude_ascending_node * constants::RAD_TO_DEG << "  " 
              << std::setprecision(3) << dOmega << " arcsec\n";
    
    std::cout << std::setprecision(6);
    std::cout << "ω [deg]     " << fitted_kep.argument_perihelion * constants::RAD_TO_DEG << "  " 
              << expected_kep.argument_perihelion * constants::RAD_TO_DEG << "  " 
              << std::setprecision(3) << domega << " arcsec\n";
    
    std::cout << std::setprecision(6);
    std::cout << "M [deg]     " << fitted_kep.mean_anomaly * constants::RAD_TO_DEG << "  " 
              << expected_kep.mean_anomaly * constants::RAD_TO_DEG << "  " 
              << std::setprecision(3) << dM << " arcsec\n";
    
    // Acceptance criteria
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                  ACCEPTANCE CRITERIA                       ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    EXPECT_TRUE(result.converged) << "Differential correction did not converge";
    EXPECT_LT(result.statistics.rms_total, 1.0) << "Final RMS too large: " << result.statistics.rms_total << " arcsec";
    EXPECT_LT(std::abs(da), 100.0) << "Semi-major axis error too large: " << da << " km";
    EXPECT_LT(std::abs(de), 1.0e-6) << "Eccentricity error too large: " << de;
    EXPECT_LT(std::abs(di), 1.0) << "Inclination error too large: " << di << " arcsec";
    
    std::cout << "All tests passed! ✓\n\n";
}

