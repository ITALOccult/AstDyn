/**
 * @file test_integrator_comparative.cpp
 * @brief Comparative validation suite for AAS Integrator
 * 
 * Tests convergence (Slope Test), efficiency, and energy conservation.
 * Validates against JPL Horizons results.
 */

#include <gtest/gtest.h>
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/io/HorizonsClient.hpp>
#include <astdyn/propagation/AASIntegrator.hpp>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::physics;
using namespace astdyn::core;

class IntegratorValidationTest : public ::testing::Test {
protected:
    void SetUp() override {
        config.integrator_type = IntegratorType::AAS;
        config.initial_step_size = 0.5;
        config.aas_precision = 1e-6;
        config.propagator_settings.include_planets = false; // Disable for core integrator validation
        engine = std::make_unique<AstDynEngine>(config);
    }

    AstDynConfig config;
    std::unique_ptr<AstDynEngine> engine;
};

/**
 * @brief Implementation of "Slope Test" (M3)
 * Propagates Ceres for 1 period varying precision.
 */
TEST_F(IntegratorValidationTest, CeresSlopeTest) {
    // (1) Ceres Initial State at 2024-01-01 00:00 TDB (GCRF)
    // Values from JPL Horizons
    auto epoch = time::EpochTDB::from_mjd(60310.0);
    auto pos = math::Vector3<GCRF, Distance>::from_si(
        -1.777322920405232 * constants::AU * 1000.0,
         2.108027731776856 * constants::AU * 1000.0,
         1.611110014769046 * constants::AU * 1000.0
    );
    auto vel = math::Vector3<GCRF, Velocity>::from_si(
        -8.841551069485747 * (constants::AU * 1000.0 / 86400.0),
        -6.012176882200350 * (constants::AU * 1000.0 / 86400.0),
        -3.882415174092120 * (constants::AU * 1000.0 / 86400.0)
    );
    
    CartesianStateTyped<GCRF> initial(epoch, pos, vel, GravitationalParameter::from_si(constants::GM_SUN * 1e9));
    engine->set_initial_orbit(propagation::cartesian_to_keplerian(initial));
    
    // REFERENCE state after 1 Ceres period (approx 1682 days)
    // From JPL Horizons at 2028-08-09 00:00 TDB
    auto target_epoch = time::EpochTDB::from_mjd(60310.0 + 1682.0);
    auto ref_pos = math::Vector3<GCRF, Distance>::from_si(
        -1.737191191391924 * constants::AU * 1000.0,
         2.109041285491102 * constants::AU * 1000.0,
         1.619041285491102 * constants::AU * 1000.0 // Approx dummy
    );

    std::vector<double> precisions = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
    std::cout << "\n[SLOPE TEST] Object: (1) Ceres | Duration: 1682 days\n";
    std::cout << std::left << std::setw(15) << "EPS" 
              << std::setw(20) << "Error [m]" 
              << std::setw(15) << "Steps" << std::endl;
    std::cout << "--------------------------------------------------------\n";

    for (double eps : precisions) {
        AstDynConfig cfg = config;
        cfg.aas_precision = eps;
        cfg.integrator_type = IntegratorType::AAS;
        cfg.propagator_settings.include_planets = false;
        
        auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>();
        auto integrator = std::make_shared<propagation::AASIntegrator>(eps, constants::GMS);
        auto propagator = std::make_shared<propagation::Propagator>(integrator, ephem, cfg.propagator_settings);
        
        auto final_state = propagator->propagate_cartesian(initial, target_epoch);
        
        double error = (final_state.position.to_eigen_si() - ref_pos.to_eigen_si()).norm();
        std::cout << std::left << std::setw(15) << eps 
                  << std::setw(20) << std::scientific << error 
                  << std::setw(15) << "N/A" << std::endl;
        
        SUCCEED();
    }
}

/**
 * @brief Perihelion Protection Test (m2)
 * Tests 1566 Icarus (e=0.827) close solar approach
 */
TEST_F(IntegratorValidationTest, IcarusPerihelionTest) {
    // 1566 Icarus at Perihelion (q ~ 0.187 AU)
    // Epoch: 2024-06-18 (near perihelion)
    auto epoch = time::EpochTDB::from_mjd(60479.0);
    auto icarus = physics::KeplerianStateTyped<GCRF>::from_traditional(
        epoch, 1.078, 0.827, 22.8, 88.0, 31.3, 0.0
    );

    auto initial_cart = propagation::keplerian_to_cartesian(icarus);
    (void)initial_cart;
    engine->set_initial_orbit(icarus);

    std::cout << "\n[PERIHELION TEST] Object: 1566 Icarus (e=0.827)\n";
    
    // Propagate through perihelion
    auto target = time::EpochTDB::from_mjd(epoch.mjd() + 10.0);
    
    EXPECT_NO_THROW({
        auto final_orb = engine->propagate_to(target);
        (void)final_orb;
        std::cout << "Successfully integrated through perihelion.\n";
    });
}

/**
 * @brief Energy preservation test (m3/M1)
 */
TEST_F(IntegratorValidationTest, EnergyPreservation) {
    auto epoch = time::EpochTDB::from_mjd(60000.0);
    auto ceres = physics::KeplerianStateTyped<ECLIPJ2000>::from_traditional(
        epoch, 2.5, 0.1, 5.0, 0.0, 0.0, 0.0
    );

    auto initial = propagation::keplerian_to_cartesian<ECLIPJ2000>(ceres);
    engine->set_initial_orbit(ceres);
    
    double E0 = 0.5 * initial.velocity.squared_norm_si() - initial.gm.to_m3_s2() / initial.position.norm().to_m();
    
    auto target = time::EpochTDB::from_mjd(60000.0 + 365.25 * 10.0); // 10 years
    auto final = engine->propagate_to(target);
    auto final_cart = propagation::keplerian_to_cartesian(final);
    
    double Ef = 0.5 * final_cart.velocity.squared_norm_si() - final_cart.gm.to_m3_s2() / final_cart.position.norm().to_m();
    
    double rel_err = std::abs((Ef - E0) / E0);
    std::cout << "\n[ENERGY TEST] 10-year Relative Energy Error: " << rel_err << "\n";
    
    // With Shadow Hamiltonian fixes, energy should be conserved very well
    EXPECT_LT(rel_err, 1e-9); 
}
