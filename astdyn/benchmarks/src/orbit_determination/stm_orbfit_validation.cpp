/**
 * @file stm_orbfit_validation.cpp
 * @brief Validation benchmark for AstDyn vs OrbFit STM mapping
 */

#include <astdyn/AstDyn.hpp>
#include <astdyn/orbit_determination/StateTransitionMatrix.hpp>
#include <astdyn/propagation/OrbFitIntegrator.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/ephemeris/DE441Provider.hpp>
#include <astdyn/core/Constants.hpp>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

using namespace astdyn;
using namespace astdyn::orbit_determination;
using namespace astdyn::physics;
using namespace astdyn::propagation;

int main() {
    std::cout << "=== AstDyn vs OrbFit STM Validation ===" << std::endl;

    // 1. Setup Propagator with OrbFit Integrator and Numerical Jacobian
    auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>();
    const std::string de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de440s.bsp";
    try {
        ephem->setProvider(std::make_shared<ephemeris::DE441Provider>(de441_path));
    } catch (const std::exception& e) {
        std::cerr << "Error loading DE441: " << e.what() << std::endl;
        return 1;
    }

    auto integrator = std::make_shared<OrbFitDPIntegrator>(0.1, 1e-12);
    
    PropagatorSettings settings;
    settings.include_planets = true;
    settings.perturb_jupiter = true;
    settings.perturb_saturn = true;
    settings.central_body_gm = constants::GMS;
    settings.integrate_in_ecliptic = true;

    auto prop = std::make_shared<Propagator>(integrator, ephem, settings);
    
    // 2. Setup STM computer
    StateTransitionMatrix<core::ECLIPJ2000> stm_engine(prop);
    stm_engine.set_use_numerical_jacobian(true);
    stm_engine.set_differentiation_step(1e-7); // OrbFit default

    // 3. Initial State (Apophis-like)
    auto epoch0 = time::EpochTDB::from_mjd(58872.0);
    Eigen::VectorXd x0_au(6);
    x0_au << -0.10303, 0.9323, 0.3440, -1.8e-2, 1.4e-3, 4.6e-5;
    auto state0 = CartesianStateTyped<core::ECLIPJ2000>::from_au_aud(epoch0, x0_au, GravitationalParameter::from_au3_d2(constants::GMS));

    // 4. Propagate 100 days
    auto target_epoch = time::EpochTDB::from_mjd(epoch0.mjd() + 100.0);
    
    std::cout << "Propagating 100 days with Numerical Jacobian (OrbFit style)..." << std::endl;
    auto result = stm_engine.compute(state0, target_epoch);

    // 5. Output STM
    std::cout << std::setprecision(12) << std::fixed;
    std::cout << "\nPropagated STM (6x6):" << std::endl;
    std::cout << result.phi << std::endl;

    // 6. Validation Checks
    double det = result.phi.determinant();
    std::cout << "\nDeterminant: " << det << " (Expected near 1.0 for Hamiltonian system)" << std::endl;
    
    if (std::abs(det - 1.0) < 1e-4) {
        std::cout << "✅ STM Symplecticity Check Passed." << std::endl;
    } else {
        std::cout << "⚠️ STM Symplecticity Warning: Determinant deviation " << std::abs(det - 1.0) << std::endl;
    }

    return 0;
}
