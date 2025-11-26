/**
 * @file test_single_residual_debug.cpp
 * @brief Debug singolo residuo per capire il problema
 */

#include <gtest/gtest.h>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <iostream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::observations;
using namespace astdyn::orbit_determination;

TEST(SingleResidualDebug, FirstObservation) {
    std::cout << "\n╔════════════════════════════════════════╗\n";
    std::cout << "║   DEBUG SINGOLO RESIDUO                ║\n";
    std::cout << "╚════════════════════════════════════════╝\n\n";
    
    // Load first observation
    auto obs_list = RWOReader::readFile("tools/203_recent100.rwo");
    ASSERT_FALSE(obs_list.empty());
    const auto& obs = obs_list[0];
    
    std::cout << "Osservazione:\n";
    std::cout << "  MJD UTC: " << std::fixed << std::setprecision(6) << obs.mjd_utc << "\n";
    std::cout << "  RA obs:  " << obs.ra * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  Dec obs: " << obs.dec * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  Obs code: " << obs.observatory_code << "\n\n";
    
    // Elements from 203_astdys.eq1 at MJD 61000.0
    KeplerianElements kep;
    kep.semi_major_axis = 2.7385249934;
    kep.eccentricity = 0.0610971811;
    kep.inclination = 3.172079 * constants::DEG_TO_RAD;
    kep.longitude_ascending_node = 347.595960 * constants::DEG_TO_RAD;
    kep.argument_perihelion = 59.961709 * constants::DEG_TO_RAD;
    kep.mean_anomaly = 64.765138 * constants::DEG_TO_RAD;
    kep.epoch_mjd_tdb = 61000.0;
    kep.gravitational_parameter = constants::GMS;
    
    // Convert to Cartesian
    auto cart = keplerian_to_cartesian(kep);
    
    std::cout << "Elementi orbitali (MJD 61000.0):\n";
    std::cout << "  a = " << kep.semi_major_axis << " AU\n";
    std::cout << "  e = " << kep.eccentricity << "\n";
    std::cout << "  Position: [" << cart.position[0] << ", " << cart.position[1] << ", " << cart.position[2] << "] AU\n\n";
    
    // Propagate to observation time
    auto ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
    auto integrator = std::make_unique<RK4Integrator>(0.1);
    PropagatorSettings settings;
    settings.include_planets = true;
    auto propagator = std::make_shared<Propagator>(std::move(integrator), ephemeris, settings);
    
    double mjd_tdb_obs = obs.mjd_utc + 69.184 / 86400.0;  // UTC to TDB approximation
    
    std::cout << "Propagazione da MJD " << kep.epoch_mjd_tdb << " a MJD " << mjd_tdb_obs << "\n";
    std::cout << "  Δt = " << (mjd_tdb_obs - kep.epoch_mjd_tdb) << " giorni\n\n";
    
    auto propagated = propagator->propagate_cartesian(cart, mjd_tdb_obs);
    
    std::cout << "Posizione propagata (MJD " << mjd_tdb_obs << "):\n";
    std::cout << "  Position: [" << propagated.position[0] << ", " << propagated.position[1] << ", " << propagated.position[2] << "] AU (eclittico)\n\n";
    
    // Compute residual
    auto residual_calc = std::make_shared<ResidualCalculator>(ephemeris, propagator);
    
    auto residual_opt = residual_calc->compute_residual(obs, propagated);
    ASSERT_TRUE(residual_opt.has_value());
    
    const auto& res = *residual_opt;
    
    std::cout << "RESIDUI:\n";
    std::cout << "  RA computed:  " << res.computed_ra * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  Dec computed: " << res.computed_dec * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  RA obs:       " << obs.ra * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  Dec obs:      " << obs.dec * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  Residual RA:  " << res.residual_ra * constants::RAD_TO_ARCSEC << " arcsec\n";
    std::cout << "  Residual Dec: " << res.residual_dec * constants::RAD_TO_ARCSEC << " arcsec\n";
    
    std::cout << "\nSE TUTTO VA BENE, i residui dovrebbero essere < 1 arcsec\n";
    std::cout << "SE C'È UN BUG, i residui saranno ~60000 arcsec (> 10 gradi)\n\n";
}
