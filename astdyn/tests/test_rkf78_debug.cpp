/**
 * Test diagnostico per investigare il problema di loop infinito in RKF78
 * Questo test cerca di capire perché RKF78 si blocca durante integrazioni lunghe
 */

#include <gtest/gtest.h>
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>
#include <chrono>
#include <thread>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::ephemeris;
using namespace astdyn::constants;

// Test con timeout per evitare che il test si blocchi indefinitamente
class TimeoutTest {
public:
    template<typename Func>
    static bool runWithTimeout(Func f, int timeout_seconds) {
        auto start = std::chrono::steady_clock::now();
        bool completed = false;
        
        std::thread t([&]() {
            f();
            completed = true;
        });
        
        t.detach();
        
        while (!completed) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
            if (elapsed > timeout_seconds) {
                return false; // Timeout
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        
        return true; // Completato
    }
};

// Test 1: Integrazione semplice con diagnostica
TEST(RKF78Debug, SimpleIntegrationWithDiagnostics) {
    std::cout << "\n=== Test 1: Simple Integration (1 day) ===" << std::endl;
    
    // Orbita semplice: Terra
    KeplerianElements kep;
    kep.epoch = time::EpochTDB::from_mjd(MJD2000);
    kep.semi_major_axis = 1.0;  // AU
    kep.eccentricity = 0.0167;
    kep.inclination = 0.0;
    kep.longitude_ascending_node = 0.0;
    kep.argument_perihelion = 0.0;
    kep.mean_anomaly = 0.0;
    kep.gravitational_parameter = GMS;
    
    auto cart = keplerian_to_cartesian(kep);
    
    // Setup propagatore
    auto integrator = std::make_unique<RKF78Integrator>(1.0, 1e-10);
    auto ephemeris = std::make_shared<PlanetaryEphemeris>();
    PropagatorSettings settings;
    settings.include_planets = false;  // Solo 2-body
    Propagator prop(std::move(integrator), ephemeris, settings);
    
    double target_mjd = MJD2000 + 1.0;  // 1 giorno
    
    std::cout << "Propagating from MJD " << MJD2000 << " to " << target_mjd << std::endl;
    
    try {
        auto result = prop.propagate_keplerian(kep, time::EpochTDB::from_mjd(target_mjd));
        std::cout << "SUCCESS: Propagation completed" << std::endl;
        std::cout << "Final a = " << result.semi_major_axis << " AU" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        FAIL() << "Propagation failed: " << e.what();
    }
}

// Test 2: Integrazione media (30 giorni)
TEST(RKF78Debug, MediumIntegration30Days) {
    std::cout << "\n=== Test 2: Medium Integration (30 days) ===" << std::endl;
    
    KeplerianElements kep;
    kep.epoch = time::EpochTDB::from_mjd(MJD2000);
    kep.semi_major_axis = 2.5;  // Asteroide nella main belt
    kep.eccentricity = 0.15;
    kep.inclination = 10.0 * DEG_TO_RAD;
    kep.longitude_ascending_node = 80.0 * DEG_TO_RAD;
    kep.argument_perihelion = 73.0 * DEG_TO_RAD;
    kep.mean_anomaly = 45.0 * DEG_TO_RAD;
    kep.gravitational_parameter = GMS;
    
    auto integrator = std::make_unique<RKF78Integrator>(1.0, 1e-10);
    auto ephemeris = std::make_shared<PlanetaryEphemeris>();
    PropagatorSettings settings;
    settings.include_planets = false;
    Propagator prop(std::move(integrator), ephemeris, settings);
    
    double target_mjd = MJD2000 + 30.0;
    
    std::cout << "Propagating asteroid orbit for 30 days" << std::endl;
    
    try {
        auto result = prop.propagate_keplerian(kep, time::EpochTDB::from_mjd(target_mjd));
        std::cout << "SUCCESS: Propagation completed" << std::endl;
        std::cout << "Initial a = " << kep.semi_major_axis << " AU" << std::endl;
        std::cout << "Final a = " << result.semi_major_axis << " AU" << std::endl;
        std::cout << "Delta a = " << (result.semi_major_axis - kep.semi_major_axis) << " AU" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        FAIL() << "Propagation failed: " << e.what();
    }
}

// Test 3: Integrazione lunga (100 giorni) con timeout
TEST(RKF78Debug, LongIntegration100DaysWithTimeout) {
    std::cout << "\n=== Test 3: Long Integration (100 days) with timeout ===" << std::endl;
    
    KeplerianElements kep;
    kep.epoch = time::EpochTDB::from_mjd(MJD2000);
    kep.semi_major_axis = 2.5;
    kep.eccentricity = 0.15;
    kep.inclination = 10.0 * DEG_TO_RAD;
    kep.longitude_ascending_node = 80.0 * DEG_TO_RAD;
    kep.argument_perihelion = 73.0 * DEG_TO_RAD;
    kep.mean_anomaly = 45.0 * DEG_TO_RAD;
    kep.gravitational_parameter = GMS;
    
    double target_mjd = MJD2000 + 100.0;
    
    std::cout << "Attempting 100-day propagation with 10-second timeout" << std::endl;
    
    bool completed = false;
    bool timeout_occurred = false;
    
    auto start = std::chrono::steady_clock::now();
    
    std::thread worker([&]() {
        try {
            auto integrator = std::make_unique<RKF78Integrator>(1.0, 1e-10);
            auto ephemeris = std::make_shared<PlanetaryEphemeris>();
            PropagatorSettings settings;
            settings.include_planets = false;
            Propagator prop(std::move(integrator), ephemeris, settings);
            
            auto result = prop.propagate_keplerian(kep, time::EpochTDB::from_mjd(target_mjd));
            completed = true;
        } catch (const std::exception& e) {
            std::cout << "ERROR in worker thread: " << e.what() << std::endl;
        }
    });
    
    // Aspetta con timeout
    const int timeout_seconds = 10;
    for (int i = 0; i < timeout_seconds * 10; ++i) {
        if (completed) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    if (!completed) {
        timeout_occurred = true;
        std::cout << "TIMEOUT after " << elapsed << " ms" << std::endl;
        std::cout << "This confirms the infinite loop issue!" << std::endl;
        
        // Detach il thread per evitare crash
        worker.detach();
    } else {
        std::cout << "SUCCESS: Completed in " << elapsed << " ms" << std::endl;
        worker.join();
    }
    
    // Il test passa se c'è timeout (conferma il bug)
    EXPECT_TRUE(timeout_occurred) << "Expected timeout to confirm infinite loop bug";
}

// Test 4: Integrazione con step size molto grande
TEST(RKF78Debug, LargeInitialStepSize) {
    std::cout << "\n=== Test 4: Large Initial Step Size ===" << std::endl;
    
    KeplerianElements kep;
    kep.epoch = time::EpochTDB::from_mjd(MJD2000);
    kep.semi_major_axis = 1.0;
    kep.eccentricity = 0.0167;
    kep.inclination = 0.0;
    kep.longitude_ascending_node = 0.0;
    kep.argument_perihelion = 0.0;
    kep.mean_anomaly = 0.0;
    kep.gravitational_parameter = GMS;
    
    // Prova con step iniziale molto grande (10 giorni)
    auto integrator = std::make_unique<RKF78Integrator>(10.0, 1e-10);
    auto ephemeris = std::make_shared<PlanetaryEphemeris>();
    PropagatorSettings settings;
    settings.include_planets = false;
    Propagator prop(std::move(integrator), ephemeris, settings);
    
    double target_mjd = MJD2000 + 10.0;
    
    std::cout << "Testing with initial step size = 10 days" << std::endl;
    
    try {
        auto result = prop.propagate_keplerian(kep, time::EpochTDB::from_mjd(target_mjd));
        std::cout << "SUCCESS with large step size" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "FAILED with large step: " << e.what() << std::endl;
    }
}

// Test 5: Confronto tolleranze
TEST(RKF78Debug, ToleranceComparison) {
    std::cout << "\n=== Test 5: Tolerance Comparison ===" << std::endl;
    
    KeplerianElements kep;
    kep.epoch = time::EpochTDB::from_mjd(MJD2000);
    kep.semi_major_axis = 2.5;
    kep.eccentricity = 0.15;
    kep.inclination = 10.0 * DEG_TO_RAD;
    kep.longitude_ascending_node = 80.0 * DEG_TO_RAD;
    kep.argument_perihelion = 73.0 * DEG_TO_RAD;
    kep.mean_anomaly = 45.0 * DEG_TO_RAD;
    kep.gravitational_parameter = GMS;
    
    double target_mjd = MJD2000 + 10.0;
    
    std::vector<double> tolerances = {1e-6, 1e-8, 1e-10, 1e-12};
    
    for (double tol : tolerances) {
        std::cout << "\nTrying tolerance = " << tol << std::endl;
        
        try {
            auto integrator = std::make_unique<RKF78Integrator>(1.0, tol);
            auto ephemeris = std::make_shared<PlanetaryEphemeris>();
            PropagatorSettings settings;
            settings.include_planets = false;
            Propagator prop(std::move(integrator), ephemeris, settings);
            
            auto start = std::chrono::steady_clock::now();
            auto result = prop.propagate_keplerian(kep, time::EpochTDB::from_mjd(target_mjd));
            auto end = std::chrono::steady_clock::now();
            
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            
            std::cout << "  SUCCESS in " << elapsed << " ms" << std::endl;
            std::cout << "  Delta a = " << (result.semi_major_axis - kep.semi_major_axis) << " AU" << std::endl;
            
        } catch (const std::exception& e) {
            std::cout << "  FAILED: " << e.what() << std::endl;
        }
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
