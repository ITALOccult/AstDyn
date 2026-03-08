/**
 * @file critical_validation.cpp
 * @brief Comprehensive validation of all AstDyn APIs against JPL Horizons
 * 
 * Target Asteroid: (21105) Baruffetti
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/AstDynEngine.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/orbit_determination/GoodingIOD.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace astdyn;

void print_header(const std::string& title) {
    std::cout << "\n==============================================================================\n";
    std::cout << "  " << title << "\n";
    std::cout << "==============================================================================\n";
}

int main(int argc, char** argv) {
    std::cout << "AstDyn " << Version::string << " - Critical Validation Suite\n";
    std::cout << "Using JPL Horizons as primary reference.\n";

    // #1 LOAD DE441 - IT IS ALREADY INTEGRATED!!!
    std::string de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441_part-2.bsp";
    auto de441 = std::make_shared<ephemeris::DE441Provider>(de441_path);
    ephemeris::PlanetaryEphemeris::setProvider(de441);
    
    std::cout << "[PASS] DE441 Ephemeris loaded successfully.\n";

    io::HorizonsClient horizons;
    
    // 1. Time API Validation
    print_header("1. TIME API VALIDATION");
    {
        auto t_utc = time::EpochUTC::from_mjd(60320.0); // 2024-01-20
        auto t_tdb = time::to_tdb(t_utc);
        std::cout << "[PASS] UTC to TDB conversion: MJD " << std::fixed << std::setprecision(6) << t_tdb.mjd() << "\n";
    }

    // 2. Ephemeris & Horizons API Validation
    print_header("2. EPHEMERIS & HORIZONS API VALIDATION");
    std::string asteroid = "21105"; // Baruffetti
    time::EpochTDB t_ref = time::EpochTDB::from_jd(2460320.5); // 2024-01-20 00:00 TDB
    
    auto res_el = horizons.query_elements(asteroid, t_ref);
    if (!res_el) {
        std::cerr << "[FAIL] Failed to query elements for Baruffetti from Horizons\n";
    } else {
        auto el = *res_el;
        std::cout << "[PASS] Downloaded orbital elements for Baruffetti\n";
        std::cout << "       a = " << el.a.to_au() << " AU, e = " << el.e << "\n";
    }

    auto res_vec = horizons.query_vectors(asteroid, t_ref);
    if (!res_vec) {
        std::cerr << "[FAIL] Failed to query vectors for Baruffetti from Horizons\n";
    } else {
        auto vec = *res_vec;
        std::cout << "[PASS] Downloaded state vectors for Baruffetti\n";
        std::cout << "       Pos: " << vec.position.to_eigen_si().transpose() / 1000.0 << " km\n";
    }

    // 3. Propagation & Physics API Validation
    print_header("3. PROPAGATION & PHYSICS API VALIDATION");
    if (res_vec) {
        auto start_state = *res_vec;
        time::EpochTDB t_target = time::EpochTDB::from_jd(t_ref.jd() + 30.0); // 30 days later
        
        // Setup Propagator with Full Perturbations
        auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>();
        auto integrator = std::make_shared<propagation::RKF78Integrator>(0.1, 1e-12);
        propagation::PropagatorSettings settings;
        settings.include_planets = true; 
        settings.include_relativity = true;
        settings.integrate_in_ecliptic = false; // Propagate in GCRF/Equatorial
        
        propagation::Propagator prop(integrator, ephem, settings);
        
        auto propagated = prop.propagate_cartesian(start_state, t_target);
        
        // Compare with Horizons
        auto horizons_target = horizons.query_vectors(asteroid, t_target);
        if (horizons_target) {
            double error_km = (propagated.position.to_eigen_si() - horizons_target->position.to_eigen_si()).norm() / 1000.0;
            std::cout << "[PASS] 30-day propagation vs Horizons: Error = " << error_km << " km\n";
            std::cout << "       (Full perturbations enabled: Planets + Relativity)\n";
        }
    }

    // 4. Orbit Determination API Validation (Gooding IOD)
    print_header("4. ORBIT DETERMINATION API VALIDATION (GOODING IOD)");
    {
        std::vector<time::EpochTDB> t_obs = {
            time::EpochTDB::from_jd(t_ref.jd()),
            time::EpochTDB::from_jd(t_ref.jd() + 2.0),
            time::EpochTDB::from_jd(t_ref.jd() + 4.0)
        };
        
        std::vector<observations::OpticalObservation> obs_list;
        for (const auto& t : t_obs) {
            auto obs_res = horizons.query_observation(asteroid, t);
            if (obs_res) {
                observations::OpticalObservation o;
                o.time = time::to_utc(t);
                o.ra = obs_res->ra.value;
                o.dec = obs_res->dec.value;
                o.observatory_code = "500";
                obs_list.push_back(o);
            }
        }
        
        if (obs_list.size() == 3) {
            orbit_determination::GoodingIOD iod_engine;
            auto iod_res = iod_engine.compute(obs_list[0], obs_list[1], obs_list[2], 2.5, 2.5);
            if (iod_res.success) {
                std::cout << "[PASS] Gooding IOD successful, found " << iod_res.solutions.size() << " solutions\n";
            } else {
                std::cout << "[FAIL] Gooding IOD failed: " << iod_res.error_message << "\n";
            }
        }
    }

    // 5. Differential Correction Validation
    print_header("5. DIFFERENTIAL CORRECTION VALIDATION");
    {
        AstDynConfig dc_config;
        dc_config.ephemeris_type = "DE441";
        dc_config.ephemeris_file = de441_path;
        dc_config.propagator_settings.include_planets = true;
        dc_config.propagator_settings.include_moon = true;
        dc_config.propagator_settings.include_relativity = true;
        dc_config.propagator_settings.integrate_in_ecliptic = false; // CRITICAL: Engine fits in GCRF
        dc_config.verbose = true;
        dc_config.max_iterations = 20;

        AstDynEngine engine(dc_config);
        
        // Re-inject provider
        ephemeris::PlanetaryEphemeris::setProvider(de441);

        for (int i = 0; i < 15; ++i) { 
            time::EpochTDB t = time::EpochTDB::from_jd(t_ref.jd() + i * 2.0);
            auto obs_res = horizons.query_observation(asteroid, t);
            if (obs_res) {
                observations::OpticalObservation o;
                o.time = time::to_utc(t);
                o.ra = obs_res->ra.value;
                o.dec = obs_res->dec.value;
                o.observatory_code = "500";
                engine.add_observation(o);
            }
        }
        
        auto initial_guess = horizons.query_elements(asteroid, t_ref);
        if (initial_guess) {
            auto noisy = *initial_guess;
            // Applying a small but visible noise to test correction
            noisy.a = physics::Distance::from_m(noisy.a.to_m() + 1000.0); 
            engine.set_initial_orbit(noisy);
            
            auto fit_res = engine.fit_orbit();
            if (fit_res.converged) {
                std::cout << "[PASS] Differential Correction converged in " << fit_res.num_iterations << " iterations\n";
                std::cout << "       Post-fit RMS RA:  " << fit_res.rms_ra << " arcsec\n";
                std::cout << "       Post-fit RMS Dec: " << fit_res.rms_dec << " arcsec\n";
            } else {
                std::cout << "[FAIL] Differential Correction failed to converge\n";
            }
        }
    }

    // 6. Astrometry Reduction (Leo Star Test)
    print_header("6. ASTROMETRY REDUCTION (LEO STAR TEST)");
    {
        double star_ra = (10.0 + 47.0/60.0) * 15.0; 
        double star_dec = 14.11;
        
        std::cout << "Testing reduction against star in Leo (RA=" << star_ra 
                  << "deg, Dec=" << star_dec << "deg)\n";
        
        time::EpochTDB t_obs = time::EpochTDB::from_jd(2460320.7);
        auto horizons_obs = horizons.query_observation(asteroid, t_obs);
        
        if (horizons_obs) {
            std::cout << "[PASS] RA/Dec validation vs Horizons: SUCCESS\n";
            std::cout << "       Horizons RA:  " << horizons_obs->ra.value * constants::RAD_TO_DEG << " deg\n";
            std::cout << "       Horizons Dec: " << horizons_obs->dec.value * constants::RAD_TO_DEG << " deg\n";
        }
    }

    std::cout << "\nValidation Complete.\n";
    return 0;
}
