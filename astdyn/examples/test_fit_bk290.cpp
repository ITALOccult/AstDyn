/**
 * @file test_fit_bk290.cpp
 * @brief Fit orbitale su BK290 con osservazioni reali AstDyS.
 *
 * Orbita iniziale = elementi AstDyS (epoca MJD 61200), osservazioni dal .rwo.
 * from_traditional vuole gli angoli in GRADI (non radianti!).
 * Filtro opzionale "ultimi N anni" via argv[2].
 * Atteso: converged=YES, RMS ~0.46 arcsec (RMSast dell'header .rwo).
 */
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/observations/ObservatoryDatabase.hpp>
#include <chrono>
#include <astdyn/observations/RWOReader.hpp>
#include <cstdlib>
#include <iostream>

using namespace astdyn;

int main(int argc, char** argv) {
    const std::string rwo_path = (argc > 1) ? argv[1] : "/tmp/820987.rwo";
    double last_years = (argc > 2) ? std::atof(argv[2]) : 0.0;

    // Engine con effemeride e perturbazioni planetarie.
    AstDynEngine engine;
    AstDynConfig cfg = engine.config();
    cfg.verbose = true;
    cfg.max_iterations = 20;
    cfg.fit_obs_years = last_years;
    cfg.fit_tolerance = (argc > 3) ? std::atof(argv[3]) : 0.0;  // 0 = tolleranza corrente  // usa il filtro della libreria, non quello locale
    cfg.tolerance = 1e-12;  // propagazione: fissa
    cfg.integrator_type = IntegratorType::RKF78;
    cfg.propagator_settings.include_planets = true;
    // Perturbatori asteroidali: AstDyS/OrbFit li includono sempre. Senza
    // asteroid_ephemeris_file, include_asteroids resta inefficace.
    cfg.asteroid_ephemeris_file = "/Users/michelebigi/.ioccultcalc/ephemerides/sb441-n16.bsp";
    cfg.propagator_settings.include_asteroids = true;
    const char* eph = std::getenv("ASTDYN_EPHEMERIS_PATH");
    if (eph) cfg.ephemeris_file = eph;
    cfg.ephemeris_type = EphemerisType::DE441;
    engine.set_config(cfg);

    // Orbita AstDyS di BK290 -- from_traditional vuole GRADI.
    auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(61200.0),
        2.68777076554469,      // a [AU]
        0.103208119159221,     // e
        11.8524777541786,      // i [deg]
        253.159103985645,      // node [deg]
        98.1398343280325,      // omega [deg]
        333.016364934386);     // M [deg]
    engine.set_initial_orbit(orbit);
    std::cout << "orbita: a=" << engine.orbit().a.to_au()
              << " e=" << engine.orbit().e
              << " i=" << engine.orbit().i.to_deg() << " deg\n";

    // Osservazioni dal .rwo (parser dedicato), con filtro opzionale.
    // Catalogo osservatori: senza, get_observer_position ricade sul geocentro
    // e si perde la parallasse topocentrica (~4" a 2.2 AU, molto di piu' sui
    // singoli casi).
    if (const char* oc = std::getenv("ASTDYN_OBSCODES")) {
        size_t n = observations::ObservatoryDatabase::getInstance().loadFromMPCFile(oc);
        std::cout << "osservatori caricati: " << n << "\n";
    } else {
        std::cout << "ATTENZIONE: ASTDYN_OBSCODES non impostata, osservatori dal geocentro\n";
    }
    auto obs = observations::RWOReader::readFile(rwo_path);
    std::cout << "osservazioni lette: " << obs.size() << "\n";
    int kept = 0;
    for (const auto& o : obs) {
        // filtro ora gestito dalla libreria via cfg.fit_obs_years
        engine.add_observation(o); ++kept;
    }
    std::cout << "osservazioni usate: " << kept << std::endl;

    // Fit.
    try {
        auto res = engine.fit_orbit();
        std::cout << "\n=== RISULTATO FIT ===\n";
        std::cout << "converged    = " << (res.converged ? "YES" : "NO") << "\n";
        std::cout << "iterations   = " << res.num_iterations << "\n";
        std::cout << "rms_total    = " << res.rms_total_arcsec << " arcsec\n";
        std::cout << "rms_ra       = " << res.rms_ra << " arcsec\n";
        std::cout << "rms_dec      = " << res.rms_dec << " arcsec\n";
        std::cout << "num_obs      = " << res.num_observations << "\n";
        std::cout << "num_rejected = " << res.num_rejected << "\n";
        std::cout << "(atteso: converged=YES, rms ~0.46 arcsec)\n";
    } catch (const std::exception& e) {
        std::cout << "FIT ECCEZIONE: " << e.what() << "\n";
    }
    return 0;
}
