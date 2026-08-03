/**
 * @file test_posizione.cpp
 * @brief Posizione apparente di un asteroide a un istante dato.
 *
 * Serve per confrontare la nostra effemeride con quella di altre predizioni
 * ALLO STESSO ISTANTE: confrontare distanze e moti apparenti a istanti diversi
 * dice solo che l'orbita e' ragionevole, non che le posizioni combacino.
 *
 * Uso:
 *   test_posizione <a> <e> <i> <node> <peri> <M> <epoca_mjd> <mjd_richiesto>
 *
 * Angoli in gradi, semiasse in AU.
 */
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/astrometry/Astrometry.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main(int argc, char** argv) {
    if (argc < 9) {
        std::cout << "uso: <a> <e> <i> <node> <peri> <M> <epoca_mjd> <mjd>\n";
        return 1;
    }
    const double a = std::atof(argv[1]), e = std::atof(argv[2]), inc = std::atof(argv[3]);
    const double node = std::atof(argv[4]), peri = std::atof(argv[5]), M = std::atof(argv[6]);
    const double ep = std::atof(argv[7]), t_req = std::atof(argv[8]);

    AstDynConfig cfg;
    cfg.verbose = false;
    cfg.integrator_type = IntegratorType::RKF78;
    cfg.tolerance = 1e-14;
    cfg.propagator_settings.include_planets = true;
    cfg.propagator_settings.include_relativity = true;
    // I perturbatori asteroidali contano: su un corpo di fascia principale
    // valgono qualche centesimo di arcsec al mese, cioe' decine di km sulla
    // traccia d'ombra.
    cfg.propagator_settings.include_asteroids = true;
    if (const char* af = std::getenv("ASTDYN_ASTEROID_FILE")) {
        cfg.asteroid_ephemeris_file = af;
        cfg.propagator_settings.asteroid_ephemeris_file = af;
    }
    if (const char* eph = std::getenv("ASTDYN_EPHEMERIS_PATH")) cfg.ephemeris_file = eph;
    cfg.ephemeris_type = EphemerisType::DE441;

    auto orbita = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(ep), a, e, inc, node, peri, M);

    const auto t = time::EpochTDB::from_mjd(t_req);

    astrometry::AstrometricSettings as;
    as.light_time_correction = true;
    as.aberrazione_differenziale = false;
    as.frame_conversion_to_equatorial = true;

    auto obs = astrometry::AstrometryReducer::compute_observation(
        orbita, orbita.epoch, t, cfg, as);

    if (!obs) { std::cout << "riduzione fallita\n"; return 1; }

    const double ra_deg = obs->ra.to_deg();
    const double dec_deg = obs->dec.to_deg();

    std::cout << std::fixed;
    std::cout << "epoca elementi MJD " << std::setprecision(4) << ep
              << "  ->  istante MJD " << t_req
              << "  (" << std::setprecision(1) << (t_req - ep) << " giorni)\n\n";
    std::cout << "RA   " << std::setprecision(7) << ra_deg << " deg"
              << "   = " << std::setprecision(8) << ra_deg / 15.0 << " h\n";
    std::cout << "Dec  " << std::setprecision(7) << dec_deg << " deg\n";
    std::cout << "dist " << std::setprecision(6) << obs->distance.to_au() << " AU\n";
    return 0;
}
