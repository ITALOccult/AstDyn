/**
 * @file test_fit_generic.cpp
 * @brief Fit orbitale su qualunque oggetto, per validare oltre il caso BK290.
 *
 * Uso:
 *   test_fit_generic <rwo> <a> <e> <i> <node> <omega> <M> <epoca_mjd> [obs_years]
 *
 * BK290 (riferimento noto: RMS 0.26 arcsec su 90 osservazioni):
 *   ... /tmp/820987.rwo 2.68777076554469 0.103208119159221 11.8524777541786 \
 *       253.159103985645 98.1398343280325 333.016364934386 61200
 *
 * Eros 433 (NEO, perielio 1.13 AU: la parallasse pesa molto di piu'):
 *   ... /tmp/433.rwo 1.4582437182 0.222877817 10.828546 \
 *       304.267954 178.918086 62.511490 61200 3
 */
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <astdyn/observations/ObservatoryDatabase.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>

using namespace astdyn;

int main(int argc, char** argv) {
    if (argc < 9) {
        std::cout << "uso: <rwo> <a> <e> <i> <node> <omega> <M> <epoca_mjd> [obs_years]\n";
        return 1;
    }
    const std::string rwo = argv[1];
    const double a = std::atof(argv[2]), ecc = std::atof(argv[3]), inc = std::atof(argv[4]);
    const double node = std::atof(argv[5]), om = std::atof(argv[6]), M = std::atof(argv[7]);
    const double epoca = std::atof(argv[8]);
    const double obs_years = (argc > 9) ? std::atof(argv[9]) : 0.0;

    // Catalogo osservatori: senza, si riduce dal geocentro e su un NEO
    // l'errore di parallasse e' enorme.
    if (const char* oc = std::getenv("ASTDYN_OBSCODES")) {
        size_t n = observations::ObservatoryDatabase::getInstance().loadFromMPCFile(oc);
        std::cout << "osservatori: " << n << "\n";
    } else {
        std::cout << "ATTENZIONE: nessun catalogo osservatori, riduzione dal geocentro\n";
    }

    AstDynEngine engine;
    AstDynConfig cfg = engine.config();
    cfg.verbose = false;
    cfg.integrator_type = IntegratorType::RKF78;
    cfg.tolerance = 1e-12;
    cfg.max_iterations = 20;
    cfg.fit_obs_years = obs_years;
    cfg.propagator_settings.include_planets = true;
    if (const char* e = std::getenv("ASTDYN_EPHEMERIS_PATH")) cfg.ephemeris_file = e;
    cfg.ephemeris_type = EphemerisType::DE441;
    engine.set_config(cfg);

    engine.set_initial_orbit(physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(epoca), a, ecc, inc, node, om, M));

    auto obs = observations::RWOReader::readFile(rwo);
    std::cout << "osservazioni nel file: " << obs.size() << "\n";

    // distribuzione per codice: quanti osservatori sono davvero in catalogo
    {
        const auto& db = observations::ObservatoryDatabase::getInstance();
        std::map<std::string,int> conteggio;
        int assenti = 0, totale = 0;
        for (const auto& o : obs) {
            ++totale;
            if (!db.getObservatory(o.observatory_code)) { ++assenti; conteggio[o.observatory_code]++; }
        }
        std::cout << "osservazioni da codici ASSENTI in catalogo: " << assenti
                  << " su " << totale;
        if (assenti > 0) {
            std::cout << "  (cadranno sul geocentro)";
            int mostrati = 0;
            std::cout << " —";
            for (const auto& [c, n] : conteggio) {
                if (mostrati++ >= 5) { std::cout << " ..."; break; }
                std::cout << " " << c << "(" << n << ")";
            }
        }
        std::cout << "\n";
    }

    for (const auto& o : obs) engine.add_observation(o);

    auto t0 = std::chrono::steady_clock::now();
    try {
        auto res = engine.fit_orbit();
        double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "\n=== RISULTATO ===\n"
                  << "converged  = " << (res.converged ? "SI" : "NO") << "\n"
                  << "motivo     = " << res.exit_reason << "\n"
                  << "iterazioni = " << res.num_iterations << "\n"
                  << "RMS        = " << res.rms_total_arcsec << " arcsec\n"
                  << "usate      = " << res.num_observations << "\n"
                  << "scartate   = " << res.num_rejected << "\n"
                  << "tempo      = " << std::fixed << std::setprecision(1) << s << " s\n";
        if (res.covariance.rows() == 6) {
            std::cout << "sigma pos  = ";
            for (int i = 0; i < 3; ++i)
                std::cout << std::scientific << std::setprecision(2)
                          << std::sqrt(std::abs(res.covariance(i,i))) * 149597870.7 << " km  ";
            std::cout << "\n";
        }
    } catch (const std::exception& e) {
        std::cout << "\nECCEZIONE: " << e.what() << "\n";
    }
    return 0;
}
