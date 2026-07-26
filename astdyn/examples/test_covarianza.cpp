/**
 * @file test_covarianza.cpp
 * @brief Confronta la covarianza calcolata dal fit con quella di AstDyS (.eq1).
 *
 * E' il confronto che decide se possiamo usare la nostra covarianza per le
 * ellissi di predizione (passo 3 della Fase 5) invece di importarla.
 *
 * Cose da guardare, in ordine:
 *  1. ORDINE DI GRANDEZZA delle deviazioni standard: se differiscono di
 *     ordini di grandezza c'e' un problema di unita' (SI contro AU/giorno).
 *  2. RAPPORTO fra le due: un fattore costante indica una normalizzazione
 *     diversa; un fattore 2 o sqrt(2) sarebbe sospetto (varianza vs sigma).
 *  3. STRUTTURA: le correlazioni dovrebbero somigliarsi anche se le scale no.
 */
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <astdyn/io/CovarianceIO.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;

int main(int argc, char** argv) {
    const std::string rwo = (argc > 1) ? argv[1] : "/tmp/820987.rwo";

    AstDynEngine engine;
    AstDynConfig cfg = engine.config();
    cfg.verbose = false;
    cfg.integrator_type = IntegratorType::RKF78;
    cfg.tolerance = 1e-12;
    cfg.propagator_settings.include_planets = true;
    const char* eph = std::getenv("ASTDYN_EPHEMERIS_PATH");
    if (eph) cfg.ephemeris_file = eph;
    cfg.ephemeris_type = EphemerisType::DE441;
    engine.set_config(cfg);

    auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(61200.0),
        2.68777076554469, 0.103208119159221,
        11.8524777541786, 253.159103985645, 98.1398343280325, 333.016364934386);
    engine.set_initial_orbit(orbit);

    auto obs = observations::RWOReader::readFile(rwo);
    for (const auto& o : obs) engine.add_observation(o);

    auto res = engine.fit_orbit();
    std::cout << "fit: converged=" << (res.converged ? "SI" : "NO")
              << "  RMS=" << res.rms_total_arcsec << " arcsec\n\n";

    if (res.covariance.rows() != 6) {
        std::cout << "covarianza non disponibile (righe: " << res.covariance.rows() << ")\n";
        return 1;
    }

    std::cout << "=== Covarianza dal fit: deviazioni standard sulla diagonale ===\n";
    for (int i = 0; i < 6; ++i) {
        double s = std::sqrt(std::abs(res.covariance(i, i)));
        const char* nome = (i < 3) ? "posizione" : "velocita'";
        std::cout << "  [" << i << "] " << nome << "  sigma = "
                  << std::scientific << std::setprecision(4) << s;
        if (i < 3) {
            std::cout << "   (se SI: " << s / 1000.0 << " km"
                      << " | se AU: " << s * 149597870.7 << " km)";
        }
        std::cout << "\n";
    }

    std::cout << "\n=== Correlazioni (adimensionali, indipendenti dalle unita') ===\n";
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            double d = std::sqrt(std::abs(res.covariance(i,i) * res.covariance(j,j)));
            double c = (d > 0) ? res.covariance(i,j) / d : 0.0;
            std::cout << std::fixed << std::setprecision(3) << std::setw(8) << c;
        }
        std::cout << "\n";
    }

    std::cout << "\nRiferimento AstDyS per BK290: l'ellisse di predizione risultante\n"
              << "e' 0.117721\" x 0.0644153\" @ PA 71.39 (cross-track 131.80 km).\n";
    return 0;
}
