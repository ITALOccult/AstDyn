/**
 * @file test_integrators_bench.cpp
 * @brief Banco di prova degli integratori: accuratezza, costo, stabilita'.
 *
 * Riferimento: RKF78 a tolleranza stretta (validato contro JPL Horizons).
 *
 * La scala di prove e' crescente, cosi' il livello a cui un metodo fallisce
 * indica dove sta il difetto:
 *   1 passo        -> coefficienti / formula dello stepper
 *   archi brevi    -> controllo del passo
 *   archi lunghi   -> stabilita', accumulo dell'errore
 *   andata-ritorno -> reversibilita'
 *   energia        -> per i simplettici, deve restare limitata
 *
 * Uso:
 *   test_integrators_bench [metodo ...]
 *   senza argomenti prova i metodi veloci (RKF78 AAS RK4 GAUSS GRKN64);
 *   i lenti (SABA4 RADAU) vanno chiesti esplicitamente perche' possono
 *   richiedere minuti.
 *
 * Ogni riga e' stampata con flush: se un metodo si pianta, i risultati
 * precedenti restano leggibili e il processo si puo' interrompere.
 */
#include <astdyn/AstDynEngine.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include <cmath>

using namespace astdyn;

static const double EPOCA = 61200.0;
static const double AU_KM = 149597870.700;

struct Metodo { const char* nome; IntegratorType tipo; };
static const std::vector<Metodo> METODI = {
    {"RKF78",  IntegratorType::RKF78},
    {"AAS",    IntegratorType::AAS},
    {"RK4",    IntegratorType::RK4},
    {"GAUSS",  IntegratorType::GAUSS},
    {"GRKN64", IntegratorType::GRKN64},
    {"SABA4",  IntegratorType::SABA4},
    {"RADAU",  IntegratorType::RADAU},
};

// archi di prova [giorni]: da un passo singolo a un decennio
static const std::vector<double> ARCHI = {0.1, 1.0, 10.0, 100.0, 1000.0, 3650.0};

static physics::KeplerianStateTyped<core::ECLIPJ2000> orbita_bk290() {
    return physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(EPOCA),
        2.68777076554469, 0.103208119159221,
        11.8524777541786, 253.159103985645, 98.1398343280325, 333.016364934386);
}

static AstDynEngine make_engine(IntegratorType tipo, double tol) {
    AstDynEngine e;
    AstDynConfig c = e.config();
    c.verbose = false;
    c.integrator_type = tipo;
    c.tolerance = tol;
    c.propagator_settings.include_planets = true;
    const char* eph = std::getenv("ASTDYN_EPHEMERIS_PATH");
    if (eph) c.ephemeris_file = eph;
    c.ephemeris_type = EphemerisType::DE441;
    e.set_config(c);
    e.set_initial_orbit(orbita_bk290());
    return e;
}

struct Esito { bool ok; double err_au; double ms; std::string nota; };

static Esito propaga(IntegratorType tipo, double tol, double mjd_target,
                     const Eigen::Vector3d& riferimento) {
    try {
        auto e = make_engine(tipo, tol);
        auto t0 = std::chrono::high_resolution_clock::now();
        auto k = e.propagate_to(time::EpochTDB::from_mjd(mjd_target));
        auto t1 = std::chrono::high_resolution_clock::now();
        auto p = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(k)
                     .position.to_eigen_si() / (AU_KM * 1000.0);
        return {true, (p - riferimento).norm(),
                std::chrono::duration<double, std::milli>(t1 - t0).count(), ""};
    } catch (const std::exception& ex) {
        return {false, 0.0, 0.0, ex.what()};
    }
}

// energia a due corpi per unita' di massa: v^2/2 - mu/r [AU^2/d^2]
static double energia(const physics::KeplerianStateTyped<core::ECLIPJ2000>& k) {
    auto c = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(k);
    double r = c.position.to_eigen_si().norm() / (AU_KM * 1000.0);
    double v = c.velocity.to_eigen_si().norm() / (AU_KM * 1000.0) * 86400.0;
    const double mu = 2.959122082855911e-4;   // GM_Sole [AU^3/d^2]
    return 0.5 * v * v - mu / r;
}

int main(int argc, char** argv) {
    std::vector<std::string> scelti;
    for (int i = 1; i < argc; ++i) scelti.push_back(argv[i]);
    if (scelti.empty()) scelti = {"RKF78", "AAS", "RK4", "GAUSS", "GRKN64"};

    std::cout << "Banco di prova integratori — BK290 dall'epoca MJD " << EPOCA << "\n"
              << "Riferimento: RKF78 @ 1e-13\n"
              << "[OK] err<1e-6 AU   [~] err<1e-4   [!!] oltre\n" << std::endl;

    // riferimenti: RKF78 stretto su ogni arco
    std::cout << "calcolo dei riferimenti..." << std::flush;
    std::vector<Eigen::Vector3d> rif;
    for (double d : ARCHI) {
        auto e = make_engine(IntegratorType::RKF78, 1e-13);
        auto k = e.propagate_to(time::EpochTDB::from_mjd(EPOCA - d));
        rif.push_back(propagation::keplerian_to_cartesian<core::ECLIPJ2000>(k)
                          .position.to_eigen_si() / (AU_KM * 1000.0));
    }
    std::cout << " fatto\n" << std::endl;

    // ---- accuratezza e costo per arco ----
    for (const auto& m : METODI) {
        bool voluto = false;
        for (const auto& s : scelti) if (s == m.nome) voluto = true;
        if (!voluto) continue;

        std::cout << "=== " << m.nome << " ===" << std::endl;
        std::cout << std::left << std::setw(12) << "arco [g]"
                  << std::setw(14) << "errore [AU]"
                  << std::setw(12) << "tempo [ms]" << "esito" << std::endl;

        for (size_t i = 0; i < ARCHI.size(); ++i) {
            auto r = propaga(m.tipo, 1e-12, EPOCA - ARCHI[i], rif[i]);
            std::cout << std::left << std::fixed << std::setprecision(1)
                      << std::setw(12) << ARCHI[i];
            if (!r.ok) {
                std::cout << "ECCEZIONE: " << r.nota << std::endl;
                break;   // inutile insistere su archi piu' lunghi
            }
            std::cout << std::scientific << std::setprecision(2) << std::setw(14) << r.err_au
                      << std::fixed << std::setprecision(1) << std::setw(12) << r.ms
                      << (r.err_au < 1e-6 ? "[OK]" : r.err_au < 1e-4 ? "[~]" : "[!!]")
                      << std::endl;
        }

        // ---- andata e ritorno su 100 giorni ----
        try {
            auto e = make_engine(m.tipo, 1e-12);
            auto part = e.orbit();
            double e0 = energia(part);
            auto avanti = e.propagate_to(time::EpochTDB::from_mjd(EPOCA - 100.0));
            e.set_initial_orbit(avanti);
            auto ritorno = e.propagate_to(time::EpochTDB::from_mjd(EPOCA));
            auto pa = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(part)
                          .position.to_eigen_si();
            auto pr = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(ritorno)
                          .position.to_eigen_si();
            double scarto = (pr - pa).norm() / (AU_KM * 1000.0);
            double de = std::abs((energia(ritorno) - e0) / e0);
            std::cout << "andata-ritorno 100 g: scarto " << std::scientific
                      << std::setprecision(2) << scarto << " AU"
                      << "   deriva energia " << de << std::endl;
        } catch (const std::exception& ex) {
            std::cout << "andata-ritorno: ECCEZIONE " << ex.what() << std::endl;
        }
        std::cout << std::endl;
    }
    return 0;
}
