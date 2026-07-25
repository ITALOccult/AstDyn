/**
 * @file test_residui_generic.cpp
 * @brief Confronta i nostri residui con quelli AstDyS del .rwo, per QUALUNQUE oggetto.
 *
 * Uso:
 *   test_residui_generic <rwo> <id> <a> <e> <i> <node> <omega> <M> <epoca_mjd> [anno_min] [nolt]
 *
 * Esempi:
 *   BK290: ... /tmp/820987.rwo 820987 2.68777076554469 0.103208119159221 \
 *              11.8524777541786 253.159103985645 98.1398343280325 333.016364934386 61200
 *   Eros:  ... /tmp/433.rwo 433 1.4582437182 0.222877817 10.828546 \
 *              304.267954 178.918086 62.511490 61200 2020
 */
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/observations/ObservatoryDatabase.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

using namespace astdyn;

int main(int argc, char** argv) {
    if (argc < 10) {
        std::cout << "uso: <rwo> <id> <a> <e> <i> <node> <omega> <M> <epoca_mjd> [anno_min] [nolt]\n";
        return 1;
    }
    const std::string rwo = argv[1];
    const std::string id  = argv[2];
    double a = std::atof(argv[3]), ecc = std::atof(argv[4]), inc = std::atof(argv[5]);
    double node = std::atof(argv[6]), omega = std::atof(argv[7]), M = std::atof(argv[8]);
    double epoch_mjd = std::atof(argv[9]);
    int anno_min = (argc > 10) ? std::atoi(argv[10]) : 0;
    bool use_lt = !(argc > 11 && std::string(argv[11]) == "nolt");

    // --- residui AstDyS dal .rwo, indicizzati per data ---
    struct Ref { double mjd; double res_ra, res_dec; };
    std::vector<Ref> ref;
    {
        std::ifstream f(rwo);
        std::string line;
        const std::string prefix = " " + id + " ";
        while (std::getline(f, line)) {
            if (line.size() < 140) continue;
            std::istringstream ss(line);
            std::vector<std::string> tok; std::string t;
            while (ss >> t) tok.push_back(t);
            if (tok.size() < 23) continue;
            if (tok[0] != id) continue;
            try {
                int anno = std::stoi(tok[3]);
                if (anno < anno_min) continue;
                int mese = std::stoi(tok[4]);
                double giorno = std::stod(tok[5]);
                // data civile -> MJD (algoritmo di Fliegel-Van Flandern)
                int y = anno, m = mese;
                if (m <= 2) { y -= 1; m += 12; }
                int A = y / 100, B = 2 - A + A / 4;
                double jd = std::floor(365.25 * (y + 4716)) + std::floor(30.6001 * (m + 1))
                          + giorno + B - 1524.5;
                ref.push_back({jd - 2400000.5, std::stod(tok[14]), std::stod(tok[22])});
            } catch (...) { continue; }
        }
    }
    {
        auto& db = observations::ObservatoryDatabase::getInstance();
        const char* obsfile = std::getenv("ASTDYN_OBSCODES");
        if (obsfile) {
            size_t k = db.loadFromMPCFile(obsfile);
            std::cout << "osservatori caricati da file: " << k << "\n";
        }
        std::cout << "osservatori nel database: " << db.size() << "\n";
    }
    std::cout << "residui AstDyS letti: " << ref.size() << " (anno >= " << anno_min << ")\n";
    if (ref.empty()) { std::cout << "nessun residuo: controllare id/formato\n"; return 1; }

    // --- nostri residui con l'orbita AstDyS, senza fit ---
    AstDynEngine engine;
    AstDynConfig cfg = engine.config();
    cfg.verbose = false;
    cfg.integrator_type = IntegratorType::RKF78;
    cfg.tolerance = 1e-10;
    cfg.propagator_settings.include_planets = true;
    const char* eph = std::getenv("ASTDYN_EPHEMERIS_PATH");
    if (eph) cfg.ephemeris_file = eph;
    cfg.ephemeris_type = EphemerisType::DE441;
    engine.set_config(cfg);

    auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(epoch_mjd), a, ecc, inc, node, omega, M);
    auto cart0 = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(orbit);

    auto all_obs = observations::RWOReader::readFile(rwo);
    std::vector<observations::OpticalObservation> obs;
    // MJD minimo corrispondente a 1 gennaio di anno_min
    double mjd_min = -1e12;
    if (anno_min > 0) {
        int y = anno_min - 1, m = 13;
        int A = y / 100, B = 2 - A + A / 4;
        mjd_min = std::floor(365.25 * (y + 4716)) + std::floor(30.6001 * (m + 1))
                + 1.0 + B - 1524.5 - 2400000.5;
    }
    for (const auto& o : all_obs) if (o.time.mjd() >= mjd_min) obs.push_back(o);
    std::cout << "nostre osservazioni: " << obs.size() << "\n";

    orbit_determination::ResidualCalculator<core::ECLIPJ2000> rc(engine.getEphemeris(),
                                                                 engine.propagator());
    rc.set_light_time_correction(use_lt);
    auto ours = rc.compute_residuals(obs, cart0);
    std::cout << "light-time: " << (use_lt ? "ON" : "OFF") << "\n\n";

    // Accoppiamento per TEMPO (tolleranza 30 secondi): robusto anche quando i
    // due filtri non selezionano esattamente le stesse righe.
    double s_dra = 0, s_ddec = 0, s_our = 0, s_ref = 0;
    size_t n = 0, mostrate = 0;
    const double tol_mjd = 30.0 / 86400.0;
    std::cout << std::left << std::setw(11) << "nostro_RA" << std::setw(11) << "astdys_RA"
              << std::setw(11) << "nostro_De" << std::setw(11) << "astdys_De" << "\n";
    std::cout << std::string(44, '-') << "\n";
    for (size_t i = 0; i < ours.size(); ++i) {
        double t_our = ours[i].time.mjd();
        const Ref* best = nullptr; double bestdt = tol_mjd;
        for (const auto& r : ref) {
            double dt = std::abs(r.mjd - t_our);
            if (dt < bestdt) { bestdt = dt; best = &r; }
        }
        if (!best) continue;
        double our_ra = ours[i].residual_ra.to_arcsec();
        double our_de = ours[i].residual_dec.to_arcsec();
        s_dra += our_ra - best->res_ra;  s_ddec += our_de - best->res_dec;
        s_our += std::abs(our_ra);       s_ref  += std::abs(best->res_ra);
        ++n;
        if (mostrate++ < 6) std::cout << std::fixed << std::setprecision(3)
                             << std::setw(11) << our_ra << std::setw(11) << best->res_ra
                             << std::setw(11) << our_de << std::setw(11) << best->res_dec << "\n";
    }
    std::cout << std::string(44, '-') << "\n";
    std::cout << "accoppiate per tempo: " << n << " osservazioni\n";
    if (n == 0) { std::cout << "NESSUN accoppiamento: controllare id/date\n"; return 1; }
    std::cout << "offset medio  RA = " << std::setprecision(4) << s_dra/n
              << "   Dec = " << s_ddec/n << " arcsec\n";
    std::cout << "|res| medio: nostro = " << s_our/n << "   AstDyS = " << s_ref/n << " arcsec\n";
    return 0;
}
