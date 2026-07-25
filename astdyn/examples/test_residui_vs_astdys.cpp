/**
 * @file test_residui_vs_astdys.cpp
 * @brief Confronta i NOSTRI residui con quelli di AstDyS, osservazione per
 *        osservazione. Il .rwo contiene i residui calcolati da AstDyS: sono
 *        l'unico riferimento onesto (l'RMS aggregato non e' confrontabile,
 *        perche' AstDyS pubblica un RMS normalizzato sui sigma).
 *
 * Lettura dei risultati:
 *   - nostri ~ loro                -> allineati, l'RMS alto era solo definizione
 *   - nostri ~ loro + costante     -> bias di catalogo non applicato
 *   - differenza crescente nel tempo -> modello dinamico che diverge
 *   - scorrelati                   -> problema nel modello di osservazione
 */
#include <astdyn/AstDynEngine.hpp>
#include <string>
#include <astdyn/observations/RWOReader.hpp>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace astdyn;

int main(int argc, char** argv) {
    const std::string rwo = (argc > 1) ? argv[1] : "/tmp/820987.rwo";

    // --- residui AstDyS dal .rwo (colonne a offset fisso, formato OrbFit v2) ---
    // RA:  rms a col 73, bias e residuo subito dopo il flag
    // Dec: idem nel blocco successivo
    struct AstDysRes { double mjd_frac; double res_ra, res_dec; };
    std::vector<AstDysRes> ref;
    {
        std::ifstream f(rwo);
        std::string line;
        while (std::getline(f, line)) {
            if (line.size() < 140) continue;
            if (line.substr(0, 8) != " 820987 ") continue;
            std::istringstream ss(line);
            std::vector<std::string> tok; std::string t;
            while (ss >> t) tok.push_back(t);
            // tok: id O C yyyy mm dd.ddddd acc hh mm ss.sss accRA rmsRA F biasRA resRA ...
            if (tok.size() < 23) continue;
            try {
                double day  = std::stod(tok[5]);
                // tok: [11]=rmsRA [12]=flag [13]=biasRA [14]=residuoRA
                //      [19]=rmsDec [20]=flag [21]=biasDec [22]=residuoDec
                double resRA  = std::stod(tok[14]);
                double resDec = std::stod(tok[22]);
                ref.push_back({day, resRA, resDec});
            } catch (...) { continue; }
        }
    }
    std::cout << "residui AstDyS letti dal .rwo: " << ref.size() << "\n";

    // --- i nostri residui, con l'orbita AstDyS (nessun fit) ---
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
        time::EpochTDB::from_mjd(61200.0),
        2.68777076554469, 0.103208119159221,
        11.8524777541786, 253.159103985645, 98.1398343280325, 333.016364934386);
    engine.set_initial_orbit(orbit);

    auto obs = observations::RWOReader::readFile(rwo);
    auto cart0 = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(orbit);
    orbit_determination::ResidualCalculator<core::ECLIPJ2000> rc(engine.getEphemeris(),
                                                                 engine.propagator());
    // argv[2] = "nolt" -> disattiva la correzione di tempo-luce (esperimento)
    bool use_lt = !(argc > 2 && std::string(argv[2]) == "nolt");
    rc.set_light_time_correction(use_lt);
    std::cout << "light-time correction: " << (use_lt ? "ON" : "OFF") << "\n";
    auto ours = rc.compute_residuals(obs, cart0);   // residui SENZA fit, orbita AstDyS

    std::cout << "nostri residui: " << ours.size() << "\n\n";
    std::cout << std::left << std::setw(12) << "MJD"
              << std::setw(11) << "nostro_RA" << std::setw(11) << "astdys_RA"
              << std::setw(11) << "nostro_De" << std::setw(11) << "astdys_De" << "\n";
    std::cout << std::string(56, '-') << "\n";

    size_t n = std::min(ours.size(), ref.size());
    double sum_dra = 0, sum_ddec = 0, sum_our_ra = 0, sum_ref_ra = 0;
    for (size_t i = 0; i < n; ++i) {
        double our_ra  = ours[i].residual_ra.to_arcsec();
        double our_dec = ours[i].residual_dec.to_arcsec();
        sum_dra  += our_ra - ref[i].res_ra;
        sum_ddec += our_dec - ref[i].res_dec;
        sum_our_ra += std::abs(our_ra);
        sum_ref_ra += std::abs(ref[i].res_ra);
        if (i < 10 || i >= n - 5) {
            std::cout << std::left << std::fixed << std::setprecision(3) << std::setw(12) << ref[i].mjd_frac
                      << std::setw(11) << our_ra  << std::setw(11) << ref[i].res_ra
                      << std::setw(11) << our_dec << std::setw(11) << ref[i].res_dec << "\n";
            if (i == 9) std::cout << "   ...\n";
        }
    }
    std::cout << std::string(56, '-') << "\n";
    std::cout << "offset medio  RA = " << std::setprecision(4) << sum_dra/n
              << "  Dec = " << sum_ddec/n << " arcsec\n";
    std::cout << "|res| medio: nostro RA = " << sum_our_ra/n
              << "   AstDyS RA = " << sum_ref_ra/n << " arcsec\n";
    return 0;
}
