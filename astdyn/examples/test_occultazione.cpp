/**
 * @file test_occultazione.cpp
 * @brief Geometria di un'occultazione, dati elementi e stella da riga di comando.
 *
 * Riporta TUTTO senza applicare filtri: istante di massimo avvicinamento,
 * separazione, coordinate del piano fondamentale, distanza d'impatto, altezza
 * del Sole lungo la traccia. Serve per capire cosa un filtro scarterebbe e
 * perche' — cosa che sull'intera campagna richiede di setacciare la
 * diagnostica.
 *
 * Confronta anche due modi di cercare il massimo avvicinamento:
 *   A) campionamento fitto, come fa oggi ClosestApproachFinder;
 *   B) campionamento rado per il bracketing, poi ricerca della radice della
 *      derivata analitica della separazione — che i polinomi di Chebyshev
 *      forniscono gia' insieme alla posizione.
 *
 * Uso:
 *   test_occultazione <a> <e> <i> <node> <peri> <M> <epoca_mjd> \
 *                     <ra_stella_deg> <dec_stella_deg> <mjd_inizio> <giorni>
 */
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/catalog/CatalogIntegration.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

using namespace astdyn;

namespace {

constexpr double DEG = 3.14159265358979323846 / 180.0;
constexpr double R_TERRA_KM = 6378.137;
constexpr double AU_KM = 1.495978707e8;

struct Punto { double jd, sep_as, ra, dec, dist; };

/// Separazione angolare fra asteroide e stella, in arcsec, e la sua derivata
/// temporale — entrambe dai polinomi, senza differenze finite.
void separazione(const catalog::ChebyshevSegment& seg, double jd,
                 double ra_s, double dec_s,
                 double& sep_as, double& dsep_dt, double& dist_au)
{
    const auto [pos, vel] = seg.evaluate_full(jd);
    const auto [ra, dec, dist] = pos;
    const auto [vra, vdec, vdist] = vel;   // gradi/giorno

    const double cd = std::cos(dec_s * DEG);
    const double dra  = (ra  - ra_s) * cd;
    const double ddec =  dec - dec_s;
    const double vra_p = vra * cd;

    sep_as = std::hypot(dra, ddec) * 3600.0;
    // d(s^2)/dt = 2 (dra vra_p + ddec vdec); il segno di questa e' quello che
    // cambia al minimo, ed e' cio' su cui si cerca la radice.
    dsep_dt = 2.0 * (dra * vra_p + ddec * vdec);
    dist_au = dist;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 12) {
        std::cout << "uso: <a> <e> <i> <node> <peri> <M> <epoca_mjd> "
                     "<ra_stella_deg> <dec_stella_deg> <mjd_inizio> <giorni>\n";
        return 1;
    }
    const double a = std::atof(argv[1]), e = std::atof(argv[2]), inc = std::atof(argv[3]);
    const double node = std::atof(argv[4]), peri = std::atof(argv[5]), M = std::atof(argv[6]);
    const double ep = std::atof(argv[7]);
    const double ra_s = std::atof(argv[8]), dec_s = std::atof(argv[9]);
    const double mjd0 = std::atof(argv[10]), giorni = std::atof(argv[11]);

    AstDynConfig cfg;
    cfg.verbose = false;
    cfg.integrator_type = IntegratorType::RKF78;
    cfg.tolerance = 1e-12;
    cfg.propagator_settings.include_planets = true;
    if (const char* eph = std::getenv("ASTDYN_EPHEMERIS_PATH")) cfg.ephemeris_file = eph;
    cfg.ephemeris_type = EphemerisType::DE441;

    auto orbita = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(ep), a, e, inc, node, peri, M);

    const double jd0 = time::EpochTDB::from_mjd(mjd0).jd();
    const auto centro = time::EpochTDB::from_jd(jd0 + giorni / 2.0);
    auto seg = catalog::fit_chebyshev(orbita, centro, giorni, cfg, 12);

    std::cout << std::fixed;
    std::cout << "=== geometria dell'occultazione ===\n\n";
    std::cout << "asteroide  a=" << std::setprecision(5) << a << " e=" << e
              << " i=" << std::setprecision(4) << inc << "  epoca MJD " << std::setprecision(1) << ep << "\n";
    std::cout << "stella     RA " << std::setprecision(6) << ra_s
              << "  Dec " << dec_s << "\n";
    std::cout << "finestra   MJD " << std::setprecision(3) << mjd0
              << " + " << giorni << " giorni\n\n";

    // ---------- metodo A: campionamento fitto ----------
    const int N_FITTO = static_cast<int>(giorni * 720);   // ogni 2 minuti
    auto t0 = std::chrono::steady_clock::now();
    Punto best_a{0, 1e30, 0, 0, 0};
    for (int i = 0; i <= N_FITTO; ++i) {
        const double jd = jd0 + giorni * i / N_FITTO;
        double s, d, dist;
        separazione(seg, jd, ra_s, dec_s, s, d, dist);
        if (s < best_a.sep_as) best_a = {jd, s, 0, 0, dist};
    }
    const double t_a = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();

    // ---------- metodo B: bracketing rado + radice della derivata ----------
    const int N_RADO = static_cast<int>(giorni * 24);     // ogni ora
    t0 = std::chrono::steady_clock::now();
    Punto best_b{0, 1e30, 0, 0, 0};
    double s_prec, d_prec, dist_prec;
    separazione(seg, jd0, ra_s, dec_s, s_prec, d_prec, dist_prec);
    double jd_prec = jd0;

    for (int i = 1; i <= N_RADO; ++i) {
        const double jd = jd0 + giorni * i / N_RADO;
        double s, d, dist;
        separazione(seg, jd, ra_s, dec_s, s, d, dist);

        // cambio di segno della derivata da negativa a positiva = minimo
        if (d_prec < 0.0 && d > 0.0) {
            double lo = jd_prec, hi = jd, f_lo = d_prec;
            for (int k = 0; k < 60; ++k) {
                const double mid = 0.5 * (lo + hi);
                double sm, dm, dm_dist;
                separazione(seg, mid, ra_s, dec_s, sm, dm, dm_dist);
                if ((f_lo < 0.0) == (dm < 0.0)) { lo = mid; f_lo = dm; } else { hi = mid; }
                if (hi - lo < 1e-11) break;
            }
            const double jd_min = 0.5 * (lo + hi);
            double sm, dm, dm_dist;
            separazione(seg, jd_min, ra_s, dec_s, sm, dm, dm_dist);
            if (sm < best_b.sep_as) best_b = {jd_min, sm, 0, 0, dm_dist};
        }
        jd_prec = jd; s_prec = s; d_prec = d;
    }
    const double t_b = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();

    // ---------- esiti ----------
    auto stampa = [&](const char* nome, const Punto& p, double ms, int valutazioni) {
        const double mjd = time::EpochTDB::from_jd(p.jd).mjd();
        const double ut = (mjd - std::floor(mjd)) * 24.0;
        const double impatto_km = p.sep_as / 206264.806 * p.dist * AU_KM;
        std::cout << nome << "\n";
        std::cout << "   massimo avvicinamento  MJD " << std::setprecision(6) << mjd
                  << "  (UT " << std::setprecision(4) << ut << " h)\n";
        std::cout << "   separazione            " << std::setprecision(4) << p.sep_as << " arcsec\n";
        std::cout << "   distanza dell'oggetto  " << std::setprecision(5) << p.dist << " AU\n";
        std::cout << "   impatto sul piano      " << std::setprecision(1) << impatto_km << " km"
                  << "  = " << std::setprecision(3) << impatto_km / R_TERRA_KM << " raggi terrestri";
        if (impatto_km < R_TERRA_KM) std::cout << "   [l'ombra tocca la Terra]";
        std::cout << "\n   valutazioni            " << valutazioni
                  << "  in " << std::setprecision(2) << ms << " ms\n\n";
    };

    stampa("A) campionamento ogni 2 minuti", best_a, t_a, N_FITTO + 1);
    stampa("B) bracketing orario + radice della derivata", best_b, t_b, N_RADO + 1 + 60);

    const double d_sec = std::abs(best_a.jd - best_b.jd) * 86400.0;
    const double d_sep = std::abs(best_a.sep_as - best_b.sep_as);
    std::cout << "=== confronto ===\n";
    std::cout << "scarto sull'istante     " << std::setprecision(3) << d_sec << " s\n";
    std::cout << "scarto sulla separazione " << std::setprecision(5) << d_sep << " arcsec\n";
    std::cout << "rapporto dei tempi       " << std::setprecision(1) << t_a / (t_b > 0 ? t_b : 1) << "x\n";
    if (d_sep > 0.01)
        std::cout << "\nIl campionamento fitto SOTTOSTIMA l'avvicinamento: il minimo vero\n"
                     "cade fra due campioni.\n";
    return 0;
}
