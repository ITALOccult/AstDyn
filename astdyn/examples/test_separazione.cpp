/**
 * @file test_separazione.cpp
 * @brief Passo 4: la traccia del satellite e' nel posto giusto?
 *
 * Il controllo e' geometrico e non richiede riferimenti esterni: la separazione
 * angolare fra primario e satellite deve corrispondere al semiasse dell'orbita
 * mutua proiettato sul cielo, e deve oscillare con il periodo orbitale.
 *
 * Per Kalliope-Linus: a = 1099 km a una distanza di ~2 AU sottende
 *     1099 / (2 x 1.496e8) x 206265 = circa 0.76 arcsec
 * con periodo 3.5951 giorni.
 *
 * Se la separazione risultasse costante, il satellite non starebbe orbitando;
 * se fosse nulla, coinciderebbe con il primario; se avesse ampiezza sbagliata,
 * il semiasse o la distanza non tornerebbero.
 */
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/astrometry/Astrometry.hpp>
#include <astdyn/astrometry/ChebyshevEphemerisManager.hpp>
#include <astdyn/astrometry/OrbitingChebyshevEphemeris.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;
using astrometry::OrbitingChebyshevEphemeris;

int main() {
    const char* eph = std::getenv("ASTDYN_EPHEMERIS_PATH");
    if (!eph) {
        std::cout << "ASTDYN_EPHEMERIS_PATH non impostata: prova saltata\n";
        return 0;
    }

    AstDynConfig cfg;
    cfg.verbose = false;
    cfg.integrator_type = IntegratorType::RKF78;
    cfg.tolerance = 1e-11;
    cfg.propagator_settings.include_planets = true;
    cfg.ephemeris_file = eph;
    cfg.ephemeris_type = EphemerisType::DE441;

    // (22) Kalliope — elementi AstDyS, epoca MJD 61200
    auto kalliope = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        time::EpochTDB::from_mjd(61200.0),
        2.9095, 0.0996, 13.702, 66.02, 356.0, 120.0);

    OrbitingChebyshevEphemeris::OrbitaMutua linus;
    linus.a_km = 1099.0;
    linus.e = 0.005;
    linus.i_deg = 103.7;
    linus.node_deg = 100.6;
    linus.peri_deg = 0.0;
    linus.M_deg = 0.0;
    linus.period_days = 3.5951;
    linus.epoch = time::EpochTDB::from_mjd(61200.0);
    linus.piano = OrbitingChebyshevEphemeris::PianoRiferimento::Equatoriale;

    const auto inizio = time::EpochTDB::from_mjd(61253.0);   // 2026-08-01
    const auto fine   = time::EpochTDB::from_mjd(61259.0);   // sei giorni: ~1.7 periodi

    astrometry::ChebyshevEphemerisManager manager(cfg);
    manager.add_asteroid("22", kalliope, inizio, fine);
    manager.add_orbiting_body("22-S1", kalliope, linus, inizio, fine);

    // diagnostica: le due posizioni cartesiane allo stesso istante
    {
        const auto t = time::EpochTDB::from_mjd(61253.0);
        AstDynEngine m; m.set_config(cfg); m.set_initial_orbit(kalliope);
        auto k = m.propagate_to(t);
        auto c = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(k);
        Eigen::Vector3d rp = c.position.to_eigen_si() / 1000.0;
        Eigen::Vector3d rel = OrbitingChebyshevEphemeris::vettore_relativo(linus, t);
        std::cout << "primario  : " << rp.transpose() << "  |r| = " << rp.norm() << " km\n";
        std::cout << "relativo  : " << rel.transpose() << "  |r| = " << rel.norm() << " km\n";
        std::cout << "rapporto  : " << rel.norm() / rp.norm() << "\n\n";
        // confronto: cosa dice il manager per il primario allo stesso istante
        const auto [pp, vv] = manager.evaluate_full("22", t);
        const auto [ra_m, dec_m, dist_m] = pp;
        std::cout << "manager 22: RA=" << ra_m << " Dec=" << dec_m
                  << " dist=" << dist_m << " AU   (diretto: " << rp.norm()/1.495978707e8 << " AU)\n";
        const auto [ps, vs] = manager.evaluate_full("22-S1", t);
        const auto [ra_x, dec_x, dist_x] = ps;
        std::cout << "manager S1: RA=" << ra_x << " Dec=" << dec_x
                  << " dist=" << dist_x << " AU\n";
        // Il primario per le due strade: polinomi (via fit_chebyshev) contro
        // fornitore diretto. Se divergono, la separazione osservata e' quella.
        astrometry::AstrometricSettings as;
        as.light_time_correction = true;
        as.aberrazione_differenziale = false;
        as.frame_conversion_to_equatorial = true;
        auto oss = astrometry::AstrometryReducer::compute_observation(
            kalliope, kalliope.epoch, t, cfg, as);
        if (oss) {
            std::cout << "reducer 22: RA=" << oss->ra.to_deg()
                      << " Dec=" << oss->dec.to_deg()
                      << " dist=" << oss->distance.to_au() << " AU\n";
        }
    }
    std::cout << "=== separazione angolare Kalliope - Linus ===\n\n";
    std::cout << std::left << std::setw(12) << "MJD"
              << std::setw(13) << "dist [AU]"
              << std::setw(15) << "separaz [\"]"
              << std::setw(15) << "attesa max [\"]" << "\n";
    std::cout << std::string(55, '-') << "\n";

    double sep_min = 1e9, sep_max = 0.0;
    double dist_media = 0.0;
    int n = 0;

    for (double mjd = 61253.0; mjd <= 61259.0; mjd += 0.25) {
        const auto t = time::EpochTDB::from_mjd(mjd);
        const auto [pos_p, vel_p] = manager.evaluate_full("22", t);
        const auto [ra_p, dec_p, dist_p] = pos_p;
        const auto [pos_s, vel_s] = manager.evaluate_full("22-S1", t);
        const auto [ra_s, dec_s, dist_s] = pos_s;

        // separazione angolare, con la proiezione in ascensione retta
        const double d2r = 3.14159265358979323846 / 180.0;
        const double dra = (ra_s - ra_p) * std::cos(dec_p * d2r);
        const double ddec = dec_s - dec_p;
        const double sep = std::sqrt(dra * dra + ddec * ddec) * 3600.0;   // arcsec

        // semiasse proiettato alla distanza corrente
        const double attesa = 1099.0 / (dist_p * 1.495978707e8) * 206264.806;

        sep_min = std::min(sep_min, sep);
        sep_max = std::max(sep_max, sep);
        dist_media += dist_p; ++n;

        std::cout << std::left << std::fixed
                  << std::setprecision(2) << std::setw(12) << mjd
                  << std::setprecision(4) << std::setw(13) << dist_p
                  << std::setprecision(4) << std::setw(15) << sep
                  << std::setprecision(4) << std::setw(15) << attesa << "\n";
    }

    dist_media /= n;
    const double attesa_max = 1099.0 / (dist_media * 1.495978707e8) * 206264.806;

    std::cout << "\nseparazione minima : " << std::setprecision(4) << sep_min << " arcsec\n";
    std::cout << "separazione massima: " << sep_max << " arcsec\n";
    std::cout << "attesa (a proiettato a " << std::setprecision(3) << dist_media
              << " AU): " << std::setprecision(4) << attesa_max << " arcsec\n\n";

    bool ok = true;
    // la separazione massima deve avvicinarsi al semiasse proiettato
    if (std::abs(sep_max - attesa_max) / attesa_max > 0.15) {
        std::cout << "[!!] separazione massima lontana dall'attesa\n"; ok = false;
    }
    // Su un'orbita quasi circolare (e=0.005) la SEPARAZIONE resta costante: e' la
    // DIREZIONE a ruotare. Il controllo giusto e' sull'angolo di posizione, che
    // deve compiere un giro completo in un periodo orbitale.
    {
        double pa_prec = 0.0, giro = 0.0;
        bool primo = true;
        const double passo = linus.period_days / 40.0;
        for (int k = 0; k <= 40; ++k) {
            const auto t = time::EpochTDB::from_mjd(61253.0 + k * passo);
            const auto [pp2, vv2] = manager.evaluate_full("22", t);
            const auto [ps2, vs2] = manager.evaluate_full("22-S1", t);
            const auto [ra_a, dec_a, d_a] = pp2;
            const auto [ra_b, dec_b, d_b] = ps2;
            const double d2r = 3.14159265358979323846 / 180.0;
            const double dra = (ra_b - ra_a) * std::cos(dec_a * d2r);
            const double ddec = dec_b - dec_a;
            const double pa = std::atan2(dra, ddec) * 180.0 / 3.14159265358979323846;
            if (!primo) {
                double d = pa - pa_prec;
                while (d >  180.0) d -= 360.0;
                while (d < -180.0) d += 360.0;
                giro += d;
            }
            pa_prec = pa; primo = false;
        }
        std::cout << "rotazione dell'angolo di posizione in un periodo: "
                  << std::setprecision(1) << std::abs(giro) << " gradi (atteso 360)\n";
        if (std::abs(std::abs(giro) - 360.0) > 30.0) {
            std::cout << "[!!] il satellite non compie un giro completo nel periodo\n";
            ok = false;
        }
    }
    // e non deve mai coincidere con il primario per l'intero intervallo
    if (sep_max < 0.01) {
        std::cout << "[!!] satellite sovrapposto al primario\n"; ok = false;
    }

    std::cout << (ok ? "Geometria coerente." : "ATTENZIONE: geometria da rivedere.") << "\n";
    return ok ? 0 : 1;
}
