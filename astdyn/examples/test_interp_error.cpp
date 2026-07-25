/**
 * @file test_interp_error.cpp
 * @brief Misura l'errore introdotto dal dense output (interpolazione Hermite).
 *
 * Per ciascun istante di test: confronta la posizione ottenuta da
 *   - propagate_to      (propagazione ESATTA fino a quell'istante, nessuna interpolazione)
 *   - propagate_ephemeris (dense output: passo naturale + interpolazione Hermite)
 * La differenza e' l'errore di interpolazione, espresso in km e in arcsec
 * come visto dalla Terra (~2 AU di distanza tipica).
 *
 * Testa sia target isolati sia un gruppo "notturno" (osservazioni a minuti di
 * distanza), che e' il caso reale del fit.
 */
#include <astdyn/AstDynEngine.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace astdyn;

static const double AU_KM = 149597870.700;

int main() {
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
    auto cart0 = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(orbit);

    // Istanti reali delle osservazioni BK290 (gruppi notturni + salti)
    std::vector<double> mjds = {
        60302.629, 60302.640, 60302.652,   // notte 2023-12-24 (3 obs a minuti)
        60328.507, 60328.532, 60328.545,   // notte 2024-01-19
        60346.429, 60346.440, 60346.450,   // notte 2024-02-06
        60820.371, 60820.378, 60820.382,   // notte 2025-05-25 (gruppo fitto)
        60848.289, 60848.302, 60848.314    // notte 2025-06-22
    };

    // 1. dense output: tutti i target in una chiamata batch
    std::vector<time::EpochTDB> targets;
    for (double m : mjds) targets.push_back(time::EpochTDB::from_mjd(m));
    auto interp = engine.propagator()->propagate_ephemeris(cart0, targets);

    std::cout << "=== Errore di interpolazione (dense output vs propagazione esatta) ===\n";
    std::cout << std::left << std::setw(12) << "MJD"
              << std::setw(14) << "err [km]"
              << std::setw(14) << "err [arcsec]" << "\n";
    std::cout << std::string(40, '-') << "\n";

    double max_km = 0.0;
    for (size_t i = 0; i < mjds.size(); ++i) {
        // 2. propagazione esatta allo stesso istante (nessuna interpolazione)
        auto exact_k = engine.propagate_to(time::EpochTDB::from_mjd(mjds[i]));
        auto exact_c = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(exact_k);

        Eigen::Vector3d pe = exact_c.position.to_eigen_si();
        Eigen::Vector3d pi = interp[i].position.to_eigen_si();
        double err_km = (pi - pe).norm() / 1000.0;
        // errore angolare visto da ~2 AU
        double err_arcsec = (err_km / (2.0 * AU_KM)) * 206264.806;
        max_km = std::max(max_km, err_km);

        std::cout << std::left << std::fixed << std::setprecision(3) << std::setw(12) << mjds[i]
                  << std::setprecision(1) << std::setw(14) << err_km
                  << std::setprecision(4) << std::setw(14) << err_arcsec << "\n";
    }

    std::cout << std::string(40, '-') << "\n";
    std::cout << "max errore = " << std::fixed << std::setprecision(1) << max_km << " km ("
              << std::setprecision(4) << (max_km / (2.0 * AU_KM)) * 206264.806 << " arcsec)\n";
    std::cout << "\nRiferimento: RMS del fit = 2.03 arcsec, atteso AstDyS = 0.46 arcsec.\n";
    std::cout << "Se l'errore di interpolazione e' >~0.5 arcsec, e' lui il problema.\n";
    return 0;
}
