/**
 * @file test_batch_vs_single.cpp
 * @brief Confronto DIRETTO propagate_to (integrate) vs propagate_ephemeris
 *        (integrate_at) sullo STESSO arco. Isola se il metodo batch e' il bug.
 */
#include <astdyn/AstDynEngine.hpp>
#include <cstdlib>
#include <iostream>
#include <chrono>

using namespace astdyn;

int main() {
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

    double target_mjd = 60302.6;  // ~2.5 anni indietro (dove il fit si blocca)

    // A. propagate_to (usa integrate, singolo target) -- atteso: liscio
    std::cout << "A. propagate_to(" << target_mjd << ") [integrate]... " << std::flush;
    auto t0 = std::chrono::high_resolution_clock::now();
    try {
        auto k = engine.propagate_to(time::EpochTDB::from_mjd(target_mjd));
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double,std::milli>(t1-t0).count();
        auto p = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(k).position.to_eigen_si();
        std::cout << "OK in " << ms << " ms, pos=(" << p.x()/1.496e11 << ","
                  << p.y()/1.496e11 << "," << p.z()/1.496e11 << ")\n";
    } catch (const std::exception& e) { std::cout << "ECCEZIONE: " << e.what() << "\n"; }

    // B. compute_ephemeris (usa propagate_ephemeris -> integrate_at, UN target)
    std::cout << "B. compute_ephemeris(" << target_mjd << ") [integrate_at, 1 target]... " << std::flush;
    t0 = std::chrono::high_resolution_clock::now();
    try {
        auto ev = engine.compute_ephemeris(time::EpochTDB::from_mjd(target_mjd),
                                           time::EpochTDB::from_mjd(target_mjd), 1.0);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double,std::milli>(t1-t0).count();
        if (!ev.empty()) {
            auto p = ev[0].position.to_eigen_si();
            std::cout << "OK in " << ms << " ms\n";
        }
    } catch (const std::exception& e) { std::cout << "ECCEZIONE: " << e.what() << "\n"; }

    // CASO C: compute_ephemeris con MOLTI target (step 30 giorni su 2.5 anni indietro).
    // Genera target multipli -> integrate_at con lista. Se esplode qui, e' l ordine/molteplicita.
    std::cout << "C. compute_ephemeris(60302..61200, step 30d) [integrate_at, ~30 target]... " << std::flush;
    t0 = std::chrono::high_resolution_clock::now();
    try {
        auto ev = engine.compute_ephemeris(time::EpochTDB::from_mjd(60302.6),
                                           time::EpochTDB::from_mjd(61200.0), 30.0);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double,std::milli>(t1-t0).count();
        std::cout << "OK in " << ms << " ms, " << ev.size() << " punti\n";
    } catch (const std::exception& e) { std::cout << "ECCEZIONE: " << e.what() << "\n"; }

    return 0;
}
