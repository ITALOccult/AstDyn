/**
 * @file test_validate_integrators.cpp
 * @brief Quadro completo degli integratori: propagazione di BK290 su piu' archi
 *        (indietro) confrontata con l'oracolo JPL Horizons, con errore e tempo.
 *
 * Oracoli JPL (eliocentrico eclittico J2000, AU) per 820987 (2015 BK290):
 *   MJD 61200.5 (epoca):  (+1.784625, -1.610900, +0.456416)
 *   MJD 60302.6 (-2.5 a): (-1.869013, +2.142570, -0.505718)
 *   MJD 59000.5 (-6.0 a): (-2.865115, -0.453251, -0.548393)
 *   MJD 57388.5 (-10 a):  (-2.867877, -0.428058, -0.550368)
 */
#include <astdyn/AstDynEngine.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>

using namespace astdyn;

struct Oracle { double mjd; double x, y, z; const char* label; };

static const std::vector<Oracle> ORACLES = {
    {61200.0,  1.781080307462242, -1.615391586772313,  0.455977411927073, "epoca"},
    {60300.0, -1.848937426673913,  2.158158630411470, -0.502632346804024, "-2.5anni"},
    {59000.0, -2.865115458643380, -0.453250895536448, -0.548392991206623, "-6anni"},
    {57388.0, -2.867876932046533, -0.428057998980700, -0.550367845691710, "-10anni"},
};

static const double AU_M = 1.495978707e11;

static void test_integrator(const char* name, IntegratorType itype, double tol) {
    for (const auto& orc : ORACLES) {
        AstDynEngine engine;
        AstDynConfig cfg = engine.config();
        cfg.verbose = false;
        cfg.integrator_type = itype;
        cfg.tolerance = tol;
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

        std::cout << std::left << std::setw(8) << name
                  << " tol=" << std::setw(7) << tol
                  << " " << std::setw(9) << orc.label << " ";
        try {
            auto t0 = std::chrono::high_resolution_clock::now();
            auto k = engine.propagate_to(time::EpochTDB::from_mjd(orc.mjd));
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            auto cart = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(k);
            auto p = cart.position.to_eigen_si();
            double x = p.x()/AU_M, y = p.y()/AU_M, z = p.z()/AU_M;
            double err = std::sqrt((x-orc.x)*(x-orc.x) + (y-orc.y)*(y-orc.y) + (z-orc.z)*(z-orc.z));
            std::cout << "err=" << std::scientific << std::setprecision(2) << err << " AU"
                      << "  t=" << std::fixed << std::setprecision(0) << std::setw(6) << ms << " ms"
                      << (err < 1e-5 ? "  [OK]" : (err < 1e-3 ? "  [~]" : "  [!!]")) << "\n";
        } catch (const std::exception& e) {
            std::cout << "ECCEZIONE: " << e.what() << "\n";
        }
    }
    std::cout << "\n";
}

int main() {
    std::cout << "=== Validazione integratori: BK290 dall'epoca MJD 61200.5 ===\n";
    std::cout << "Oracolo JPL. [OK]<1e-5 AU  [~]<1e-3  [!!]diverge\n\n";
    test_integrator("RKF78", IntegratorType::RKF78, 1e-12);
    test_integrator("RKF78", IntegratorType::RKF78, 1e-10);
    test_integrator("AAS",   IntegratorType::AAS,   1e-12);
    test_integrator("SABA4", IntegratorType::SABA4, 1e-12);
    test_integrator("RADAU", IntegratorType::RADAU, 1e-12);
    test_integrator("GAUSS", IntegratorType::GAUSS, 1e-12);
    test_integrator("GRKN64",IntegratorType::GRKN64,1e-12);
    test_integrator("RK4",   IntegratorType::RK4,   1e-12);
    return 0;
}
