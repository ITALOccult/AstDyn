/**
 * @file test_jacobian.cpp
 * @brief Valida la jacobiana ANALITICA vs NUMERICA confrontando Phi della STM
 *        su un arco breve di BK290. Se coincidono (~1e-4 rel), l'analitica e'
 *        validata e puo' diventare il default (risolve il collo di bottiglia fit).
 */
#include <astdyn/AstDynEngine.hpp>
#include <astdyn/orbit_determination/StateTransitionMatrix.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace astdyn;

int main() {
    AstDynEngine engine;
    AstDynConfig cfg = engine.config();
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
    auto cart0 = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(orbit);

    auto ephem = engine.getEphemeris();
    auto prop  = engine.propagator();

    // Arco breve: 60 giorni indietro (numerico completa in fretta)
    std::vector<time::EpochTDB> targets = { time::EpochTDB::from_mjd(60302.0) };
    // observer positions fittizie (non usate per Phi)
    std::vector<math::Vector3<core::GCRF, physics::Distance>> obs_pos = {
        math::Vector3<core::GCRF, physics::Distance>::from_si(1.5e11, 0.0, 0.0)
    };

    auto run = [&](bool numerical) {
        orbit_determination::StateTransitionMatrix<core::ECLIPJ2000> stm(prop, ephem);
        stm.set_use_numerical_jacobian(numerical);
        auto t0 = std::chrono::high_resolution_clock::now();
        auto res = stm.compute_batch(cart0, targets, obs_pos);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double,std::milli>(t1-t0).count();
        std::cout << (numerical ? "NUMERICA " : "ANALITICA") << " tempo="
                  << std::fixed << std::setprecision(1) << ms << " ms\n";
        return res[0].phi;
    };

    std::cout << "=== Validazione jacobiana STM: Phi(61200->61140), BK290 ===\n";
    auto phi_num = run(true);
    auto phi_ana = run(false);

    astdyn::Matrix6d diff = phi_ana - phi_num;
    double max_abs = diff.cwiseAbs().maxCoeff();
    double rel = max_abs / phi_num.cwiseAbs().maxCoeff();

    std::cout << "\nmax |Phi_ana - Phi_num| = " << std::scientific << std::setprecision(3) << max_abs << "\n";
    std::cout << "relativo               = " << rel << "\n";
    std::cout << (rel < 1e-4 ? "[OK] jacobiana analitica VALIDATA\n"
                            : "[!!] discrepanza: analitica da rivedere\n");
    return 0;
}
