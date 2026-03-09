/**
 * @file benchmark_multi_integrator.cpp
 * @brief Confronto AAS vs RKF78 su conservazione dell'energia
 *        per 4 asteroidi su 10 e 100 anni.
 *
 * Output: astdyn_energy.csv con colonne:
 *   asteroid, integrator, time_years, delta_H, delta_H_shadow, nfe
 *
 * Aggiungilo in tools/ e nel CMakeLists.txt:
 *   add_executable(benchmark_multi_integrator
 *       tools/benchmark_multi_integrator.cpp)
 *   target_link_libraries(benchmark_multi_integrator astdyn)
 *
 * Uso:
 *   ./benchmark_multi_integrator
 */

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <limits>

using namespace astdyn;

// ---------------------------------------------------------------------------
// Dati asteroide (elementi kepleriani + stato cartesiano iniziale da Horizons)
// ---------------------------------------------------------------------------
struct AsteroidData {
    std::string name;
    // Elementi kepleriani non usati direttamente — usiamo cartesiano da Horizons
    // Stato cartesiano eliocentrico J2000 GCRF da JPL Horizons
    double epoch_mjd;
    double x_au, y_au, z_au;       // posizione [AU]
    double vx_aud, vy_aud, vz_aud; // velocità [AU/day]
};

// Coordinate da JPL Horizons, epoca MJD 60310.0 (2024-Jan-01 TDB)
// Eliocentrico, J2000 GCRF
static const std::vector<AsteroidData> ASTEROIDS = {
    { 
        "Ceres", 60310.0,
        -1.0959125546760351, -2.3700844523978368, -0.8947830136316418,
        0.0089737783853470, -0.0037458210685192, -0.0035936597554197
    },
    { 
        "Icarus", 60310.0,
        1.0872131858370557, -0.8732191349734646, -0.8953721605377515,
        -0.0005195618936165, 0.0081797267817210, 0.0039326179149579
    },
    { 
        "Baruffetti", 60310.0,
        -2.8487026231700034, 0.6666974742659076, 0.8662121456383369,
        -0.0036002743896638, -0.0081998385956099, -0.0041380272752745
    },
    { 
        "Pallas", 60310.0,
        -0.7152424144751038, -0.0862557980736527, 0.0064391698273529,
        0.0020291931062381, -0.0183773350136283, -0.0083975508881355
    },
};

static const std::vector<double> DURATIONS_YR = { 10.0, 100.0 };
static const double SAMPLE_YR = 0.5;

// ---------------------------------------------------------------------------
// Costruzione stato iniziale dalla struttura dati
// ---------------------------------------------------------------------------
physics::CartesianStateTyped<core::GCRF> make_initial_state(const AsteroidData& ast)
{
    auto epoch = time::EpochTDB::from_mjd(ast.epoch_mjd);

    auto pos = math::Vector3<core::GCRF, physics::Distance>::from_si(
        physics::Distance::from_au(ast.x_au).to_m(),
        physics::Distance::from_au(ast.y_au).to_m(),
        physics::Distance::from_au(ast.z_au).to_m()
    );
    auto vel = math::Vector3<core::GCRF, physics::Velocity>::from_si(
        physics::Velocity::from_au_d(ast.vx_aud).to_ms(),
        physics::Velocity::from_au_d(ast.vy_aud).to_ms(),
        physics::Velocity::from_au_d(ast.vz_aud).to_ms()
    );

    return physics::CartesianStateTyped<core::GCRF>(
        epoch, pos, vel, physics::GravitationalParameter::sun());
}

// ---------------------------------------------------------------------------
// Energia specifica [m^2/s^2] — usa GM dal GravitationalParameter dello stato
// ---------------------------------------------------------------------------
template <typename Frame>
double specific_energy_si(const physics::CartesianStateTyped<Frame>& state)
{
    double mu  = state.gm.to_m3_s2();
    auto   r   = state.position.to_eigen_si();
    auto   v   = state.velocity.to_eigen_si();
    return 0.5 * v.squaredNorm() - mu / r.norm();
}

// ---------------------------------------------------------------------------
// Record energia
// ---------------------------------------------------------------------------
struct EnergyRecord {
    double time_years;
    double delta_H;
    double delta_H_shadow;
    long   nfe;
};

// ---------------------------------------------------------------------------
// Propagazione con campionamento energetico
// ---------------------------------------------------------------------------
std::vector<EnergyRecord> run_benchmark(
    const AsteroidData& ast,
    IntegratorType      integ_type,
    double              precision,
    double              duration_yr,
    bool                get_shadow,
    const AstDynConfig& base_config)
{
    AstDynConfig config    = base_config;
    config.integrator_type = integ_type;

    if (integ_type == IntegratorType::AAS)
        config.aas_precision = precision;
    else if (integ_type == IntegratorType::RKF78)
        config.tolerance = precision;

    auto cart0 = make_initial_state(ast);
    double E0  = specific_energy_si(cart0);

    AstDynEngine engine;
    engine.set_config(config);
    engine.set_initial_orbit(propagation::cartesian_to_keplerian(cart0));

    int n_steps = static_cast<int>(duration_yr / SAMPLE_YR);
    std::vector<EnergyRecord> records;
    records.reserve(n_steps + 1);
    records.push_back({ 0.0, 0.0, 0.0, 0L });

    long total_nfe = 0;

    for (int k = 1; k <= n_steps; ++k) {
        double target_mjd   = ast.epoch_mjd + k * SAMPLE_YR * 365.25;
        auto   target_epoch = time::EpochTDB::from_mjd(target_mjd);

        auto final_kep  = engine.propagate_to(target_epoch);
        auto final_cart = propagation::keplerian_to_cartesian(final_kep);

        double E  = specific_energy_si(final_cart);
        double dH = std::abs((E - E0) / std::abs(E0));

        double dHs = std::numeric_limits<double>::quiet_NaN();
        if (get_shadow) {
            try { dHs = engine.shadow_hamiltonian_drift(); } catch (...) {}
        }
        try { total_nfe = engine.total_force_evaluations(); } catch (...) {}

        records.push_back({ k * SAMPLE_YR, dH, dHs, total_nfe });
    }

    return records;
}

// ---------------------------------------------------------------------------
// CSV
// ---------------------------------------------------------------------------
void write_csv(std::ofstream& csv, const std::string& asteroid,
               const std::string& integ, const std::vector<EnergyRecord>& recs)
{
    for (const auto& r : recs) {
        csv << asteroid << "," << integ << ","
            << std::fixed      << std::setprecision(6) << r.time_years    << ","
            << std::scientific << std::setprecision(6) << r.delta_H       << ","
            << r.delta_H_shadow << "," << r.nfe << "\n";
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main()
{
    std::cout << "=== AstDyn Multi-Integrator Energy Benchmark ===\n\n";

    // Sanity check energia su Ceres
    {
        auto cart0 = make_initial_state(ASTEROIDS[0]);
        double E0  = specific_energy_si(cart0);
        std::cout << "Sanity check Ceres E0 = " << std::scientific
                  << E0 << " m^2/s^2\n";
        std::cout << "Atteso: ~ -4.43e+08 m^2/s^2  (negativo = orbita ellittica)\n\n";
        if (E0 > 0) {
            std::cerr << "ERRORE: energia positiva — verifica coordinate!\n";
            return 1;
        }
    }

    // Force model: due corpi puri
    AstDynConfig base_config;
    base_config.propagator_settings.include_planets    = false;
    base_config.propagator_settings.include_relativity = false;
    base_config.propagator_settings.include_moon       = false;
    base_config.propagator_settings.include_asteroids  = false;

    struct IntegConfig { std::string tag; IntegratorType type; double prec; bool shadow; };
    const std::vector<IntegConfig> INTEGRATORS = {
        { "AAS_1e-4", IntegratorType::AAS,   1e-4,  true  },
        { "AAS_1e-6", IntegratorType::AAS,   1e-6,  true  },
        { "AAS_1e-8", IntegratorType::AAS,   1e-8,  true  },
        { "RKF78",    IntegratorType::RKF78, 1e-10, false },
    };

    std::ofstream csv("astdyn_energy.csv");
    csv << "asteroid,integrator,time_years,delta_H,delta_H_shadow,nfe\n";

    auto t0_total = std::chrono::steady_clock::now();

    for (const auto& ast : ASTEROIDS) {
        for (double dur : DURATIONS_YR) {
            std::cout << "----------------------------------------------\n"
                      << ast.name << " — " << dur << " yr\n";

            for (const auto& ic : INTEGRATORS) {
                std::cout << "  [" << ic.tag << "] ..." << std::flush;
                auto t0 = std::chrono::steady_clock::now();
                try {
                    auto recs = run_benchmark(ast, ic.type, ic.prec,
                                              dur, ic.shadow, base_config);
                    write_csv(csv, ast.name, ic.tag, recs);
                    double elapsed = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t0).count();
                    std::cout << " OK  dH=" << std::scientific << std::setprecision(2)
                              << recs.back().delta_H
                              << "  t=" << std::fixed << std::setprecision(1)
                              << elapsed << "s\n";
                } catch (const std::exception& ex) {
                    std::cout << " ERRORE: " << ex.what() << "\n";
                }
            }
        }
    }

    csv.close();
    double total_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0_total).count();

    std::cout << "\nBenchmark completato in " << total_s << " s\n"
              << "Output: astdyn_energy.csv\n\n"
              << "Passo successivo:\n"
              << "  python benchmark_aas_vs_ias15.py --astdyn-csv astdyn_energy.csv\n";
    return 0;
}
