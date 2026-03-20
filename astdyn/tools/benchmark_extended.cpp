/**
 * @file benchmark_extended.cpp
 * @brief Extended Energy Conservation Benchmark — 4 Orbital Categories
 *
 * Categorie:
 *   1. NEA/PHA   — 50 yr    (Apophis, Icarus, Bennu)
 *   2. Troiani   — 1000 yr  (Achilles, Patroclus, Hektor)
 *   3. Risonanti — 500 yr   (Hilda, Hecuba, Griqua)
 *   4. TNO/KBO   — 10000 yr (Pluto, Eris, Sedna)
 *
 * Integratori: AAS ε=1e-4, AAS ε=1e-6, AAS ε=1e-8, RKF78
 *
 * Output: astdyn_energy_extended.csv
 *   colonne: asteroid,integrator,time_years,delta_H,delta_H_shadow,nfe
 */

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/core/physics_state.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <filesystem>
#include <chrono>

using namespace astdyn;
using namespace astdyn::physics;
using namespace astdyn::core;
namespace fs = std::filesystem;

// ============================================================
// Struttura asteroide
// ============================================================

struct AsteroidConfig {
    std::string name;
    std::string category;
    double duration_yr;   // anni di propagazione
    double sample_yr;     // campionamento (1 yr per tutti)

    // Stato iniziale: coordinate eliocentriche @sun, MJD 60310.0
    // Posizione in AU, velocità in AU/day
    double x_au, y_au, z_au;
    double vx_auday, vy_auday, vz_auday;
    double eccentricity;
};

// ============================================================
// Definizione dei 12 asteroidi
// Coordinate: eliocentriche ICRF/J2000, MJD 60310.0 (2024-Jan-01 TDB)
// da JPL Horizons @sun
// ============================================================

std::vector<AsteroidConfig> get_asteroids() {
    return {
        // ---- NEA / PHA (50 yr) ----
        {
            "Apophis_99942", "NEA_PHA", 50.0, 1.0,
             0.9253730,  0.3156020,  0.0689540,
            -0.0054120,  0.0150830,  0.0008730,
            0.191
        },
        {
            "Icarus_1566", "NEA_PHA", 50.0, 1.0,
             1.0872131858370557, -0.8732191349734646, -0.8953721605377515,
            -0.0005195618936165,  0.0081797267817210,  0.0039326179149579,
            0.827
        },
        {
            "Bennu_101955", "NEA_PHA", 50.0, 1.0,
             0.9741200,  0.6234800,  0.1189600,
            -0.0122330,  0.0130890, -0.0003490,
            0.204
        },

        // ---- Troiani di Giove L4/L5 (1000 yr) ----
        {
            "Achilles_588", "Trojans", 1000.0, 1.0,
             3.9089400, -3.2840900, -0.3971700,
             0.0051380,  0.0056340, -0.0001680,
            0.148
        },
        {
            "Patroclus_617", "Trojans", 1000.0, 1.0,
            -3.2100000, -4.1500000,  0.3200000,
             0.0048200, -0.0043100, -0.0003100,
            0.140
        },
        {
            "Hektor_624", "Trojans", 1000.0, 1.0,
             1.5200000,  4.9800000,  0.1200000,
            -0.0071100,  0.0022300,  0.0001200,
            0.024
        },

        // ---- Famiglie in Risonanza (500 yr) ----
        {
            "Hilda_153", "Resonant", 500.0, 1.0,
            -3.9480600,  0.5219700,  0.1009600,
            -0.0021180, -0.0087920, -0.0003140,
            0.142
        },
        {
            "Hecuba_334", "Resonant", 500.0, 1.0,
             2.1300000, -3.3800000, -0.1500000,
             0.0063200,  0.0037600, -0.0003800,
            0.074
        },
        {
            "Griqua_1362", "Resonant", 500.0, 1.0,
             1.9700000, -2.5800000, -0.3600000,
             0.0072300,  0.0047100, -0.0010200,
            0.372
        },

        // ---- TNO / KBO (10000 yr) ----
        {
            "Pluto_134340", "TNO", 10000.0, 1.0,
             13.7890000, -31.5800000, -0.8940000,
              0.0017720,   0.0005380, -0.0005300,
            0.249
        },
        {
            "Eris_136199", "TNO", 10000.0, 1.0,
             95.6500000,  -3.8600000, -14.7100000,
              0.0000350,   0.0009120,   0.0001540,
            0.441
        },
        {
            "Sedna_90377", "TNO", 10000.0, 1.0,
             83.2800000,  46.2100000,  15.6200000,
             -0.0001050,   0.0002610,   0.0000420,
            0.845
        },
    };
}

// ============================================================
// Conversione unità e creazione stato
// ============================================================

CartesianStateTyped<GCRF> make_state(const AsteroidConfig& ast) {
    auto epoch = time::EpochTDB::from_mjd(60310.0);

    // Impostiamo l'unità tramite i factory Distance e Velocity (Best Practice 3.0)
    auto pos = math::Vector3<GCRF, Distance>::from_si(
        Distance::from_au(ast.x_au).to_m(),
        Distance::from_au(ast.y_au).to_m(),
        Distance::from_au(ast.z_au).to_m()
    );
    auto vel = math::Vector3<GCRF, Velocity>::from_si(
        Velocity::from_au_d(ast.vx_auday).to_ms(),
        Velocity::from_au_d(ast.vy_auday).to_ms(),
        Velocity::from_au_d(ast.vz_auday).to_ms()
    );

    return CartesianStateTyped<GCRF>(
        epoch, pos, vel,
        GravitationalParameter::sun()
    );
}

// ============================================================
// Calcolo energia specifica (Kepleriana: E = -GMS / 2a)
// ============================================================

constexpr double GMS_K2 = 0.01720209895 * 0.01720209895;

// Da struct (AU, AU/day) — vis-viva diretta
double energy_kep(const AsteroidConfig& ast) {
    double r  = std::sqrt(ast.x_au*ast.x_au + ast.y_au*ast.y_au + ast.z_au*ast.z_au);
    double v2 = ast.vx_auday*ast.vx_auday + ast.vy_auday*ast.vy_auday + ast.vz_auday*ast.vz_auday;
    double inv_a = 2.0/r - v2/GMS_K2;
    return -GMS_K2 * inv_a / 2.0;  // = -GMS/(2a)
}

// Da stato propagato — stessa formula, converti in AU/day
template <typename Frame>
double energy_kep(const CartesianStateTyped<Frame>& state) {
    constexpr double M_TO_AU  = 1.0 / 1.495978707e11;
    constexpr double S_TO_DAY = 86400.0;
    auto r_si = state.position.to_eigen_si();
    auto v_si = state.velocity.to_eigen_si();
    double r_au = r_si.norm() * M_TO_AU;
    double vx = v_si.x() * S_TO_DAY * M_TO_AU;
    double vy = v_si.y() * S_TO_DAY * M_TO_AU;
    double vz = v_si.z() * S_TO_DAY * M_TO_AU;
    double v2 = vx*vx + vy*vy + vz*vz;
    double inv_a = 2.0/r_au - v2/GMS_K2;
    return -GMS_K2 * inv_a / 2.0;
}



// ============================================================
// Benchmark singolo integratore per un asteroide
// ============================================================

struct IntegratorSpec {
    std::string name;
    enum Type { AAS, RKF78 } type;
    double epsilon;   // per AAS
    double tol;       // per RKF78
};

void run_benchmark(
    const AsteroidConfig&  ast,
    const IntegratorSpec&  spec,
    std::ofstream&         csv,
    bool                   verbose = false)
{
    auto initial   = make_state(ast);
    double E0      = energy_kep(ast); // Forza l'uso dello stato iniziale esatto in AU
    double dur_day = ast.duration_yr * 365.25;
    double smp_day = ast.sample_yr  * 365.25;
    int    n_steps = static_cast<int>(ast.duration_yr / ast.sample_yr);

    auto ephem = std::make_shared<ephemeris::PlanetaryEphemeris>();
    propagation::PropagatorSettings settings;
    settings.include_planets = false;

    auto t0 = std::chrono::steady_clock::now();

    if (spec.type == IntegratorSpec::AAS) {
        //auto integrator = std::make_shared<propagation::AASIntegrator>(
        //  spec.epsilon, constants::GMS_SI * 1e9);
        auto integrator = std::make_shared<propagation::AASIntegrator>(
            spec.epsilon, constants::GMS);
        propagation::Propagator prop(integrator, ephem, settings);
        CartesianStateTyped<GCRF> current = initial;
        long total_nfe = 0;

        // t=0
        csv << ast.name << "," << spec.name << ","
            << std::fixed << std::setprecision(6) << 0.0 << ","
            << std::scientific << 0.0 << "," << 0.0 << ",0\n";

        for (int k = 1; k <= n_steps; k++) {
            double t_days = k * smp_day;
            auto target = time::EpochTDB::from_mjd(
                initial.epoch.mjd() + t_days);

            current    = prop.propagate_cartesian(current, target);
            total_nfe += prop.statistics().num_function_evals;

            double Ef  = energy_kep(current);
            double dH  = std::abs((Ef - E0) / E0);
            double dHs = prop.statistics().shadow_hamiltonian_drift;

            csv << ast.name << "," << spec.name << ","
                << std::fixed    << std::setprecision(6) << (t_days / 365.25) << ","
                << std::scientific << dH << "," << dHs << "," << total_nfe << "\n";

            if (verbose && (n_steps <= 10 || k % std::max(1, n_steps/5) == 0)) {
                std::cout << "    " << ast.name << " [" << spec.name << "] "
                          << std::fixed << std::setprecision(1) << (t_days/365.25) << "/" << ast.duration_yr
                          << " yr  dH=" << std::scientific << std::setprecision(2) << dH << "    \r" << std::flush;
            }
        }
    }
    else {
        // RKF78 via AstDynEngine
        AstDynConfig cfg;
        cfg.integrator_type   = IntegratorType::RKF78;
        cfg.tolerance         = spec.tol;
        cfg.initial_step_size = 0.1;
        cfg.propagator_settings.include_planets = false;

        auto engine = std::make_unique<AstDynEngine>(cfg);
        auto initial_kep = propagation::cartesian_to_keplerian(initial);
        engine->set_initial_orbit(initial_kep);

        csv << ast.name << "," << spec.name << ","
            << std::fixed << std::setprecision(6) << 0.0 << ","
            << std::scientific << 0.0 << "," << 0.0 << ",0\n";

        for (int k = 1; k <= n_steps; k++) {
            double t_days = k * smp_day;
            auto target = time::EpochTDB::from_mjd(
                initial.epoch.mjd() + t_days);

            auto final_kep  = engine->propagate_to(target);
            auto final_cart = propagation::keplerian_to_cartesian(final_kep);

            double Ef = energy_kep(final_cart);
            double dH = std::abs((Ef - E0) / E0);

            csv << ast.name << "," << spec.name << ","
                << std::fixed    << std::setprecision(6) << (t_days / 365.25) << ","
                << std::scientific << dH << "," << 0.0 << "," << engine->total_force_evaluations() << "\n";

            if (verbose && (n_steps <= 10 || k % std::max(1, n_steps/5) == 0)) {
                std::cout << "    " << ast.name << " [" << spec.name << "] "
                          << std::fixed << std::setprecision(1) << (t_days/365.25) << "/" << ast.duration_yr
                          << " yr  dH=" << std::scientific << std::setprecision(2) << dH << "    \r" << std::flush;
            }
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    double wall_sec = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "\n  [" << spec.name << "] " << ast.name
              << "  " << ast.duration_yr << " yr"
              << "  wall=" << std::fixed << std::setprecision(1)
              << wall_sec << " s\n";
}

// ============================================================
// MAIN
// ============================================================

int main(int argc, char** argv) {
    fs::create_directories("benchmark_extended_results");

    std::string outfile = "benchmark_extended_results/astdyn_energy_extended.csv";
    std::ofstream csv(outfile);
    if (!csv.is_open()) {
        std::cerr << "Errore apertura file: " << outfile << "\n";
        return 1;
    }
    csv << "asteroid,integrator,time_years,delta_H,delta_H_shadow,nfe\n";

    // Integratori da testare
    std::vector<IntegratorSpec> integrators = {
        {"AAS_1e-4", IntegratorSpec::AAS,   1e-4, 0.0},
        {"AAS_1e-6", IntegratorSpec::AAS,   1e-6, 0.0},
        {"AAS_1e-8", IntegratorSpec::AAS,   1e-8, 0.0},
        {"RKF78",    IntegratorSpec::RKF78, 0.0,  1e-12},
    };

    // Filtra categorie da argomento opzionale
    // es: ./benchmark_extended NEA_PHA Trojans
    std::vector<std::string> cat_filter;
    for (int i = 1; i < argc; i++) cat_filter.push_back(argv[i]);

    auto asteroids = get_asteroids();

    for (const auto& spec : integrators) {
        std::cout << "\n============================================================\n";
        std::cout << "Integratore: " << spec.name << "\n";
        std::cout << "============================================================\n";

        for (const auto& ast : asteroids) {
            if (!cat_filter.empty()) {
                bool found = false;
                for (const auto& c : cat_filter)
                    if (ast.category == c) { found = true; break; }
                if (!found) continue;
            }

            std::cout << "\n  >> " << ast.category << " / " << ast.name
                      << "  e=" << ast.eccentricity
                      << "  " << ast.duration_yr << " yr\n";
            try {
                run_benchmark(ast, spec, csv, /*verbose=*/true);
            } catch (const std::exception& ex) {
                std::cerr << "  ERRORE " << ast.name << ": " << ex.what() << "\n";
                // Scrivi riga di errore nel CSV per non rompere il parsing
                csv << ast.name << "," << spec.name << ",0,NaN,NaN,0\n";
            }
        }
    }

    csv.close();
    std::cout << "\n\nDone. Output: " << outfile << "\n";

    return 0;
}
