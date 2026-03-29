/**
 * @file benchmark_integrators.cpp
 * @brief AAS vs SABA4 vs RKF78 — riproduce tutti i CSV benchmark + A6 + A7
 *
 * Stato iniziale: JPL Horizons J2000.0 (JD 2451545.0),
 * frame ICRF, unità AU e AU/day.
 */

// Low-level integrators
#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/propagation/saba4_integrator.hpp"
#include "astdyn/propagation/Integrator.hpp"
// High-level covariance (A7)
#include "astdyn/propagation/CovariancePropagator.hpp"
#include "astdyn/propagation/PropagatorSettings.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
// Core
#include "astdyn/core/Constants.hpp"
#include "astdyn/core/physics_state_au.hpp"
// STL
#include <Eigen/Dense>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <curl/curl.h>

namespace {

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::constants;

// ─── Costanti ─────────────────────────────────────────────────────────────────

static constexpr double MU_AU3D2         = GMS;
static constexpr double T_1YR_D          = DAYS_PER_YEAR;
static constexpr double T_2YR_D          = 2.0 * DAYS_PER_YEAR;
static constexpr double T_30D_D          = 30.0;
static constexpr double DELTA_R_TARGET_AU = 1.67e-11; // 2.5 μm (marker)

static const std::filesystem::path OUTPUT_DIR =
    "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn"
    "/examples/benchmark_results";

// ─── Strutture dati ───────────────────────────────────────────────────────────

struct AsteroidState {
    std::string           name;
    std::string           horizons_cmd;
    double                e;
    std::array<double, 6> xv_au_d; // {x,y,z,vx,vy,vz} AU, AU/day
};

struct BenchmarkConfig {
    std::string integrator_name;
    double      precision_or_step;
    int         n_steps_max;
};

struct BenchmarkResult {
    double delta_E;
    int    n_steps;
    double cpu_ms;
};

struct EnergyRecord {
    std::string integrator_name, asteroid;
    double      precision_or_step, delta_E;
    int         n_steps;
};

// ─── Stato iniziale da CSV ────────────────────────────────────────────────────

static const std::array<AsteroidState, 4> FALLBACK_STATES = {{
    {"Ceres",      "'Ceres'", 0.078368,
     {-2.379327705915647e+00,  5.456711318631370e-01,  7.412254807065680e-01,
      -3.584228273182221e-03, -9.845217307745824e-03, -3.904543826033526e-03}},
    {"Apophis",    "99942",  0.191393,
     {-1.037925696095382e+00, -1.268092611155217e-01, -7.404282940865300e-02,
       4.227374301983514e-03, -1.412107207092049e-02, -5.145341154826345e-03}},
    {"Phaethon",   "3200",   0.890158,
     { 1.603652907608655e+00,  1.218350965043925e+00,  1.181715803389789e+00,
      -6.194883762149516e-04,  4.278066650583167e-03,  1.425200387721508e-03}},
    {"Baruffetti", "309704", 0.050315,
     { 2.828732608224157e+00, -7.765821462552001e-01, -9.299867446588169e-01,
       3.121362656476206e-03,  8.244939831915332e-03,  4.224843180355973e-03}},
}};

static AsteroidState enrich_from_fallback(AsteroidState ast) {
    for (const auto& fb : FALLBACK_STATES)
        if (fb.name == ast.name) { ast.horizons_cmd = fb.horizons_cmd; ast.e = fb.e; }
    return ast;
}

std::vector<AsteroidState> load_initial_states(const std::string& csv_path) {
    std::ifstream file(csv_path);
    if (!file) {
        std::cout << "  [warn] " << csv_path << " not found — using fallback\n";
        return {FALLBACK_STATES.begin(), FALLBACK_STATES.end()};
    }
    std::vector<AsteroidState> states;
    std::string line;
    std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> fields;
        while (std::getline(ss, token, ',')) fields.push_back(token);
        if (fields.size() < 7) continue;
        AsteroidState ast;
        ast.name = fields[0];
        for (int i = 0; i < 6; ++i) ast.xv_au_d[i] = std::stod(fields[i + 1]);
        states.push_back(enrich_from_fallback(std::move(ast)));
    }
    return states.empty() ? std::vector<AsteroidState>(FALLBACK_STATES.begin(), FALLBACK_STATES.end())
                          : states;
}

void write_initial_states(const std::vector<AsteroidState>& states,
                          const std::string& out_dir) {
    std::ofstream out(out_dir + "/initial_states.csv");
    out << "asteroid,x,y,z,vx,vy,vz,epoch_jd\n";
    for (const auto& ast : states) {
        const auto& xv = ast.xv_au_d;
        out << std::setprecision(15) << ast.name
            << "," << xv[0] << "," << xv[1] << "," << xv[2]
            << "," << xv[3] << "," << xv[4] << "," << xv[5]
            << std::fixed << std::setprecision(1) << "," << JD2000 << "\n"
            << std::defaultfloat;
    }
}

// ─── Helpers fisici ───────────────────────────────────────────────────────────

Eigen::VectorXd to_eigen(const std::array<double, 6>& xv_au_d) {
    Eigen::VectorXd y(6);
    for (int i = 0; i < 6; ++i) y[i] = xv_au_d[i];
    return y;
}

double energy_au_d(const Eigen::VectorXd& y, double mu_au3d2) {
    return 0.5 * y.tail<3>().squaredNorm() - mu_au3d2 / y.head<3>().norm();
}

// Laplace-Runge-Lenz: e_vec = v×h/μ − r̂
double eccentricity(const Eigen::VectorXd& y, double mu_au3d2) {
    const Eigen::Vector3d r     = y.head<3>();
    const Eigen::Vector3d v     = y.tail<3>();
    const Eigen::Vector3d e_vec = v.cross(r.cross(v)) / mu_au3d2 - r.normalized();
    return e_vec.norm();
}

DerivativeFunction make_dynamics() {
    return [](double, const Eigen::VectorXd& y) -> Eigen::VectorXd {
        Eigen::VectorXd dy(6);
        const double r3 = std::pow(y.head<3>().norm(), 3);
        dy.head<3>() = y.tail<3>();
        dy.tail<3>() = (-MU_AU3D2 / r3) * y.head<3>();
        return dy;
    };
}

Eigen::Matrix3d keplerian_hessian(const Eigen::Vector3d& q) {
    const double r  = q.norm();
    const double r3 = r * r * r;
    const double r5 = r3 * r * r;
    return (-MU_AU3D2/r3)*Eigen::Matrix3d::Identity()
         + (3.0*MU_AU3D2/r5)*(q*q.transpose());
}

Eigen::Vector3d rotate_icrf_to_eclip(const Eigen::Vector3d& v_icrf) {
    const double eps = OBLIQUITY_J2000;
    const double ce = std::cos(eps);
    const double se = std::sin(eps);
    return { v_icrf.x(),
             v_icrf.y() * ce + v_icrf.z() * se,
            -v_icrf.y() * se + v_icrf.z() * ce };
}

// Shadow Hamiltonian O(h^4) — Yoshida-4 (Preto & Tremaine 1999)
double compute_shadow_H(const Eigen::VectorXd& y, double dt_d) {
    const Eigen::Vector3d q = y.head<3>(), v = y.tail<3>();
    const double r    = q.norm();
    const double H0   = 0.5*v.squaredNorm() - MU_AU3D2/r;
    const double kk   = std::pow(2.0, 1.0/3.0);
    const double w1   = 1.0 / (2.0 - kk);
    const double w0   = 1.0 - 2.0*w1;
    const double s4   = (2.0*std::pow(w1,4) + std::pow(w0,4)) / 24.0;
    const double Phi2 = MU_AU3D2 / (r*r*r);
    const double gP2  = std::pow(MU_AU3D2/(r*r), 2);
    return H0 + std::pow(dt_d, 4) * s4 * (gP2*Phi2 + v.dot(keplerian_hessian(q)*v)*Phi2);
}

// ─── B1: Energy conservation ──────────────────────────────────────────────────

BenchmarkResult run_one_aas(const AsteroidState& ast, double prec) {
    auto f = make_dynamics();
    AASIntegrator integrator(prec, {MU_AU3D2});
    const auto    y0  = to_eigen(ast.xv_au_d);
    const double  e0  = energy_au_d(y0, MU_AU3D2);
    const auto    t0  = std::chrono::high_resolution_clock::now();
    const auto    yf  = integrator.integrate(f, y0, 0.0, T_1YR_D);
    const double  cpu = std::chrono::duration<double, std::milli>(
                            std::chrono::high_resolution_clock::now()-t0).count();
    return {std::abs((energy_au_d(yf,MU_AU3D2)-e0)/e0), integrator.statistics().num_steps, cpu};
}

BenchmarkResult run_one_saba4(const AsteroidState& ast, double step_d) {
    auto f = make_dynamics();
    SABA4Integrator integrator(step_d, 1e-30, step_d, step_d);
    const auto      y0  = to_eigen(ast.xv_au_d);
    const double    e0  = energy_au_d(y0, MU_AU3D2);
    const auto      t0  = std::chrono::high_resolution_clock::now();
    const auto      yf  = integrator.integrate(f, y0, 0.0, T_1YR_D);
    const double    cpu = std::chrono::duration<double, std::milli>(
                              std::chrono::high_resolution_clock::now()-t0).count();
    return {std::abs((energy_au_d(yf,MU_AU3D2)-e0)/e0), integrator.statistics().num_steps, cpu};
}

BenchmarkResult run_one_rkf78(const AsteroidState& ast, double tol) {
    auto f = make_dynamics();
    RKF78Integrator integrator(1.0, tol);
    const auto      y0  = to_eigen(ast.xv_au_d);
    const double    e0  = energy_au_d(y0, MU_AU3D2);
    const auto      t0  = std::chrono::high_resolution_clock::now();
    const auto      yf  = integrator.integrate(f, y0, 0.0, T_1YR_D);
    const double    cpu = std::chrono::duration<double, std::milli>(
                              std::chrono::high_resolution_clock::now()-t0).count();
    return {std::abs((energy_au_d(yf,MU_AU3D2)-e0)/e0), integrator.statistics().num_steps, cpu};
}

void emit_energy_row(std::ofstream& out, std::vector<EnergyRecord>& recs,
                     const std::string& label, const AsteroidState& ast,
                     double prec, const BenchmarkResult& r) {
    out << label << "," << ast.name << "," << prec << ","
        << r.delta_E << "," << r.n_steps << "," << r.cpu_ms << "\n";
    recs.push_back({label, ast.name, prec, r.delta_E, r.n_steps});
}

std::vector<EnergyRecord> run_energy_benchmark(
    const std::vector<AsteroidState>& states, const std::string& out_dir) {
    std::ofstream out(out_dir + "/energy_vs_precision.csv");
    out << "integrator,asteroid,precision,delta_E_over_E,n_steps,cpu_ms\n";
    const std::array aas_precs   = {1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 1e-6};
    const std::array saba4_steps = {2.0, 1.0, 0.5, 0.25, 0.1, 0.05};
    const std::array rkf78_tols  = {1e-6, 1e-8, 1e-10, 1e-12};
    std::vector<EnergyRecord> records;
    for (const auto& ast : states) {
        for (double p : aas_precs)
            emit_energy_row(out, records, "AAS",   ast, p, run_one_aas(ast, p));
        for (double s : saba4_steps)
            emit_energy_row(out, records, "SABA4", ast, s, run_one_saba4(ast, s));
        for (double t : rkf78_tols)
            emit_energy_row(out, records, "RKF78", ast, t, run_one_rkf78(ast, t));
        std::cout << "  [B1] " << ast.name << "\n";
    }
    return records;
}

// ─── B2: Eccentricity divergence ─────────────────────────────────────────────

double find_divergence_eccentricity(const std::vector<Eigen::VectorXd>& nom,
                                    const std::vector<Eigen::VectorXd>& pert,
                                    const std::vector<double>& targets,
                                    double e0) {
    static constexpr double DIVERGENCE_THRESHOLD = 0.1;
    for (size_t i = 0; i < nom.size(); ++i) {
        if (std::abs(eccentricity(pert[i], MU_AU3D2) - e0) > DIVERGENCE_THRESHOLD * e0)
            return targets[i];
    }
    return -1.0;
}

void run_divergence_for_asteroid(std::ofstream& out, const AsteroidState& ast,
                                 const DerivativeFunction& f,
                                 const std::vector<double>& targets) {
    const auto y0 = to_eigen(ast.xv_au_d);
    AASIntegrator   aas_nom(1e-4, {MU_AU3D2});
    SABA4Integrator s4_nom(0.5, 1e-30, 0.5, 0.5);
    const auto nom_aas = aas_nom.integrate_at(f, y0, 0.0, targets);
    const auto nom_s4  = s4_nom.integrate_at(f, y0, 0.0, targets);
    for (int k = 1; k <= 8; ++k) {
        Eigen::VectorXd y0p = y0;
        y0p[0] += static_cast<double>(k) * 1e-12;
        AASIntegrator   aas_p(1e-4, {MU_AU3D2});
        SABA4Integrator s4_p(0.5, 1e-30, 0.5, 0.5);
        const auto pert_aas = aas_p.integrate_at(f, y0p, 0.0, targets);
        const auto pert_s4  = s4_p.integrate_at(f, y0p, 0.0, targets);
        out << "AAS,"  << ast.name << ",1e-4," << k << ","
            << find_divergence_eccentricity(nom_aas, pert_aas, targets, ast.e) << "\n";
        out << "SABA4," << ast.name << ",0.5," << k << ","
            << find_divergence_eccentricity(nom_s4,  pert_s4,  targets, ast.e) << "\n";
    }
    std::cout << "  [B2] " << ast.name << "\n";
}

void run_divergence_benchmark(const std::vector<AsteroidState>& states,
                              const std::string& out_dir) {
    std::ofstream out(out_dir + "/divergence_time.csv");
    out << "integrator,asteroid,precision,perturb_index,t_divergence_days\n";
    auto f = make_dynamics();
    std::vector<double> targets;
    targets.reserve(static_cast<size_t>(T_2YR_D) + 1);
    for (int d = 1; d <= static_cast<int>(T_2YR_D); ++d)
        targets.push_back(static_cast<double>(d));
    for (int ai = 0; ai < 2 && ai < static_cast<int>(states.size()); ++ai)
        run_divergence_for_asteroid(out, states[ai], f, targets);
}

// ─── B3: Step distribution ────────────────────────────────────────────────────

void run_step_distribution_benchmark(const AsteroidState& apophis,
                                     const std::string& out_dir) {
    std::ofstream out(out_dir + "/step_distribution.csv");
    out << "t_days,dt_days,r_au,integrator\n";
    auto f = make_dynamics();
    const auto y0 = to_eigen(apophis.xv_au_d);
    std::vector<double> t_out;
    std::vector<Eigen::VectorXd> y_out;
    AASIntegrator integrator(1e-4, {MU_AU3D2});
    integrator.integrate_steps(f, y0, 0.0, T_1YR_D, t_out, y_out);
    for (size_t i = 1; i < t_out.size(); ++i) {
        if ((i - 1) % 10 != 0) continue;
        const double dt_d = t_out[i] - t_out[i - 1];
        const double r_au = y_out[i].head<3>().norm();
        out << t_out[i] << "," << dt_d << "," << r_au << ",AAS\n";
    }
    std::cout << "  [B3] " << t_out.size() << " steps\n";
}

// ─── B4: Shadow Hamiltonian ───────────────────────────────────────────────────

void write_shadow_rows(std::ofstream& out,
                       const std::vector<double>& t_out,
                       const std::vector<Eigen::VectorXd>& y_out,
                       double H0, double Hs0,
                       double step0_d, const std::string& label, bool is_fixed_step) {
    for (size_t i = 0; i < t_out.size(); i += 10) {
        const double dt_d = is_fixed_step ? step0_d
                          : (i > 0 ? t_out[i] - t_out[i-1] : step0_d);
        const double H    = energy_au_d(y_out[i], MU_AU3D2);
        const double Hs   = compute_shadow_H(y_out[i], dt_d);
        const double dH   = std::abs(H0)  > 0.0 ? std::abs((H  - H0 ) / H0)  : 0.0;
        const double dHs  = std::abs(Hs0) > 0.0 ? std::abs((Hs - Hs0) / Hs0) : 0.0;
        out << t_out[i] << "," << dH << "," << dHs << "," << label << "\n";
    }
}

void run_shadow_benchmark(const AsteroidState& ceres, const std::string& out_dir) {
    std::ofstream out(out_dir + "/shadow_hamiltonian.csv");
    out << "t_days,delta_H_physical,delta_H_shadow,integrator\n";
    auto f = make_dynamics();
    const auto y0 = to_eigen(ceres.xv_au_d);
    {
        std::vector<double> t_out; std::vector<Eigen::VectorXd> y_out;
        AASIntegrator integrator(1e-4, {MU_AU3D2});
        integrator.integrate_steps(f, y0, 0.0, T_2YR_D, t_out, y_out);
        const double dt0_d = t_out.size() > 1 ? t_out[1] - t_out[0] : 1.0;
        write_shadow_rows(out, t_out, y_out, energy_au_d(y0, MU_AU3D2),
                          compute_shadow_H(y0, dt0_d), dt0_d, "AAS", false);
        std::cout << "  [B4] AAS\n";
    }
    {
        const double step_d = 0.5;
        std::vector<double> t_out; std::vector<Eigen::VectorXd> y_out;
        SABA4Integrator integrator(step_d, 1e-30, step_d, step_d);
        integrator.integrate_steps(f, y0, 0.0, T_2YR_D, t_out, y_out);
        write_shadow_rows(out, t_out, y_out, energy_au_d(y0, MU_AU3D2),
                          compute_shadow_H(y0, step_d), step_d, "SABA4", true);
        std::cout << "  [B4] SABA4\n";
    }
}

// ─── B5: STM accuracy ────────────────────────────────────────────────────────

Eigen::VectorXd make_stm_y0(const Eigen::VectorXd& y0) {
    Eigen::VectorXd y(42);
    y.head(6) = y0;
    y.tail(36).setZero();
    for (int i = 0; i < 6; ++i) y[6 + i*6 + i] = 1.0; // row-major identity
    return y;
}

using PerturbedStates = std::array<std::vector<Eigen::VectorXd>, 6>;

PerturbedStates perturbed_trajectories(const DerivativeFunction& f,
                                       const Eigen::VectorXd& y0,
                                       const std::array<double, 6>& deltas,
                                       const std::vector<double>& targets,
                                       bool is_plus) {
    PerturbedStates result;
    for (int k = 0; k < 6; ++k) {
        Eigen::VectorXd yp = y0;
        yp[k] += is_plus ? deltas[k] : -deltas[k];
        AASIntegrator ai(1e-4, {MU_AU3D2});
        result[k] = ai.integrate_at(f, yp, 0.0, targets);
    }
    return result;
}

void run_stm_benchmark(const AsteroidState& apophis, const std::string& out_dir) {
    std::ofstream out(out_dir + "/stm_accuracy.csv");
    out << "t_days,stm_error_frobenius,method\n";
    auto f = make_dynamics();
    const auto y0 = to_eigen(apophis.xv_au_d);
    std::vector<double> targets;
    for (int d = 0; d <= 30; ++d) targets.push_back(static_cast<double>(d));
    AASIntegrator stm_int(1e-4, {MU_AU3D2});
    const auto stm_states = stm_int.integrate_at(f, make_stm_y0(y0), 0.0, targets);
    const std::array<double, 6> deltas = {1e-6, 1e-6, 1e-6, 1e-9, 1e-9, 1e-9};
    const auto plus_traj  = perturbed_trajectories(f, y0, deltas, targets, true);
    const auto minus_traj = perturbed_trajectories(f, y0, deltas, targets, false);
    for (size_t ti = 0; ti < targets.size(); ++ti) {
        const Eigen::MatrixXd phi_aas =
            Eigen::Map<const Eigen::MatrixXd>(stm_states[ti].data() + 6, 6, 6);
        Eigen::MatrixXd phi_num(6, 6);
        for (int k = 0; k < 6; ++k)
            phi_num.col(k) = (plus_traj[k][ti] - minus_traj[k][ti]) / (2.0 * deltas[k]);
        out << targets[ti] << "," << (phi_aas - phi_num).norm() << ",AAS_analytic\n";
    }
    std::cout << "  [B5] done\n";
}

// ─── Curl helpers condivisi (B6, A6) ─────────────────────────────────────────

static size_t curl_write_cb(void* ptr, size_t size, size_t nmemb, void* userdata) {
    auto* data = static_cast<std::string*>(userdata);
    data->append(static_cast<char*>(ptr), size * nmemb);
    return size * nmemb;
}

std::string fetch_horizons_raw(const std::string& cmd, double jd_start, double jd_stop,
                                const std::string& step_size) {
    std::ostringstream url;
    url << "https://ssd.jpl.nasa.gov/api/horizons.api?"
        << "format=text&COMMAND=" << cmd
        << "&OBJ_DATA=NO&MAKE_EPHEM=YES&TABLE_TYPE=VECTORS"
        << "&CENTER=500%4010&REF_PLANE=FRAME"
        << "&START_TIME=JD" << std::fixed << std::setprecision(1) << jd_start
        << "&STOP_TIME=JD"  << jd_stop
        << "&STEP_SIZE=" << step_size << "&VEC_TABLE=2&REF_SYSTEM=ICRF"
        << "&VEC_CORR=NONE&VEC_LABELS=NO&VEC_DELTA_T=NO"
        << "&CSV_FORMAT=YES&OUT_UNITS=AU-D";
    std::string response;
    CURL* curl = curl_easy_init();
    if (!curl) return "";
    curl_easy_setopt(curl, CURLOPT_URL,           url.str().c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curl_write_cb);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA,     &response);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT,       30L);
    const CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    return (res == CURLE_OK) ? response : "";
}

std::vector<Eigen::Vector3d> parse_horizons_positions(const std::string& raw) {
    std::vector<Eigen::Vector3d> positions;
    bool in_ephem = false;
    std::istringstream stream(raw);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.find("$$SOE") != std::string::npos) { in_ephem = true;  continue; }
        if (line.find("$$EOE") != std::string::npos) break;
        if (!in_ephem) continue;
        std::istringstream lss(line);
        std::string token;
        std::vector<std::string> fields;
        while (std::getline(lss, token, ',')) fields.push_back(token);
        if (fields.size() < 5) continue;
        try {
            positions.emplace_back(std::stod(fields[2]),
                                   std::stod(fields[3]),
                                   std::stod(fields[4]));
        } catch (...) {}
    }
    return positions;
}

// ─── B6: Validazione Horizons (step 30d, 1 anno) ─────────────────────────────

template<typename Integ>
void emit_horizons_rows(std::ofstream& out, Integ& integ, const std::string& label,
                        const AsteroidState& ast,
                        const std::vector<Eigen::Vector3d>& hzn_pos,
                        const std::vector<double>& targets) {
    auto f = make_dynamics();
    const auto states = integ.integrate_at(f, to_eigen(ast.xv_au_d), 0.0, targets);
    for (size_t i = 0; i < std::min(states.size(), hzn_pos.size()); ++i) {
        const Eigen::VectorXd& sv = states[i];
        out << label << "," << ast.name << "," << targets[i] << ","
            << (sv.head<3>() - hzn_pos[i]).norm() * AU << "\n";
    }
}

void run_horizons_validation_benchmark(const std::vector<AsteroidState>& states,
                                       const std::string& out_dir) {
    std::ofstream out(out_dir + "/horizons_validation.csv");
    out << "integrator,asteroid,t_days,dr_km\n";
    std::vector<double> targets;
    for (double t = T_30D_D; t <= T_1YR_D; t += T_30D_D)
        targets.push_back(t);
    for (const auto& ast : states) {
        const auto raw = fetch_horizons_raw(ast.horizons_cmd, JD2000, JD2000+T_1YR_D, "30d");
        if (raw.empty()) { std::cout << "  [B6] SKIP " << ast.name << "\n"; continue; }
        const auto hzn_all = parse_horizons_positions(raw);
        if (hzn_all.empty()) { std::cout << "  [B6] PARSE FAIL " << ast.name << "\n"; continue; }
        // Skip t=0 element so hzn_pos[i] matches targets[i]
        const std::vector<Eigen::Vector3d> hzn_pos(hzn_all.begin() + 1, hzn_all.end());
        const double r_hzn_au   = hzn_all[0].norm();
        const double r_state_au = to_eigen(ast.xv_au_d).head<3>().norm();
        std::cout << "    t=0 check: Horizons r=" << r_hzn_au
                  << " AU, state r=" << r_state_au
                  << " AU, Δ=" << std::abs(r_hzn_au - r_state_au) << "\n";
        AASIntegrator   aas(1e-6, {MU_AU3D2});
        SABA4Integrator saba4(0.1, 1e-30, 0.1, 0.1);
        RKF78Integrator rkf78(1.0, 1e-12);
        emit_horizons_rows(out, aas,   "AAS",   ast, hzn_pos, targets);
        emit_horizons_rows(out, saba4, "SABA4", ast, hzn_pos, targets);
        emit_horizons_rows(out, rkf78, "RKF78", ast, hzn_pos, targets);
        std::cout << "  [B6] " << ast.name << " (" << hzn_pos.size() << " epoche)\n";
    }
}

// ─── B7: Efficienza ───────────────────────────────────────────────────────────

static int n_evals_per_step(const std::string& integ_name) {
    if (integ_name == "AAS")   return 7;
    if (integ_name == "SABA4") return 4;
    return 13;
}

void run_efficiency_benchmark(const std::vector<EnergyRecord>& records,
                               const std::string& out_dir) {
    std::ofstream out(out_dir + "/efficiency.csv");
    out << "integrator,asteroid,precision,n_func_evals,delta_E_over_E\n";
    for (const auto& rec : records) {
        if (rec.asteroid != "Phaethon" && rec.asteroid != "Apophis") continue;
        const long long n_evals = static_cast<long long>(rec.n_steps)
                                * n_evals_per_step(rec.integrator_name);
        out << rec.integrator_name << "," << rec.asteroid << ","
            << rec.precision_or_step << "," << n_evals << "," << rec.delta_E << "\n";
    }
    std::cout << "  [B7] done\n";
}

// ─── A6: Residui Horizons short-term (step 1d, 30 giorni) ────────────────────

template<typename Integ>
void emit_horizons_short_rows(std::ofstream& out, Integ& integ, const std::string& label,
                               const AsteroidState& ast,
                               const std::vector<Eigen::Vector3d>& hzn_pos_ecl,
                               const std::vector<double>& targets) {
    auto f = make_dynamics();
    Eigen::VectorXd y0_icrf = to_eigen(ast.xv_au_d);
    Eigen::VectorXd y0_ecl(6);
    y0_ecl.head<3>() = rotate_icrf_to_eclip(y0_icrf.head<3>());
    y0_ecl.tail<3>() = rotate_icrf_to_eclip(y0_icrf.tail<3>());

    const auto states = integ.integrate_at(f, y0_ecl, 0.0, targets);
    for (size_t i = 0; i < std::min(states.size(), hzn_pos_ecl.size()); ++i) {
        const Eigen::VectorXd& sv = states[i];
        const double delta_r_au = (sv.head<3>() - hzn_pos_ecl[i]).norm();
        out << label << "," << ast.name << "," << targets[i] << "," << delta_r_au << "\n";
    }
    if (!states.empty() && !hzn_pos_ecl.empty()) {
        const Eigen::VectorXd& sv_back = states.back();
        const double dr30_au = (sv_back.head<3>() - hzn_pos_ecl.back()).norm();
        std::cout << "    " << std::setw(5) << label << "  t=30d: Δr=" << std::scientific << std::setprecision(4) << dr30_au << " AU\n" << std::defaultfloat;
    }
}

void run_horizons_short_benchmark(const std::vector<AsteroidState>& states,
                                   const std::string& out_dir) {
    std::ofstream out(out_dir + "/horizons_short.csv");
    out << "integrator,asteroid,t_days,delta_r_au\n";
    std::vector<double> targets;
    for (int d = 1; d <= static_cast<int>(T_30D_D); ++d)
        targets.push_back(static_cast<double>(d));
    for (const auto& ast : states) {
        if (ast.name == "Apophis") {
             std::cout << "    [1a] Propagator t_start = 0.0 (JD " << std::fixed << std::setprecision(1) << JD2000 << ")\n";
             std::cout << "    [1a] Horizons START_TIME = JD " << JD2000 << "\n";
        }
        const auto raw = fetch_horizons_raw(ast.horizons_cmd, JD2000, JD2000+T_30D_D, "1d");
        if (raw.empty()) { std::cout << "  [A6] SKIP " << ast.name << "\n"; continue; }
        
        if (ast.name == "Apophis") {
             std::cout << "    [1c] Horizons RAW extract (first lines):\n";
             size_t soe = raw.find("$$SOE");
             if (soe != std::string::npos) std::cout << raw.substr(soe, 200) << "\n";
        }

        const auto hzn_all_icrf = parse_horizons_positions(raw);
        if (hzn_all_icrf.empty()) { std::cout << "  [A6] PARSE FAIL " << ast.name << "\n"; continue; }
        
        std::vector<Eigen::Vector3d> hzn_all_ecl;
        for (const auto& p : hzn_all_icrf) hzn_all_ecl.push_back(rotate_icrf_to_eclip(p));

        // Skip t=0 element so hzn_pos[i] matches targets[i]
        const std::vector<Eigen::Vector3d> hzn_pos_ecl(hzn_all_ecl.begin() + 1, hzn_all_ecl.end());
        const double r_hzn_au   = hzn_all_ecl[0].norm();
        const double r_state_au = to_eigen(ast.xv_au_d).head<3>().norm();
        std::cout << "    t=0 check: Horizons (Ecl) r=" << std::setprecision(12) << r_hzn_au
                  << " AU, state r=" << r_state_au
                  << " AU, Δ=" << std::scientific << std::abs(r_hzn_au - r_state_au) << " AU\n" << std::defaultfloat;
        
        AASIntegrator   aas(1e-4, {MU_AU3D2});
        SABA4Integrator saba4(0.5, 1e-30, 0.5, 0.5);
        RKF78Integrator rkf78(1.0, 1e-10);
        std::cout << "  [A6] " << ast.name << "\n";
        emit_horizons_short_rows(out, aas,   "AAS",   ast, hzn_pos_ecl, targets);
        emit_horizons_short_rows(out, saba4, "SABA4", ast, hzn_pos_ecl, targets);
        emit_horizons_short_rows(out, rkf78, "RKF78", ast, hzn_pos_ecl, targets);
    }
}

// ─── A7: Propagazione dell'incertezza ────────────────────────────────────────

static std::shared_ptr<Propagator> make_two_body_propagator_aas() {
    PropagatorSettings settings;
    settings.include_planets    = false;
    settings.include_moon       = false;
    settings.include_asteroids  = false;
    settings.include_relativity = false;
    settings.include_earth_j2   = false;
    settings.include_sun_j2     = false;
    return std::make_shared<Propagator>(
        std::make_shared<AASIntegrator>(1e-4, std::vector<double>{MU_AU3D2}),
        std::make_shared<ephemeris::PlanetaryEphemeris>(),
        settings
    );
}

static physics::CartesianStateTyped<core::ECLIPJ2000> to_si_state(
    const AsteroidState& ast, double t0_mjd) {
    const double au_m   = AU * 1000.0;
    const double aud_ms = au_m / SECONDS_PER_DAY;
    const auto& xv = ast.xv_au_d;
    return physics::CartesianStateTyped<core::ECLIPJ2000>::from_si(
        time::EpochTDB::from_mjd(t0_mjd),
        xv[0]*au_m,   xv[1]*au_m,   xv[2]*au_m,
        xv[3]*aud_ms, xv[4]*aud_ms, xv[5]*aud_ms,
        GM_SUN * 1e9
    );
}

static physics::CovarianceMatrixAU<core::ECLIPJ2000> make_diagonal_P0_au() {
    static constexpr double SIGMA_POS_AU  = 1e-8;
    static constexpr double SIGMA_VEL_AUD = 1e-11;
    Eigen::Matrix<double,6,6> p0 = Eigen::Matrix<double,6,6>::Zero();
    p0.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * (SIGMA_POS_AU  * SIGMA_POS_AU);
    p0.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * (SIGMA_VEL_AUD * SIGMA_VEL_AUD);
    return physics::CovarianceMatrixAU<core::ECLIPJ2000>(p0);
}

// Estrae sigma_R, sigma_AT, sigma_CT dalla covarianza 3×3 in AU² nel frame RTN
static std::tuple<double,double,double> rtn_sigmas_au(
    const Eigen::Matrix3d& P_pos_au2,
    const Eigen::Vector3d& r_au,
    const Eigen::Vector3d& v_aud) {
    const Eigen::Vector3d R_hat = r_au.normalized();
    const Eigen::Vector3d N_hat = r_au.cross(v_aud).normalized();
    const Eigen::Vector3d T_hat = N_hat.cross(R_hat);
    Eigen::Matrix3d rot; rot.row(0) = R_hat; rot.row(1) = T_hat; rot.row(2) = N_hat;
    const Eigen::Matrix3d P_rtn = rot * P_pos_au2 * rot.transpose();
    return {std::sqrt(std::max(0.0, P_rtn(0,0))),
            std::sqrt(std::max(0.0, P_rtn(1,1))),
            std::sqrt(std::max(0.0, P_rtn(2,2)))};
}

static std::tuple<double,double,double> cov_prop_rtn_one_day(
    CovariancePropagator& cp,
    const physics::CartesianStateTyped<core::ECLIPJ2000>& state0,
    const physics::CovarianceMatrixAU<core::ECLIPJ2000>& P0,
    const physics::CovarianceMatrixAU<core::ECLIPJ2000>& Q0,
    int day) {
    const auto target = time::EpochTDB::from_mjd(MJD2000 + day);
    const auto result = cp.propagate_with_covariance(state0, P0, target, Q0);
    const auto mu_au  = physics::CartesianStateAU<core::ECLIPJ2000>::from_si(result.mean_state);
    const Eigen::Vector3d r_au(mu_au.x,  mu_au.y,  mu_au.z);
    const Eigen::Vector3d v_aud(mu_au.vx, mu_au.vy, mu_au.vz);
    return rtn_sigmas_au(result.covariance.matrix().block<3,3>(0,0), r_au, v_aud);
}

void run_cov_prop_method(std::ofstream& out, const AsteroidState& apophis) {
    const auto state0 = to_si_state(apophis, MJD2000);
    const auto P0     = make_diagonal_P0_au();
    const auto Q_zero = physics::CovarianceMatrixAU<core::ECLIPJ2000>{};
    auto cp           = CovariancePropagator::make(make_two_body_propagator_aas());
    for (int d = 0; d <= static_cast<int>(T_30D_D); ++d) {
        const auto [sig_R, sig_AT, sig_CT] = cov_prop_rtn_one_day(cp, state0, P0, Q_zero, d);
        out << d << "," << sig_AT << "," << sig_CT << "," << sig_R << ",CovProp\n";
    }
    std::cout << "  [A7] CovProp done\n";
}

// Genera N_SAMPLES stati perturbati da N(y0, P0) usando Cholesky
static std::vector<Eigen::VectorXd> generate_mc_samples_au(
    const Eigen::VectorXd& y0_au, int n_samples) {
    static constexpr double SIGMA_POS_AU  = 1e-8;
    static constexpr double SIGMA_VEL_AUD = 1e-11;
    Eigen::Matrix<double,6,6> P0 = Eigen::Matrix<double,6,6>::Zero();
    P0.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * (SIGMA_POS_AU  * SIGMA_POS_AU);
    P0.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * (SIGMA_VEL_AUD * SIGMA_VEL_AUD);
    const Eigen::Matrix<double,6,6> L = P0.llt().matrixL();
    std::mt19937_64 rng(42);
    std::normal_distribution<double> nd;
    std::vector<Eigen::VectorXd> samples(n_samples);
    for (int s = 0; s < n_samples; ++s) {
        Eigen::Matrix<double,6,1> z;
        for (int i = 0; i < 6; ++i) z[i] = nd(rng);
        samples[s] = y0_au + L * z;
    }
    return samples;
}

// Sigma RTN per un giorno dall'ensemble di posizioni
static std::tuple<double,double,double> mc_rtn_sigma_one_day_au(
    const std::vector<std::vector<Eigen::VectorXd>>& traj, size_t day_idx) {
    const int n = static_cast<int>(traj.size());
    Eigen::Vector3d r_mean = Eigen::Vector3d::Zero(), v_mean = Eigen::Vector3d::Zero();
    for (const auto& t : traj) { r_mean += t[day_idx].head<3>(); v_mean += t[day_idx].tail<3>(); }
    r_mean /= n; v_mean /= n;
    const Eigen::Vector3d R = r_mean.normalized();
    const Eigen::Vector3d N = r_mean.cross(v_mean).normalized();
    const Eigen::Vector3d T = N.cross(R);
    double var_R = 0, var_AT = 0, var_CT = 0;
    for (const auto& t : traj) {
        const Eigen::Vector3d d = t[day_idx].head<3>() - r_mean;
        var_R += d.dot(R)*d.dot(R); var_AT += d.dot(T)*d.dot(T); var_CT += d.dot(N)*d.dot(N);
    }
    return {std::sqrt(var_R/n), std::sqrt(var_AT/n), std::sqrt(var_CT/n)};
}

void run_monte_carlo_method(std::ofstream& out, const AsteroidState& apophis) {
    static constexpr int N_SAMPLES = 500;
    const auto y0 = to_eigen(apophis.xv_au_d);
    const auto samples = generate_mc_samples_au(y0, N_SAMPLES);
    std::vector<double> targets;
    for (int d = 0; d <= static_cast<int>(T_30D_D); ++d) targets.push_back(d);
    auto f = make_dynamics();
    std::vector<std::vector<Eigen::VectorXd>> traj(N_SAMPLES);
    for (int s = 0; s < N_SAMPLES; ++s) {
        AASIntegrator ai(1e-4, {MU_AU3D2});
        traj[s] = ai.integrate_at(f, samples[s], 0.0, targets);
    }
    for (size_t di = 0; di < targets.size(); ++di) {
        const auto [sig_R, sig_AT, sig_CT] = mc_rtn_sigma_one_day_au(traj, di);
        out << targets[di] << "," << sig_AT << "," << sig_CT << "," << sig_R << ",MonteCarlo\n";
    }
    std::cout << "  [A7] MonteCarlo done\n";
}

void run_uncertainty_benchmark(const AsteroidState& apophis, const std::string& out_dir) {
    std::ofstream out(out_dir + "/uncertainty.csv");
    out << "t_days,sigma_AT_au,sigma_CT_au,sigma_R_au,method\n";
    run_cov_prop_method(out, apophis);
    run_monte_carlo_method(out, apophis);
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────────

int main() {
    const std::string out_dir = OUTPUT_DIR.string();
    std::filesystem::create_directories(out_dir);
    std::cout << "Output: " << out_dir << "\n";

    const auto states = load_initial_states(out_dir + "/initial_states.csv");
    write_initial_states(states, out_dir);
    std::cout << "[A2] initial_states.csv\n";

    std::cout << "[B1] Energy conservation...\n";
    const auto energy_records = run_energy_benchmark(states, out_dir);

    std::cout << "[B2] Eccentricity divergence...\n";
    run_divergence_benchmark(states, out_dir);

    std::cout << "[B3] Step distribution...\n";
    run_step_distribution_benchmark(states[1], out_dir); // Apophis

    std::cout << "[B4] Shadow Hamiltonian...\n";
    run_shadow_benchmark(states[0], out_dir); // Ceres

    std::cout << "[B5] STM accuracy...\n";
    run_stm_benchmark(states[1], out_dir); // Apophis

    std::cout << "[B6] Horizons validation...\n";
    run_horizons_validation_benchmark(states, out_dir);

    std::cout << "[B7] Efficiency...\n";
    run_efficiency_benchmark(energy_records, out_dir);

    std::cout << "[A6] Short-term Horizons residuals...\n";
    run_horizons_short_benchmark(states, out_dir);

    std::cout << "[A7] Uncertainty propagation...\n";
    run_uncertainty_benchmark(states[1], out_dir); // Apophis

    std::cout << "All benchmarks complete.\n";
    return 0;
}
