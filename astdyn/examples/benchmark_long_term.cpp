/**
 * @file benchmark_long_term.cpp
 * @brief AAS vs SABA4 vs RKF78 — long-term benchmarks, 12 asteroids, 6 tests.
 *
 * State vectors: JPL Horizons, JD 2460310.5 (MJD 60310.0), ICRF, AU-D.
 * Output: examples/benchmark_results/long_term_tests.csv
 *         examples/benchmark_results/lyapunov_series.csv
 */
#include "astdyn/propagation/AASIntegrator.hpp"
#include "astdyn/propagation/saba4_integrator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {
using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::constants;

// ─── Physical constants ───────────────────────────────────────────────────────
static constexpr double GMS_AU3D2          = GMS;
static double A_JUPITER_AU       = 5.2044; // Initialized from Horizons
static double GM_JUPITER_AU3D2   = GMS_AU3D2 * JUPITER_MASS_RATIO;
static constexpr double D0_LYAP_AU         = 1e-6;
static constexpr double TAU_LYAP_D         = 10.0 * DAYS_PER_YEAR;
static double THETA0_JUPITER_RAD = 0.0; // Calculated at runtime
static double N_JUPITER_DAV      = 0.0; // Calculated at runtime

Eigen::Vector3d rotate_icrf_to_eclip(const Eigen::Vector3d& v_icrf) {
    const double ce = std::cos(OBLIQUITY_J2000);
    const double se = std::sin(OBLIQUITY_J2000);
    return { v_icrf.x(),
             v_icrf.y() * ce + v_icrf.z() * se,
            -v_icrf.y() * se + v_icrf.z() * ce };
}

static const std::filesystem::path OUT_DIR =
    "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn"
    "/examples/benchmark_results";

// ─── Data structures ─────────────────────────────────────────────────────────
struct AsteroidSpec {
    std::string          name;
    std::string          category;
    double               T_yr;
    int                  lyap_N;
    std::array<double,6> xv_au_d;
};

struct IntegratorSpec {
    std::string name;
    double      precision;
};

struct CsvRow {
    std::string asteroid, integrator, test;
    double      T_yr, value;
    std::string unit, note;
};

struct LyapRow {
    std::string asteroid, integrator;
    int         interval_index;
    double      lambda_i_yr;
};

struct RunResult {
    std::vector<double>      times_yr;
    std::vector<double>      dH_vals;
    std::vector<double>      dL_vals;
    Eigen::VectorXd          y_final;
};

// ─── Asteroid catalog at JD 2460310.5 in ICRF ────────────────────────────────
static const std::vector<AsteroidSpec> ASTEROIDS = {
    // NEA/PHA — 50 yr, lyap_N=5
    {"Apophis",  "NEA",       50.0,    5,
     {-7.9651378683175167e-01, -5.1147558537755800e-01, -2.1018833794803629e-01,
       1.2113921469451630e-02, -1.1302603935421080e-02, -3.8947578481756361e-03}},
    {"Icarus",   "NEA",       50.0,    5,
     { 1.0872131858370560e+00, -8.7321913497346448e-01, -8.9537216053775170e-01,
      -5.1956189361648244e-04,  8.1797267817210334e-03,  3.9326179149579013e-03}},
    {"Phaethon", "NEA",       50.0,    5,
     { 1.2389289867467810e+00,  3.7146229370786499e-01,  6.8935792068876300e-01,
       7.5285852629341302e-03,  8.4084561586449322e-03,  6.6062188359328781e-03}},
    // Trojans — 1000 yr, lyap_N=100
    {"Achilles",  "Trojan", 1000.0,  100,
     {-2.0872570925811651e+00,  3.5453971326231488e+00,  1.8253263113571241e+00,
      -7.8079076461476027e-03, -2.5687698399377260e-03, -2.6720907075462822e-03}},
    {"Patroclus", "Trojan", 1000.0,  100,
     { 3.8001930000000000e+00, -1.0375360000000000e+00, -2.2006130000000000e+00,
       3.7519500000000000e-03,  6.5799700000000000e-03,  4.1182390000000000e-03}},
    {"Hektor",    "Trojan", 1000.0,  100,
     {-2.2101480000000000e+00,  3.6819150000000000e+00,  2.9112470000000000e+00,
      -6.7950240000000000e-03, -2.1420500000000000e-03, -2.6820210000000000e-03}},
    // Resonant — 500 yr, lyap_N=50
    {"Hilda",    "Resonant",  500.0,   50,
     { 1.9129084628411279e+00, -2.8414207759505330e+00, -7.2857399649781429e-01,
       8.4506808527998519e-03,  4.2218122471573459e-03,  2.2938716619912980e-03}},
    {"Thule",    "Resonant",  500.0,   50,
     {-2.0325425212207922e+00,  3.2027991108574470e+00,  1.5234228841835140e+00,
      -7.5880309277105594e-03, -3.9622180680430634e-03, -1.4555970738063470e-03}},
    {"Griqua",   "Resonant",  500.0,   50,
     { 8.3656836631057996e-01,  1.9811874592418870e+00,  4.0817601809689928e-02,
      -1.0364934774173450e-02,  6.9660257510013339e-03,  5.2069220065296294e-03}},
    // TNO — 10000 yr, lyap_N=1000
    {"Pluto",    "TNO",    10000.0, 1000,
     { 1.7209973448800749e+01, -2.7147271035327631e+01, -1.3654654732267741e+01,
       2.7979116842486770e-03,  1.1459876132036480e-03, -4.9067466069460349e-04}},
    {"Eris",     "TNO",    10000.0, 1000,
     { 8.5555681957690126e+01,  4.2859102087585242e+01, -1.2370402920363570e+00,
      -4.5487245258154862e-04,  4.3114848079874659e-04,  1.2037461565282980e-03}},
    {"Sedna",    "TNO",    10000.0, 1000,
     { 4.1210526013738388e+01,  7.1631143957832890e+01,  1.2312543261197961e+01,
      -2.4698263380263448e-03,  5.0478923254402858e-04,  4.3174001133092959e-04}},
};

// ─── Helpers ─────────────────────────────────────────────────────────────────
static Eigen::VectorXd to_eigen(const std::array<double,6>& xv_au_d) {
    Eigen::VectorXd y(6);
    for (int i = 0; i < 6; ++i) y[i] = xv_au_d[i];
    return y;
}

static double energy_au2d2(const Eigen::VectorXd& y) {
    return 0.5 * y.tail<3>().squaredNorm() - GMS_AU3D2 / y.head<3>().norm();
}

static double ang_mom_au2d(const Eigen::VectorXd& y) {
    return y.head<3>().cross(y.tail<3>()).norm();
}

static DerivativeFunction make_dynamics() {
    return [](double /*t*/, const Eigen::VectorXd& y) {
        Eigen::VectorXd dy(6);
        const double r3_au3 = std::pow(y.head<3>().norm(), 3);
        dy.head<3>() = y.tail<3>();
        dy.tail<3>() = -(GMS_AU3D2 / r3_au3) * y.head<3>();
        return dy;
    };
}

static DerivativeFunction make_dynamics_crtbp() {
    return [](double t_d, const Eigen::VectorXd& y) {
        Eigen::VectorXd dy(6);
        const double theta   = N_JUPITER_DAV * t_d + THETA0_JUPITER_RAD;
        const Eigen::Vector3d r_J{A_JUPITER_AU * std::cos(theta),
                                   A_JUPITER_AU * std::sin(theta), 0.0};
        const Eigen::Vector3d r_au  = y.head<3>();
        const double r_sun3 = std::pow(r_au.norm(), 3);
        const Eigen::Vector3d dr_J  = r_au - r_J;
        const double rJ3    = std::pow(dr_J.norm(), 3);
        dy.head<3>() = y.tail<3>();
        dy.tail<3>() = -(GMS_AU3D2 / r_sun3) * r_au
                      - (GM_JUPITER_AU3D2 / rJ3) * dr_J;
        return dy;
    };
}

// ─── Integrator factory ───────────────────────────────────────────────────────
static double saba4_step_d(const std::string& category) {
    if (category == "TNO")    return 5.0;
    if (category == "Trojan") return 1.0;
    return 0.5;
}

static double rkf78_maxstep_d(const std::string& category) {
    if (category == "TNO")    return 1000.0;
    if (category == "Trojan") return 100.0;
    return 10.0;
}

static std::unique_ptr<Integrator> make_integrator(const IntegratorSpec& spec,
                                                    const std::string& category) {
    if (spec.name == "AAS")
        return std::make_unique<AASIntegrator>(spec.precision, std::vector<double>{GMS_AU3D2});
    if (spec.name == "SABA4") {
        const double h_d = saba4_step_d(category);
        return std::make_unique<SABA4Integrator>(h_d, 1e-30, h_d, h_d);
    }
    const double max_h_d = rkf78_maxstep_d(category);
    return std::make_unique<RKF78Integrator>(1.0, spec.precision, 1e-6, max_h_d);
}

// ─── Jacobi constant (restricted 3-body: Sun + Jupiter) ──────────────────────
static Eigen::Vector3d rotating_pos(const Eigen::VectorXd& y, double ct, double st) {
    return { y[0]*ct + y[1]*st, -y[0]*st + y[1]*ct, y[2] };
}

static Eigen::Vector3d rotating_vel(const Eigen::VectorXd& y, double ct, double st,
                                     double n_J_dav, const Eigen::Vector3d& rr_au) {
    return { y[3]*ct + y[4]*st + n_J_dav*rr_au[1],
            -y[3]*st + y[4]*ct - n_J_dav*rr_au[0],
             y[5] };
}

static double jacobi_constant(const Eigen::VectorXd& y, double t_d) {
    const double n_J_dav = std::sqrt(GMS_AU3D2 / std::pow(A_JUPITER_AU, 3));
    const double theta   = n_J_dav * t_d + THETA0_JUPITER_RAD;
    const double ct = std::cos(theta), st = std::sin(theta);
    const auto rr_au = rotating_pos(y, ct, st);
    const auto vr_aud = rotating_vel(y, ct, st, n_J_dav, rr_au);
    const double r1_au = rr_au.norm();
    const double r2_au = (rr_au - Eigen::Vector3d{A_JUPITER_AU, 0.0, 0.0}).norm();
    return n_J_dav*n_J_dav*(rr_au[0]*rr_au[0]+rr_au[1]*rr_au[1])
         + 2.0*(GMS_AU3D2/r1_au + GM_JUPITER_AU3D2/r2_au)
         - vr_aud.squaredNorm();
}

// ─── Trajectory integration (year-by-year to stay within step limits) ─────────
static void update_run_result(RunResult& res, const Eigen::VectorXd& y,
                               double H0, double L0_au2d, double yr) {
    res.times_yr.push_back(yr);
    res.dH_vals.push_back(std::abs(energy_au2d2(y) - H0) / std::abs(H0));
    res.dL_vals.push_back(std::abs(ang_mom_au2d(y) - L0_au2d) / L0_au2d);
}

static RunResult run_trajectory(Integrator& integ, const AsteroidSpec& ast,
                                const DerivativeFunction& f) {
    const Eigen::VectorXd y0 = to_eigen(ast.xv_au_d);
    const double H0      = energy_au2d2(y0);
    const double L0_au2d = ang_mom_au2d(y0);
    RunResult res;
    Eigen::VectorXd y = y0;
    for (int yr = 1; yr <= static_cast<int>(ast.T_yr); ++yr) {
        y = integ.integrate(f, y, 0.0, DAYS_PER_YEAR);
        update_run_result(res, y, H0, L0_au2d, static_cast<double>(yr));
    }
    res.y_final = y;
    return res;
}

// ─── Test 1 — Energy conservation ────────────────────────────────────────────
static void run_energy_test(const RunResult& run, const AsteroidSpec& ast,
                             const IntegratorSpec& spec, std::vector<CsvRow>& rows) {
    rows.push_back({ast.name, spec.name, "energy",
                    ast.T_yr, run.dH_vals.back(), "dimensionless", ""});
}

// ─── Test 3 — Secular drift (uses energy series from Test 1) ─────────────────
static void linear_fit_yr(const std::vector<double>& xs, const std::vector<double>& ys,
                           double& slope_out, double& r2_out) {
    const double n_d = static_cast<double>(xs.size());
    double s_x = 0.0, s_y = 0.0, s_xx = 0.0, s_xy = 0.0;
    for (size_t i = 0; i < xs.size(); ++i) {
        s_x += xs[i]; s_y += ys[i];
        s_xx += xs[i]*xs[i]; s_xy += xs[i]*ys[i];
    }
    slope_out = (n_d*s_xy - s_x*s_y) / (n_d*s_xx - s_x*s_x);
    const double b_val  = (s_y - slope_out*s_x) / n_d;
    const double mean_y = s_y / n_d;
    double ss_res = 0.0, ss_tot = 0.0;
    for (size_t i = 0; i < xs.size(); ++i) {
        const double yhat = slope_out*xs[i] + b_val;
        ss_res += (ys[i]-yhat)*(ys[i]-yhat);
        ss_tot += (ys[i]-mean_y)*(ys[i]-mean_y);
    }
    r2_out = (ss_tot > 0.0) ? 1.0 - ss_res/ss_tot : 1.0;
}

static void run_secular_drift_test(const RunResult& run, const AsteroidSpec& ast,
                                    const IntegratorSpec& spec, std::vector<CsvRow>& rows) {
    double slope_1_yr = 0.0, r2 = 0.0;
    linear_fit_yr(run.times_yr, run.dH_vals, slope_1_yr, r2);
    rows.push_back({ast.name, spec.name, "secular_slope", ast.T_yr, slope_1_yr, "1/yr", ""});
    rows.push_back({ast.name, spec.name, "secular_R2",    ast.T_yr, r2, "dimensionless", ""});
}

// ─── Test 4 — Angular momentum ────────────────────────────────────────────────
static void run_angular_momentum_test(const RunResult& run, const AsteroidSpec& ast,
                                       const IntegratorSpec& spec, std::vector<CsvRow>& rows) {
    rows.push_back({ast.name, spec.name, "angular_momentum",
                    ast.T_yr, run.dL_vals.back(), "dimensionless", ""});
}

// ─── Test 6 — Jacobi constant via CRTBP (Trojans only) ───────────────────────
static Eigen::VectorXd propagate_crtbp_years(Integrator& integ, int n_years,
                                               const Eigen::VectorXd& y0) {
    const auto f_crtbp = make_dynamics_crtbp();
    Eigen::VectorXd y = y0;
    for (int yr = 1; yr <= n_years; ++yr) {
        const double t0_d = (yr - 1) * DAYS_PER_YEAR;
        y = integ.integrate(f_crtbp, y, t0_d, t0_d + DAYS_PER_YEAR);
    }
    return y;
}

static void run_jacobi_crtbp_test(const AsteroidSpec& ast, const IntegratorSpec& spec,
                                    std::vector<CsvRow>& rows) {
    if (ast.category != "Trojan") return;
    auto integ = make_integrator(spec, ast.category);
    const Eigen::VectorXd y0 = to_eigen(ast.xv_au_d);
    const double CJ0 = jacobi_constant(y0, 0.0);
    std::cout << "    CJ(0)=" << CJ0 << " AU²/d²";
    const int n_yr = static_cast<int>(ast.T_yr);
    const Eigen::VectorXd yf = propagate_crtbp_years(*integ, n_yr, y0);
    const double t_final_d = ast.T_yr * DAYS_PER_YEAR;
    const double dCJ = std::abs(jacobi_constant(yf, t_final_d) - CJ0) / std::abs(CJ0);
    rows.push_back({ast.name, spec.name, "jacobi", ast.T_yr, dCJ, "dimensionless", ""});
}

// ─── Test 2 — Reversibility ───────────────────────────────────────────────────
static Eigen::VectorXd propagate_years(Integrator& integ, const DerivativeFunction& f,
                                        const Eigen::VectorXd& y0, double T_yr) {
    Eigen::VectorXd y = y0;
    for (int yr = 0; yr < static_cast<int>(T_yr); ++yr)
        y = integ.integrate(f, y, 0.0, DAYS_PER_YEAR);
    return y;
}

static std::pair<double,double> reversibility_errors(Integrator& integ,
                                                       const AsteroidSpec& ast,
                                                       const DerivativeFunction& f) {
    const Eigen::VectorXd y0   = to_eigen(ast.xv_au_d);
    const Eigen::VectorXd yT   = propagate_years(integ, f, y0, ast.T_yr);
    Eigen::VectorXd yrev = yT; yrev.tail<3>() *= -1.0;
    const Eigen::VectorXd yf   = propagate_years(integ, f, yrev, ast.T_yr);
    Eigen::VectorXd yfinal = yf; yfinal.tail<3>() *= -1.0;
    return {(yfinal.head<3>() - y0.head<3>()).norm() / y0.head<3>().norm(),
            (yfinal.tail<3>() - y0.tail<3>()).norm() / y0.tail<3>().norm()};
}

static bool needs_reversibility(const std::string& name) {
    return name == "Apophis" || name == "Achilles" || name == "Pluto";
}

static void run_reversibility_test(Integrator& integ, const AsteroidSpec& ast,
                                    const DerivativeFunction& f, const IntegratorSpec& spec,
                                    std::vector<CsvRow>& rows) {
    const auto [eps_r, eps_v] = reversibility_errors(integ, ast, f);
    rows.push_back({ast.name, spec.name, "reversibility_r", ast.T_yr, eps_r, "dimensionless", ""});
    rows.push_back({ast.name, spec.name, "reversibility_v", ast.T_yr, eps_v, "dimensionless", ""});
}

// ─── Test 5 — Lyapunov exponent (mLCE, Benettin 1980) ────────────────────────
static double lyapunov_interval(Integrator& integ, const DerivativeFunction& f,
                                 Eigen::VectorXd& y_nom, Eigen::VectorXd& y_pert) {
    const Eigen::VectorXd yf_nom  = integ.integrate(f, y_nom,  0.0, TAU_LYAP_D);
    const Eigen::VectorXd yf_pert = integ.integrate(f, y_pert, 0.0, TAU_LYAP_D);
    const double d_tau_au   = (yf_pert.head<3>() - yf_nom.head<3>()).norm();
    const double lambda_i_yr = std::log(d_tau_au / D0_LYAP_AU) * DAYS_PER_YEAR / TAU_LYAP_D;
    const double scale = D0_LYAP_AU / std::max(d_tau_au, 1e-30);
    y_pert = yf_nom + (yf_pert - yf_nom) * scale;
    y_nom  = yf_nom;
    return lambda_i_yr;
}

static void run_lyapunov_test(Integrator& integ, const AsteroidSpec& ast,
                               const DerivativeFunction& f, const IntegratorSpec& spec,
                               std::vector<CsvRow>& rows, std::vector<LyapRow>& lyap_rows) {
    Eigen::VectorXd y_nom  = to_eigen(ast.xv_au_d);
    Eigen::VectorXd y_pert = y_nom; y_pert[0] += D0_LYAP_AU;
    double lambda_sum_yr = 0.0;
    for (int i = 0; i < ast.lyap_N; ++i) {
        const double lambda_i_yr = lyapunov_interval(integ, f, y_nom, y_pert);
        lambda_sum_yr += lambda_i_yr;
        lyap_rows.push_back({ast.name, spec.name, i + 1, lambda_i_yr});
    }
    const double mLCE_yr = lambda_sum_yr / static_cast<double>(ast.lyap_N);
    rows.push_back({ast.name, spec.name, "lyapunov", ast.T_yr, mLCE_yr, "1/yr", ""});
}

// ─── Dispatcher ───────────────────────────────────────────────────────────────
static void run_all_tests(const AsteroidSpec& ast, const IntegratorSpec& spec,
                          std::vector<CsvRow>& rows, std::vector<LyapRow>& lyap_rows) {
    auto integ       = make_integrator(spec, ast.category);
    const auto f     = make_dynamics();
    const auto run   = run_trajectory(*integ, ast, f);
    run_energy_test(run, ast, spec, rows);
    run_angular_momentum_test(run, ast, spec, rows);
    run_secular_drift_test(run, ast, spec, rows);
    run_jacobi_crtbp_test(ast, spec, rows);
    if (needs_reversibility(ast.name))
        run_reversibility_test(*integ, ast, f, spec, rows);
    run_lyapunov_test(*integ, ast, f, spec, rows, lyap_rows);
}

// ─── CSV output ───────────────────────────────────────────────────────────────
static void write_csv(const std::vector<CsvRow>& rows, const std::string& path_out) {
    std::ofstream out(path_out);
    out << "asteroid,integrator,test,T_yr,value,unit,note\n"
        << std::setprecision(15) << std::scientific;
    for (const auto& r : rows)
        out << r.asteroid << "," << r.integrator << "," << r.test << ","
            << r.T_yr << "," << r.value << "," << r.unit << "," << r.note << "\n";
    std::cout << "  wrote " << path_out << "  (" << rows.size() << " rows)\n";
}

static void write_lyap_csv(const std::vector<LyapRow>& rows, const std::string& path_out) {
    std::ofstream out(path_out);
    out << "asteroid,integrator,interval_index,lambda_i\n"
        << std::setprecision(10) << std::scientific;
    for (const auto& r : rows)
        out << r.asteroid << "," << r.integrator << ","
            << r.interval_index << "," << r.lambda_i_yr << "\n";
    std::cout << "  wrote " << path_out << "  (" << rows.size() << " rows)\n";
}

} // namespace

int main() {
    std::filesystem::create_directories(OUT_DIR);

    // [2b] Get Jupiter theta_0 from Horizons at MJD 60310.0
    {
        io::HorizonsClient hzn;
        auto res = hzn.query_vectors("599", time::EpochTDB::from_mjd(60310.0), "500");
        if (res) {
            Eigen::Vector3d r_icrf = res->to_eigen_au_aud().head<3>();
            Eigen::Vector3d r_ecl = rotate_icrf_to_eclip(r_icrf);
            THETA0_JUPITER_RAD = std::atan2(r_ecl.y(), r_ecl.x());
            A_JUPITER_AU       = r_ecl.norm();
            N_JUPITER_DAV      = std::sqrt(GMS_AU3D2 / std::pow(A_JUPITER_AU, 3));
            std::cout << "[2b] Jupiter at MJD 60310.0: theta=" << THETA0_JUPITER_RAD 
                      << " rad, a=" << A_JUPITER_AU << " AU, n=" << N_JUPITER_DAV << " rad/d\n";
        } else {
            std::cerr << "[2b] ERROR: Failed to get Jupiter state from Horizons. Using fallback constants.\n";
            THETA0_JUPITER_RAD = 0.75676560; 
            N_JUPITER_DAV      = std::sqrt(GMS_AU3D2 / (A_JUPITER_AU * A_JUPITER_AU * A_JUPITER_AU));
        }
    }

    const std::vector<IntegratorSpec> INT_SPECS = {
        {"AAS",   1e-4},
        {"SABA4", 0.0},
        {"RKF78", 1e-10},
    };
    std::vector<CsvRow>  rows;
    std::vector<LyapRow> lyap_rows;
    for (const auto& ast : ASTEROIDS) {
        std::cout << "[LT] " << ast.name << " (" << ast.category
                  << ", T=" << ast.T_yr << " yr)\n";
        for (const auto& spec : INT_SPECS) {
            std::cout << "  [" << spec.name << "] ..."; std::cout.flush();
            run_all_tests(ast, spec, rows, lyap_rows);
            std::cout << " done\n";
        }
    }
    write_csv(rows,      (OUT_DIR / "long_term_tests.csv").string());
    write_lyap_csv(lyap_rows, (OUT_DIR / "lyapunov_series.csv").string());
    std::cout << "All done.\n";
    return 0;
}
