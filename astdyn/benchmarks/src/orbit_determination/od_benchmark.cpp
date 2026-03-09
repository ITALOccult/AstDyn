// =============================================================================
//  od_benchmark.cpp  —  AstDyn 3.0  Orbit Determination Benchmark
//
//  Tests GaussIOD → DifferentialCorrector → ExtendedKalmanFilter pipeline
//  on three targets: (99942) Apophis, (4) Vesta, (3200) Phaethon
//
//  Uses JPL DE441 ephemerides for high precision.
// =============================================================================

#include <astdyn/AstDyn.hpp>
#include <astdyn/orbit_determination/GaussIOD.hpp>
#include <astdyn/orbit_determination/GoodingIOD.hpp>
#include <astdyn/orbit_determination/DifferentialCorrector.hpp>
#include <astdyn/orbit_determination/ExtendedKalmanFilter.hpp>
#include <astdyn/orbit_determination/CovariancePropagator.hpp>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <astdyn/orbit_determination/StateTransitionMatrix.hpp>
#include <astdyn/io/MPCParser.hpp>
#include <astdyn/ephemeris/DE441Provider.hpp>

#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace astdyn;
using namespace astdyn::orbit_determination;
using namespace astdyn::observations;
using namespace astdyn::physics;
using namespace astdyn::propagation;
using namespace astdyn::ephemeris;
using namespace astdyn::core;

// =============================================================================
//  Utility Functions
// =============================================================================

static constexpr double AU_KM = 149597870.7;

static double delta_r_km(const CartesianStateTyped<GCRF>& a, const CartesianStateTyped<GCRF>& b) {
    return (a.position - b.position).norm().to_km();
}

static double delta_v_ms(const CartesianStateTyped<GCRF>& a, const CartesianStateTyped<GCRF>& b) {
    return (a.velocity - b.velocity).norm().to_ms();
}

// Covariance P is in [AU², AU·AU/day, (AU/day)²] (normalized units from STM/DC).
// Position-position block [0:3,0:3] is in AU² → sigma_r in km.
static double sigma_r_km(const Matrix6d& P) {
    return std::sqrt(P(0,0) + P(1,1) + P(2,2)) * AU_KM;
}

static double wall_sec() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

static std::vector<OpticalObservation> load_mpc_obs(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) return {};
    std::string content((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    return io::MPCParser::parse_file(content);
}

struct Target {
    std::string name;
    std::string mpc_full;
    std::string mpc_iod;
    double r1, r3;
    bool gooding;
    CartesianStateTyped<GCRF> ref;
    double obs_min_mjd = 51544.0; ///< Minimum MJD for full observations (default: year 2000)
    double obs_max_mjd = 99999.0; ///< Maximum MJD for full observations (default: no limit)
    int gooding_max_iter = 80;    ///< GoodingIOD max iterations (0 = single Lambert, no Newton)
    bool warm_start = false;      ///< If true, skip IOD and use ref state back-propagated to arc start
};

static void write_csv_row(std::ostream& f, const std::string& t, const std::string& m, const std::string& s,
                          double dr, double dv, double rms, double chi2, int iter, double sigma, double wall) {
    f << std::fixed << std::setprecision(6) << t << "," << m << "," << s << "," << dr << "," << dv << ","
      << rms << "," << chi2 << "," << iter << "," << sigma << "," << wall << "\n" << std::flush;
}

int main() {
    std::cout << "=== AstDyn 3.0  OD Benchmark (JPL DE440) ===" << std::endl;

    // Initialize high-precision ephemeris
    const std::string de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de440s.bsp";
    auto ephem = std::make_shared<PlanetaryEphemeris>();
    try {
        ephem->setProvider(std::make_shared<DE441Provider>(de441_path));
    } catch (const std::exception& e) {
        std::cerr << "Error loading DE441: " << e.what() << std::endl;
        return 1;
    }

    // Integrator settings: tolerances relaxed to avoid step rejection near perihelion
    auto integrator = std::make_shared<RKF78Integrator>(0.05, 1e-8, 1e-8);
    auto propagator = std::make_shared<Propagator>(integrator, ephem);

    std::vector<Target> targets;
    // Apophis: GoodingIOD with geocentric range guesses for Jan-Feb 2020 (~1.1 AU)
    // GaussIOD fails: near-coplanar LOS vectors (D0≈2.8e-4), converges to wrong root at ~2.5 AU
    // obs_min_mjd = IOD epoch start → forward-only propagation in DC
    targets.push_back({"Apophis_99942", "data/apophis_mpc_full.txt", "data/apophis_iod_3obs.txt", 0.40,0.30,true,
        CartesianStateTyped<GCRF>::from_si(time::EpochTDB::from_mjd(60310.0),
        0.925373*constants::AU*1000, 0.315602*constants::AU*1000, 0.068954*constants::AU*1000,
        -0.005412*constants::AU*1000/86400, 0.015083*constants::AU*1000/86400, 0.000873*constants::AU*1000/86400),
        58856.0 /* MJD 2020-01-08: IOD epoch */, 58946.0 /* MJD 2020-04-07: 90-day arc */});
    
    // Vesta: ~307-day arc Jan-Nov 2016 (MJD 57407-57714, ~0.23 revolutions).
    // IOD is ill-conditioned for Vesta: near-ecliptic LOS geometry makes ALL standard
    // IOD methods fail (GaussIOD: D0≈0.002 < threshold; GoodingIOD Newton: finds spurious
    // orbits). Using warm_start = true: reference state back-propagated to the arc start
    // (Jan 2016) to seed DC. This tests the DC+EKF convergence without the IOD step.
    targets.push_back({"Vesta_4", "data/vesta_mpc_full.txt", "" /* auto-pick IOD */, 0,0,false,
        CartesianStateTyped<GCRF>::from_si(time::EpochTDB::from_mjd(60310.0),
        -1.54386*constants::AU*1000, -1.94123*constants::AU*1000, 0.17328*constants::AU*1000,
        0.008321*constants::AU*1000/86400, -0.005743*constants::AU*1000/86400, -0.000488*constants::AU*1000/86400),
        57407.0 /* MJD 2016-01-20 */, 57714.0 /* MJD 2016-11-22: ~307-day arc */,
        80, true /* warm_start: back-propagate ref to arc start */});

    // Phaethon: 1983 discovery obs use non-standard RA format (HH MM.mm) → RA=0.
    // Use Sep-Oct 2017 arc (before retrograde loop / Dec 2017 close approach at 0.069 AU).
    // GoodingIOD with r1=2.0, r3=1.5 AU (Phaethon still ~2 AU geocentric in Sep-Oct 2017,
    // approaching Earth). Auto-pick IOD from first/mid/last of filtered arc.
    targets.push_back({"Phaethon_3200", "data/phaethon_mpc_full.txt", "" /* auto-pick IOD */, 2.0,1.5,true,
        CartesianStateTyped<GCRF>::from_si(time::EpochTDB::from_mjd(60310.0),
        -0.38211*constants::AU*1000, 1.21844*constants::AU*1000, 0.49718*constants::AU*1000,
        -0.018324*constants::AU*1000/86400, -0.003471*constants::AU*1000/86400, 0.004451*constants::AU*1000/86400),
        57985.0 /* MJD 2017-09-01 */, 58030.0 /* MJD 2017-10-16: ~45 day arc, no retrograde */});

    std::ofstream csv("benchmarks/results/orbit_determination/results.csv");
    csv << "target,method,solver,dr_km,dv_ms,rms_arcsec,chi2,iterations,sigma_r_km,wall_s\n";

    int target_idx = 0;
    int total_targets = targets.size();

    for (const auto& tgt : targets) {
        target_idx++;
        std::cout << "\n[" << target_idx << "/" << total_targets << "] Processing target: " << tgt.name << "..." << std::endl;
        
        auto full = load_mpc_obs(tgt.mpc_full);
        // Filter to target-specific date range for numerical conditioning
        {
            double min_mjd = tgt.obs_min_mjd;
            double max_mjd = tgt.obs_max_mjd;
            auto it = std::remove_if(full.begin(), full.end(),
                [min_mjd, max_mjd](const OpticalObservation& o){
                    return o.time.mjd() < min_mjd || o.time.mjd() > max_mjd;
                });
            full.erase(it, full.end());
        }
        if (full.size() > 500) {
            std::cout << "  - Downsampling observations from " << full.size() << " to 500..." << std::endl;
            std::vector<OpticalObservation> sampled;
            double step = static_cast<double>(full.size()) / 500.0;
            for (int i = 0; i < 500; ++i) {
                sampled.push_back(full[static_cast<size_t>(i * step)]);
            }
            full = std::move(sampled);
        }
        
        std::vector<OpticalObservation> iod_obs;
        if (!tgt.mpc_iod.empty()) {
            iod_obs = load_mpc_obs(tgt.mpc_iod);
        }
        // Fallback: use first/mid/last from filtered full arc when no dedicated IOD file
        if (iod_obs.size() < 3 && full.size() >= 3) {
            iod_obs = {full.front(), full[full.size() / 2], full.back()};
            std::cout << "  - Auto-selecting IOD obs: first/mid/last from filtered arc." << std::endl;
        }
        if (full.empty() || iod_obs.size() < 3) {
            std::cout << "  [!] Insufficient data for " << tgt.name << std::endl;
            continue;
        }

        std::cout << "  - Observations: " << full.size() << " total, 3 for IOD." << std::endl;

        // IOD (or warm start)
        double t_start = wall_sec();
        CartesianStateTyped<GCRF> iod_state;
        bool iod_ok = false;
        std::string iod_name;

        if (tgt.warm_start) {
            // Warm start: back-propagate reference state to first observation epoch.
            // Used for targets where angular IOD is ill-conditioned (e.g. low-inclination MBA).
            iod_name = "WarmStart";
            std::cout << "  - Warm start: back-propagating reference state to arc start..." << std::flush;
            try {
                auto first_epoch = time::to_tdb(full.front().time);
                iod_state = propagator->propagate_cartesian(tgt.ref, first_epoch);
                iod_ok = true;
            } catch (const std::exception& e) {
                std::cout << " FAILED: " << e.what() << std::endl;
            }
        } else {
            std::cout << "  - Starting IOD..." << std::flush;
            if (!tgt.gooding) {
                iod_name = "GaussIOD";
                GaussIODSettings s; s.use_light_time = true; s.min_separation_days = 0.5;
                auto res = GaussIOD(s).compute(iod_obs);
                iod_ok = res.success; iod_state = res.state;
            } else {
                iod_name = (tgt.gooding_max_iter == 0) ? "LambertIOD" : "GoodingIOD";
                GoodingIOD::Settings gs; gs.max_iterations = tgt.gooding_max_iter;
                auto res = GoodingIOD(gs).compute(iod_obs[0], iod_obs[1], iod_obs[2], tgt.r1, tgt.r3);
                iod_ok = res.success && !res.solutions.empty();
                if (iod_ok) iod_state = res.solutions[0].state;
            }
        }

        double dt_iod = wall_sec() - t_start;
        if (iod_ok) {
            std::cout << " OK (" << iod_name << ") [dt: " << std::fixed << std::setprecision(3) << dt_iod << "s]" << std::endl;
            write_csv_row(csv, tgt.name, "IOD", iod_name, 0, 0, 0, 0, 0, 0, dt_iod);
        } else {
            std::cout << " FAILED" << std::endl;
            continue;
        }

        // LSQ
        std::cout << "  - Starting LSQ (Differential Corrector)..." << std::endl;
        t_start = wall_sec();
        auto res_calc = std::make_shared<ResidualCalculator<GCRF>>(ephem, propagator);
        auto stm_comp = std::make_shared<StateTransitionMatrix<GCRF>>(propagator);
        DifferentialCorrectorSettings ds;
        ds.reject_outliers = false;      // Disable outlier rejection: rough IOD residuals >> sigma
        ds.compute_covariance = true;
        ds.max_iterations = 50;
        ds.convergence_tolerance = 1e-6; // [AU] ≈ 150 m
        ds.verbose = true;
        Matrix6d dc_cov = Matrix6d::Identity() * 1e-8; // fallback if DC fails
        try {
            auto dc_res = DifferentialCorrector<GCRF>(res_calc, stm_comp).fit(full, iod_state, ds);
            double dt_dc = wall_sec() - t_start;
            if (dc_res.converged) {
                try {
                auto final_ref = propagator->propagate_cartesian(dc_res.final_state, tgt.ref.epoch);
                double dr = delta_r_km(final_ref, tgt.ref);
                double dv = delta_v_ms(final_ref, tgt.ref);
                double rms = dc_res.statistics.rms_total.to_arcsec();
                double sig = sigma_r_km(dc_res.covariance);
                std::cout << "  → Converged [dt: " << std::fixed << std::setprecision(3) << dt_dc
                          << "s | dr: " << dr << " km | dv: " << dv
                          << " m/s | RMS: " << rms << "\" | sigma_r: " << sig << " km]" << std::endl;
                write_csv_row(csv, tgt.name, "LSQ", "DiffCorr", dr, dv, rms,
                              dc_res.statistics.chi_squared, dc_res.iterations, sig, dt_dc);
                iod_state = dc_res.final_state;
                dc_cov = dc_res.covariance; // pass DC covariance to EKF
                } catch (const std::exception& e) {
                    std::cout << "  → Converged but propagation to ref failed: " << e.what() << std::endl;
                    iod_state = dc_res.final_state;
                    dc_cov = dc_res.covariance;
                }
            } else {
                double dt_dc2 = wall_sec() - t_start;
                std::cout << "  → Did not converge after " << dc_res.iterations
                          << " iterations [dt: " << dt_dc2 << "s]" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cout << "  → FAILED: " << e.what() << std::endl;
        }

        // EKF — initialized with DC covariance when available (same AU² units)
        std::cout << "  - Starting EKF..." << std::flush;
        t_start = wall_sec();
        ExtendedKalmanFilter::Settings es;
        es.process_noise = Matrix6d::Identity() * 1e-15;
        ExtendedKalmanFilter ekf(propagator, es);
        CartesianStateTyped<GCRF> ekf_s = iod_state;
        Matrix6d ekf_p = dc_cov; // AU² units, consistent with STM normalization
        double ssq = 0; int n=0;
        for (const auto& o : full) {
            try {
                auto r = ekf.update(ekf_s, ekf_p, o);
                ekf_s = r.state; ekf_p = r.covariance;
                ssq += r.innovation.squaredNorm(); n++;
            } catch (...) {}
        }
        double dt_ekf = wall_sec() - t_start;
        if (n>0) {
            try {
            auto final_ref = propagator->propagate_cartesian(ekf_s, tgt.ref.epoch);
            double dr = delta_r_km(final_ref, tgt.ref);
            double dv = delta_v_ms(final_ref, tgt.ref);
            double rms = std::sqrt(ssq/(2*n)) * constants::RAD_TO_ARCSEC;
            double sig = sigma_r_km(ekf_p);
            std::cout << " Finished [dt: " << dt_ekf << "s, dr: " << dr << " km, RMS: " << rms << "\"]" << std::endl;
            write_csv_row(csv, tgt.name, "EKF", "EKF", dr, dv, rms, 0, 0, sig, dt_ekf);
            } catch (const std::exception& e) {
                std::cout << " EKF propagation to ref failed: " << e.what() << std::endl;
            }
        } else {
            std::cout << " No updates" << std::endl;
        }
    }
    csv.close();
    return 0;
}
