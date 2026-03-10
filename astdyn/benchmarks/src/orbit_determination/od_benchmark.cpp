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
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/propagation/GaussIntegrator.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/orbit_determination/ODSmartPolicy.hpp>

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

// Utility functions

static double delta_r_km(const CartesianStateTyped<GCRF>& a, const CartesianStateTyped<GCRF>& b) {
    return (a.position - b.position).norm().to_km();
}

static double delta_v_ms(const CartesianStateTyped<GCRF>& a, const CartesianStateTyped<GCRF>& b) {
    return (a.velocity - b.velocity).norm().to_ms();
}

// Covariance P is in [m², m·m/s, (m/s)²] (SI units from DC).
// Position-position block [0:3,0:3] is in m² → sigma_r in km.
static double sigma_r_km(const Matrix6d& P) {
    return std::sqrt(P(0,0) + P(1,1) + P(2,2)) / 1000.0;
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
    double yarkovsky_a2 = 0.0;   ///< Yarkovsky transverse parameter [AU/d²], 0 = disabled
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

    // --- Base propagator settings (full N-body + GR + Moon) ---
    PropagatorSettings base_settings;
    base_settings.include_planets  = true;
    base_settings.perturb_mercury  = true;
    base_settings.perturb_venus    = true;
    base_settings.perturb_earth    = true;
    base_settings.perturb_mars     = true;
    base_settings.perturb_jupiter  = true;
    base_settings.perturb_saturn   = true;
    base_settings.perturb_uranus   = true;
    base_settings.perturb_neptune  = true;
    base_settings.include_relativity = true;
    base_settings.include_moon     = true;

    // Factory: creates an RKF78-based propagator, optionally with Yarkovsky.
    // RKF78 (explicit 7/8th-order adaptive) is the standard integrator for OD:
    // fast, accurate, and used by OrbFit/JPL for arcs up to ~10 years.
    // GaussIntegrator would be ~10x slower here with no meaningful accuracy gain
    // for multi-year OD arcs.
    // The StateTransitionMatrix uses its own separate RKF78 for the 42-component ODE.
    auto make_propagator = [&](double yarkovsky_a2 = 0.0) {
        auto integr = std::make_shared<RKF78Integrator>(
            0.5,    // initial step [days]
            1e-12,  // relative tolerance
            1e-12   // absolute tolerance
        );
        PropagatorSettings s = base_settings;
        s.include_yarkovsky = (std::abs(yarkovsky_a2) > 0.0);
        s.yarkovsky_a2      = yarkovsky_a2;
        return std::make_shared<Propagator>(integr, ephem, s);
    };

    std::vector<Target> targets;
    // All reference states from JPL Horizons (heliocentric ICRF/J2000, Sun-centred).
    // Epochs are chosen to be within or near the observation arc to minimise warm-start
    // back-propagation time (<50 days vs. the previous 4-year back-propagation that gave
    // completely wrong positions due to incorrect placeholder state vectors).

    // Apophis — arc MJD 58856–58946 (2020-Jan-08 to 2020-Apr-07, 90 days)
    // Reference: JPL Horizons heliocentric ICRF, 2020-Jan-24 00:00 TDB (MJD 58872.0)
    {
        // Position [km], velocity [km/s] directly from Horizons (no frame conversion needed)
        Eigen::Vector3d r_km(-1.541414301798325E+07,  1.394708232005084E+08,  5.145546099858084E+07);
        Eigen::Vector3d v_kms(-2.845748009596083E+01,  2.125072123515800E+00,  6.954818026860536E-02);

        targets.push_back({"Apophis_99942",
            "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/apophis_mpc_full.txt",
            "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/apophis_iod_3obs.txt",
            1.1, 1.0, true,
            CartesianStateTyped<GCRF>::from_si(time::EpochTDB::from_mjd(58872.0),
                r_km.x() * 1000.0, r_km.y() * 1000.0, r_km.z() * 1000.0,
                v_kms.x() * 1000.0, v_kms.y() * 1000.0, v_kms.z() * 1000.0),
            58856.0 /* MJD 2020-01-08: arc start */, 58946.0 /* MJD 2020-04-07: arc end */,
            80, true, -2.901e-13 /* Yarkovsky A2 [AU/d²] from JPL solution */});
    }

    // Vesta — arc MJD 57407–57714 (2016-Jan-20 to 2016-Nov-22, 307 days)
    // Reference: JPL Horizons heliocentric ICRF, 2016-Mar-01 00:00 TDB (MJD 57448.0)
    {
        Eigen::Vector3d r_km( 2.737877963027807E+08,  2.540494479148592E+08,  6.538351637650730E+07);
        Eigen::Vector3d v_kms(-1.159539810770006E+01,  1.221289828462967E+01,  6.383069804659519E+00);

        targets.push_back({"Vesta_4",
            "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/vesta_mpc_full.txt",
            "" /* auto-pick IOD */,
            2.4, 2.3, true,
            CartesianStateTyped<GCRF>::from_si(time::EpochTDB::from_mjd(57448.0),
                r_km.x() * 1000.0, r_km.y() * 1000.0, r_km.z() * 1000.0,
                v_kms.x() * 1000.0, v_kms.y() * 1000.0, v_kms.z() * 1000.0),
            57407.0 /* MJD 2016-01-20: arc start */, 57714.0 /* MJD 2016-11-22: arc end */,
            80, true});
    }

    // Phaethon — arc MJD 57985–58085 (2017-Aug-20 to 2017-Nov-28, 100 days)
    // Reference: JPL Horizons heliocentric ICRF, 2017-Sep-15 00:00 TDB (MJD 58011.0)
    {
        Eigen::Vector3d r_km( 1.390964499994716E+08,  2.198526146466129E+08,  1.474819321622483E+08);
        Eigen::Vector3d v_kms(-1.145826169909928E+01, -3.917012727493627E+00, -6.560664426417666E+00);

        targets.push_back({"Phaethon_3200",
            "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/phaethon_mpc_full.txt",
            "" /* auto-pick IOD */,
            1.2, 1.1, true,
            CartesianStateTyped<GCRF>::from_si(time::EpochTDB::from_mjd(58011.0),
                r_km.x() * 1000.0, r_km.y() * 1000.0, r_km.z() * 1000.0,
                v_kms.x() * 1000.0, v_kms.y() * 1000.0, v_kms.z() * 1000.0),
            57985.0 /* MJD 2017-08-20: arc start */, 58085.0 /* MJD 2017-Nov-28: arc end */,
            80, true, -6.142920483398e-15 /* Yarkovsky A2 [AU/d²] from JPL, g(r)=(1AU/r)^2 */});
    }

    std::ofstream csv("/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/benchmarks/results/orbit_determination/results.csv");
    csv << "target,method,solver,dr_km,dv_ms,rms_arcsec,chi2,iterations,sigma_r_km,wall_s\n";

    int target_idx = 0;
    int total_targets = targets.size();

    for (const auto& tgt : targets) {
        target_idx++;
        std::cout << "\n[" << target_idx << "/" << total_targets << "] Processing target: " << tgt.name << "..." << std::endl;

        // Build a per-target propagator (includes Yarkovsky if a2 != 0)
        auto prop = make_propagator(tgt.yarkovsky_a2);

        // DIAGNOSTIC CORE PHYSICS CHECK
        {
            auto t_ref = tgt.ref.epoch;
            auto sun_bary = ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(t_ref);
            auto earth_helio = ephemeris::PlanetaryEphemeris::getPosition(CelestialBody::EARTH, t_ref);
            std::cout << "  [USER-DIAG] Sun Barycentric:    " << (sun_bary.to_eigen_si() / (constants::AU*1000.0)).transpose() << " AU" << std::endl;
            std::cout << "  [USER-DIAG] Earth Heliocentric:  " << (earth_helio.to_eigen_si() / (constants::AU*1000.0)).transpose() << " AU | dist: " << earth_helio.norm().to_au() << " AU" << std::endl;
        }
        
        auto full = load_mpc_obs(tgt.mpc_full);
        if (full.empty()) {
            std::cout << "  [!] Error: Could not load " << tgt.mpc_full << std::endl;
        } else {
            std::cout << "  - Loaded " << full.size() << " observations from " << tgt.mpc_full << std::endl;
        }

        // Filter to target-specific date range for numerical conditioning
        {
            double min_mjd = tgt.obs_min_mjd;
            double max_mjd = tgt.obs_max_mjd;
            auto it = std::remove_if(full.begin(), full.end(),
                [min_mjd, max_mjd](const OpticalObservation& o){
                    return o.time.mjd() < min_mjd || o.time.mjd() > max_mjd;
                });
            full.erase(it, full.end());
            std::cout << "  - After filtering (" << min_mjd << " to " << max_mjd << "): " << full.size() << " obs remaining." << std::endl;
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
                iod_state = prop->propagate_cartesian(tgt.ref, first_epoch);
                iod_ok = true;
            } catch (const std::exception& e) {
                std::cout << " FAILED: " << e.what() << std::endl;
            }
        } else {
            std::cout << "  - Starting IOD..." << std::flush;
            if (!tgt.gooding) {
                iod_name = "GaussIOD";
                GaussIODSettings s; s.use_light_time = true; s.min_separation_days = 0.5; s.verbose = true;
                auto res = GaussIOD(s).compute(iod_obs);
                iod_ok = res.success; iod_state = res.state;
            } else {
                iod_name = (tgt.gooding_max_iter == 0) ? "LambertIOD" : "GoodingIOD";
                GoodingIOD::Settings gs; gs.max_iterations = tgt.gooding_max_iter; gs.verbose = true;
                auto res = GoodingIOD(gs).compute(iod_obs[0], iod_obs[1], iod_obs[2], tgt.r1, tgt.r3);
                iod_ok = res.success && !res.solutions.empty();
                if (iod_ok) iod_state = res.solutions[0].state;
            }
        }

        double dt_iod = wall_sec() - t_start;
        if (iod_ok) {
            std::cout << " OK (" << iod_name << ") [dt: " << std::fixed << std::setprecision(3) << dt_iod << "s]" << std::endl;

            // Rule 3 — Apply automatic MPC-style weights to all observations
            ODPolicyEngine::apply_auto_weights(full);

            // Rule 1+2 — Analyze orbit and arc, then rebuild propagator if needed
            auto policy = ODPolicyEngine::analyze(
                physics::CartesianStateTyped<core::GCRF>::from_si(
                    iod_state.epoch,
                    iod_state.position.x_si(), iod_state.position.y_si(), iod_state.position.z_si(),
                    iod_state.velocity.x_si(), iod_state.velocity.y_si(), iod_state.velocity.z_si(),
                    iod_state.gm.to_m3_s2()),
                full, base_settings);
            policy.print();

            if (!policy.arc.is_admissible) {
                std::cout << "  [!] Arc inadmissible: " << policy.arc.rejection_reason
                          << " — DC may produce unreliable results.\n";
            }

            // Rebuild propagator if force model changed (e.g. PPN toggled on)
            if (policy.propagator_settings.include_relativity != base_settings.include_relativity) {
                std::cout << "  [ODPolicy] Rebuilding propagator with updated force model.\n";
                auto integr2 = std::make_shared<RKF78Integrator>(
                    policy.force.suggested_step_days, 1e-12, 1e-12);
                PropagatorSettings new_s = policy.propagator_settings;
                new_s.include_yarkovsky = (std::abs(tgt.yarkovsky_a2) > 0.0);
                new_s.yarkovsky_a2      = tgt.yarkovsky_a2;
                prop = std::make_shared<Propagator>(integr2, ephem, new_s);
            }
            
            // Back-propagate or forward-propagate reference to check IOD quality
            try {
                auto iod_ref_comp = prop->propagate_cartesian(iod_state, tgt.ref.epoch);
                double dr_iod = delta_r_km(iod_ref_comp, tgt.ref);
                double dv_iod = delta_v_ms(iod_ref_comp, tgt.ref);
                write_csv_row(csv, tgt.name, "IOD", iod_name, dr_iod, dv_iod, 0, 0, 0, 0, dt_iod);
                std::cout << "  - IOD Error @ Ref Epoch: dr=" << dr_iod << " km, dv=" << dv_iod << " m/s" << std::endl;
            } catch (...) {
                write_csv_row(csv, tgt.name, "IOD", iod_name, 0, 0, 0, 0, 0, 0, dt_iod);
            }
        } else {
            std::cout << " FAILED" << std::endl;
            continue;
        }

        // Residual diagnostic: print first observation vs computed for debugging
        {
            auto diag_calc = ResidualCalculator<GCRF>(ephem, prop);
            auto sorted_diag = full;
            std::sort(sorted_diag.begin(), sorted_diag.end(), [](const auto& a, const auto& b){ return a.time < b.time; });
            auto first_res = diag_calc.compute_residuals({sorted_diag.front()}, iod_state);
            if (!first_res.empty()) {
                const auto& fr = first_res[0];
                std::cout << "  [DIAG] 1st obs:"
                          << " obs_ra="   << std::fixed << std::setprecision(4) << sorted_diag.front().ra.to_deg()  << "°"
                          << " obs_dec="  << sorted_diag.front().dec.to_deg()  << "°"
                          << " cmp_ra="   << fr.computed_ra.to_deg()   << "°"
                          << " cmp_dec="  << fr.computed_dec.to_deg()  << "°"
                          << " res_ra="   << std::setprecision(2) << fr.residual_ra.to_arcsec()  << "\""
                          << " res_dec="  << fr.residual_dec.to_arcsec() << "\"" << std::endl;
                std::cout << "  [DIAG] iod_state epoch MJD=" << iod_state.epoch.mjd()
                          << " pos_au=(" << iod_state.position.to_eigen_si().x() / constants::AU / 1000.0
                          << "," << iod_state.position.to_eigen_si().y() / constants::AU / 1000.0
                          << "," << iod_state.position.to_eigen_si().z() / constants::AU / 1000.0 << ")" << std::endl;
            }
        }

        // LSQ
        std::cout << "  - Starting LSQ (Differential Corrector)..." << std::endl;
        t_start = wall_sec();
        auto res_calc = std::make_shared<ResidualCalculator<GCRF>>(ephem, prop);
        auto stm_comp = std::make_shared<StateTransitionMatrix<GCRF>>(prop);
        DifferentialCorrectorSettings ds;
        ds.reject_outliers = false;      // Disable outlier rejection: rough IOD residuals >> sigma
        ds.compute_covariance = true;
        ds.max_iterations = 50;
        ds.convergence_tolerance = physics::Distance::from_au(1e-6); 
        ds.verbose = true;
        Matrix6d dc_cov = Matrix6d::Identity() * 1e-8; // fallback if DC fails
        try {
            auto dc_res = DifferentialCorrector<GCRF>(res_calc, stm_comp).fit(full, iod_state, ds);
            double dt_dc = wall_sec() - t_start;
            if (dc_res.converged) {
                try {
                auto final_ref = prop->propagate_cartesian(dc_res.final_state, tgt.ref.epoch);
                double dr = delta_r_km(final_ref, tgt.ref);
                double dv = delta_v_ms(final_ref, tgt.ref);
                double rms = dc_res.statistics.rms_total.to_arcsec();
                double sig = sigma_r_km(dc_res.covariance);
                
                std::cout << "  → Converged [dt: " << std::fixed << std::setprecision(3) << dt_dc
                          << "s | dr: " << dr << " km | dv: " << dv
                          << " m/s | RMS: " << rms << "\" | sigma_r: " << sig << " km]" << std::endl;
                 
                std::cout << "    [USER-DIAG] final_ref (AU): " << (final_ref.position.to_eigen_si() / (constants::AU*1000.0)).transpose() << " | dist=" << (final_ref.position.norm().to_au()) << " AU" << std::endl;
                std::cout << "    [USER-DIAG] horizons  (AU): " << (tgt.ref.position.to_eigen_si() / (constants::AU*1000.0)).transpose() << " | dist=" << (tgt.ref.position.norm().to_au()) << " AU" << std::endl;
                
                // Kep elements check
                auto mu_sun = physics::GravitationalParameter::sun().to_m3_s2();
                double r_mag = final_ref.position.norm().to_m();
                double v_mag2 = final_ref.velocity.to_eigen_si().squaredNorm();
                double a_inv = 2.0 / r_mag - v_mag2 / mu_sun;
                double a_au = (1.0 / a_inv) / (constants::AU * 1000.0);
                std::cout << "    [USER-DIAG] Semi-major axis (a): " << a_au << " AU" << std::endl;
                
                std::cout << "    [USER-DIAG] Epoch: " << final_ref.epoch.mjd() << " MJD vs Ref: " << tgt.ref.epoch.mjd() << " MJD" << std::endl;
                
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

        // EKF
        std::cout << "  - Starting EKF..." << std::flush;
        t_start = wall_sec();
        ExtendedKalmanFilter::Settings es;
        es.process_noise = Matrix6d::Identity() * 1e-18;
        ExtendedKalmanFilter ekf(prop, es);

        CartesianStateTyped<GCRF> ekf_state = iod_state;
        Matrix6d ekf_cov = dc_cov; 
        double ssq_rad = 0;
        int n_updates = 0;

        for (const auto& obs : full) {
            try {
                auto update_res = ekf.update(ekf_state, ekf_cov, obs);
                ekf_state = update_res.state;
                ekf_cov = update_res.covariance;
                ssq_rad += update_res.innovation.squaredNorm();
                n_updates++;
            } catch (...) {
                // Skip failed updates
            }
        }
        
        double dt_ekf = wall_sec() - t_start;

        if (n_updates > 0) {
            try {
                auto final_ref = prop->propagate_cartesian(ekf_state, tgt.ref.epoch);
                double dr = delta_r_km(final_ref, tgt.ref);
                double dv = delta_v_ms(final_ref, tgt.ref);
                double rms_arcsec = std::sqrt(ssq_rad / (2.0 * n_updates)) * constants::RAD_TO_ARCSEC;
                
                std::cout << " Finished [dt: " << std::fixed << std::setprecision(3) << dt_ekf
                          << "s, dr: " << dr << " km, RMS: " << rms_arcsec << "\"]" << std::endl;
                
                write_csv_row(csv, tgt.name, "EKF", "EKF", dr, dv, rms_arcsec, 0, 0, sigma_r_km(ekf_cov), dt_ekf);
            } catch (const std::exception& e) {
                std::cout << " Finished but propagation to ref failed: " << e.what() << std::endl;
            }
        } else {
            std::cout << " FAILED (no updates)" << std::endl;
        }
    }
    csv.close();
    return 0;
}
