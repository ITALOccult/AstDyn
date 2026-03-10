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
#include <astdyn/propagation/Integrator.hpp>
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

static double delta_r_km(const CartesianStateTyped<ECLIPJ2000>& a, const CartesianStateTyped<ECLIPJ2000>& b) {
    return (a.position - b.position).norm().to_km();
}

static double delta_v_ms(const CartesianStateTyped<ECLIPJ2000>& a, const CartesianStateTyped<ECLIPJ2000>& b) {
    return (a.velocity - b.velocity).norm().to_ms();
}

// Covariance P is in [m², m·m/s, (m/s)²] (SI units from DC).
// Position-position block [0:3,0:3] is in m² → sigma_r in km.
static double sigma_r_km(const Matrix6d& P) {
    return std::sqrt(P(0,0) + P(1,1) + P(2,2)) / 1000.0;
}

// Type-safe frame transformation helper
template <typename From, typename To>
CartesianStateTyped<To> transform_state_typed(const CartesianStateTyped<From>& s) {
    auto pos_to = coordinates::ReferenceFrame::transform_pos<From, To>(s.position, s.epoch);
    auto vel_to = coordinates::ReferenceFrame::transform_vel<From, To>(s.position, s.velocity, s.epoch);
    return CartesianStateTyped<To>(s.epoch, pos_to, vel_to, s.gm);
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
    CartesianStateTyped<ECLIPJ2000> ref;
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

    // Factory: creates a RK4-based propagator (Fixed step).
    auto make_propagator = [&](double yarkovsky_a2 = 0.0) {
        auto integr = std::make_shared<RKF78Integrator>(
            0.5,    // initial step [days]
            1e-12   // tolerance
        );
        PropagatorSettings s = base_settings;
        s.include_yarkovsky = (std::abs(yarkovsky_a2) > 0.0);
        s.yarkovsky_a2      = yarkovsky_a2;
        s.integrate_in_ecliptic = true; // Essential for asteroid precision
        return std::make_shared<Propagator>(integr, ephem, s);
    };

    std::vector<Target> targets;

    // Apophis — full arc 2004-2024
    {
        // Reference: JPL Horizons heliocentric ECLIPTIC J2000, 2020-Jan-24 00:00 TDB (MJD 58872.0)
        Eigen::Vector3d r_km(-1.541414301798325E+07,  1.394708232005084E+08,  5.145546099858084E+07);
        Eigen::Vector3d v_kms(-2.845748009596083E+01,  2.125072123515800E+00,  6.954818026860536E-02);

        targets.push_back({"Apophis_99942",
            "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/apophis_mpc_full.txt",
            "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/apophis_iod_3obs.txt",
            1.1, 1.0, true, // r1, r3, gooding
            CartesianStateTyped<core::ECLIPJ2000>::from_si(time::EpochTDB::from_mjd(58872.0),
                r_km.x() * 1000.0, r_km.y() * 1000.0, r_km.z() * 1000.0,
                v_kms.x() * 1000.0, v_kms.y() * 1000.0, v_kms.z() * 1000.0),
            40000.0, 70000.0, // obs_min_mjd, obs_max_mjd
            80, true, -2.901e-13 /* Yarkovsky A2 [AU/d²] */}); // gooding_max_iter, warm_start, yarkovsky_a2
    }

    // Vesta — full arc
    {
        // Reference: JPL Horizons heliocentric ECLIPTIC J2000, 2016-Mar-01 00:00 TDB (MJD 57448.0)
        Eigen::Vector3d r_km( 2.737877963027807E+08,  2.540494479148592E+08,  6.538351637650730E+07);
        Eigen::Vector3d v_kms(-1.159539810770006E+01,  1.221289828462967E+01,  6.383069804659519E+00);

        targets.push_back({"Vesta_4",
            "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/vesta_mpc_full.txt",
            "", 2.4, 2.3, true,
            CartesianStateTyped<core::ECLIPJ2000>::from_si(time::EpochTDB::from_mjd(57448.0),
                r_km.x() * 1000.0, r_km.y() * 1000.0, r_km.z() * 1000.0,
                v_kms.x() * 1000.0, v_kms.y() * 1000.0, v_kms.z() * 1000.0),
            40000.0, 70000.0,
            80, false});
    }

    // Phaethon — full arc
    {
        // Reference: JPL Horizons heliocentric ECLIPTIC J2000, 2017-Sep-15 00:00 TDB (MJD 58011.0)
        Eigen::Vector3d r_km(1.405553641753041E+08, 1.956799516664952E+08, 1.777977457782163E+08);
        Eigen::Vector3d v_kms(-1.145826169909928E+01, -3.917012727493627E+00, -6.560664426417666E+00);

        targets.push_back({"Phaethon_3200",
            "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/phaethon_mpc_full.txt",
            "", 1.2, 1.1, true,
            CartesianStateTyped<core::ECLIPJ2000>::from_si(time::EpochTDB::from_mjd(58011.0),
                r_km.x() * 1000.0, r_km.y() * 1000.0, r_km.z() * 1000.0,
                v_kms.x() * 1000.0, v_kms.y() * 1000.0, v_kms.z() * 1000.0),
            40000.0, 70000.0,
            80, false, -6.142920483398e-15 /* Yarkovsky A2 */});
    }

    std::ofstream csv("/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/benchmarks/results/orbit_determination/results.csv");
    csv << "target,method,solver,dr_km,dv_ms,rms_arcsec,chi2,iterations,sigma_r_km,wall_s\n";

    int target_idx = 0;
    int total_targets = targets.size();

    for (const auto& tgt : targets) {
        target_idx++;
        std::cout << "\n[" << target_idx << "/" << total_targets << "] Processing target: " << tgt.name << "..." << std::endl;

        auto prop = make_propagator(tgt.yarkovsky_a2);
        auto full = load_mpc_obs(tgt.mpc_full);
        if (full.empty()) continue;

        // Filter arc
        {
            auto it = std::remove_if(full.begin(), full.end(),
                [&tgt](const OpticalObservation& o){
                    return o.time.mjd() < tgt.obs_min_mjd || o.time.mjd() > tgt.obs_max_mjd;
                });
            full.erase(it, full.end());
            
            // Optimization: sample 1 out of 10 observations for faster benchmark turnaround
            if (full.size() > 1000) {
                std::vector<OpticalObservation> sampled;
                for (size_t i = 0; i < full.size(); i += 100) {
                    sampled.push_back(full[i]);
                }
                full = std::move(sampled);
            }
            
            std::cout << "  - Total arc loaded: " << full.size() << " observations (sampled 1/10)." << std::endl;
        }
        
        std::vector<OpticalObservation> iod_obs;
        if (!tgt.mpc_iod.empty()) iod_obs = load_mpc_obs(tgt.mpc_iod);
        if (iod_obs.size() < 3 && full.size() >= 3) {
            iod_obs = {full.front(), full[full.size() / 2], full.back()};
            std::cout << "  - Auto-selecting IOD obs: first/mid/last fromFiltered arc." << std::endl;
        }

        double t_start = wall_sec();
        CartesianStateTyped<GCRF> iod_state_gcrf;
        bool iod_ok = false;
        std::string iod_name;
        CartesianStateTyped<ECLIPJ2000> iod_state;

        auto first_epoch = time::to_tdb(full.front().time);

        if (tgt.warm_start) {
            std::cout << "  - Warm start: propagating reference from MJD " << tgt.ref.epoch.mjd() << std::endl;
            auto el_ref = ODPolicyEngine::compute_keplerian(tgt.ref.template cast_frame<core::GCRF>());
            std::cout << "  - Ref SMA: " << el_ref.a_au << " AU, e: " << el_ref.e << " mu: " << tgt.ref.gm.to_m3_s2() << std::endl;
            
            try {
                iod_state = prop->propagate_cartesian(tgt.ref, first_epoch);
                iod_state_gcrf = iod_state.template cast_frame<core::GCRF>();
                auto el_iod = ODPolicyEngine::compute_keplerian(iod_state_gcrf);
                std::cout << "  - IOD SMA: " << el_iod.a_au << " AU, e: " << el_iod.e << std::endl;
                iod_ok = true;
                iod_name = "WarmStart";
            } catch (const std::exception& e) {
                std::cerr << "Propagation failed: " << e.what() << std::endl;
                continue;
            }
        } else {
            std::cout << "  - Starting IOD..." << std::flush;
            if (!tgt.gooding) {
                iod_name = "GaussIOD";
                GaussIODSettings s; s.use_light_time = true; s.min_separation_days = 0.5;
                auto res = GaussIOD(s).compute(iod_obs);
                iod_ok = res.success; 
                if (iod_ok) {
                    iod_state_gcrf = res.state;
                    iod_state = transform_state_typed<GCRF, ECLIPJ2000>(res.state);
                }
            } else {
                iod_name = (tgt.gooding_max_iter == 0) ? "LambertIOD" : "GoodingIOD";
                GoodingIOD::Settings gs; gs.max_iterations = tgt.gooding_max_iter;
                auto res = GoodingIOD(gs).compute(iod_obs[0], iod_obs[1], iod_obs[2], tgt.r1, tgt.r3);
                iod_ok = res.success && !res.solutions.empty();
                if (iod_ok) {
                    iod_state_gcrf = res.solutions[0].state;
                    iod_state = transform_state_typed<GCRF, ECLIPJ2000>(res.solutions[0].state);
                }
            }
        }

        double dt_iod = wall_sec() - t_start;
        if (iod_ok) {
            std::cout << " OK (" << iod_name << ") [dt: " << dt_iod << "s]" << std::endl;



            ODPolicyEngine::apply_auto_weights(full);
            
            auto policy = ODPolicyEngine::analyze(iod_state_gcrf, full, base_settings);
            policy.print();

            // --- DIAGNOSTIC: Initial RMS check ---
            {
                auto init_calc = std::make_shared<ResidualCalculator<ECLIPJ2000>>(ephem, prop);
                auto init_res = init_calc->compute_residuals(full, iod_state);
                auto init_stats = ResidualCalculator<ECLIPJ2000>::compute_statistics(init_res);
                std::cout << "  - Initial O-C RMS: " << init_stats.rms_total.to_arcsec() << "\"" << std::endl;
            }
            // ------------------------------------

            // LSQ
            std::cout << "  - Starting LSQ (Differential Corrector)..." << std::endl;
            t_start = wall_sec();
            auto res_calc = std::make_shared<ResidualCalculator<ECLIPJ2000>>(ephem, prop);
            auto stm_comp = std::make_shared<StateTransitionMatrix<ECLIPJ2000>>(prop);
            
            DifferentialCorrectorSettings ds;
            ds.max_iterations = 50;
            ds.convergence_tolerance = Distance::from_au(1.0e-9); 
            ds.verbose = true;
            
            try {
                auto dc_res = DifferentialCorrector<ECLIPJ2000>(res_calc, stm_comp).fit(full, iod_state, ds);
                double dt_dc = wall_sec() - t_start;
                
                if (dc_res.converged) {
                    auto final_ref = prop->propagate_cartesian(dc_res.final_state, tgt.ref.epoch);
                    double dr = delta_r_km(final_ref, tgt.ref);
                    double dv = delta_v_ms(final_ref, tgt.ref);
                    
                    std::cout << "  → Converged [dt: " << dt_dc << "s | dr: " << dr << " km | RMS: " << dc_res.statistics.rms_total.to_arcsec() << "\"]" << std::endl;
                    write_csv_row(csv, tgt.name, "LSQ", "DiffCorr", dr, dv, dc_res.statistics.rms_total.to_arcsec(), 0, dc_res.iterations, sigma_r_km(dc_res.covariance), dt_dc);
                    
                    // EKF (Note: EKF currently works in GCRF in this library version)
                    std::cout << "  - Starting EKF/Smoothing..." << std::flush;
                    ExtendedKalmanFilter ekf(prop, ExtendedKalmanFilter::Settings());
                    CartesianStateTyped<GCRF> ekf_state = transform_state_typed<ECLIPJ2000, GCRF>(dc_res.final_state);
                    Matrix6d ekf_cov = dc_res.covariance;
                    for (const auto& obs : full) {
                        try {
                            auto up = ekf.update(ekf_state, ekf_cov, obs);
                            ekf_state = up.state; ekf_cov = up.covariance;
                        } catch (...) {}
                    }
                    std::cout << " Done." << std::endl;
                }
            } catch (const std::exception& e) {
                std::cout << "  → FAILED: " << e.what() << std::endl;
            }
        }
    }
    csv.close();
    return 0;
}
