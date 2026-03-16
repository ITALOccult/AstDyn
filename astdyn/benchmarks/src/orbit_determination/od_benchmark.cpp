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
#include <astdyn/orbit_determination/GaussIOD.hpp>
#include <astdyn/orbit_determination/OrbFitIOD.hpp>
#include <astdyn/propagation/OrbFitIntegrator.hpp>
#include <astdyn/core/Configurator.hpp>

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

template <typename F1, typename F2>
static double delta_r_km(const CartesianStateTyped<F1>& a, const CartesianStateTyped<F2>& b) {
    auto b_cast = b.template cast_frame<F1>();
    return (a.position.to_eigen_si() - b_cast.position.to_eigen_si()).norm() / 1000.0;
}

template <typename F1, typename F2>
static double delta_v_ms(const CartesianStateTyped<F1>& a, const CartesianStateTyped<F2>& b) {
    auto b_cast = b.template cast_frame<F1>();
    return (a.velocity.to_eigen_si() - b_cast.velocity.to_eigen_si()).norm();
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
    CartesianStateTyped<GCRF> ref; // Ref states from Horizons are usually GCRF (J2000 Equatorial)
    double obs_min_mjd = 51544.0;
    double obs_max_mjd = 99999.0;
    int gooding_max_iter = 80;
    bool warm_start = false;
    double yarkovsky_a2 = 0.0;
};

static void write_csv_row(std::ostream& f, const std::string& t, const std::string& m, const std::string& s,
                          double dr, double dv, double rms, double chi2, int iter, double sigma, double wall) {
    f << std::fixed << std::setprecision(6) << t << "," << m << "," << s << "," << dr << "," << dv << ","
      << rms << "," << chi2 << "," << iter << "," << sigma << "," << wall << "\n" << std::flush;
}

int main(int argc, char** argv) {
    std::string filter = "";
    std::string config_file = "";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--target=") == 0) {
            filter = arg.substr(9);
        } else if (arg == "-t" && i + 1 < argc) {
            filter = argv[++i];
        } else if (arg.find("--config=") == 0) {
            config_file = arg.substr(9);
        } else if (arg == "-c" && i + 1 < argc) {
            config_file = argv[++i];
        }
    }
    
    std::cout << "=== AstDyn 3.0  OD Benchmark (Build 1130) ===" << std::endl;
    if (!filter.empty()) {
        std::cout << "Target Filter Active: [" << filter << "]" << std::endl;
    }

    // Initialize high-precision ephemeris
    const std::string de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de440s.bsp";
    auto ephem = std::make_shared<PlanetaryEphemeris>();
    try {
        ephem->setProvider(std::make_shared<DE441Provider>(de441_path));
    } catch (const std::exception& e) {
        std::cerr << "Error loading DE441: " << e.what() << std::endl;
        return 1;
    }

    // --- Base propagator settings ---
    PropagatorSettings base_settings;
    DifferentialCorrectorSettings lsq_settings;
    GaussIODSettings iod_settings;
    STMSettings stm_settings;

    if (!config_file.empty()) {
        std::cout << "Loading configuration from: " << config_file << std::endl;
        std::ifstream ifs(config_file);
        if (ifs.is_open()) {
            Configurator::loadFromStream(ifs, base_settings, lsq_settings, iod_settings, stm_settings);
        } else {
            std::cerr << "Warning: Could not open config file, using defaults." << std::endl;
        }
    } else {
        // Fallback hardcoded defaults (same as before)
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

        lsq_settings.max_iterations = 30;
        lsq_settings.outlier_sigma = 5.0;      
        lsq_settings.outlier_max_sigma = 100.0;
        lsq_settings.outlier_min_sigma = 3.0;
        lsq_settings.reject_outliers = true;
        lsq_settings.convergence_tolerance = Distance::from_au(1.0e-9);
        lsq_settings.verbose = true;
    }

    // Factory: creates a propagator.
    auto make_propagator = [&](double yarkovsky_a2 = 0.0) {
        if (!config_file.empty()) {
            std::ifstream ifs(config_file);
            nlohmann::json j;
            ifs >> j;
            auto s = base_settings;
            s.include_yarkovsky = (std::abs(yarkovsky_a2) > 0.0);
            s.yarkovsky_a2      = yarkovsky_a2;
            return Configurator::createPropagator(j, ephem);
        } else {
            auto integr = std::make_shared<OrbFitDPIntegrator>(0.1, 1e-12, 1e-6, 5.0);
            PropagatorSettings s = base_settings;
            s.include_yarkovsky = (std::abs(yarkovsky_a2) > 0.0);
            s.yarkovsky_a2      = yarkovsky_a2;
            s.integrate_in_ecliptic = true;
            return std::make_shared<Propagator>(integr, ephem, s);
        }
    };

    std::vector<Target> targets;

    // Apophis — full arc 2004-2024
    {
        // Reference: JPL Horizons heliocentric J2000 EQUATORIAL, 2020-Jan-24 00:00 TDB (MJD 58872.0)
        Eigen::Vector3d r_km(-1.541414301798325E+07,  1.394708232005084E+08,  5.145546099858084E+07);
        Eigen::Vector3d v_kms(-2.845748009596083E+01,  2.125072123515800E+00,  6.954818026860536E-02);

        Target t;
        t.name = "Apophis_99942";
        t.mpc_full = "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/apophis_mpc_full.txt";
        t.r1 = 1.1; t.r3 = 1.0; t.gooding = true;
        t.ref = CartesianStateTyped<GCRF>::from_si(time::EpochTDB::from_mjd(58872.0),
                r_km.x() * 1000.0, r_km.y() * 1000.0, r_km.z() * 1000.0,
                v_kms.x() * 1000.0, v_kms.y() * 1000.0, v_kms.z() * 1000.0);
        t.obs_min_mjd = 40000.0; t.obs_max_mjd = 70000.0;
        t.gooding_max_iter = 80; t.warm_start = true;
        t.yarkovsky_a2 = -2.89e-14;
        targets.push_back(t);
    }

    // Vesta — full arc
    {
        // Reference: JPL Horizons heliocentric J2000 EQUATORIAL, 2016-Mar-01 00:00 TDB (MJD 57448.0)
        Eigen::Vector3d r_km( 2.737877963027807E+08,  2.540494479148592E+08,  6.538351637650730E+07);
        Eigen::Vector3d v_kms(-1.159539810770006E+01,  1.221289828462967E+01,  6.383069804659519E+00);

        Target t;
        t.name = "Vesta_4";
        t.mpc_full = "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/vesta_mpc_full.txt";
        t.r1 = 2.4; t.r3 = 2.3; t.gooding = true;
        t.ref = CartesianStateTyped<GCRF>::from_si(time::EpochTDB::from_mjd(57448.0),
                r_km.x() * 1000.0, r_km.y() * 1000.0, r_km.z() * 1000.0,
                v_kms.x() * 1000.0, v_kms.y() * 1000.0, v_kms.z() * 1000.0);
        t.obs_min_mjd = 40000.0; t.obs_max_mjd = 70000.0;
        t.warm_start = true;
        targets.push_back(t);
    }

    // Phaethon — full arc
    {
        // Reference: JPL Horizons heliocentric J2000 EQUATORIAL, 2017-Sep-15 00:00 TDB (MJD 58011.0)
        Eigen::Vector3d r_km(1.405553641753041E+08, 1.956799516664952E+08, 1.777977457782163E+08);
        Eigen::Vector3d v_kms(-1.145826169909928E+01, -3.917012727493627E+00, -6.560664426417666E+00);

        Target t;
        t.name = "Phaethon_3200";
        t.mpc_full = "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/phaethon_mpc_full.txt";
        t.r1 = 1.2; t.r3 = 1.1; t.gooding = true;
        t.ref = CartesianStateTyped<GCRF>::from_si(time::EpochTDB::from_mjd(58011.0),
                r_km.x() * 1000.0, r_km.y() * 1000.0, r_km.z() * 1000.0,
                v_kms.x() * 1000.0, v_kms.y() * 1000.0, v_kms.z() * 1000.0);
        t.obs_min_mjd = 40000.0; t.obs_max_mjd = 70000.0;
        t.warm_start = true;
        t.yarkovsky_a2 = -1.1e-14;
        targets.push_back(t);
    }

    std::ofstream csv("/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/benchmarks/results/orbit_determination/results.csv");
    csv << "target,method,solver,dr_km,dv_ms,rms_arcsec,chi2,iterations,sigma_r_km,wall_s\n";

    int target_idx = 0;
    int total_targets = targets.size();

    for (const auto& tgt : targets) {
        if (!filter.empty()) {
            // Case-insensitive/fuzzy match
            std::string name_lower = tgt.name;
            std::string filter_lower = filter;
            std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);
            std::transform(filter_lower.begin(), filter_lower.end(), filter_lower.begin(), ::tolower);
            if (name_lower.find(filter_lower) == std::string::npos) continue;
        }
        target_idx++;
        std::cout << "\n[" << target_idx << "/" << total_targets << "] Processing target: " << tgt.name << "..." << std::endl;

        auto prop = make_propagator(tgt.yarkovsky_a2);
        auto full = load_mpc_obs(tgt.mpc_full);
        if (full.empty()) continue;

        // Chronological sort is essential for EKF
        std::sort(full.begin(), full.end(), [](const auto& a, const auto& b) {
            return a.time < b.time;
        });

        // Filter arc
        {
            auto it = std::remove_if(full.begin(), full.end(),
                [&tgt](const OpticalObservation& o){
                    return o.time.mjd() < tgt.obs_min_mjd || o.time.mjd() > tgt.obs_max_mjd;
                });
            full.erase(it, full.end());
            
            // Optimization: sample 1 out of 3 observations for faster benchmark turnaround
            if (full.size() > 3) {
                std::vector<OpticalObservation> sampled;
                for (size_t i = 0; i < full.size(); i += 3) {
                    sampled.push_back(full[i]);
                }
                full = std::move(sampled);
            }
            
            std::cout << "  - Total arc loaded: " << full.size() << " observations (1/3 sampling)." << std::endl;
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
                iod_state = prop->propagate_cartesian(tgt.ref, first_epoch).template cast_frame<core::ECLIPJ2000>();
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
                iod_name = "OrbFitIOD";
                GaussIODSettings s; s.use_light_time = true; s.min_separation_days = 0.5;
                auto res = OrbFitIOD(ephem, s).compute(iod_obs);
                iod_ok = res.success; 
                if (iod_ok) {
                    iod_state_gcrf = res.state;
                    iod_state = transform_state_typed<GCRF, ECLIPJ2000>(res.state);
                }
            } else {
                iod_name = (tgt.gooding_max_iter == 0) ? "LambertIOD" : "GoodingIOD";
                GoodingIOD::Settings gs; gs.max_iterations = tgt.gooding_max_iter;
                auto res = GoodingIOD(ephem, gs).compute(iod_obs[0], iod_obs[1], iod_obs[2], tgt.r1, tgt.r3);
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

            // --- PRE-CONDITIONING: Shift epoch to arc start ---
            std::cout << "  - Shifting guess epoch to arc start (" << full[0].time.mjd() << ")..." << std::flush;
            auto start_epoch = time::EpochTDB::from_mjd(full[0].time.mjd());
            auto start_state = prop->propagate_cartesian(iod_state, start_epoch);
            std::cout << " Done." << std::endl;

            // --- DIAGNOSTIC: Initial RMS check ---
            {
                auto init_calc = std::make_shared<ResidualCalculator<ECLIPJ2000>>(ephem, prop);
                auto init_res = init_calc->compute_residuals(full, start_state);
                auto init_stats = ResidualCalculator<ECLIPJ2000>::compute_statistics(init_res);
                std::cout << "  - Initial O-C RMS at start epoch: " << init_stats.rms_total.to_arcsec() << "\"" << std::endl;
                if (!init_res.empty()) {
                    std::cout << "  - First obs innovation: RA=" << init_res[0].residual_ra.to_arcsec() 
                              << "\" Dec=" << init_res[0].residual_dec.to_arcsec() << "\"" << std::endl;
                }
            }
            // ------------------------------------

            auto res_calc = std::make_shared<ResidualCalculator<GCRF>>(ephem, prop);
            auto stm_comp = std::make_shared<StateTransitionMatrix<GCRF>>(prop);
            stm_comp->apply_settings(stm_settings);

            CartesianStateTyped<GCRF> pre_ekf_state = transform_state_typed<ECLIPJ2000, GCRF>(start_state);
            pre_ekf_state.gm = physics::GravitationalParameter::from_au3_d2(constants::GMS);

            if (tgt.warm_start) {
                std::vector<OpticalObservation> short_arc;
                size_t req_obs = std::min(size_t(30), full.size());
                for (size_t i = 0; i < req_obs; ++i) {
                    short_arc.push_back(full[i]);
                }

                double timespan = short_arc.empty() ? 0 : (short_arc.back().time.mjd() - short_arc.front().time.mjd());
                std::cout << "  - [Mini-LSQ] Refining short arc (" << short_arc.size() << " obs spanning " 
                          << timespan << " days) for stable EKF init..." << std::endl;
                DifferentialCorrectorSettings short_lsq_settings;
                short_lsq_settings.verbose = true;
                short_lsq_settings.max_iterations = 30;
                short_lsq_settings.use_line_search = true;
                short_lsq_settings.outlier_sigma = 10.0;
                short_lsq_settings.outlier_max_sigma = 1000.0;
                short_lsq_settings.outlier_min_sigma = 5.0;
                auto dc_short = DifferentialCorrector<GCRF>(res_calc, stm_comp).fit(short_arc, pre_ekf_state, short_lsq_settings);
                if (dc_short.converged) {
                    pre_ekf_state = dc_short.final_state;
                    std::cout << "    ...Mini-LSQ Converged: RMS " << dc_short.statistics.rms_total.to_arcsec() << "\"" << std::endl;
                } else {
                    std::cout << "    ...Mini-LSQ Failed to converge. Fallback to raw state. (Final RMS: " << dc_short.statistics.rms_total.to_arcsec() << "\")" << std::endl;
                }
            }

            // --- STAGE 1: EKF Recursive Acquisition ---
            CartesianStateTyped<GCRF> stage2_input_state;
            std::cout << "  - [Stage 1] EKF recursive acquisition (" << full.size() << " obs)..." << std::endl;
            ExtendedKalmanFilter::Settings ekf_settings;
            ekf_settings.process_noise = Matrix6d::Identity() * 1e-12; 
            Matrix6d ekf_cov = Matrix6d::Identity() * 1e10;  // reduced covariance because state is already good
            
            ExtendedKalmanFilter ekf(prop, ekf_settings);
            CartesianStateTyped<GCRF> ekf_state = pre_ekf_state;
            
            int ekf_processed = 0;
            for (const auto& obs : full) {
                try {
                    auto res = ekf.update(ekf_state, ekf_cov, obs);
                    ekf_state = res.state;
                    ekf_cov = res.covariance;
                    ekf_processed++;
                    if (ekf_processed % 100 == 0 || ekf_processed < 20) {
                        std::cout << "    ...processed " << ekf_processed << " obs, current innovation: " 
                                  << res.innovation.norm() * constants::RAD_TO_ARCSEC << "\"" << std::endl;
                    }
                } catch (...) { break; }
            }
            stage2_input_state = ekf_state;

            // --- STAGE 2: Final LSQ refinement ---
            std::cout << "  - [Stage 2] Final LSQ refinement..." << std::endl;
            
            t_start = wall_sec();
            try {
                auto dc_res = DifferentialCorrector<GCRF>(res_calc, stm_comp).fit(full, stage2_input_state, lsq_settings);
                double dt_dc = wall_sec() - t_start;
                
                if (dc_res.converged) {
                    auto final_ref = prop->propagate_cartesian(dc_res.final_state, tgt.ref.epoch);
                    double dr = delta_r_km(final_ref, tgt.ref);
                    double dv = delta_v_ms(final_ref, tgt.ref);
                    
                    std::cout << "  → Converged [dt: " << dt_dc << "s | dr: " << dr << " km | RMS: " << dc_res.statistics.rms_total.to_arcsec() << "\"]" << std::endl;
                    write_csv_row(csv, tgt.name, "Hybrid", "EKF+LSQ", dr, dv, dc_res.statistics.rms_total.to_arcsec(), 0, dc_res.iterations, sigma_r_km(dc_res.covariance), dt_dc);
                } else {
                    std::cout << "  → FAILED to converge in Stage 2." << std::endl;
                }
            } catch (const std::exception& e) {
                std::cout << "  → ERROR: " << e.what() << std::endl;
            }
        }
    }
    csv.close();
    return 0;
}
