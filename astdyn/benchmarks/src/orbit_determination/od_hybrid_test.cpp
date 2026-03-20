/**
 * @file od_hybrid_test.cpp
 * @brief Test for OrbFit IOD + AstDyn LSQ/EKF
 */

#include <astdyn/AstDyn.hpp>
#include <astdyn/orbit_determination/GaussIOD.hpp>
#include <astdyn/orbit_determination/OrbFitIOD.hpp>
#include <astdyn/orbit_determination/DifferentialCorrector.hpp>
#include <astdyn/orbit_determination/ExtendedKalmanFilter.hpp>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <astdyn/orbit_determination/StateTransitionMatrix.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/io/MPCParser.hpp>
#include <astdyn/ephemeris/DE441Provider.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/core/Constants.hpp>
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

int main() {
    std::cout << "=== AstDyn Hybrid OD Test (OrbFit IOD Seed) ===" << std::endl;

    // 1. Setup Environment
    const std::string de441_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de440s.bsp";
    auto ephem = std::make_shared<PlanetaryEphemeris>();
    ephem->setProvider(std::make_shared<DE441Provider>(de441_path));

    // Classical AstDyn Integrator: RKF78
    auto integr = std::make_shared<RKF78Integrator>(0.1, 1e-12);
    PropagatorSettings settings;
    settings.include_planets = true;
    settings.perturb_mercury = true;
    settings.perturb_venus = true;
    settings.perturb_earth = true;
    settings.perturb_mars = true;
    settings.perturb_jupiter = true;
    settings.perturb_saturn = true;
    settings.perturb_uranus = true;
    settings.perturb_neptune = true;
    settings.include_relativity = true;
    settings.include_moon = true;
    settings.integrate_in_ecliptic = true;

    auto prop = std::make_shared<Propagator>(integr, ephem, settings);

    // 2. Load Data (Apophis)
    std::string mpc_path = "/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/data/apophis_mpc_full.txt";
    auto full_data = load_mpc_obs(mpc_path);
    if (full_data.empty()) {
        std::cerr << "Error: Could not load Apophis data." << std::endl;
        return 1;
    }
    std::sort(full_data.begin(), full_data.end(), [](const auto& a, const auto& b) { return a.time < b.time; });
    std::cout << "Loaded " << full_data.size() << " total observations." << std::endl;

    // 3. Stage 0: OrbFit IOD
    std::cout << "Step 1: Running OrbFit IOD (selecting 60-day mini-arc from unsampled data)..." << std::endl;
    std::vector<OpticalObservation> iod_mini_arc;
    double t0 = full_data.front().time.mjd();
    for (const auto& o : full_data) {
        if (o.time.mjd() - t0 <= 90.0) { // Increased to 90 days for better stability
            iod_mini_arc.push_back(o);
        } else {
            break;
        }
    }
    
    if (iod_mini_arc.size() < 3) {
        std::cerr << "Error: Not enough observations for IOD." << std::endl;
        return 1;
    }

    GaussIODSettings iod_settings;
    iod_settings.use_light_time = true;
    iod_settings.verbose = true;
    
    // Explicitly pick three observations with good separation
    std::vector<OpticalObservation> selected_iod;
    selected_iod.push_back(full_data[0]);
    
    // Find second obs ~15 days later
    for (const auto& o : full_data) {
        if (o.time.mjd() - selected_iod[0].time.mjd() >= 15.0) {
            selected_iod.push_back(o);
            break;
        }
    }
    // Find third obs ~15 days AFTER second
    if (selected_iod.size() == 2) {
        for (const auto& o : full_data) {
            if (o.time.mjd() - selected_iod[1].time.mjd() >= 15.0) {
                selected_iod.push_back(o);
                break;
            }
        }
    }

    if (selected_iod.size() < 3) {
        std::cerr << "Error: Not enough observations for 30-day IOD baseline." << std::endl;
        return 1;
    }

    auto iod_res = OrbFitIOD(iod_settings).compute(selected_iod);
    
    if (!iod_res.success) {
        std::cerr << " OrbFit IOD Failed: " << iod_res.error_message << std::endl;
        return 1;
    }
    std::cout << " Success. Epoch: " << iod_res.epoch.mjd() << std::endl;
    auto iod_state = iod_res.state.template cast_frame<GCRF>();
    std::cout << " IOD State (AU): " << std::endl;
    std::cout << iod_state.to_eigen_au_aud() << std::endl;

    // 4. Sample observations for refinement
    std::vector<OpticalObservation> full;
    for (size_t i = 0; i < full_data.size(); i += 50) full.push_back(full_data[i]);
    std::cout << "Refining on " << full.size() << " observations (sampled 1/50)." << std::endl;

    // 5. Classical AstDyn Pipeline (EKF + LSQ)
    std::cout << "Step 2: Classical AstDyn Refinement (Analytical Jacobian, RKF78)..." << std::endl;
    
    auto res_calc = std::make_shared<ResidualCalculator<GCRF>>(ephem, prop);
    auto stm_comp = std::make_shared<StateTransitionMatrix<GCRF>>(prop);
    stm_comp->set_use_numerical_jacobian(false); // ANALYTICAL (Classico)

    // Stage 1: EKF
    std::cout << "  - Initializing EKF at first observation epoch..." << std::endl;
    ExtendedKalmanFilter::Settings ekf_settings;
    ExtendedKalmanFilter ekf(prop, ekf_settings);
    
    // Propagate IOD state to the first observation of the full set
    auto ekf_state = prop->propagate_cartesian(iod_state, time::to_tdb(full.front().time));
    Matrix6d ekf_cov = Matrix6d::Identity() * 1e12; // Broad uncertainty

    int count = 0;
    int skipped = 0;
    for (const auto& obs : full) {
        try {
            auto res = ekf.update(ekf_state, ekf_cov, obs);
            if (std::isnan(res.state.position.to_eigen_si()[0])) {
                skipped++;
                continue;
            }
            ekf_state = res.state;
            ekf_cov = res.covariance;
            if (++count % 50 == 0) std::cout << "  - Processed " << count << " obs..." << std::endl;
        } catch (const std::exception& e) { 
            skipped++;
            continue; 
        }
    }
    std::cout << "  - EKF processing complete. Processed: " << count << " Skipped: " << skipped << std::endl;

    // Stage 2: LSQ
    DifferentialCorrectorSettings dc_settings;
    dc_settings.verbose = true;
    dc_settings.max_iterations = 20;
    dc_settings.reject_outliers = true;
    
    auto dc_res = DifferentialCorrector<GCRF>(res_calc, stm_comp).fit(full, ekf_state, dc_settings);
    
    // 5. Output Results
    std::cout << "\n=== FINAL RESULTS (Hybrid) ===" << std::endl;
    dc_res.print();
    
    std::ofstream csv("/Users/michelebigi/Documents/Develop/ASTDYN/IOccultLibrary/astdyn/benchmarks/results/orbit_determination/results_hybrid.csv");
    csv << "target,method,iod,rms_arcsec,iterations,converged\n";
    csv << "Apophis,Hybrid,OrbFitIOD," << dc_res.statistics.rms_total.to_arcsec() << "," << dc_res.iterations << "," << (dc_res.converged ? "YES" : "NO") << "\n";
    csv.close();
    
    std::cout << "Results saved to results_hybrid.csv" << std::endl;

    return 0;
}
