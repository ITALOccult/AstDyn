#include "astdyn/orbit_determination/OrbitDetermination.hpp"
#include "astdyn/orbit_determination/ResidualCalculator.hpp"
#include "astdyn/propagation/STMPropagator.hpp"
#include "astdyn/propagation/AnalyticalJacobian.hpp"
#include "astdyn/ephemeris/VSOP87Provider.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn::orbit_determination;
using namespace astdyn::propagation;

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║            ORBIT DETERMINATION - SIMULATION TEST                 ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";

    // 1. Setup Components
    constexpr double GM_SUN = 2.959122082855911e-4;
    
    // Explicit std::function wrappers to help constructor
    STMPropagator::ForceFunction force_func = [](double t, const Eigen::Vector<double, 6>& x) {
        Eigen::Vector3d r = x.head<3>();
        Eigen::Vector3d v = x.tail<3>();
        double r_norm = r.norm();
        Eigen::Vector3d a = -GM_SUN * r / std::pow(r_norm, 3);
        
        Eigen::Vector<double, 6> dxdt;
        dxdt.head<3>() = v;
        dxdt.tail<3>() = a;
        return dxdt;
    };
    
    STMPropagator::JacobianFunction jac_func = [](double t, const Eigen::Vector<double, 6>& x) {
        return AnalyticalJacobian::two_body(x, GM_SUN);
    };
    
    auto stm_prop = std::make_unique<STMPropagator>(
        std::make_unique<RKF78Integrator>(0.1, 1e-12),
        force_func,
        jac_func
    );
    
    ResidualCalculator res_calc;
    auto ephem = std::make_shared<astdyn::ephemeris::VSOP87Provider>();
    res_calc.set_ephemeris_provider(ephem);
    
    LeastSquaresFitter fitter;
    
    // 2. Define True State (Asteroid) @ MJD 60000
    double epoch0 = 60000.0;
    Eigen::Vector<double, 6> state_true;
    // a=2.76, e=0.08, i=5 deg
    state_true << 2.55, 0.45, 0.1, -0.003, 0.009, 0.0005; 
    
    std::cout << "TRUE State:\n" << state_true.transpose() << "\n\n";

    // 3. Generate Synthetic Observations
    std::vector<Observation> synthetic_obs;
    std::vector<double> obs_times = {60010, 60020, 60030, 60040, 60050}; 
    
    std::cout << "Generating " << obs_times.size() << " synthetic observations...\n";
    
    for (double t : obs_times) {
        auto prop_res = stm_prop->propagate(state_true, epoch0, t);
        
        // Exact same calculation as ResidualCalculator to ensure consistency
        Eigen::Vector3d r_earth = ephem->getPosition(astdyn::ephemeris::CelestialBody::EARTH, t);
        Eigen::Vector3d r_ast = prop_res.state.template head<3>();
        Eigen::Vector3d r_topo_ecl = r_ast - r_earth;
        
        // Light time
        double lt = r_topo_ecl.norm() / 173.1446;
        auto prop_res_ret = stm_prop->propagate(state_true, epoch0, t - lt);
        r_ast = prop_res_ret.state.template head<3>();
        r_topo_ecl = r_ast - r_earth; 
        
        // Rotate Ecliptic -> Equatorial (+epsilon rotation for OBSERVER, -epsilon for Vector? Wait.)
        // We want coords in Equatorial.
        // Vector is in Ecliptic. Rotation to Equatorial is Rx(-eps).
        double eps = 23.4392911 * M_PI / 180.0;
        double ce = std::cos(eps), se = std::sin(eps);
        Eigen::Vector3d r_topo_eq;
        r_topo_eq.x() = r_topo_ecl.x();
        r_topo_eq.y() = r_topo_ecl.y() * ce - r_topo_ecl.z() * se;
        r_topo_eq.z() = r_topo_ecl.y() * se + r_topo_ecl.z() * ce;
        
        double ra_rad = std::atan2(r_topo_eq.y(), r_topo_eq.x());
        if (ra_rad < 0) ra_rad += 2*M_PI;
        double dec_rad = std::asin(r_topo_eq.z() / r_topo_eq.norm());
        
        Observation obs;
        obs.epoch_mjd = t;
        obs.observatory_code = "500"; // Geocentric to match our simple generation
        obs.ra_deg = ra_rad * 180.0 / M_PI;
        obs.dec_deg = dec_rad * 180.0 / M_PI;
        obs.weight = 1.0;
        obs.rejected = false;
        
        synthetic_obs.push_back(obs);
        std::cout << "  Obs @" << t << ": RA=" << obs.ra_deg << " deg, Dec=" << obs.dec_deg << " deg\n";
    }
    
    // 4. Perturb Initial State
    Eigen::Vector<double, 6> state_guess = state_true;
    state_guess(0) += 0.05; 
    state_guess(4) += 0.0001; 
    
    std::cout << "\nINITIAL GUESS State:\n" << state_guess.transpose() << "\n";
    std::cout << "Delta norm: " << (state_guess - state_true).norm() << "\n\n";
    
    // 5. Run Fit
    auto residual_func = [&](const Eigen::Vector<double, 6>& state, double epoch) {
        std::vector<ObservationResidual> residuals;
        for (const auto& obs : synthetic_obs) {
            auto res = stm_prop->propagate(state, epoch, obs.epoch_mjd);
            auto calc = res_calc.compute_residual(obs, res.state, obs.epoch_mjd);
            
            ObservationResidual or_;
            or_.epoch_mjd = calc.epoch_mjd;
            or_.ra_residual_arcsec = calc.ra_residual_arcsec;
            or_.dec_residual_arcsec = calc.dec_residual_arcsec;
            or_.weight_ra = 1.0;
            or_.weight_dec = 1.0;
            or_.rejected = false;
            residuals.push_back(or_);
        }
        return residuals;
    };
    
    auto stm_func = [&](const Eigen::Vector<double, 6>& state, double t0, double tf) {
        auto res = stm_prop->propagate(state, t0, tf);
        return std::make_pair(res.state, res.stm);
    };
    
    fitter.set_max_iterations(10);
    fitter.set_tolerance(1e-7);
    
    std::cout << "Starting Least Squares Fit...\n";
    auto result = fitter.fit(state_guess, epoch0, residual_func, stm_func);
    
    std::cout << "\nCONVERGED: " << (result.converged ? "YES" : "NO") << "\n";
    std::cout << "Iterations: " << result.num_iterations << "\n";
    std::cout << "Final RMS: " << result.rms_total_arcsec << " arcsec\n";
    std::cout << "Final State:\n" << result.state.transpose() << "\n\n";
    
    double error_norm = (result.state - state_true).norm();
    std::cout << "Error vs Truth: " << error_norm << "\n";
    
    if (result.converged && error_norm < 1e-4) {
        std::cout << "✅ TEST PASSED: Orbit recovered successfully.\n";
        return 0;
    } else {
        std::cout << "❌ TEST FAILED: Did not recover orbit correctly.\n";
        return 1;
    }
}
