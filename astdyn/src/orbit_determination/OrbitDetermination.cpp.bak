/**
 * @file OrbitDetermination.cpp
 * @brief Implementation of complete orbit determination system
 */

#include "astdyn/orbit_determination/OrbitDetermination.hpp"
#include "astdyn/propagation/AnalyticalJacobian.hpp"
#include "astdyn/ephemeris/VSOP87Provider.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace astdyn::orbit_determination {

// Gravitational parameter (Sun) in AU^3/day^2
// k = 0.01720209895 -> k^2 = 0.0002959122082855911
constexpr double GM_SUN_OD = 2.959122082855911e-4;

// 2-body force function
static Eigen::Vector<double, 6> kepler_force(double t, const Eigen::Vector<double, 6>& x) {
    Eigen::Vector3d r = x.head<3>();
    Eigen::Vector3d v = x.tail<3>();
    
    double r3 = std::pow(r.norm(), 3);
    Eigen::Vector3d a = -GM_SUN_OD * r / r3;
    
    Eigen::Vector<double, 6> dxdt;
    dxdt.head<3>() = v;
    dxdt.tail<3>() = a;
    return dxdt;
}

// Jacobian function
static Eigen::Matrix<double, 6, 6> kepler_jacobian(double t, const Eigen::Vector<double, 6>& x) {
    return propagation::AnalyticalJacobian::two_body(x, GM_SUN_OD);
}

OrbitDetermination::OrbitDetermination() {
    // Create components
    // RKF78 with strictly tight tolerance for OD
    auto integrator = std::make_unique<propagation::RKF78Integrator>(0.1, 1e-13);
    stm_propagator_ = std::make_unique<propagation::STMPropagator>(
        std::move(integrator), kepler_force, kepler_jacobian
    );
    
    residual_calculator_ = std::make_unique<ResidualCalculator>();
    
    // Use VSOP87 for Earth ephemeris
    // Note: VSOP87Provider must be fully defined
    auto ephemeris = std::make_shared<ephemeris::VSOP87Provider>();
    residual_calculator_->set_ephemeris_provider(ephemeris);
    
    fitter_ = std::make_unique<LeastSquaresFitter>();
}

void OrbitDetermination::load_elements(const std::string& eq1_file) {
    io::parsers::OrbFitEQ1Parser parser;
    auto elements = parser.parse(eq1_file);
    
    epoch_mjd_ = elements.epoch_mjd_tdb;
    state_ = elements_to_cartesian(elements);
    
    std::cout << "✓ Loaded elements from " << eq1_file << "\n";
    std::cout << "  Epoch: MJD " << epoch_mjd_ << "\n";
    std::cout << "  a = " << elements.semi_major_axis << " AU\n";
    std::cout << "  e = " << elements.eccentricity << "\n";
}

void OrbitDetermination::load_observations(const std::string& rwo_file, size_t max_obs) {
    io::parsers::AstDysRWOParser parser;
    auto opt_obs = parser.parse(rwo_file, max_obs);
    
    observations_.clear();
    observations_.reserve(opt_obs.size());
    
    for (const auto& obs : opt_obs) {
        observations_.push_back(convert_observation(obs));
    }
    
    std::cout << "✓ Loaded " << observations_.size() << " observations from " << rwo_file << "\n";
}

void OrbitDetermination::set_initial_state(const Eigen::Vector<double, 6>& state, double epoch_mjd) {
    state_ = state;
    epoch_mjd_ = epoch_mjd;
}

FitResult OrbitDetermination::fit() {
    std::cout << "\n╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║              ORBIT DETERMINATION - LEAST SQUARES FIT             ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";
    
    std::cout << "Initial state:\n";
    std::cout << "  Position: " << state_.head<3>().transpose() << " AU\n";
    std::cout << "  Velocity: " << state_.tail<3>().transpose() << " AU/day\n";
    std::cout << "  Epoch: MJD " << epoch_mjd_ << "\n\n";
    
    std::cout << "Observations: " << observations_.size() << "\n";
    std::cout << "Max iterations: " << max_iterations_ << "\n";
    std::cout << "Tolerance: " << tolerance_ << "\n\n";
    
    // Setup fitter
    fitter_->set_max_iterations(max_iterations_);
    fitter_->set_tolerance(tolerance_);
    fitter_->set_outlier_threshold(outlier_threshold_);
    
    // Residual function - Propagate to observation epoch!
    auto residual_func = [this](const Eigen::Vector<double, 6>& state, double epoch) {
        std::vector<ObservationResidual> residuals;
        for (const auto& obs : observations_) {
            // Propagate from epoch to observation time
            auto prop_result = stm_propagator_->propagate(state, epoch, obs.epoch_mjd);
            
            // Compute residual with propagated state (Topocentric + Light-time in ResidualCalculator)
            auto res_calc_result = residual_calculator_->compute_residual(
                obs, prop_result.state, obs.epoch_mjd
            );
            
            ObservationResidual obs_res;
            obs_res.epoch_mjd = res_calc_result.epoch_mjd;
            obs_res.ra_obs_deg = obs.ra_deg;
            obs_res.dec_obs_deg = obs.dec_deg;
            obs_res.ra_computed_deg = res_calc_result.ra_computed_deg;
            obs_res.dec_computed_deg = res_calc_result.dec_computed_deg;
            obs_res.ra_residual_arcsec = res_calc_result.ra_residual_arcsec;
            obs_res.dec_residual_arcsec = res_calc_result.dec_residual_arcsec;
            obs_res.weight_ra = obs.weight;
            obs_res.weight_dec = obs.weight;
            obs_res.rejected = obs.rejected;
            
            residuals.push_back(obs_res);
        }
        return residuals;
    };
    
    // STM function
    auto stm_func = [this](const Eigen::Vector<double, 6>& state, double t0, double tf) {
        auto result = stm_propagator_->propagate(state, t0, tf);
        return std::make_pair(result.state, result.stm);
    };
    
    // Perform fit
    std::cout << "Starting fit...\n";
    auto result = fitter_->fit(state_, epoch_mjd_, residual_func, stm_func);
    
    // Update state
    state_ = result.state;
    
    // Print results
    std::cout << "\n╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                         FIT RESULTS                              ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n\n";
    
    std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
    std::cout << "Iterations: " << result.num_iterations << "\n";
    std::cout << "Observations used: " << (result.num_observations - result.num_rejected) 
              << " / " << result.num_observations << "\n";
    std::cout << "Rejected: " << result.num_rejected << "\n\n";
    
    std::cout << "RMS Residuals:\n";
    std::cout << "  RA:    " << std::fixed << std::setprecision(3) 
              << result.rms_ra_arcsec << " arcsec\n";
    std::cout << "  Dec:   " << result.rms_dec_arcsec << " arcsec\n";
    std::cout << "  Total: " << result.rms_total_arcsec << " arcsec\n\n";
    
    std::cout << "Final state:\n";
    std::cout << "  Position: " << std::setprecision(8) 
              << result.state.head<3>().transpose() << " AU\n";
    std::cout << "  Velocity: " << result.state.tail<3>().transpose() << " AU/day\n\n";
    
    return result;
}

Eigen::Vector<double, 6> OrbitDetermination::elements_to_cartesian(
    const io::IOrbitParser::OrbitalElements& elem
) {
    // Keplerian to Cartesian conversion
    double a = elem.semi_major_axis;
    double e = elem.eccentricity;
    double i = elem.inclination;
    double Omega = elem.longitude_asc_node;
    double omega = elem.argument_perihelion;
    double M = elem.mean_anomaly;
    
    // Solve Kepler's equation for E
    double E = M;
    for (int iter = 0; iter < 10; ++iter) {
        E = M + e * std::sin(E);
    }
    
    // True anomaly
    double nu = 2.0 * std::atan2(
        std::sqrt(1 + e) * std::sin(E / 2.0),
        std::sqrt(1 - e) * std::cos(E / 2.0)
    );
    
    // Distance
    double r = a * (1 - e * std::cos(E));
    
    // Position and velocity in orbital plane
    double cos_nu = std::cos(nu);
    double sin_nu = std::sin(nu);
    
    Eigen::Vector3d r_orb(r * cos_nu, r * sin_nu, 0.0);
    
    double h = std::sqrt(GM_SUN_OD * a * (1 - e * e));
    Eigen::Vector3d v_orb(-GM_SUN_OD / h * sin_nu, GM_SUN_OD / h * (e + cos_nu), 0.0);
    
    // Rotation matrices
    double cos_Omega = std::cos(Omega);
    double sin_Omega = std::sin(Omega);
    double cos_omega = std::cos(omega);
    double sin_omega = std::sin(omega);
    double cos_i = std::cos(i);
    double sin_i = std::sin(i);
    
    Eigen::Matrix3d R;
    R(0, 0) = cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i;
    R(0, 1) = -cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i;
    R(0, 2) = sin_Omega * sin_i;
    R(1, 0) = sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i;
    R(1, 1) = -sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i;
    R(1, 2) = -cos_Omega * sin_i;
    R(2, 0) = sin_omega * sin_i;
    R(2, 1) = cos_omega * sin_i;
    R(2, 2) = cos_i;
    
    Eigen::Vector<double, 6> state;
    state.head<3>() = R * r_orb;
    state.tail<3>() = R * v_orb;
    
    return state;
}

Observation OrbitDetermination::convert_observation(
    const io::IObservationParser::OpticalObservation& opt_obs
) {
    Observation obs;
    obs.epoch_mjd = opt_obs.mjd_utc;
    obs.ra_deg = opt_obs.ra * 180.0 / M_PI;
    obs.dec_deg = opt_obs.dec * 180.0 / M_PI;
    obs.ra_sigma_arcsec = opt_obs.sigma_ra;
    obs.dec_sigma_arcsec = opt_obs.sigma_dec;
    obs.observatory_code = opt_obs.obs_code;
    obs.weight = 1.0 / (opt_obs.sigma_ra * opt_obs.sigma_ra);
    obs.rejected = false;
    return obs;
}

} // namespace astdyn::orbit_determination
