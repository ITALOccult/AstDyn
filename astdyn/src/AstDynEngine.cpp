/**
 * @file AstDynEngine.cpp
 * @brief Implementation of main AstDyn engine
 */

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/orbit_determination/GaussIOD.hpp"
#include "astdyn/observations/ObservatoryDatabase.hpp"
#include "astdyn/observations/MPCReader.hpp"
#include "astdyn/coordinates/KeplerianElements.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/io/AstDynConfig.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <stdexcept>

namespace astdyn {

using namespace propagation;
using namespace observations;
using namespace orbit_determination;
using namespace close_approach;

// ============================================================================
// Construction and Initialization
// ============================================================================

AstDynEngine::AstDynEngine()
    : config_()
{
    ephemeris_ = std::make_shared<ephemeris::PlanetaryEphemeris>();
    update_propagator();
    
    // CloseApproachDetector requires propagator, not just ephemeris
    ca_detector_ = std::make_unique<CloseApproachDetector>(propagator_, config_.ca_settings);
}

AstDynEngine::AstDynEngine(const AstDynConfig& config)
    : config_(config)
{
    ephemeris_ = std::make_shared<ephemeris::PlanetaryEphemeris>();
    update_propagator();
    
    // CloseApproachDetector requires propagator, not just ephemeris
    ca_detector_ = std::make_unique<CloseApproachDetector>(propagator_, config_.ca_settings);
}

void AstDynEngine::update_propagator() {
    // Create integrator based on configuration
    std::unique_ptr<Integrator> integrator;
    
    if (config_.integrator_type == "RKF78") {
        integrator = std::make_unique<RKF78Integrator>(
            config_.initial_step_size,
            config_.tolerance);
    } else if (config_.integrator_type == "RK4") {
        integrator = std::make_unique<RK4Integrator>(
            config_.initial_step_size);
    } else {
        // Default to RK4
        if (config_.verbose) {
            std::cerr << "WARNING: Unknown integrator type '" << config_.integrator_type 
                      << "', using RK4\n";
        }
        integrator = std::make_unique<RK4Integrator>(
            config_.initial_step_size);
    }
    
    
    // Update Ephemeris Provider based on config
    if (config_.ephemeris_type == "DE441" && !config_.ephemeris_file.empty()) {
        static std::string loaded_file = "";
        static bool is_loaded = false;
        
        // Avoid reloading if same file
        if (!is_loaded || loaded_file != config_.ephemeris_file) {
            try {
                if (config_.verbose) std::cout << "Loading DE441 Ephemeris: " << config_.ephemeris_file << "...\n";
                auto provider = std::make_shared<ephemeris::DE441Provider>(config_.ephemeris_file);
                ephemeris::PlanetaryEphemeris::setProvider(provider);
                loaded_file = config_.ephemeris_file;
                is_loaded = true;
                if (config_.verbose) std::cout << "DE441 Loaded successfully.\n";
            } catch (const std::exception& e) {
                std::cerr << "Error loading DE441: " << e.what() << ". Using analytical ephemeris.\n";
                // Fallback
                ephemeris::PlanetaryEphemeris::setProvider(nullptr);
                is_loaded = false;
            }
        }
    } else {
        // Default / Analytical
        // ephemeris::PlanetaryEphemeris::setProvider(nullptr); 
        // Note: we don't reset to nullptr automatically to preserve "manual" setting if user did it elsewhere? 
        // No, config should rule.
        if (config_.ephemeris_type == "Analytical") {
             ephemeris::PlanetaryEphemeris::setProvider(nullptr);
        }
    }

        // Ephemeris Provider based on config
    // ...
    
    // Sync asteroid file to settings
    config_.propagator_settings.asteroid_ephemeris_file = config_.asteroid_ephemeris_file;

    propagator_ = std::make_shared<Propagator>(
        std::move(integrator),
        ephemeris_,
        config_.propagator_settings);
}

void AstDynEngine::load_config(const std::string& oop_file) {
    config::OptionFileParser parser;
    
    if (!parser.load(oop_file)) {
        if (config_.verbose) {
            std::cerr << "WARNING: Could not load config file " << oop_file 
                      << ", using defaults\n";
        }
        return;
    }
    
    // Load integrator configuration
    config_.integrator_type = parser.getString("integrator.type", "RK4");
    config_.initial_step_size = parser.getDouble("integrator.step_size", 0.1);
    config_.tolerance = parser.getDouble("integrator.tolerance", 1e-12);
    
    // Load perturbation settings
    config_.propagator_settings.include_planets = parser.getBool("perturb.planets", true);
    config_.propagator_settings.include_asteroids = parser.getBool("perturb.asteroids", false);
    config_.propagator_settings.include_relativity = parser.getBool("perturb.relativity", false);
    config_.propagator_settings.include_moon = parser.getBool("perturb.moon", true);
    
    // Detailed planetary settings (only effective if perturb.planets is true)
    config_.propagator_settings.perturb_mercury = parser.getBool("perturb.mercury", true);
    config_.propagator_settings.perturb_venus = parser.getBool("perturb.venus", true);
    config_.propagator_settings.perturb_earth = parser.getBool("perturb.earth", true);
    config_.propagator_settings.perturb_mars = parser.getBool("perturb.mars", true);
    config_.propagator_settings.perturb_jupiter = parser.getBool("perturb.jupiter", true);
    config_.propagator_settings.perturb_saturn = parser.getBool("perturb.saturn", true);
    config_.propagator_settings.perturb_uranus = parser.getBool("perturb.uranus", true);
    config_.propagator_settings.perturb_neptune = parser.getBool("perturb.neptune", true);

    // Non-Gravitational
    config_.propagator_settings.include_yarkovsky = parser.getBool("perturb.yarkovsky", false);
    config_.propagator_settings.yarkovsky_a2 = parser.getDouble("perturb.yarkovsky_a2", 0.0);
    
    // Load Ephemeris settings
    config_.ephemeris_type = parser.getString("ephemeris.type", "Analytical");
    config_.ephemeris_file = parser.getString("ephemeris.file", "");
    config_.asteroid_ephemeris_file = parser.getString("ephemeris.asteroid_file", "");
    
    // Load EOP settings
    config_.eop_file = parser.getString("eop.file", "");
    if (!config_.eop_file.empty()) {
        if (config_.verbose) std::cout << "Loading EOP file: " << config_.eop_file << "...\n";
        if (time::load_dut1_data(config_.eop_file)) {
             if (config_.verbose) std::cout << "EOP data loaded successfully.\n";
        } else {
             if (config_.verbose) std::cerr << "WARNING: Failed to load EOP file: " << config_.eop_file << "\n";
        }
    }
    
    // Load differential correction settings
    config_.max_iterations = parser.getInt("diffcorr.max_iter", 20);
    config_.convergence_threshold = parser.getDouble("diffcorr.convergence", 1e-6);
    config_.outlier_sigma = parser.getDouble("diffcorr.outlier_threshold", 3.0);
    config_.outlier_max_sigma = parser.getDouble("diffcorr.outlier_max", 10.0);
    config_.outlier_min_sigma = parser.getDouble("diffcorr.outlier_min", 3.0);
    
    // Load residual settings
    config_.aberration_correction = parser.getBool("residuals.aberration", true);
    config_.light_time_correction = parser.getBool("residuals.light_time", true);
    
    // Load output settings
    config_.verbose = parser.getBool("verbose", true);
    
    // Update propagator with new settings
    update_propagator();
    
    // Always print configuration details
    std::cout << "   Configuration loaded:\n";
    std::cout << "   • Integrator: " << config_.integrator_type << "\n";
    std::cout << "   • Step size: " << config_.initial_step_size << " days\n";
    std::cout << "   • Tolerance: " << config_.tolerance << "\n";
    std::cout << "   • Max iterations: " << config_.max_iterations << "\n";
    if (!config_.eop_file.empty()) {
        std::cout << "   • EOP File: " << config_.eop_file << "\n";
    }
}

// ============================================================================
// Observation Management
// ============================================================================

int AstDynEngine::load_observations(const std::string& filename) {
    if (config_.verbose) {
        std::cout << "Loading observations from: " << filename << "\n";
    }
    
    observations_ = MPCReader::readFile(filename);
    
    validate_observations();
    
    if (config_.verbose) {
        std::cout << "Loaded " << observations_.size() << " observations\n";
        if (!observations_.empty()) {
            std::cout << "Time span: " 
                     << observations_.front().mjd_utc << " - "
                     << observations_.back().mjd_utc << " MJD UTC\n";
        }
    }
    
    return observations_.size();
}

void AstDynEngine::load_catalog_biases(const std::string& filename) {
    catalog_biases_.clear();
    std::ifstream file(filename);
    if (!file.is_open()) {
        if (config_.verbose) std::cerr << "Warning: Could not open catalog bias file " << filename << std::endl;
        return;
    }
    
    std::string line;
    // Skip header if present
    std::getline(file, line);
    if (line.find("Code") == std::string::npos && line.find("code") == std::string::npos) {
        // No header? Reset stream
        file.clear();
        file.seekg(0);
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string code, ra_bias_str, dec_bias_str;
        
        if (std::getline(ss, code, ',') && 
            std::getline(ss, ra_bias_str, ',') && 
            std::getline(ss, dec_bias_str, ',')) { // Simple CSV: Code,RA_Bias,Dec_Bias
            
            try {
                // Trim code
                code.erase(0, code.find_first_not_of(" \t\r\n"));
                code.erase(code.find_last_not_of(" \t\r\n") + 1);
                
                double ra_bias = std::stod(ra_bias_str); // arcsec
                double dec_bias = std::stod(dec_bias_str); // arcsec
                catalog_biases_[code] = {ra_bias, dec_bias};
            } catch (...) {
                continue;
            }
        }
    }
    
    if (config_.verbose) {
        std::cout << "Loaded biases for " << catalog_biases_.size() << " observatories from " << filename << std::endl;
    }
}

void AstDynEngine::add_observation(const observations::OpticalObservation& obs) {
    observations::OpticalObservation corrected_obs = obs;
    
    // Apply catalog bias if available
    auto it = catalog_biases_.find(obs.observatory_code);
    if (it != catalog_biases_.end()) {
        // Biases are usually defined as (Observed - True), so we subtract bias to get "True"
        // Wait, bias correction means: Corrected = Observed - Bias?
        // Usually catalogs have systematic errors. 
        // If star catalog is shifted by +0.1", obs is +0.1". 
        // to correct, we subtract 0.1".
        // So Corrected = Observed - Bias.
        // Bias units: arcsec.
        // RA/Dec units: degrees (in OpticalObservation structure).
        
        double ra_bias_deg = it->second.first / 3600.0;
        double dec_bias_deg = it->second.second / 3600.0;
        
        corrected_obs.ra -= ra_bias_deg;
        corrected_obs.dec -= dec_bias_deg;
    }

    observations_.push_back(corrected_obs);
}

void AstDynEngine::validate_observations() {
    if (observations_.empty()) {
        return;
    }
    
    // Sort by time
    std::sort(observations_.begin(), observations_.end(),
        [](const OpticalObservation& a, const OpticalObservation& b) {
            return a.mjd_utc < b.mjd_utc;
        });
    
    // Check for valid observatory codes
    ObservatoryDatabase& db = ObservatoryDatabase::getInstance();
    int unknown_count = 0;
    
    for (const auto& obs : observations_) {
        if (!db.hasObservatory(obs.observatory_code)) {
            unknown_count++;
        }
    }
    
    if (unknown_count > 0 && config_.verbose) {
        std::cout << "Warning: " << unknown_count 
                 << " observations from unknown observatories\n";
    }
}

// ============================================================================
// Orbit Determination
// ============================================================================

void AstDynEngine::set_initial_orbit(const KeplerianElements& elements) {
    current_orbit_ = elements;
    has_orbit_ = true;
    
    if (config_.verbose) {
        std::cout << "\nInitial orbit set:\n";
        std::cout << "  Epoch: " << std::fixed << std::setprecision(6) 
                 << elements.epoch_mjd_tdb << " MJD TDB\n";
        std::cout << "  a = " << elements.semi_major_axis << " AU\n";
        std::cout << "  e = " << elements.eccentricity << "\n";
        std::cout << "  i = " << (elements.inclination * constants::RAD_TO_DEG) << " deg\n";
    }
}

KeplerianElements AstDynEngine::initial_orbit_determination() {
    if (observations_.size() < 3) {
        throw std::runtime_error(
            "Insufficient observations for IOD (need at least 3, have " +
            std::to_string(observations_.size()) + ")");
    }
    
    if (config_.verbose) {
        std::cout << "\n=== Initial Orbit Determination ===\n";
        std::cout << "Using Gauss method with " << observations_.size() 
                 << " observations\n";
    }
    
    // Convert observations to optical format
    std::vector<observations::OpticalObservation> optical_obs;
    for (const auto& obs : observations_) {
        observations::OpticalObservation opt_obs;
        opt_obs.mjd_utc = obs.mjd_utc;
        opt_obs.ra = obs.ra;
        opt_obs.dec = obs.dec;
        opt_obs.sigma_ra = 1.0 / 206265.0;  // ~1 arcsec default
        opt_obs.sigma_dec = 1.0 / 206265.0;
        // Copy observatory code
        optical_obs.push_back(opt_obs);
    }
    
    // Setup Gauss IOD
    orbit_determination::GaussIODSettings gauss_settings;
    gauss_settings.verbose = config_.verbose;
    gauss_settings.max_iterations = 50;
    gauss_settings.tolerance = 1e-8;
    
    orbit_determination::GaussIOD gauss(gauss_settings);
    
    // Compute preliminary orbit
    auto result = gauss.compute(optical_obs);
    
    if (!result.success) {
        throw std::runtime_error("Gauss IOD failed: " + result.error_message);
    }
    
    if (config_.verbose) {
        result.print_summary();
    }
    
    // Store the Cartesian state and compute basic Keplerian elements
    has_orbit_ = true;
    
    // Simple conversion to Keplerian elements for initial orbit
    const auto& r = result.state.position;
    const auto& v = result.state.velocity;
    
    double r_mag = r.norm();
    double v_mag = v.norm();
    
    // Specific angular momentum
    Vector3d h = r.cross(v);
    double h_mag = h.norm();
    
    // Inclination
    double i = std::acos(h[2] / h_mag);
    
    // Node vector
    Vector3d n = Vector3d(0, 0, 1).cross(h);
    double n_mag = n.norm();
    
    // Longitude of ascending node
    double Omega = (n_mag > 1e-10) ? std::atan2(n[1], n[0]) : 0.0;
    if (Omega < 0) Omega += 2.0 * constants::PI;
    
    // Eccentricity vector
    Vector3d e_vec = ((v_mag*v_mag - constants::GMS/r_mag) * r - r.dot(v) * v) / constants::GMS;
    double e = e_vec.norm();
    
    // Argument of periapsis
    double omega = (n_mag > 1e-10 && e > 1e-10) ? std::acos(n.dot(e_vec) / (n_mag * e)) : 0.0;
    if (e_vec[2] < 0) omega = 2.0 * constants::PI - omega;
    
    // True anomaly
    double nu = (e > 1e-10) ? std::acos(e_vec.dot(r) / (e * r_mag)) : 0.0;
    if (r.dot(v) < 0) nu = 2.0 * constants::PI - nu;
    
    // Semi-major axis
    double a = 1.0 / (2.0/r_mag - v_mag*v_mag/constants::GMS);
    
    // Eccentric anomaly and mean anomaly
    double E = 2.0 * std::atan(std::sqrt((1-e)/(1+e)) * std::tan(nu/2.0));
    double M = E - e * std::sin(E);
    
    // Create Keplerian elements
    KeplerianElements kep_result;
    kep_result.semi_major_axis = a;
    kep_result.eccentricity = e;
    kep_result.inclination = i;
    kep_result.longitude_ascending_node = Omega;
    kep_result.argument_perihelion = omega;
    kep_result.mean_anomaly = M;
    kep_result.epoch_mjd_tdb = result.epoch_mjd_tdb;
    kep_result.gravitational_parameter = constants::GMS;
    
    current_orbit_ = kep_result;
    
    return kep_result;
}

OrbitDeterminationResult AstDynEngine::fit_orbit() {
    if (!has_orbit_) {
        throw std::runtime_error(
            "No initial orbit available. Call set_initial_orbit() or "
            "initial_orbit_determination() first.");
    }
    
    if (observations_.empty()) {
        throw std::runtime_error("No observations loaded");
    }
    
    if (config_.verbose) {
        std::cout << "\n=== Differential Correction ===\n";
        std::cout << "Observations: " << observations_.size() << "\n";
        std::cout << "Max iterations: " << config_.max_iterations << "\n";
        std::cout << "Convergence threshold: " << config_.convergence_threshold << "\n";
    }
    
    // Create residual calculator and STM computer (both need propagator)
    auto residual_calc = std::make_shared<orbit_determination::ResidualCalculator>(ephemeris_, propagator_);
    residual_calc->set_aberration_correction(config_.aberration_correction);
    residual_calc->set_light_time_correction(config_.light_time_correction);
    
    auto stm_computer = std::make_shared<orbit_determination::StateTransitionMatrix>(propagator_);
    
    // Create differential corrector
    auto corrector = std::make_unique<orbit_determination::DifferentialCorrector>(
        residual_calc, stm_computer);
    
    // Convert initial orbit to Cartesian at reference epoch
    propagation::CartesianElements initial_state = keplerian_to_cartesian(current_orbit_);
    
    // Setup differential corrector settings
    orbit_determination::DifferentialCorrectorSettings dc_settings;
    dc_settings.max_iterations = config_.max_iterations;
    dc_settings.convergence_tolerance = config_.convergence_threshold;
    // Relax outlier threshold to avoid rejecting all data initially
    // Standard OrbFit practice: Initial large sigma, then tighten.
    // Ideally DC handles this strategy. For now, set to 10.0 if default 3.0 is likely too tight for initial guess.
    dc_settings.outlier_sigma = config_.outlier_sigma;
    dc_settings.outlier_max_sigma = config_.outlier_max_sigma;
    dc_settings.outlier_min_sigma = config_.outlier_min_sigma;
    
    dc_settings.verbose = config_.verbose;
    
    // Perform differential correction using fit() method
    auto result_dc = corrector->fit(observations_, initial_state, dc_settings);
    
    // Convert final state back to Keplerian
    current_orbit_ = cartesian_to_keplerian(result_dc.final_state);
    
    // Build result structure
    OrbitDeterminationResult result;
    result.orbit = current_orbit_;
    result.covariance = result_dc.covariance;
    result.rms_ra = result_dc.statistics.rms_ra;  
    result.rms_dec = result_dc.statistics.rms_dec;
    result.chi_squared = result_dc.statistics.chi_squared;
    result.num_observations = result_dc.statistics.num_observations;
    result.num_rejected = result_dc.statistics.num_outliers;
    result.num_iterations = result_dc.iterations;
    result.converged = result_dc.converged;
    
    // Extract residuals (convert rad to arcsec)
    for (const auto& res : result_dc.residuals) {
        if (!res.outlier) {
            result.residuals_ra.push_back(res.residual_ra * 3600.0 * constants::RAD_TO_DEG);
            result.residuals_dec.push_back(res.residual_dec * 3600.0 * constants::RAD_TO_DEG);
        }
    }
    
    last_result_ = result;
    
    if (config_.verbose) {
        std::cout << "\n=== Final Orbit ===\n";
        print_orbit_summary();
        print_residuals_summary();
    }
    
    return result;
}

// ============================================================================
// Ephemeris Generation
// ============================================================================

std::vector<CartesianElements> AstDynEngine::compute_ephemeris(
    double start_mjd,
    double end_mjd,
    double step_days)
{
    if (!has_orbit_) {
        throw std::runtime_error("No orbit available");
    }
    
    if (config_.verbose) {
        std::cout << "\n=== Computing Ephemeris ===\n";
        std::cout << "Start: " << start_mjd << " MJD\n";
        std::cout << "End:   " << end_mjd << " MJD\n";
        std::cout << "Step:  " << step_days << " days\n";
    }
    
    // Generate epoch list
    std::vector<double> epochs;
    for (double mjd = start_mjd; mjd <= end_mjd; mjd += step_days) {
        epochs.push_back(mjd);
    }
    if (epochs.back() < end_mjd) {
        epochs.push_back(end_mjd);
    }
    
    // Convert orbit to Cartesian
    CartesianElements initial = keplerian_to_cartesian(current_orbit_);
    
    // Propagate to all epochs
    auto ephemeris = propagator_->propagate_ephemeris(initial, epochs);
    
    if (config_.verbose) {
        std::cout << "Generated " << ephemeris.size() << " ephemeris points\n";
    }
    
    return ephemeris;
}

KeplerianElements AstDynEngine::propagate_to(double target_mjd) {
    if (!has_orbit_) {
        throw std::runtime_error("No orbit available");
    }
    
    return propagator_->propagate_keplerian(current_orbit_, target_mjd);
}

// ============================================================================
// Close Approach Analysis
// ============================================================================

std::vector<CloseApproach> AstDynEngine::find_close_approaches(
    double start_mjd,
    double end_mjd)
{
    if (!has_orbit_) {
        throw std::runtime_error("No orbit available");
    }
    
    if (config_.verbose) {
        std::cout << "\n=== Close Approach Search ===\n";
        std::cout << "Searching from " << start_mjd << " to " << end_mjd << " MJD\n";
    }
    
    // Convert orbit to Cartesian
    CartesianElements initial = keplerian_to_cartesian(current_orbit_);
    
    // Search for close approaches
    auto approaches = ca_detector_->find_approaches(
        initial,
        start_mjd,
        end_mjd);
    
    if (config_.verbose) {
        std::cout << "Found " << approaches.size() << " close approaches\n";
        for (const auto& ca : approaches) {
            std::cout << "  " << static_cast<int>(ca.body) 
                     << " at " << ca.mjd_tdb << " MJD: "
                     << ca.distance << " AU ("
                     << ca.distance_in_radii(0.0001) << " radii)\n";  // Need planet radius
        }
    }
    
    return approaches;
}

double AstDynEngine::compute_moid(ephemeris::CelestialBody planet) {
    if (!has_orbit_) {
        throw std::runtime_error("No orbit available");
    }
    
    if (config_.verbose) {
        std::cout << "\n=== MOID Computation ===\n";
        std::cout << "Computing MOID with planet " 
                 << static_cast<int>(planet) << "\n";
    }
    
    // TODO: Implement MOID computation from original Fortran software
    // For now, use a simplified approach: search minimum distance over orbit period
    double moid = 999.0;  // Placeholder
    
    if (config_.verbose) {
        std::cout << "MOID = " << moid << " AU (placeholder)\n";
    }
    
    return moid;
}

// ============================================================================
// Output and Reporting
// ============================================================================

void AstDynEngine::print_orbit_summary(std::ostream& os) const {
    if (!has_orbit_) {
        os << "No orbit available\n";
        return;
    }
    
    os << std::fixed << std::setprecision(6);
    os << "Epoch:    " << current_orbit_.epoch_mjd_tdb << " MJD TDB\n";
    os << "a:        " << std::setprecision(9) << current_orbit_.semi_major_axis << " AU\n";
    os << "e:        " << std::setprecision(9) << current_orbit_.eccentricity << "\n";
    os << "i:        " << std::setprecision(6) 
       << (current_orbit_.inclination * constants::RAD_TO_DEG) << " deg\n";
    os << "Ω:        " << (current_orbit_.longitude_ascending_node * constants::RAD_TO_DEG) << " deg\n";
    os << "ω:        " << (current_orbit_.argument_perihelion * constants::RAD_TO_DEG) << " deg\n";
    os << "M:        " << (current_orbit_.mean_anomaly * constants::RAD_TO_DEG) << " deg\n";
    os << "Period:   " << std::setprecision(3) << (current_orbit_.period() / constants::DAY) << " days\n";
    os << "q:        " << std::setprecision(6) << current_orbit_.perihelion_distance() << " AU\n";
    os << "Q:        " << current_orbit_.aphelion_distance() << " AU\n";
}

void AstDynEngine::print_residuals_summary(std::ostream& os) const {
    if (last_result_.num_observations == 0) {
        os << "No residuals available\n";
        return;
    }
    
    os << "\n=== Residuals Summary ===\n";
    os << "Observations:  " << last_result_.num_observations << "\n";
    os << "Rejected:      " << last_result_.num_rejected << "\n";
    os << "RMS (RA):      " << std::fixed << std::setprecision(3) 
       << last_result_.rms_ra << " arcsec\n";
    os << "RMS (Dec):     " << last_result_.rms_dec << " arcsec\n";
    os << "Chi²:          " << std::setprecision(2) << last_result_.chi_squared << "\n";
    os << "Iterations:    " << last_result_.num_iterations << "\n";
    os << "Converged:     " << (last_result_.converged ? "Yes" : "No") << "\n";
}

void AstDynEngine::export_orbit(const std::string& filename, 
                                const std::string& format)
{
    if (!has_orbit_) {
        throw std::runtime_error("No orbit to export");
    }
    
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    if (format == "oef" || format == "OEF") {
        // AstDyn Orbit Element Format
        file << std::fixed << std::setprecision(9);
        file << "! AstDyn orbit elements\n";
        file << "! Epoch (MJD TDB): " << current_orbit_.epoch_mjd_tdb << "\n";
        file << current_orbit_.semi_major_axis << "  ! a (AU)\n";
        file << current_orbit_.eccentricity << "  ! e\n";
        file << (current_orbit_.inclination * constants::RAD_TO_DEG) << "  ! i (deg)\n";
        file << (current_orbit_.longitude_ascending_node * constants::RAD_TO_DEG) << "  ! Omega (deg)\n";
        file << (current_orbit_.argument_perihelion * constants::RAD_TO_DEG) << "  ! omega (deg)\n";
        file << (current_orbit_.mean_anomaly * constants::RAD_TO_DEG) << "  ! M (deg)\n";
    }
    else if (format == "mpc" || format == "MPC") {
        // MPC format (simplified)
        file << "! MPC orbital elements\n";
        file << "! Epoch: " << current_orbit_.epoch_mjd_tdb << " MJD\n";
        file << "! a=" << current_orbit_.semi_major_axis 
             << " e=" << current_orbit_.eccentricity 
             << " i=" << (current_orbit_.inclination * constants::RAD_TO_DEG) << "\n";
    }
    else {
        throw std::runtime_error("Unknown format: " + format);
    }
    
    if (config_.verbose) {
        std::cout << "Orbit exported to: " << filename << " (" << format << " format)\n";
    }
}

} // namespace astdyn
