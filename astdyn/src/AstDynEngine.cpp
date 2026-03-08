/**
 * @file AstDynEngine.cpp
 * @brief Implementation of main AstDyn engine
 */

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/orbit_determination/StateTransitionMatrix.hpp"
#include "astdyn/orbit_determination/DifferentialCorrector.hpp"
#include "astdyn/orbit_determination/Residuals.hpp"
#include "astdyn/orbit_determination/GaussIOD.hpp"
#include "astdyn/observations/ObservatoryDatabase.hpp"
#include "astdyn/observations/MPCReader.hpp"
#include "astdyn/coordinates/KeplerianElements.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/io/AstDynConfig.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/propagation/GaussIntegrator.hpp"
#include "astdyn/propagation/RadauIntegrator.hpp"
#include "astdyn/propagation/saba4_integrator.hpp"
#include "astdyn/propagation/AASIntegrator.hpp"
#include <nlohmann/json.hpp>
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
}

AstDynEngine::AstDynEngine(const AstDynConfig& config)
    : config_(config)
{
    ephemeris_ = std::make_shared<ephemeris::PlanetaryEphemeris>();
    update_propagator();
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
    } else if (config_.integrator_type == "SABA4") {
        // Force minimum step to 0.5 for SABA4 to ensure symplectic speed-up (prevent 7s stall on 0.2s steps)
        double saba_step = std::max(0.5, config_.initial_step_size);
        integrator = std::make_unique<SABA4Integrator>(
            saba_step,
            config_.tolerance);
    } else if (config_.integrator_type == "GAUSS") {
        integrator = std::make_unique<GaussIntegrator>(
            config_.initial_step_size,
            config_.tolerance);
    } else if (config_.integrator_type == "RADAU") {
        integrator = std::make_unique<RadauIntegrator>(
            config_.initial_step_size,
            config_.tolerance);
    } else if (config_.integrator_type == "AAS") {
        // AstDyn-Adaptive Symplectic (Precision-based)
        // Ensure MU matches propagator units (AU/day)
        double mu_val = config_.propagator_settings.central_body_gm;
        integrator = std::make_unique<AASIntegrator>(
            config_.aas_precision, 
            mu_val);
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
        // Avoid reloading if same file
        if (!ephemeris_loaded_ || loaded_ephemeris_file_ != config_.ephemeris_file) {
            try {
                if (config_.verbose) std::cout << "Loading DE441 Ephemeris: " << config_.ephemeris_file << "...\n";
                auto provider = std::make_shared<ephemeris::DE441Provider>(config_.ephemeris_file);
                ephemeris::PlanetaryEphemeris::setProvider(provider);
                loaded_ephemeris_file_ = config_.ephemeris_file;
                ephemeris_loaded_ = true;
                if (config_.verbose) std::cout << "DE441 Loaded successfully.\n";
            } catch (const std::exception& e) {
                std::cerr << "Error loading DE441: " << e.what() << ". Using analytical ephemeris.\n";
                // Fallback
                ephemeris::PlanetaryEphemeris::setProvider(nullptr);
                ephemeris_loaded_ = false;
            }
        }
    } else {
        // Default / Analytical
        if (config_.ephemeris_type == "Analytical") {
             ephemeris::PlanetaryEphemeris::setProvider(nullptr);
             ephemeris_loaded_ = false;
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
    
    // Update CloseApproachDetector with the new propagator and settings
    ca_detector_ = std::make_unique<CloseApproachDetector>(propagator_, config_.ca_settings);
}

void AstDynEngine::load_config(const std::string& config_file) {
    if (config_.verbose) {
        std::cout << "Loading configuration from: " << config_file << "\n";
    }

    std::ifstream f(config_file);
    if (!f.is_open()) {
        if (config_.verbose) {
            std::cerr << "WARNING: Could not open config file " << config_file 
                      << ", using defaults\n";
        }
        return;
    }

    nlohmann::json j;
    try {
        f >> j;
    } catch (const nlohmann::json::parse_error& e) {
        if (config_.verbose) {
            std::cerr << "WARNING: JSON parse error in " << config_file << ": " 
                      << e.what() << ", using defaults\n";
        }
        return;
    }
    
    // Helper lambda for safe value retrieval
    auto get_val = [&j](const std::string& key, auto default_val) {
        return j.value(key, default_val);
    };

    // Load Integrator settings
    if (j.contains("integrator")) {
        auto& ji = j["integrator"];
        config_.integrator_type = ji.value("type", "RK4");
        config_.initial_step_size = ji.value("step_size", 0.1);
        config_.tolerance = ji.value("tolerance", 1e-12);
        config_.aas_precision = ji.value("aas_precision", 1e-4);
    } else {
        // Flat fallback or EngineConfig structure (AsteroidFitConfig style)
        config_.integrator_type = j.value("integrator_type", "RK4");
        config_.initial_step_size = j.value("initial_step_size", 0.1);
        config_.tolerance = j.value("tolerance", 1e-12);
        config_.aas_precision = j.value("aas_precision", 1e-4);
    }
    
    // Load Perturbation settings
    // Support both flattened (OptionFileParser style keys converted to json?) or structured
    // We will stick to the structure defined in AsteroidFitConfig for consistency
    
    // Check for "perturb" object first
    if (j.contains("perturb")) {
        auto& jp = j["perturb"];
        config_.propagator_settings.include_planets = jp.value("planets", true);
        config_.propagator_settings.include_asteroids = jp.value("asteroids", false);
        config_.propagator_settings.include_relativity = jp.value("relativity", false);
        config_.propagator_settings.include_moon = jp.value("moon", true);
        
        config_.propagator_settings.perturb_mercury = jp.value("mercury", true);
        config_.propagator_settings.perturb_venus = jp.value("venus", true);
        config_.propagator_settings.perturb_earth = jp.value("earth", true);
        config_.propagator_settings.perturb_mars = jp.value("mars", true);
        config_.propagator_settings.perturb_jupiter = jp.value("jupiter", true);
        config_.propagator_settings.perturb_saturn = jp.value("saturn", true);
        config_.propagator_settings.perturb_uranus = jp.value("uranus", true);
        config_.propagator_settings.perturb_neptune = jp.value("neptune", true);
        
        config_.propagator_settings.include_yarkovsky = jp.value("yarkovsky", false);
        config_.propagator_settings.yarkovsky_a2 = jp.value("yarkovsky_a2", 0.0);
    } else {
        // Fallback to top-level keys (matching AsteroidFitConfig::EngineConfig)
        config_.propagator_settings.include_planets = j.value("include_planets", true);
        config_.propagator_settings.include_asteroids = j.value("include_asteroids", false);
        config_.propagator_settings.include_relativity = j.value("include_relativity", false);
        config_.propagator_settings.include_moon = j.value("include_moon", true);
        
        config_.propagator_settings.perturb_mercury = j.value("perturb_mercury", true);
        config_.propagator_settings.perturb_venus = j.value("perturb_venus", true);
        config_.propagator_settings.perturb_earth = j.value("perturb_earth", true);
        config_.propagator_settings.perturb_mars = j.value("perturb_mars", true);
        config_.propagator_settings.perturb_jupiter = j.value("perturb_jupiter", true);
        config_.propagator_settings.perturb_saturn = j.value("perturb_saturn", true);
        config_.propagator_settings.perturb_uranus = j.value("perturb_uranus", true);
        config_.propagator_settings.perturb_neptune = j.value("perturb_neptune", true);
        
        config_.propagator_settings.include_yarkovsky = j.value("include_yarkovsky", false);
        config_.propagator_settings.yarkovsky_a2 = j.value("yarkovsky_a2", 0.0);
    }

    // Ephemeris Settings
    if (j.contains("ephemeris")) {
        auto& je = j["ephemeris"];
        config_.ephemeris_type = je.value("type", "Analytical");
        config_.ephemeris_file = je.value("file", "");
        config_.asteroid_ephemeris_file = je.value("asteroid_file", "");
    } else {
        // Top level fallback
        // infer type from file existence or specific key
        config_.ephemeris_file = j.value("ephemeris_file", "");
        config_.asteroid_ephemeris_file = j.value("asteroid_ephemeris_file", "");
        if (!config_.ephemeris_file.empty()) {
             config_.ephemeris_type = "DE441";
        }
    }
    
    // EOP
    if (j.contains("eop")) {
        config_.eop_file = j["eop"].value("file", "");
    } else {
        config_.eop_file = j.value("eop_file", "");
    }
    
    if (!config_.eop_file.empty()) {
        if (config_.verbose) std::cout << "Loading EOP file: " << config_.eop_file << "...\n";
        if (time::load_dut1_data(config_.eop_file)) {
             if (config_.verbose) std::cout << "EOP data loaded successfully.\n";
        } else {
             if (config_.verbose) std::cerr << "WARNING: Failed to load EOP file: " << config_.eop_file << "\n";
        }
    }

    // Catalog Bias
    if (j.contains("catalog")) {
        config_.catalog_bias_file = j["catalog"].value("bias_file", "");
    } else {
        config_.catalog_bias_file = j.value("catalog_bias_file", "");
    }
    
    if (!config_.catalog_bias_file.empty()) {
        load_catalog_biases(config_.catalog_bias_file);
    }
    
    // Differential Correction
    if (j.contains("diffcorr")) {
        auto& jd = j["diffcorr"];
        config_.max_iterations = jd.value("max_iter", 20);
        config_.convergence_threshold = jd.value("convergence", 1e-6);
        config_.outlier_sigma = jd.value("outlier_threshold", 3.0);
        config_.outlier_max_sigma = jd.value("outlier_max", 10.0);
        config_.outlier_min_sigma = jd.value("outlier_min", 3.0);
    } else {
        config_.max_iterations = j.value("max_iterations", 20); // Note key difference potentially
        // Fallback for flat structure
        config_.outlier_max_sigma = j.value("outlier_max_sigma", 10.0);
        config_.outlier_min_sigma = j.value("outlier_min_sigma", 3.0);
    }
    
    // Residuals settings
    if (j.contains("residuals")) {
        auto& jr = j["residuals"];
        config_.aberration_correction = jr.value("aberration", true);
        config_.light_time_correction = jr.value("light_time", true);
    }
    
    // Close Approach settings
    if (j.contains("close_approach")) {
        auto& jc = j["close_approach"];
        config_.ca_settings.detection_distance = jc.value("detection_distance", 0.05 * constants::AU);
        config_.ca_settings.min_distance = jc.value("min_distance", 1e-6 * constants::AU);
        config_.ca_settings.compute_b_plane = jc.value("compute_b_plane", true);
        config_.ca_settings.refine_time = jc.value("refine_time", true);
        config_.ca_settings.time_tolerance = jc.value("time_tolerance", 1e-6);
        config_.ca_settings.max_refinement_iter = jc.value("max_iter", 10);
    } else {
        // Flat fallback
        config_.ca_settings.detection_distance = j.value("ca_detection_distance", 0.05 * constants::AU);
        config_.ca_settings.min_distance = j.value("ca_min_distance", 1e-6 * constants::AU);
        config_.ca_settings.compute_b_plane = j.value("ca_compute_b_plane", true);
        config_.ca_settings.refine_time = j.value("ca_refine_time", true);
        config_.ca_settings.time_tolerance = j.value("ca_time_tolerance", 1e-6);
        config_.ca_settings.max_refinement_iter = j.value("ca_max_iter", 10);
    }

    // General
    config_.verbose = j.value("verbose", true);

    // Update propagator with new settings
    update_propagator();
    
    // Always print configuration details
    std::cout << "   Configuration loaded (JSON):\n";
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
                     << observations_.front().time.mjd() << " - "
                     << observations_.back().time.mjd() << " MJD UTC\n";
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
        // Corrected = Observed - Bias.
        // Bias units: arcsec.
        // RA/Dec units: radians (in OpticalObservation structure).
        
        double ra_bias_rad = it->second.first * constants::ARCSEC_TO_RAD;
        double dec_bias_rad = it->second.second * constants::ARCSEC_TO_RAD;
        
        corrected_obs.ra -= ra_bias_rad;
        corrected_obs.dec -= dec_bias_rad;
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
            return a.time.mjd() < b.time.mjd();
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

void AstDynEngine::set_initial_orbit_ecl(const physics::KeplerianStateTyped<core::ECLIPJ2000>& elements) {
    current_orbit_ = elements;
    has_orbit_ = true;
    
    if (config_.verbose) {
        std::cout << "\nInitial orbit set:\n";
        std::cout << "  Epoch: " << std::fixed << std::setprecision(6) 
                 << elements.epoch.mjd() << " MJD TDB\n";
        std::cout << "  a = " << elements.a.to_au() << " AU\n";
        std::cout << "  e = " << elements.e << "\n";
        std::cout << "  i = " << elements.i.to_deg() << " deg\n";
    }
}

physics::KeplerianStateTyped<core::ECLIPJ2000> AstDynEngine::initial_orbit_determination() {
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
        opt_obs.time = obs.time;
        opt_obs.ra = obs.ra;
        opt_obs.dec = obs.dec;
        opt_obs.sigma_ra = obs.sigma_ra;
        opt_obs.sigma_dec = obs.sigma_dec;
        // Copy observatory code
        opt_obs.observatory_code = obs.observatory_code;
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
    
    // 4. Bring into AstDyn 3.0 Type Safety (Ecliptic J2000)
    // Convert GCRF result to Ecliptic J2000
    auto cart_gcrf = result.state;
    auto pos_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(cart_gcrf.position, cart_gcrf.epoch);
    auto vel_ecl = coordinates::ReferenceFrame::transform_vel<core::GCRF, core::ECLIPJ2000>(cart_gcrf.position, cart_gcrf.velocity, cart_gcrf.epoch);
    
    auto cart_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>(
        cart_gcrf.epoch, pos_ecl, vel_ecl, cart_gcrf.gm
    );
    
    current_orbit_ = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(cart_ecl);
    has_orbit_ = true;
    return current_orbit_;
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
    }
    
    bool in_ecl = propagator_->settings().integrate_in_ecliptic;
    
    // Setup differential corrector settings
    orbit_determination::DifferentialCorrectorSettings dc_settings;
    dc_settings.max_iterations = config_.max_iterations;
    dc_settings.convergence_tolerance = config_.convergence_threshold;
    dc_settings.outlier_sigma = config_.outlier_sigma;
    dc_settings.outlier_max_sigma = std::max(100.0, config_.outlier_sigma * 10.0); 
    dc_settings.outlier_min_sigma = config_.outlier_sigma;
    dc_settings.reject_outliers = true;
    dc_settings.verbose = config_.verbose;

    OrbitDeterminationResult result;
    
    if (in_ecl) {
        auto residual_calc = std::make_shared<orbit_determination::ResidualCalculator<core::ECLIPJ2000>>(ephemeris_, propagator_);
        residual_calc->set_aberration_correction(config_.aberration_correction);
        residual_calc->set_light_time_correction(config_.light_time_correction);
        
        auto stm_computer = std::make_shared<orbit_determination::StateTransitionMatrix<core::ECLIPJ2000>>(propagator_);
        auto corrector = std::make_unique<orbit_determination::DifferentialCorrector<core::ECLIPJ2000>>(residual_calc, stm_computer);
        
        auto initial_state_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(current_orbit_);
        auto result_dc = corrector->fit(observations_, initial_state_ecl, dc_settings);
        
        current_orbit_ = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(result_dc.final_state);
        
        result.orbit = current_orbit_;
        result.converged = result_dc.converged;
        result.rms_ra = result_dc.statistics.rms_ra;
        result.rms_dec = result_dc.statistics.rms_dec;
        result.rms_total_arcsec = result_dc.statistics.rms_total;
        result.num_observations = result_dc.statistics.num_observations;
        result.num_iterations = result_dc.iterations;
        result.covariance = result_dc.covariance;
        result.chi_squared = result_dc.statistics.chi_squared;
        result.num_rejected = result_dc.statistics.num_outliers;
        
        for (const auto& res : result_dc.residuals) {
            if (!res.outlier) {
                result.residuals_ra.push_back(res.residual_ra * 3600.0 * constants::RAD_TO_DEG);
                result.residuals_dec.push_back(res.residual_dec * 3600.0 * constants::RAD_TO_DEG);
            }
        }
    } else {
        auto residual_calc = std::make_shared<orbit_determination::ResidualCalculator<core::GCRF>>(ephemeris_, propagator_);
        residual_calc->set_aberration_correction(config_.aberration_correction);
        residual_calc->set_light_time_correction(config_.light_time_correction);
        
        auto stm_computer = std::make_shared<orbit_determination::StateTransitionMatrix<core::GCRF>>(propagator_);
        auto corrector = std::make_unique<orbit_determination::DifferentialCorrector<core::GCRF>>(residual_calc, stm_computer);
        
        auto initial_state_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(current_orbit_);
        auto pos_g = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(initial_state_ecl.position, initial_state_ecl.epoch);
        auto vel_g = coordinates::ReferenceFrame::transform_vel<core::ECLIPJ2000, core::GCRF>(initial_state_ecl.position, initial_state_ecl.velocity, initial_state_ecl.epoch);
        auto initial_gcrf = physics::CartesianStateTyped<core::GCRF>(initial_state_ecl.epoch, pos_g, vel_g, initial_state_ecl.gm);
        
        auto result_dc = corrector->fit(observations_, initial_gcrf, dc_settings);
        
        auto final_gcrf = result_dc.final_state;
        auto pos_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(final_gcrf.position, final_gcrf.epoch);
        auto vel_ecl = coordinates::ReferenceFrame::transform_vel<core::GCRF, core::ECLIPJ2000>(final_gcrf.position, final_gcrf.velocity, final_gcrf.epoch);
        auto final_ecl = physics::CartesianStateTyped<core::ECLIPJ2000>(final_gcrf.epoch, pos_ecl, vel_ecl, final_gcrf.gm);
        
        current_orbit_ = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(final_ecl);
        
        result.orbit = current_orbit_;
        result.converged = result_dc.converged;
        result.rms_ra = result_dc.statistics.rms_ra;
        result.rms_dec = result_dc.statistics.rms_dec;
        result.rms_total_arcsec = result_dc.statistics.rms_total;
        result.num_observations = result_dc.statistics.num_observations;
        result.num_iterations = result_dc.iterations;
        result.covariance = result_dc.covariance;
        result.chi_squared = result_dc.statistics.chi_squared;
        result.num_rejected = result_dc.statistics.num_outliers;
        
        for (const auto& res : result_dc.residuals) {
            if (!res.outlier) {
                result.residuals_ra.push_back(res.residual_ra * 3600.0 * constants::RAD_TO_DEG);
                result.residuals_dec.push_back(res.residual_dec * 3600.0 * constants::RAD_TO_DEG);
            }
        }
    }
    
    last_result_ = result;
    has_orbit_ = true;
    
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

std::vector<physics::CartesianStateTyped<core::GCRF>> AstDynEngine::compute_ephemeris(
    time::EpochTDB start_time,
    time::EpochTDB end_time,
    double step_days)
{
    if (!has_orbit_) {
        throw std::runtime_error("No orbit available");
    }
    
    if (config_.verbose) {
        std::cout << "\n=== Computing Ephemeris ===\n";
        std::cout << "Start: " << start_time.mjd() << " MJD\n";
        std::cout << "End:   " << end_time.mjd() << " MJD\n";
        std::cout << "Step:  " << step_days << " days\n";
    }
    
    // Generate epoch list
    std::vector<time::EpochTDB> epochs;
    if (step_days <= 0) {
        epochs.push_back(start_time);
    } else {
        for (double mjd = start_time.mjd(); mjd <= end_time.mjd(); mjd += step_days) {
            epochs.push_back(time::EpochTDB::from_mjd(mjd));
        }
        if (epochs.empty() || std::abs(epochs.back().mjd() - end_time.mjd()) > 1e-8) {
            epochs.push_back(end_time);
        }
    }
    
    bool in_ecl = propagator_->settings().integrate_in_ecliptic;
    std::vector<physics::CartesianStateTyped<core::GCRF>> results_gcrf;
    results_gcrf.reserve(epochs.size());

    if (in_ecl) {
        // Correct bridge: Start and stay in Ecliptic J2000
        auto initial_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(current_orbit_);
        auto ephemeris_ecl = propagator_->propagate_ephemeris(initial_ecl, epochs);
        
        // Final bridge: Convert results to GCRF for return type
        for (const auto& state : ephemeris_ecl) {
            auto pos_g = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(state.position, state.epoch);
            auto vel_g = coordinates::ReferenceFrame::transform_vel<core::ECLIPJ2000, core::GCRF>(state.position, state.velocity, state.epoch);
            results_gcrf.push_back(physics::CartesianStateTyped<core::GCRF>(state.epoch, pos_g, vel_g, state.gm));
        }
    } else {
        // Classical GCRF propagation
        auto initial_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(current_orbit_);
        auto pos_g = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(initial_ecl.position, initial_ecl.epoch);
        auto vel_g = coordinates::ReferenceFrame::transform_vel<core::ECLIPJ2000, core::GCRF>(initial_ecl.position, initial_ecl.velocity, initial_ecl.epoch);
        auto initial_gcrf = physics::CartesianStateTyped<core::GCRF>(initial_ecl.epoch, pos_g, vel_g, initial_ecl.gm);
        
        results_gcrf = propagator_->propagate_ephemeris(initial_gcrf, epochs);
    }
    
    if (config_.verbose) {
        std::cout << "Generated " << results_gcrf.size() << " ephemeris points\n";
    }
    
    return results_gcrf;
}

physics::KeplerianStateTyped<core::ECLIPJ2000> AstDynEngine::propagate_to(time::EpochTDB target_time) {
    if (!has_orbit_) {
        throw std::runtime_error("No orbit available");
    }
    
    return propagator_->propagate_keplerian(current_orbit_, target_time);
}

// ============================================================================
// Close Approach Analysis
// ============================================================================

std::vector<CloseApproach> AstDynEngine::find_close_approaches(
    time::EpochTDB start_time,
    time::EpochTDB end_time)
{
    if (!has_orbit_) {
        throw std::runtime_error("No orbit available");
    }
    
    if (config_.verbose) {
        std::cout << "\n=== Close Approach Search ===\n";
        std::cout << "Searching from " << start_time.mjd() << " to " << end_time.mjd() << " MJD\n";
    }
    
    // Convert orbit to Cartesian J2000 Equatorial
    auto initial_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(current_orbit_);
    auto pos_gcrf = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(initial_ecl.position, initial_ecl.epoch);
    auto vel_gcrf = coordinates::ReferenceFrame::transform_vel<core::ECLIPJ2000, core::GCRF>(initial_ecl.position, initial_ecl.velocity, initial_ecl.epoch);
    
    auto initial = physics::CartesianStateTyped<core::GCRF>(
        initial_ecl.epoch, pos_gcrf, vel_gcrf, initial_ecl.gm
    );
    
    // Search for close approaches
    auto approaches = ca_detector_->find_approaches(
        initial,
        start_time,
        end_time);
    
    if (config_.verbose) {
        std::cout << "Found " << approaches.size() << " close approaches\n";
        for (const auto& ca : approaches) {
            std::cout << "  " << static_cast<int>(ca.body) 
                     << " at " << ca.time.mjd() << " MJD: "
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

ApparentPlace AstDynEngine::compute_asteroid_apparent_place(time::EpochTDB t_occult, const std::string& observatory_code) {
    if (!has_orbit_) throw std::runtime_error("Orbit must be fitted or set before computing apparent place.");
    
    // 1. Propagate orbit to occultation time
    auto state_ecl = propagator_->propagate_cartesian(propagation::keplerian_to_cartesian(current_orbit_), t_occult);
    
    // 2. Setup residual calculator for precise geometry
    auto calculator = orbit_determination::ResidualCalculator<core::ECLIPJ2000>(ephemeris_, propagator_);
    calculator.set_aberration_correction(config_.aberration_correction);
    calculator.set_light_time_correction(config_.light_time_correction);
    
    observations::OpticalObservation obs;
    obs.time = time::to_utc(t_occult);
    obs.observatory_code = observatory_code;
    
    auto res_opt = calculator.compute_residual(obs, state_ecl);
    if (!res_opt) throw std::runtime_error("Failed to compute asteroid apparent place.");
    
    ApparentPlace result;
    result.ra = res_opt->computed_ra;
    result.dec = res_opt->computed_dec;
    result.distance = res_opt->range;
    result.light_time = res_opt->range / (constants::C_LIGHT * 86400.0);
    
    return result;
}

ApparentPlace AstDynEngine::compute_star_apparent_place(
    double ra_j2000, double dec_j2000,
    double pm_ra, double pm_dec,
    time::EpochTDB t_occult, const std::string& observatory_code) {
    
    // 1. Proper motion to occultation epoch
    double dt_yr = (t_occult.mjd() - constants::MJD2000) / constants::DAYS_PER_YEAR; // Julian years
    double ra_epoch = ra_j2000 + pm_ra * dt_yr;
    double dec_epoch = dec_j2000 + pm_dec * dt_yr;
    
    // 2. Build direction vector in GCRF (Ecliptic logic handles GCRF transformation)
    Eigen::Vector3d u_gcrf;
    u_gcrf << std::cos(dec_epoch) * std::cos(ra_epoch),
              std::cos(dec_epoch) * std::sin(ra_epoch),
              std::sin(dec_epoch);
               
    // 3. Setup calculator
    auto calculator = orbit_determination::ResidualCalculator<core::ECLIPJ2000>(ephemeris_);
    calculator.set_aberration_correction(config_.aberration_correction);
    
    observations::OpticalObservation obs;
    obs.time = time::to_utc(t_occult);
    obs.observatory_code = observatory_code;
    
    // 4. Place star dummy state far away (100,000 AU) in GCRF
    auto obs_pos_opt = calculator.get_observer_position(obs);
    if (!obs_pos_opt) throw std::runtime_error("Failed to get observer state for star apparent place.");
    
    math::Vector3<core::GCRF, physics::Distance> u_vec = math::Vector3<core::GCRF, physics::Distance>::from_si(u_gcrf[0], u_gcrf[1], u_gcrf[2]);
    math::Vector3<core::GCRF, physics::Distance> rho_star = u_vec * (100000.0 * constants::AU * 1000.0);
    
    physics::CartesianStateTyped<core::GCRF> star_state;
    star_state.epoch = t_occult;
    star_state.position = *obs_pos_opt + rho_star;
    star_state.velocity = math::Vector3<core::GCRF, physics::Velocity>::from_si(0.0, 0.0, 0.0);
    star_state.gm = physics::GravitationalParameter::from_si(1.32712440018e20); // Dummy SUN GM
    
    // 5. Compute 'residual' to get topocentric apparent RA/Dec (includes aberration, light bending)
    auto pos_ecl = coordinates::ReferenceFrame::transform_pos<core::GCRF, core::ECLIPJ2000>(star_state.position, star_state.epoch);
    auto vel_ecl = coordinates::ReferenceFrame::transform_vel<core::GCRF, core::ECLIPJ2000>(star_state.position, star_state.velocity, star_state.epoch);
    physics::CartesianStateTyped<core::ECLIPJ2000> star_state_ecl(star_state.epoch, pos_ecl, vel_ecl, star_state.gm);
    
    auto star_res_opt = calculator.compute_residual(obs, star_state_ecl);
    
    if (!star_res_opt) throw std::runtime_error("Failed to compute star apparent place.");
    
    ApparentPlace result;
    result.ra = star_res_opt->computed_ra;
    result.dec = star_res_opt->computed_dec;
    result.distance = star_res_opt->range;
    result.light_time = 0.0;
    
    return result;
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
    os << "Epoch:    " << current_orbit_.epoch.mjd() << " MJD TDB\n";
    os << "a:        " << std::setprecision(9) << current_orbit_.a.to_au() << " AU\n";
    os << "e:        " << std::setprecision(9) << current_orbit_.e << "\n";
    os << "i:        " << std::setprecision(6) 
       << current_orbit_.i.to_deg() << " deg\n";
    os << "Ω:        " << current_orbit_.node.to_deg() << " deg\n";
    os << "ω:        " << current_orbit_.omega.to_deg() << " deg\n";
    os << "M:        " << current_orbit_.M.to_deg() << " deg\n";
    os << "Period:   " << std::setprecision(3) << current_orbit_.period_days() << " days\n";
    os << "q:        " << std::setprecision(6) << current_orbit_.perihelion_au() << " AU\n";
    os << "Q:        " << current_orbit_.aphelion_au() << " AU\n";
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
        file << "! Epoch (MJD TDB): " << current_orbit_.epoch.mjd() << "\n";
        file << current_orbit_.a.to_au() << "  ! a (AU)\n";
        file << current_orbit_.e << "  ! e\n";
        file << current_orbit_.i.to_deg() << "  ! i (deg)\n";
        file << current_orbit_.node.to_deg() << "  ! Omega (deg)\n";
        file << current_orbit_.omega.to_deg() << "  ! omega (deg)\n";
        file << current_orbit_.M.to_deg() << "  ! M (deg)\n";
    }
    else if (format == "mpc" || format == "MPC") {
        // MPC format (simplified)
        file << "! MPC orbital elements\n";
        file << "! Epoch: " << current_orbit_.epoch.mjd() << " MJD\n";
        file << "! a=" << current_orbit_.a.to_au() 
             << " e=" << current_orbit_.e 
             << " i=" << current_orbit_.i.to_deg() << "\n";
    }
    else {
        throw std::runtime_error("Unknown format: " + format);
    }
    
    if (config_.verbose) {
        std::cout << "Orbit exported to: " << filename << " (" << format << " format)\n";
    }
}

} // namespace astdyn
