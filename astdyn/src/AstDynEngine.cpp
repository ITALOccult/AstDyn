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
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/io/AstDynConfig.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/propagation/GaussIntegrator.hpp"
#include "astdyn/propagation/RadauIntegrator.hpp"
#include "astdyn/propagation/saba4_integrator.hpp"
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
    } else {
        // Flat fallback or EngineConfig structure (AsteroidFitConfig style)
        config_.integrator_type = j.value("integrator_type", "RK4");
        config_.initial_step_size = j.value("initial_step_size", 0.1);
        config_.tolerance = j.value("tolerance", 1e-12);
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

void AstDynEngine::set_initial_orbit(const physics::KeplerianStateTyped<core::ECLIPJ2000>& elements) {
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
        opt_obs.sigma_ra = 1.0 / 206265.0;  // ~1 arcsec default
        opt_obs.sigma_dec = 1.0 / 206265.0;
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
    
    // 4. Bring into AstDyn 3.0 Type Safety
    auto cart_typed = result.state; // GaussIODResult::state is already CartesianStateTyped<core::GCRF>
    
    // Convert to Keplerian (Ecliptic J2000)
    // Bridge back to un-typed for returning to Keplerian
    CartesianElements cart_old;
    cart_old.epoch = cart_typed.epoch;
    
    // Rotate back to Ecliptic
    auto pos_eq = cart_typed.position.to_eigen_si();
    auto vel_eq = cart_typed.velocity.to_eigen_si();
    auto mat_ecl_to_eq = coordinates::ReferenceFrame::ecliptic_to_j2000();
    auto pos_ecl = mat_ecl_to_eq.transpose() * pos_eq;
    auto vel_ecl = mat_ecl_to_eq.transpose() * vel_eq;
    
    cart_old.position = types::Vector3<core::GCRF, core::Meter>(pos_ecl);
    cart_old.velocity = types::Vector3<core::GCRF, core::Meter>(vel_ecl);
    cart_old.gravitational_parameter = constants::GM_SUN * 1e9;
    
    KeplerianElements kep_old = cartesian_to_keplerian(cart_old);
    
    current_orbit_ = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        cart_typed.epoch,
        kep_old.semi_major_axis, kep_old.eccentricity,
        kep_old.inclination * constants::RAD_TO_DEG,
        kep_old.longitude_ascending_node * constants::RAD_TO_DEG,
        kep_old.argument_perihelion * constants::RAD_TO_DEG,
        kep_old.mean_anomaly * constants::RAD_TO_DEG,
        physics::GravitationalParameter::from_si(cart_old.gravitational_parameter)
    );
    
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
    // 1. Bridge to un-typed format for conversion math
    KeplerianElements kep_start_old;
    kep_start_old.epoch = current_orbit_.epoch;
    kep_start_old.semi_major_axis = current_orbit_.a.to_au();
    kep_start_old.eccentricity = current_orbit_.e;
    kep_start_old.inclination = current_orbit_.i.to_rad();
    kep_start_old.longitude_ascending_node = current_orbit_.node.to_rad();
    kep_start_old.argument_perihelion = current_orbit_.omega.to_rad();
    kep_start_old.mean_anomaly = current_orbit_.M.to_rad();
    kep_start_old.gravitational_parameter = current_orbit_.gm.to_au3_d2();
    
    // 2. Convert to un-typed Cartesian (SI internally)
    CartesianElements cart_start_old = keplerian_to_cartesian(kep_start_old);
    
    // 3. Rotate to ICRF
    auto pos_ecl_start = cart_start_old.position.to_eigen();
    auto vel_ecl_start = cart_start_old.velocity.to_eigen();
    auto pos_eq_start = coordinates::ReferenceFrame::ecliptic_to_j2000() * pos_ecl_start;
    auto vel_eq_start = coordinates::ReferenceFrame::ecliptic_to_j2000() * vel_ecl_start;
    
    // 4. Bring into AstDyn 3.0 Type Safety
    auto initial_state = physics::CartesianStateTyped<core::GCRF>::from_si(
        current_orbit_.epoch,
        pos_eq_start.x(), pos_eq_start.y(), pos_eq_start.z(),
        vel_eq_start.x(), vel_eq_start.y(), vel_eq_start.z(),
        current_orbit_.gm.to_m3_s2()
    );
    
    // Setup differential corrector settings
    orbit_determination::DifferentialCorrectorSettings dc_settings;
    dc_settings.max_iterations = config_.max_iterations;
    dc_settings.convergence_tolerance = config_.convergence_threshold;
    
    // Outlier strategy: if initial guess is poor, use a very large threshold for the first step
    dc_settings.outlier_sigma = config_.outlier_sigma;
    dc_settings.outlier_max_sigma = std::max(100.0, config_.outlier_sigma * 10.0); 
    dc_settings.outlier_min_sigma = config_.outlier_sigma;
    dc_settings.reject_outliers = true;
    
    dc_settings.verbose = config_.verbose;
    
    // Perform differential correction using fit() method
    auto result_dc = corrector->fit(observations_, initial_state, dc_settings);
    
    // Convert final state back to Keplerian (Ecliptic J2000)
    auto cart_final = result_dc.final_state;
    auto pos_eq_final = cart_final.position.to_eigen_si();
    auto vel_eq_final = cart_final.velocity.to_eigen_si();
    auto mat_ecl_to_eq_final = coordinates::ReferenceFrame::ecliptic_to_j2000();
    auto pos_ecl_final = mat_ecl_to_eq_final.transpose() * pos_eq_final;
    auto vel_ecl_final = mat_ecl_to_eq_final.transpose() * vel_eq_final;
    
    CartesianElements cart_final_old;
    cart_final_old.epoch = cart_final.epoch;
    cart_final_old.position = types::Vector3<core::GCRF, core::Meter>(pos_ecl_final);
    cart_final_old.velocity = types::Vector3<core::GCRF, core::Meter>(vel_ecl_final);
    cart_final_old.gravitational_parameter = cart_final.gm.to_m3_s2();
    
    KeplerianElements kep_final_old = cartesian_to_keplerian(cart_final_old);
    
    current_orbit_ = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        cart_final.epoch,
        kep_final_old.semi_major_axis, kep_final_old.eccentricity,
        kep_final_old.inclination * constants::RAD_TO_DEG,
        kep_final_old.longitude_ascending_node * constants::RAD_TO_DEG,
        kep_final_old.argument_perihelion * constants::RAD_TO_DEG,
        kep_final_old.mean_anomaly * constants::RAD_TO_DEG,
        cart_final.gm
    );
    
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
    for (double mjd = start_time.mjd(); mjd <= end_time.mjd(); mjd += step_days) {
        epochs.push_back(time::EpochTDB::from_mjd(mjd));
    }
    if (epochs.back().mjd() < end_time.mjd()) {
        epochs.push_back(end_time);
    }
    
    // Convert orbit to Cartesian J2000 Equatorial
    // 1. Bridge to un-typed format
    KeplerianElements kep_old;
    kep_old.epoch = current_orbit_.epoch;
    kep_old.semi_major_axis = current_orbit_.a.to_au();
    kep_old.eccentricity = current_orbit_.e;
    kep_old.inclination = current_orbit_.i.to_rad();
    kep_old.longitude_ascending_node = current_orbit_.node.to_rad();
    kep_old.argument_perihelion = current_orbit_.omega.to_rad();
    kep_old.mean_anomaly = current_orbit_.M.to_rad();
    kep_old.gravitational_parameter = current_orbit_.gm.to_au3_d2();
    
    CartesianElements cart_old = keplerian_to_cartesian(kep_old);
    
    auto pos_ecl = cart_old.position.to_eigen();
    auto vel_ecl = cart_old.velocity.to_eigen();
    auto pos_eq = astdyn::coordinates::ReferenceFrame::ecliptic_to_j2000() * pos_ecl;
    auto vel_eq = astdyn::coordinates::ReferenceFrame::ecliptic_to_j2000() * vel_ecl;
    
    auto initial = physics::CartesianStateTyped<core::GCRF>::from_si(
        current_orbit_.epoch,
        pos_eq.x(), pos_eq.y(), pos_eq.z(),
        vel_eq.x(), vel_eq.y(), vel_eq.z(),
        current_orbit_.gm.to_m3_s2()
    );
    
    // Propagate to all epochs
    auto ephemeris = propagator_->propagate_ephemeris(initial, epochs);
    
    if (config_.verbose) {
        std::cout << "Generated " << ephemeris.size() << " ephemeris points\n";
    }
    
    return ephemeris;
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
    // 1. Bridge to un-typed format
    KeplerianElements kep_old;
    kep_old.epoch = current_orbit_.epoch;
    kep_old.semi_major_axis = current_orbit_.a.to_au();
    kep_old.eccentricity = current_orbit_.e;
    kep_old.inclination = current_orbit_.i.to_rad();
    kep_old.longitude_ascending_node = current_orbit_.node.to_rad();
    kep_old.argument_perihelion = current_orbit_.omega.to_rad();
    kep_old.mean_anomaly = current_orbit_.M.to_rad();
    kep_old.gravitational_parameter = current_orbit_.gm.to_au3_d2();
    
    CartesianElements cart_old = keplerian_to_cartesian(kep_old);
    
    auto pos_ecl = cart_old.position.to_eigen();
    auto vel_ecl = cart_old.velocity.to_eigen();
    auto pos_eq = astdyn::coordinates::ReferenceFrame::ecliptic_to_j2000() * pos_ecl;
    auto vel_eq = astdyn::coordinates::ReferenceFrame::ecliptic_to_j2000() * vel_ecl;
    
    auto initial = physics::CartesianStateTyped<core::GCRF>::from_si(
        current_orbit_.epoch,
        pos_eq.x(), pos_eq.y(), pos_eq.z(),
        vel_eq.x(), vel_eq.y(), vel_eq.z(),
        current_orbit_.gm.to_m3_s2()
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
