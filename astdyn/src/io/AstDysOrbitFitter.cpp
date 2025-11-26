/**
 * @file AstDysOrbitFitter.cpp
 * @brief Implementation of AstDyS orbit fitting utility
 */

#include <astdyn/io/AstDysOrbitFitter.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

namespace astdyn {
namespace io {

AstDysOrbitFitter::AstDysOrbitFitter() 
    : verbose_(false) {
}

void AstDysOrbitFitter::set_observations_file(const std::string& filename, 
                                               const std::string& format) {
    std::string actual_format = format;
    if (format == "auto") {
        actual_format = detect_format(filename);
    }
    
    if (actual_format == "rwo") {
        load_rwo_file(filename);
    } else if (actual_format == "mpc") {
        load_mpc_file(filename);
    } else {
        throw std::runtime_error("Unknown observation format: " + actual_format);
    }
    
    if (verbose_) {
        std::cout << "Loaded " << observations_.size() << " observations from " 
                  << filename << " (format: " << actual_format << ")\n";
    }
}

void AstDysOrbitFitter::set_observations(
    const std::vector<observations::OpticalObservation>& observations) {
    observations_ = observations;
    
    if (verbose_) {
        std::cout << "Set " << observations_.size() << " observations directly\n";
    }
}

void AstDysOrbitFitter::set_elements_file(const std::string& filename,
                                          const std::string& format) {
    std::string actual_format = format;
    if (format == "auto") {
        actual_format = detect_format(filename);
    }
    
    if (actual_format == "eq1") {
        load_eq1_file(filename);
    } else if (actual_format == "oel") {
        load_oel_file(filename);
    } else {
        throw std::runtime_error("Unknown elements format: " + actual_format);
    }
    
    if (verbose_) {
        std::cout << "Loaded orbital elements from " << filename 
                  << " (format: " << actual_format << ")\n";
    }
}

void AstDysOrbitFitter::set_elements(const propagation::KeplerianElements& elements) {
    initial_elements_ = elements;
    
    if (verbose_) {
        std::cout << "Set orbital elements directly\n";
    }
}

void AstDysOrbitFitter::set_config_file(const std::string& filename) {
    load_oop_file(filename);
    
    if (verbose_) {
        std::cout << "Loaded configuration from " << filename << "\n";
    }
}

void AstDysOrbitFitter::set_reference_orbit(
    const propagation::KeplerianElements& elements) {
    reference_orbit_ = elements;
}

AstDysFitResult AstDysOrbitFitter::fit() {
    AstDysFitResult result;
    
    // Validate inputs
    if (observations_.empty()) {
        throw std::runtime_error("No observations loaded");
    }
    if (!initial_elements_.has_value()) {
        throw std::runtime_error("No initial orbital elements provided");
    }
    
    result.num_observations_loaded = observations_.size();
    
    if (verbose_) {
        std::cout << "\n=== Starting Orbit Fit ===\n";
        std::cout << "Observations: " << observations_.size() << "\n";
        std::cout << "Initial epoch: MJD " << initial_elements_->epoch_mjd_tdb << "\n";
        std::cout << "Initial a: " << initial_elements_->semi_major_axis << " AU\n\n";
    }
    
    // Create AstDyn engine
    AstDynEngine engine;
    
    // Apply configuration if provided
    if (config_.has_value()) {
        engine.set_config(*config_);
    }
    
    // Set verbose mode
    engine.set_verbose(verbose_);
    
    // Load observations
    for (const auto& obs : observations_) {
        engine.add_observation(obs);
    }
    
    // Set initial orbit
    engine.set_initial_orbit(*initial_elements_);
    
    // Run differential correction
    auto fit_result = engine.fit_orbit();
    
    // Fill result structure
    result.fitted_orbit = fit_result.orbit;
    result.converged = fit_result.converged;
    result.num_iterations = fit_result.num_iterations;
    result.rms_ra = fit_result.rms_ra;
    result.rms_dec = fit_result.rms_dec;
    result.chi_squared = fit_result.chi_squared;
    result.num_observations_used = fit_result.num_observations;
    result.num_outliers = fit_result.num_rejected;
    
    // Compare with reference orbit if provided
    if (reference_orbit_.has_value()) {
        result.reference_orbit = *reference_orbit_;
        
        double da = (result.fitted_orbit.semi_major_axis - 
                    reference_orbit_->semi_major_axis) * constants::AU;
        double de = result.fitted_orbit.eccentricity - reference_orbit_->eccentricity;
        double di = (result.fitted_orbit.inclination - 
                    reference_orbit_->inclination) * 180.0 / constants::PI;
        
        result.delta_a_km = da / 1000.0;
        result.delta_e = de;
        result.delta_i_arcsec = di * 3600.0;
        
        if (verbose_) {
            std::cout << "\n=== Comparison with Reference ===\n";
            std::cout << "Δa = " << *result.delta_a_km << " km\n";
            std::cout << "Δe = " << std::scientific << *result.delta_e << std::fixed << "\n";
            std::cout << "Δi = " << *result.delta_i_arcsec << " arcsec\n";
        }
    }
    
    if (verbose_) {
        std::cout << "\n=== Fit Complete ===\n";
        std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
        std::cout << "RMS RA: " << result.rms_ra << " arcsec\n";
        std::cout << "RMS Dec: " << result.rms_dec << " arcsec\n";
        std::cout << "Observations used: " << result.num_observations_used 
                  << " / " << result.num_observations_loaded << "\n";
    }
    
    return result;
}

// Private helper methods

void AstDysOrbitFitter::load_rwo_file(const std::string& filename) {
    // RWO format: First 80 columns are MPC format + extended fields
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open RWO file: " + filename);
    }
    
    // Extract MPC format (first 80 columns) to temporary file
    std::string temp_file = "/tmp/astdyn_temp_mpc.txt";
    std::ofstream temp(temp_file);
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.length() >= 80) {
            temp << line.substr(0, 80) << "\n";
        }
    }
    temp.close();
    
    // Parse with MPCReader
    observations_ = observations::MPCReader::readFile(temp_file);
}

void AstDysOrbitFitter::load_mpc_file(const std::string& filename) {
    observations_ = observations::MPCReader::readFile(filename);
}

void AstDysOrbitFitter::load_eq1_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open EQ1 file: " + filename);
    }
    
    // Parse equinoctial elements
    double a, h, k, p, q, lambda, mjd;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '!' || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string key;
        double value;
        
        if (iss >> key >> value) {
            if (key == "a") a = value;
            else if (key == "h") h = value;
            else if (key == "k") k = value;
            else if (key == "p") p = value;
            else if (key == "q") q = value;
            else if (key == "lambda") lambda = value;
            else if (key == "MJD") mjd = value;
        }
    }
    
    initial_elements_ = equinoctial_to_keplerian(a, h, k, p, q, lambda, mjd);
}

void AstDysOrbitFitter::load_oel_file(const std::string& filename) {
    // OEL format is the same as EQ1 (equinoctial elements)
    load_eq1_file(filename);
}

void AstDysOrbitFitter::load_oop_file(const std::string& filename) {
    // Parse OrbFit options file
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open OOP file: " + filename);
    }
    
    // For now, create default config
    // Full implementation would parse .oop and set all options
    config_ = AstDynConfig();
    
    // Parse basic options
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '!' || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string key, equals, value;
        
        if (iss >> key >> equals >> value) {
            if (key == "max_iterations") {
                config_->max_iterations = std::stoi(value);
            } else if (key == "convergence_threshold") {
                config_->convergence_threshold = std::stod(value);
            }
        }
    }
}

propagation::KeplerianElements AstDysOrbitFitter::equinoctial_to_keplerian(
    double a, double h, double k, double p, double q, double lambda, double mjd) {
    
    propagation::KeplerianElements kep;
    
    kep.semi_major_axis = a;
    
    // Eccentricity
    double e = std::sqrt(h*h + k*k);
    kep.eccentricity = e;
    
    // Inclination
    double tan_half_i = std::sqrt(p*p + q*q);
    kep.inclination = 2.0 * std::atan(tan_half_i);
    
    // Longitude of ascending node
    double Omega = std::atan2(p, q);
    if (Omega < 0) Omega += 2.0 * constants::PI;
    kep.longitude_ascending_node = Omega;
    
    // Argument of perihelion
    double omega_plus_Omega = std::atan2(h, k);
    if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * constants::PI;
    double omega = omega_plus_Omega - Omega;
    if (omega < 0) omega += 2.0 * constants::PI;
    kep.argument_perihelion = omega;
    
    // Mean anomaly
    double M = lambda - omega_plus_Omega;
    while (M < 0) M += 2.0 * constants::PI;
    while (M >= 2.0*constants::PI) M -= 2.0 * constants::PI;
    kep.mean_anomaly = M;
    
    kep.epoch_mjd_tdb = mjd;
    kep.gravitational_parameter = constants::GMS;
    
    return kep;
}

std::string AstDysOrbitFitter::detect_format(const std::string& filename) {
    // Detect format from file extension
    if (filename.find(".rwo") != std::string::npos) {
        return "rwo";
    } else if (filename.find(".eq1") != std::string::npos) {
        return "eq1";
    } else if (filename.find(".oel") != std::string::npos) {
        return "oel";
    } else if (filename.find(".oop") != std::string::npos) {
        return "oop";
    } else if (filename.find(".txt") != std::string::npos || 
               filename.find(".obs") != std::string::npos) {
        return "mpc";
    }
    
    // Default to MPC for observations, eq1 for elements
    return "mpc";
}

AstDysFitResult AstDysOrbitFitter::fit_from_astdys(
    const std::string& object_name, const std::string& download_dir) {
    
    // URLs for AstDyS files
    // Format: https://newton.spacedys.com/~astdys2/propsynth/{numbered|unnumbered}/
    // For numbered: https://newton.spacedys.com/~astdys2/propsynth/numbered/0/203.oel
    // For obs: https://newton.spacedys.com/~astdys2/mpcobs/numbered/0/203.rwo
    
    std::string base_url_elements = "https://newton.spacedys.com/~astdys2/propsynth/numbered/";
    std::string base_url_obs = "https://newton.spacedys.com/~astdys2/mpcobs/numbered/";
    
    // Parse object name to get number
    std::string obj_num = object_name;
    // Remove non-numeric characters
    obj_num.erase(std::remove_if(obj_num.begin(), obj_num.end(), 
        [](char c) { return !std::isdigit(c); }), obj_num.end());
    
    if (obj_num.empty()) {
        throw std::runtime_error("fit_from_astdys: Cannot parse object number from '" + 
                                object_name + "'");
    }
    
    // Determine subdirectory (first digit for numbered asteroids)
    std::string subdir = obj_num.substr(0, 1);
    
    // Construct file paths
    std::string oel_file = download_dir + "/" + obj_num + ".oel";
    std::string rwo_file = download_dir + "/" + obj_num + ".rwo";
    
    // Construct URLs
    std::string oel_url = base_url_elements + subdir + "/" + obj_num + ".oel";
    std::string rwo_url = base_url_obs + subdir + "/" + obj_num + ".rwo";
    
    std::cout << "\n=== Downloading from AstDyS ===\n";
    std::cout << "Object: " << object_name << " (" << obj_num << ")\n";
    std::cout << "Elements URL: " << oel_url << "\n";
    std::cout << "Observations URL: " << rwo_url << "\n";
    
    // Download files using curl or wget
    // For now, provide instructions and throw
    std::string instructions = 
        "\nTo use fit_from_astdys(), please download files manually:\n"
        "  wget " + oel_url + " -O " + oel_file + "\n"
        "  wget " + rwo_url + " -O " + rwo_file + "\n"
        "\nThen use:\n"
        "  fitter.set_observations_file(\"" + rwo_file + "\", \"rwo\");\n"
        "  fitter.set_elements_file(\"" + oel_file + "\", \"oel\");\n"
        "  auto result = fitter.fit();\n";
    
    throw std::runtime_error("fit_from_astdys: Automatic download not yet implemented.\n" + 
                            instructions);
    
    // Future implementation would use libcurl:
    // 1. Download .oel file
    // 2. Download .rwo file  
    // 3. Call set_observations_file() and set_elements_file()
    // 4. Call fit()
    // 5. Return result
}

} // namespace io
} // namespace astdyn
