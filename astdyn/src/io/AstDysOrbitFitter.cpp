/**
 * @file AstDysOrbitFitter.cpp
 * @brief Implementation of AstDyS orbit fitting utility
 */

#include <astdyn/io/AstDysOrbitFitter.hpp>
#include <astdyn/observations/MPCReader.hpp>
#include <curl/curl.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <nlohmann/json.hpp>

namespace {
    size_t write_data(void* ptr, size_t size, size_t nmemb, FILE* stream) {
        return fwrite(ptr, size, nmemb, stream);
    }

    bool download_file(const std::string& url, const std::string& path) {
        CURL* curl = curl_easy_init();
        if (!curl) return false;
        FILE* fp = fopen(path.c_str(), "wb");
        if (!fp) { curl_easy_cleanup(curl); return false; }
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        CURLcode res = curl_easy_perform(curl);
        fclose(fp); curl_easy_cleanup(curl);
        return (res == CURLE_OK);
    }
}

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

void AstDysOrbitFitter::set_elements(const physics::KeplerianStateTyped<core::ECLIPJ2000>& elements) {
    initial_elements_ = elements;
    
    if (verbose_) {
        std::cout << "Set orbital elements directly\n";
    }
}

void AstDysOrbitFitter::set_config_file(const std::string& filename) {
    load_json_file(filename);
    
    if (verbose_) {
        std::cout << "Loaded configuration from " << filename << "\n";
    }
}

void AstDysOrbitFitter::set_reference_orbit(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& elements) {
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
        std::cout << "Initial epoch: MJD " << initial_elements_->epoch.mjd() << " MJD TDB\n";
        std::cout << "Initial a: " << initial_elements_->a.to_au() << " AU\n\n";
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
    result.object_name = "Unknown Object"; // Default or extracted
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
        
        double da = (result.fitted_orbit.a.to_au() - 
                    reference_orbit_->a.to_au()) * constants::AU;
        double de = result.fitted_orbit.e - reference_orbit_->e;
        double di = (result.fitted_orbit.i.to_rad() - 
                    reference_orbit_->i.to_rad()) * 180.0 / constants::PI;
        
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
        // RWO lines often start with a space. The MPC record (80 chars) starts after that.
        if (line.empty()) continue;
        size_t start = 0;
        if (line[0] == ' ') start = 1;
        
        if (line.length() >= start + 80) {
            temp << line.substr(start, 80) << "\n";
        }
    }
    temp.close();
    
    // Parse with MPCReader and filter for post-1970 observations to avoid EOP issues
    auto all_obs = observations::MPCReader::readFile(temp_file);
    observations_.clear();
    for (const auto& obs : all_obs) {
        auto [y, m, d, f] = time::mjd_to_calendar(obs.time.mjd());
        if (y >= 1970) {
            observations_.push_back(obs);
        }
    }
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
        
        if (iss >> key) {
            if (key == "EQU") {
                iss >> a >> h >> k >> p >> q >> lambda;
            } else if (key == "MJD") {
                double val;
                std::string tdt;
                iss >> val >> tdt;
                mjd = val;
            }
        }
    }
    
    initial_elements_ = equinoctial_to_keplerian(a, h, k, p, q, lambda, time::EpochTDB::from_mjd(mjd));
}

void AstDysOrbitFitter::load_oel_file(const std::string& filename) {
    // OEL format is the same as EQ1 (equinoctial elements)
    load_eq1_file(filename);
}

void AstDysOrbitFitter::load_json_file(const std::string& filename) {
    // Parse JSON options file
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JSON config file: " + filename);
    }
    
    config_ = AstDynConfig();
    nlohmann::json j;
    try {
        file >> j;
    } catch (...) {
         throw std::runtime_error("Failed to parse JSON config: " + filename);
    }
    
    if (j.contains("diffcorr")) {
        auto& jd = j["diffcorr"];
        if (jd.contains("max_iter") || jd.contains("max_iterations")) 
             config_->max_iterations = jd.value("max_iter", jd.value("max_iterations", 20));
             
        if (jd.contains("convergence") || jd.contains("convergence_threshold"))
             config_->convergence_threshold = jd.value("convergence", jd.value("convergence_threshold", 1e-6));
    } else {
        // top level fallback
        if (j.contains("max_iterations")) config_->max_iterations = j["max_iterations"];
        if (j.contains("convergence_threshold")) config_->convergence_threshold = j["convergence_threshold"];
    }
}

physics::KeplerianStateTyped<core::ECLIPJ2000> AstDysOrbitFitter::equinoctial_to_keplerian(
    double a, double h, double k, double p, double q, double lambda, time::EpochTDB epoch) {
    
    // Eccentricity
    double e = std::sqrt(h*h + k*k);
    
    // Inclination
    double tan_half_i = std::sqrt(p*p + q*q);
    double i_rad = 2.0 * std::atan(tan_half_i);
    
    // Longitude of ascending node
    double Omega = std::atan2(p, q);
    if (Omega < 0) Omega += 2.0 * constants::PI;
    
    // Argument of perihelion
    double omega_plus_Omega = std::atan2(h, k);
    if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * constants::PI;
    double omega = omega_plus_Omega - Omega;
    if (omega < 0) omega += 2.0 * constants::PI;
    
    // Mean anomaly
    double M = lambda - omega_plus_Omega;
    while (M < 0) M += 2.0 * constants::PI;
    while (M >= 2.0*constants::PI) M -= 2.0 * constants::PI;
    
    return physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
        epoch,
        a * constants::AU, e, i_rad * constants::RAD_TO_DEG,
        Omega * constants::RAD_TO_DEG,
        omega * constants::RAD_TO_DEG,
        M * constants::RAD_TO_DEG,
        physics::GravitationalParameter::sun()
    );
}

std::string AstDysOrbitFitter::detect_format(const std::string& filename) {
    // Detect format from file extension
    if (filename.find(".rwo") != std::string::npos) {
        return "rwo";
    } else if (filename.find(".eq1") != std::string::npos) {
        return "eq1";
    } else if (filename.find(".oel") != std::string::npos) {
        return "oel";
    } else if (filename.find(".json") != std::string::npos) {
        return "json";
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
    
    if (!download_file(oel_url, oel_file)) throw std::runtime_error("Failed to download elements from: " + oel_url);
    if (!download_file(rwo_url, rwo_file)) throw std::runtime_error("Failed to download observations from: " + rwo_url);

    AstDysOrbitFitter instance;
    instance.set_observations_file(rwo_file, "rwo");
    instance.set_elements_file(oel_file, "oel");
    return instance.fit();
}

} // namespace io
} // namespace astdyn
