/**
 * @file AstDysOrbitFitter.hpp
 * @brief Utility class to perform orbit fitting using AstDyS data
 * 
 * This class downloads or accepts observations, orbital elements, and configuration
 * from AstDyS/OrbFit format files and performs differential correction.
 */

#pragma once

#include <astdyn/AstDynEngine.hpp>
#include <astdyn/io/AstDynConfig.hpp>
#include <astdyn/observations/MPCReader.hpp>
#include <astdyn/core/Constants.hpp>
#include <string>
#include <vector>
#include <memory>
#include <optional>

namespace astdyn {
namespace io {

/**
 * @brief Results from AstDyS orbit fitting
 */
struct AstDysFitResult {
    // Input information
    std::string object_name;
    int num_observations_loaded;
    int num_observations_used;
    
    // Fitted orbit
    propagation::KeplerianElements fitted_orbit;
    
    // Fit quality
    bool converged;
    int num_iterations;
    double rms_ra;       ///< arcsec
    double rms_dec;      ///< arcsec
    double chi_squared;
    int num_outliers;
    
    // Comparison with reference (if provided)
    std::optional<propagation::KeplerianElements> reference_orbit;
    std::optional<double> delta_a_km;
    std::optional<double> delta_e;
    std::optional<double> delta_i_arcsec;
};

/**
 * @brief Orbit fitter using AstDyS/OrbFit format files
 * 
 * This class simplifies the workflow:
 * 1. Load observations from .rwo (OrbFit with residuals) or .txt (MPC format)
 * 2. Load initial orbital elements from .eq1 (equinoctial) or .oel (OrbFit)
 * 3. Load configuration from .oop (OrbFit options)
 * 4. Run differential correction
 * 5. Compare with reference orbit
 * 
 * Example usage:
 * @code
 *   AstDysOrbitFitter fitter;
 *   fitter.set_observations_file("203.rwo");  // or .txt for MPC
 *   fitter.set_elements_file("203.eq1");      // or .oel
 *   fitter.set_config_file("203.oop");        // optional
 *   
 *   auto result = fitter.fit();
 *   if (result.converged) {
 *       std::cout << "RMS: " << result.rms_ra << " arcsec\n";
 *   }
 * @endcode
 */
class AstDysOrbitFitter {
public:
    /**
     * @brief Constructor
     */
    AstDysOrbitFitter();
    
    /**
     * @brief Set observations file (.rwo or .txt/.obs)
     * @param filename Path to observations file
     * @param format "rwo" (OrbFit with residuals) or "mpc" (standard MPC format)
     */
    void set_observations_file(const std::string& filename, 
                              const std::string& format = "auto");
    
    /**
     * @brief Set observations directly
     * @param observations Vector of optical observations
     */
    void set_observations(const std::vector<observations::OpticalObservation>& observations);
    
    /**
     * @brief Set orbital elements file (.eq1 equinoctial or .oel OrbFit)
     * @param filename Path to elements file
     * @param format "eq1" (equinoctial) or "oel" (OrbFit)
     */
    void set_elements_file(const std::string& filename,
                          const std::string& format = "auto");
    
    /**
     * @brief Set orbital elements directly
     * @param elements Keplerian elements
     */
    void set_elements(const propagation::KeplerianElements& elements);
    
    /**
     * @brief Set configuration file (.oop OrbFit options)
     * @param filename Path to configuration file
     */
    void set_config_file(const std::string& filename);
    
    /**
     * @brief Set reference orbit for comparison
     * @param elements Reference Keplerian elements
     */
    void set_reference_orbit(const propagation::KeplerianElements& elements);
    
    /**
     * @brief Enable/disable verbose output
     */
    void set_verbose(bool verbose) { verbose_ = verbose; }
    
    /**
     * @brief Perform orbit fitting
     * @return Fit results including fitted orbit and statistics
     */
    AstDysFitResult fit();
    
    /**
     * @brief Download files from AstDyS and perform fit
     * @param object_name Object designation (e.g., "203" or "Pompeja")
     * @param download_dir Directory to save downloaded files
     * @return Fit results
     */
    static AstDysFitResult fit_from_astdys(const std::string& object_name,
                                           const std::string& download_dir = "/tmp");

private:
    // Input data
    std::vector<observations::OpticalObservation> observations_;
    std::optional<propagation::KeplerianElements> initial_elements_;
    std::optional<propagation::KeplerianElements> reference_orbit_;
    
    // Configuration
    std::optional<AstDynConfig> config_;
    bool verbose_;
    
    // Helper methods
    void load_rwo_file(const std::string& filename);
    void load_mpc_file(const std::string& filename);
    void load_eq1_file(const std::string& filename);
    void load_oel_file(const std::string& filename);
    void load_oop_file(const std::string& filename);
    
    propagation::KeplerianElements equinoctial_to_keplerian(
        double a, double h, double k, double p, double q, double lambda, double mjd);
    
    std::string detect_format(const std::string& filename);
};

} // namespace io
} // namespace astdyn
