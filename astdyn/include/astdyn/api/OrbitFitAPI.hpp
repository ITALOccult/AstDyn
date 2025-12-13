/**
 * @file OrbitFitAPI.hpp
 * @brief Simple API for OrbFit-style orbit fitting workflow
 * @author ITALOccult AstDyn Team
 * @date 2025-12-13
 */

#ifndef ASTDYN_API_ORBITFITAPI_HPP
#define ASTDYN_API_ORBITFITAPI_HPP

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include <string>
#include <optional>

namespace astdyn {
namespace api {

/**
 * @brief Result structure for OrbitFitAPI
 */
struct OrbitFitResult {
    bool success;
    std::string message;
    
    // Fitted Orbit (Equatorial J2000)
    propagation::KeplerianElements fitted_orbit;
    
    // Statistics
    double rms_ra;        // [arcsec]
    double rms_dec;       // [arcsec]
    int iterations;
    int num_observations;
    int num_outliers;
    bool converged;
    
    // Comparison with initial orbit (if provided)
    double delta_a_km;    // [km]
    double delta_e;       // dimensionless
    double delta_i_arcsec;// [arcsec]
};

/**
 * @brief Simplified API for performing orbit fits using OrbFit file formats
 */
class OrbitFitAPI {
public:
    /**
     * @brief Run a complete orbit fit using .eq1 and .rwo files
     * 
     * This function encapsulates the entire workflow:
     * 1. Parses the .eq1 file for initial Equinoctial elements (Ecliptic J2000).
     * 2. Loads observations from the .rwo file.
     * 3. Transforms valid initial elements to Equatorial J2000.
     * 4. Configures the engine (loading .oop if provided, or defaults).
     * 5. Runs the differential correction.
     * 
     * @param eq1_file Path to OrbFit .eq1 file (initial orbit)
     * @param rwo_file Path to OrbFit .rwo file (observations)
     * @param oop_file Path to OrbFit .oop file (configuration), optional
     * @param verbose  Enable verbose output to stdout
     * @return OrbitFitResult structure containing fit statistics and final orbit
     */
    static OrbitFitResult run_fit(
        const std::string& eq1_file,
        const std::string& rwo_file,
        const std::string& oop_file = "",
        bool verbose = true
    );

    /**
     * @brief Parse an OrbFit .eq1 file to proper orbital elements
     * 
     * Handles specific OrbFit format quirks and units.
     */
    static propagation::EquinoctialElements parse_eq1(const std::string& filepath);

private:
     // Helper to trim strings
    static std::string ltrim(const std::string& s);
};

} // namespace api
} // namespace astdyn

#endif // ASTDYN_API_ORBITFITAPI_HPP
