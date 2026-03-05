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
#include "src/core/frame_tags.hpp"
#include "src/core/units.hpp"
#include "src/types/orbital_state.hpp"
#include <string>
#include <optional>

namespace astdyn {
namespace api {

/**
 * @brief Result structure for OrbitFitAPI
 */
struct OrbitFitResult {
    bool success = false;
    std::string message;
    
    // Fitted Orbit (GCRF / Equatorial J2000)
    std::optional<types::OrbitalState<core::GCRF, types::KeplerianTag>> fitted_state;
    
    // Statistics
    core::MilliArcSecond rms_ra = core::MilliArcSecond(0.0);
    core::MilliArcSecond rms_dec = core::MilliArcSecond(0.0);
    int iterations = 0;
    int num_observations = 0;
    int num_outliers = 0;
    bool converged = false;
    
    // Comparison metrics
    double delta_a_km = 0.0;
    double delta_e = 0.0;

    OrbitFitResult() = default;
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
     * @param config_file Path to configuration file (.json), optional
     * @param verbose  Enable verbose output to stdout
     * @return OrbitFitResult structure containing fit statistics and final orbit
     */
    static OrbitFitResult run_fit(
        const std::string& eq1_file,
        const std::string& rwo_file,
        const std::string& config_file = "",
        bool verbose = true,
        const std::string& de441_path = ""
    );

    /**
     * @brief Parse an OrbFit .eq1 file to proper orbital elements
     * 
     * Handles specific OrbFit format quirks and units.
     */
    static propagation::EquinoctialElements parse_eq1(const std::string& filepath);
    
    /**
     * @brief Convert Mean Equinoctial elements (Ecliptic) to Osculating Keplerian (Equatorial)
     */
    static types::OrbitalState<core::GCRF, types::KeplerianTag> 
    prepare_initial_state(const propagation::EquinoctialElements& mean_equ);

private:
     // Helper to trim strings
    static std::string ltrim(const std::string& s);
};

} // namespace api
} // namespace astdyn

#endif // ASTDYN_API_ORBITFITAPI_HPP
