#ifndef ASTDYN_IO_HORIZONS_CLIENT_HPP
#define ASTDYN_IO_HORIZONS_CLIENT_HPP

#include "src/utils/time_types.hpp"
#include "astdyn/core/physics_state.hpp"
#include "astdyn/astrometry/AstrometricTypes.hpp"
#include <string>
#include <expected>
#include <vector>

namespace astdyn::io {

enum class HorizonsError {
    NetworkError,
    InvalidResponse,
    TargetNotFound,
    ParsingError,
    Timeout
};

/**
 * @brief JPL Horizons REST API Client
 * 
 * Provides high-level access to JPL Horizons data (Vectors, Elements, Observer).
 * Adheres to CTFYH standard: Type-safe, functional, immutable results.
 */
class HorizonsClient {
public:
    struct Config {
        std::string base_url;
        int timeout_seconds;
        
        Config() : base_url("https://ssd.jpl.nasa.gov/api/horizons.api"), timeout_seconds(30) {}
    };

    explicit HorizonsClient(const Config& config = Config());
    ~HorizonsClient();

    /**
     * @brief Query orbital elements for a specific target and epoch
     * @return Strongly-typed KeplerianState in ECLIPJ2000 frame
     */
    std::expected<physics::KeplerianStateTyped<core::ECLIPJ2000>, HorizonsError> 
    query_elements(const std::string& target, const utils::Instant& epoch);

    /**
     * @brief Query state vectors for a specific target and epoch.
     * @return Strongly-typed CartesianState in GCRF/ICRF frame.
     */
    std::expected<physics::CartesianStateTyped<core::GCRF>, HorizonsError> 
    query_vectors(const std::string& target, const utils::Instant& epoch, const std::string& center = "500@0");

    /**
     * @brief Query observed position (RA/Dec) for a specific target and epoch
     * @return AstrometricObservation
     */
    std::expected<astrometry::AstrometricObservation, HorizonsError> 
    query_observation(const std::string& target, const utils::Instant& epoch, const std::string& observer = "500@399");

private:
    Config config_;
    
    // Internal helper to perform HTTP GET
    std::expected<std::string, HorizonsError> fetch_url(const std::string& url);
    
    // URL Builders
    std::string build_elements_url(const std::string& target, const utils::Instant& epoch);
    std::string build_vectors_url(const std::string& target, const utils::Instant& epoch, const std::string& center);
    std::string build_observer_url(const std::string& target, const utils::Instant& epoch, const std::string& center);
};

} // namespace astdyn::io

#endif // ASTDYN_IO_HORIZONS_CLIENT_HPP
