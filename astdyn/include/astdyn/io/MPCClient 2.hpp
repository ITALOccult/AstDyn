#ifndef ASTDYN_IO_MPCCLIENT_HPP
#define ASTDYN_IO_MPCCLIENT_HPP

#include "astdyn/observations/ObservationManager.hpp"
#include <string>
#include <vector>
#include <expected>

namespace astdyn::io {

enum class MPCError {
    NetworkError,
    TargetNotFound,
    ParsingError
};

/**
 * @brief Client for downloading asteroid observations from the Minor Planet Center (MPC).
 */
class MPCClient {
public:
    struct Config {
        std::string base_url;
        long timeout_seconds;
        
        Config() 
            : base_url("https://www.minorplanetcenter.net/db_search/show_object"),
              timeout_seconds(30) {}
    };

    explicit MPCClient(const Config& config = Config());
    ~MPCClient();

    /**
     * @brief Downloads observations for a given object.
     * @param target Object designation (e.g. "Vesta" or "4").
     * @return List of observations or error.
     */
    std::expected<std::vector<observations::OpticalObservation>, MPCError> 
    fetch_observations(const std::string& target);

private:
    Config config_;
    void* curl_handle_; // Internal CURL handle
};

} // namespace astdyn::io

#endif // ASTDYN_IO_MPCCLIENT_HPP
