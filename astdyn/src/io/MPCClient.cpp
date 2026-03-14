#include "astdyn/io/MPCClient.hpp"
#include "astdyn/observations/MPCReader.hpp"
#include <curl/curl.h>
#include <sstream>
#include <iostream>
#include <regex>

namespace astdyn::io {

// --- CURL Helpers ---
static size_t mpc_write_callback(void* contents, size_t size, size_t nmemb, void* userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

MPCClient::MPCClient(const Config& config) : config_(config) {
    curl_global_init(CURL_GLOBAL_DEFAULT);
    curl_handle_ = curl_easy_init();
}

MPCClient::~MPCClient() {
    if (curl_handle_) curl_easy_cleanup(curl_handle_);
    curl_global_cleanup();
}

std::expected<std::vector<observations::OpticalObservation>, MPCError> 
MPCClient::fetch_observations(const std::string& target) {
    if (!curl_handle_) return std::unexpected(MPCError::NetworkError);

    // MPC Query URL (simplified for example, real implementation would use MPC's export API)
    // Here we use a hypothetical endpoint or a known one like 'https://www.minorplanetcenter.net/db_search/show_object?object_id=...'
    // For this task, we assume we want to download the MPC format file.
    
    std::string url = "https://www.minorplanetcenter.net/db_search/show_object?object_id=" + target + "&format=mpc";
    
    std::string readBuffer;
    curl_easy_setopt(curl_handle_, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl_handle_, CURLOPT_WRITEFUNCTION, mpc_write_callback);
    curl_easy_setopt(curl_handle_, CURLOPT_WRITEDATA, &readBuffer);
    curl_easy_setopt(curl_handle_, CURLOPT_TIMEOUT, config_.timeout_seconds);
    curl_easy_setopt(curl_handle_, CURLOPT_FOLLOWLOCATION, 1L);

    CURLcode res = curl_easy_perform(curl_handle_);
    if (res != CURLE_OK) return std::unexpected(MPCError::NetworkError);

    // If result contains MPC format lines, they are usually in a <pre> block or direct text
    // The MPCReader::fromString handles a multi-line string in MPC 80-column format.
    try {
        // Simple heuristic: search for the 80-char lines
        // A more robust way would be to parse the HTML if it's not a direct format.
        // For now, we assume 'readBuffer' contains the MPC data.
        return observations::MPCReader::parseString(readBuffer);
    } catch (...) {
        return std::unexpected(MPCError::ParsingError);
    }
}

} // namespace astdyn::io
