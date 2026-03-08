#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <sstream>
#include <iomanip>
#include <regex>

namespace astdyn::io {

using json = nlohmann::json;

// --- CURL Helpers ---
static size_t write_callback(void* contents, size_t size, size_t nmemb, void* userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

static std::string url_encode(const std::string& value) {
    CURL* curl = curl_easy_init();
    if (!curl) return value;
    char* output = curl_easy_escape(curl, value.c_str(), value.length());
    if (!output) {
        curl_easy_cleanup(curl);
        return value;
    }
    std::string res(output);
    curl_free(output);
    curl_easy_cleanup(curl);
    return res;
}

HorizonsClient::HorizonsClient(const Config& config) : config_(config) {
    curl_global_init(CURL_GLOBAL_DEFAULT);
}

HorizonsClient::~HorizonsClient() {
    curl_global_cleanup();
}

std::expected<std::string, HorizonsError> HorizonsClient::fetch_url(const std::string& url) {
    CURL* curl = curl_easy_init();
    if (!curl) return std::unexpected(HorizonsError::NetworkError);

    std::string readBuffer;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, config_.timeout_seconds);

    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);

    if (res != CURLE_OK) return std::unexpected(HorizonsError::NetworkError);
    return readBuffer;
}

// --- Query Implementations ---

std::expected<physics::KeplerianStateTyped<core::ECLIPJ2000>, HorizonsError> 
HorizonsClient::query_elements(const std::string& target, const time::EpochTDB& epoch) {
    auto url = build_elements_url(target, epoch);
    auto raw = fetch_url(url);
    if (!raw) return std::unexpected(raw.error());

    try {
        auto data = json::parse(*raw);
        std::string result = data["result"];
        
        auto find_val = [&](const std::string& token) -> double {
            size_t pos = result.find(token);
            if (pos == std::string::npos) throw std::runtime_error("Token not found: " + token);
            size_t val_pos = result.find_first_of("0123456789+-.", pos + token.length());
            if (val_pos == std::string::npos) throw std::runtime_error("Value not found for: " + token);
            return std::stod(result.substr(val_pos));
        };

        double a_au = find_val("A =");
        double e = find_val("EC=");
        double i = find_val("IN=");
        double om = find_val("OM=");
        double w = find_val("W =");
        double ma = find_val("MA=");

        return physics::KeplerianStateTyped<core::ECLIPJ2000>{
            epoch,
            physics::Distance::from_au(a_au),
            e,
            astrometry::Angle::from_deg(i),
            astrometry::Angle::from_deg(om),
            astrometry::Angle::from_deg(w),
            astrometry::Angle::from_deg(ma),
            physics::GravitationalParameter::sun()
        };
    } catch (...) {
        return std::unexpected(HorizonsError::ParsingError);
    }
}

std::expected<physics::CartesianStateTyped<core::GCRF>, HorizonsError> 
HorizonsClient::query_vectors(const std::string& target, const time::EpochTDB& epoch, const std::string& center) {
    auto url = build_vectors_url(target, epoch, center);
    auto raw = fetch_url(url);
    if (!raw) return std::unexpected(raw.error());

    try {
        auto data = json::parse(*raw);
        std::string result = data["result"];
        
        // Search for the Start of the Data block $$SOE
        size_t start = result.find("$$SOE");
        if (start == std::string::npos) return std::unexpected(HorizonsError::TargetNotFound);
        
        std::string block = result.substr(start);
        auto find_coord = [&](const std::string& token) -> double {
            size_t pos = block.find(token);
            if (pos == std::string::npos) throw std::runtime_error("Coord not found: " + token);
            size_t val_pos = block.find_first_of("0123456789+-.", pos + token.length());
            return std::stod(block.substr(val_pos));
        };

        // Horizons VECTOR ephemeris returns positions in km, velocities in km/s.
        // Standardize to meters and m/s to match convention. Return type is MeterTag.
        // Determine GM based on center (simplified)
        double gm = constants::GM_SUN;
        if (center.find("399") != std::string::npos || center.find("Earth") != std::string::npos) {
            gm = constants::GM_EARTH;
        }

        constexpr double KM_TO_M = 1000.0;
        return physics::CartesianStateTyped<core::GCRF>::from_si(
            epoch,
            find_coord(" X =") * KM_TO_M, find_coord(" Y =") * KM_TO_M, find_coord(" Z =") * KM_TO_M,
            find_coord(" VX=") * KM_TO_M, find_coord(" VY=") * KM_TO_M, find_coord(" VZ=") * KM_TO_M,
            gm * 1e9 // GM in m^3/s^2
        );
    } catch (...) {
        return std::unexpected(HorizonsError::ParsingError);
    }
}

std::expected<astrometry::AstrometricObservation, HorizonsError> 
HorizonsClient::query_observation(const std::string& target, const time::EpochTDB& epoch, const std::string& observer) {
    auto url = build_observer_url(target, epoch, observer);
    auto raw = fetch_url(url);
    if (!raw) return std::unexpected(raw.error());

    try {
        auto data = json::parse(*raw);
        std::string result = data["result"];
        
        size_t start = result.find("$$SOE");
        if (start == std::string::npos) return std::unexpected(HorizonsError::TargetNotFound);
        
        size_t line_start = result.find('\n', start);
        std::string line = result.substr(line_start + 1);
        std::istringstream iss(line);
        std::string date, time;
        iss >> date >> time;
        
        std::string ra_h, ra_m, ra_s;
        std::string dec_d, dec_m, dec_s;
        iss >> ra_h >> ra_m >> ra_s >> dec_d >> dec_m >> dec_s;
        
        double ra = (std::stod(ra_h) + std::stod(ra_m)/60.0 + std::stod(ra_s)/3600.0) * 15.0 * constants::DEG_TO_RAD;
        
        bool negative = (dec_d.find('-') != std::string::npos);
        double d_deg = std::abs(std::stod(dec_d));
        double dec = (d_deg + std::stod(dec_m)/60.0 + std::stod(dec_s)/3600.0) * constants::DEG_TO_RAD;
        if (negative) dec = -dec;

        double delta;
        iss >> delta;

        return astrometry::AstrometricObservation{
            core::Radian(ra), core::Radian(dec), core::Meter(delta * constants::AU * 1000.0)
        };
    } catch (...) {
        return std::unexpected(HorizonsError::ParsingError);
    }
}

// --- URL Builders ---

std::string HorizonsClient::build_elements_url(const std::string& target, const time::EpochTDB& epoch) {
    std::stringstream ss;
    ss << config_.base_url << "?format=json&COMMAND='" << url_encode(target) << "'"
       << "&OBJ_DATA='NO'&MAKE_EPHEM='YES'&EPHEM_TYPE='ELEMENTS'&CENTER='500@10'"
       << "&START_TIME='JD" << std::fixed << std::setprecision(6) << (epoch.mjd() + 2400000.5) << "'"
       << "&STOP_TIME='JD" << (epoch.mjd() + 2400000.6) << "'"
       << "&STEP_SIZE='1d'&OUT_UNITS='AU-D'";
    return ss.str();
}

std::string HorizonsClient::build_vectors_url(const std::string& target, const time::EpochTDB& epoch, const std::string& center) {
    std::stringstream ss;
    ss << config_.base_url << "?format=json&COMMAND='" << url_encode(target) << "'"
       << "&OBJ_DATA='NO'&MAKE_EPHEM='YES'&EPHEM_TYPE='VECTORS'&CENTER='" << center << "'"
       << "&START_TIME='JD" << std::fixed << std::setprecision(6) << (epoch.mjd() + 2400000.5) << "'"
       << "&STOP_TIME='JD" << (epoch.mjd() + 2400000.6) << "'"
       << "&STEP_SIZE='1d'&CSV_FORMAT='NO'&REF_PLANE='FRAME'&OUT_UNITS='KM-S'&VEC_TABLE='3'";
    return ss.str();
}

std::string HorizonsClient::build_observer_url(const std::string& target, const time::EpochTDB& epoch, const std::string& observer) {
    std::stringstream ss;
    auto epoch_utc = time::to_utc(epoch);
    ss << config_.base_url << "?format=json&COMMAND='" << url_encode(target) << "'"
       << "&OBJ_DATA='NO'&MAKE_EPHEM='YES'&EPHEM_TYPE='OBSERVER'&CENTER='" << observer << "'"
       << "&START_TIME='JD" << std::fixed << std::setprecision(6) << (epoch_utc.mjd() + 2400000.5) << "'"
       << "&STOP_TIME='JD" << (epoch_utc.mjd() + 2400000.6) << "'"
       << "&STEP_SIZE='1d'&QUANTITIES='1,20'&EXTRA_PREC='YES'";
    return ss.str();
}

} // namespace astdyn::io
