#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <mutex>

namespace astdyn::io {

static std::once_flag curl_init_flag;
static void global_curl_init() {
    curl_global_init(CURL_GLOBAL_DEFAULT);
}

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

static std::string format_command(const std::string& target) {
    return target;
}

HorizonsClient::HorizonsClient(const Config& config) : config_(config) {
    std::call_once(curl_init_flag, global_curl_init);
}

HorizonsClient::~HorizonsClient() {
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

std::expected<physics::KeplerianStateTyped<core::ECLIPJ2000>, HorizonsError> HorizonsClient::query_elements(const std::string& target, const time::EpochTDB& epoch) {
    auto raw = fetch_url(build_elements_url(target, epoch));
    if (!raw) return std::unexpected(raw.error());
    try {
        auto res = json::parse(*raw)["result"].get<std::string>();
        size_t start = res.find("$$SOE");
        if (start == std::string::npos) return std::unexpected(HorizonsError::TargetNotFound);
        std::string block = res.substr(start);
        
        return physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
            epoch, 
            parse_token(block, "A ="), parse_token(block, "EC="),
            parse_token(block, "IN="), parse_token(block, "OM="),
            parse_token(block, "W ="), parse_token(block, "MA="),
            physics::GravitationalParameter::sun()
        );
    } catch (...) { return std::unexpected(HorizonsError::ParsingError); }
}

double HorizonsClient::parse_token(const std::string& text, const std::string& token) {
    size_t pos = text.find(token);
    if (pos == std::string::npos) throw std::runtime_error("Horizons token not found");
    return std::stod(text.substr(text.find_first_of("0123456789+-.", pos + token.length())));
}

std::expected<physics::CartesianStateTyped<core::GCRF>, HorizonsError> HorizonsClient::query_vectors(const std::string& target, const time::EpochTDB& epoch, const std::string& center) {
    auto raw = fetch_url(build_vectors_url(target, epoch, center));
    if (!raw) return std::unexpected(raw.error());
    try {
        auto res = json::parse(*raw)["result"].get<std::string>();
        size_t start = res.find("$$SOE");
        if (start == std::string::npos) return std::unexpected(HorizonsError::TargetNotFound);
        std::string block = res.substr(start);
        double gm = (center.find("399") != std::string::npos) ? constants::GM_EARTH : constants::GM_SUN;
        return physics::CartesianStateTyped<core::GCRF>::from_si(epoch,
            parse_token(block, " X =") * 1000.0, parse_token(block, " Y =") * 1000.0, parse_token(block, " Z =") * 1000.0,
            parse_token(block, " VX=") * 1000.0, parse_token(block, " VY=") * 1000.0, parse_token(block, " VZ=") * 1000.0,
            gm * 1e9);
    } catch (...) { return std::unexpected(HorizonsError::ParsingError); }
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

        return astrometry::AstrometricObservation(
            astrometry::RightAscension::from_rad(ra),
            astrometry::Declination::from_rad(dec),
            physics::Distance::from_au(delta)
        );
    } catch (...) {
        return std::unexpected(HorizonsError::ParsingError);
    }
}

// --- URL Builders ---

std::string HorizonsClient::build_elements_url(const std::string& target, const time::EpochTDB& epoch) {
    std::stringstream ss;
    ss << config_.base_url << "?format=json&COMMAND='" << url_encode(format_command(target)) << "'"
       << "&OBJ_DATA='NO'&MAKE_EPHEM='YES'&EPHEM_TYPE='ELEMENTS'&CENTER='500@10'"
       << "&START_TIME='JD" << std::fixed << std::setprecision(6) << (epoch.mjd() + 2400000.5) << "'"
       << "&STOP_TIME='JD" << (epoch.mjd() + 2400000.6) << "'"
       << "&STEP_SIZE='1d'&OUT_UNITS='AU-D'";
    return ss.str();
}

std::string HorizonsClient::build_vectors_url(const std::string& target, const time::EpochTDB& epoch, const std::string& center) {
    std::stringstream ss;
    ss << config_.base_url << "?format=json&COMMAND='" << url_encode(format_command(target)) << "'"
       << "&OBJ_DATA='YES'&MAKE_EPHEM='YES'&EPHEM_TYPE='VECTORS'"
       << "&CENTER='" << url_encode(center) << "'"
       << "&START_TIME='JD" << std::fixed << std::setprecision(6) << (epoch.mjd() + 2400000.5) << "'"
       << "&STOP_TIME='JD" << std::fixed << std::setprecision(6) << (epoch.mjd() + 2400000.6) << "'"
       << "&STEP_SIZE='1d'&REF_PLANE='FRAME'&REF_SYSTEM='ICRF'&CSV_FORMAT='NO'&VEC_TABLE='2'";
    
    std::cout << "[Horizons] Query URL: " << ss.str() << "\n";
    return ss.str();
}

std::string HorizonsClient::build_observer_url(const std::string& target, const time::EpochTDB& epoch, const std::string& observer) {
    std::stringstream ss;
    auto epoch_utc = time::to_utc(epoch);
    ss << config_.base_url << "?format=json&COMMAND='" << url_encode(format_command(target)) << "'"
       << "&OBJ_DATA='NO'&MAKE_EPHEM='YES'&EPHEM_TYPE='OBSERVER'&CENTER='" << observer << "'"
       << "&START_TIME='JD" << std::fixed << std::setprecision(6) << (epoch_utc.mjd() + 2400000.5) << "'"
       << "&STOP_TIME='JD" << (epoch_utc.mjd() + 2400000.6) << "'"
       << "&STEP_SIZE='1d'&QUANTITIES='1,20'&EXTRA_PREC='YES'";
    return ss.str();
}

} // namespace astdyn::io

namespace astdyn::io {

std::expected<PhysicalProperties, HorizonsError> 
HorizonsClient::query_physical_properties(const std::string& target) {
    std::stringstream ss;
    ss << config_.base_url << "?format=json&COMMAND='" << url_encode(format_command(target)) << "'"
       << "&OBJ_DATA='YES'&MAKE_EPHEM='NO'";
    
    auto raw = fetch_url(ss.str());
    if (!raw) return std::unexpected(raw.error());

    try {
        auto data = json::parse(*raw);
        std::string result = data["result"];
        PhysicalProperties props;

        // Use regex for robust finding
        std::regex h_regex(R"(H\s*=\s*([0-9.]+))");
        std::regex diam_regex(R"(Radius\s*\(km\)\s*=\s*([0-9.]+))");
        std::regex diam_regex2(R"(Diameter\s*\(km\)\s*=\s*([0-9.]+))");

        std::smatch match;
        if (std::regex_search(result, match, h_regex)) {
            props.h_mag = std::stod(match[1]);
        }
        if (std::regex_search(result, match, diam_regex)) {
            props.diameter_km = std::stod(match[1]) * 2.0; // Radius to Diameter
        } else if (std::regex_search(result, match, diam_regex2)) {
            props.diameter_km = std::stod(match[1]);
        }

        return props;
    } catch (...) {
        return std::unexpected(HorizonsError::ParsingError);
    }
}

} // namespace astdyn::io
