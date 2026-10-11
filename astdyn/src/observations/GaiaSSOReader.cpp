/**
 * @file GaiaSSOReader.cpp
 * @brief Gaia DR3/FPR sso_observation CSV reader.
 */
#include "astdyn/observations/GaiaSSOReader.hpp"
#include "astdyn/core/Constants.hpp"
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace astdyn::observations {

namespace {

/// epoch_utc counts days from J2010.0 UTC, 2010-01-01T00:00.
constexpr double JD_J2010_UTC = 2455197.5;
constexpr double MAS_TO_ARCSEC = 1.0e-3;
const std::string GAIA_OBSERVATORY_CODE = "258";

using Record = std::map<std::string, std::string>;

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    for (std::string f; std::getline(ss, f, ',');) fields.push_back(f);
    return fields;
}

double number(const Record& r, const std::string& column) {
    const auto it = r.find(column);
    if (it == r.end()) throw std::runtime_error("Gaia SSO: missing column " + column);
    return std::stod(it->second);
}

/// Covariance of (alpha*, delta) for one error component [arcsec^2].
Eigen::Matrix2d component_covariance(const Record& r, const std::string& kind) {
    const double sigma_ra = number(r, "ra_error_" + kind) * MAS_TO_ARCSEC;
    const double sigma_dec = number(r, "dec_error_" + kind) * MAS_TO_ARCSEC;
    const double rho = number(r, "ra_dec_correlation_" + kind);
    Eigen::Matrix2d c;
    c << sigma_ra * sigma_ra, rho * sigma_ra * sigma_dec,
         rho * sigma_ra * sigma_dec, sigma_dec * sigma_dec;
    return c;
}

OpticalObservation to_observation(const Record& r, int observations_in_transit) {
    const double au_m = constants::AU * 1000.0;
    OpticalObservation obs;
    obs.observatory_code = GAIA_OBSERVATORY_CODE;
    obs.object_designation = r.count("number_mp") ? r.at("number_mp") : std::string();
    obs.time = time::EpochUTC::from_jd(JD_J2010_UTC + number(r, "epoch_utc"));
    obs.ra = astrometry::RightAscension::from_deg(number(r, "ra"));
    obs.dec = astrometry::Declination::from_deg(number(r, "dec"));
    obs.space_observer_gcrf = math::Vector3<core::GCRF, physics::Distance>::from_si(
        number(r, "x_gaia_geocentric") * au_m, number(r, "y_gaia_geocentric") * au_m,
        number(r, "z_gaia_geocentric") * au_m);
    // The systematic error is common to the ~8 CCD observations of a transit:
    // counted once per observation it was trusted N times over. Inflating it by
    // N keeps the transit mean's systematic uncertainty whole (as Veres 2017
    // does for same-night batches).
    const Eigen::Matrix2d cov = observations_in_transit * component_covariance(r, "systematic")
                              + component_covariance(r, "random");
    obs.radec_covariance_arcsec2 = cov;
    obs.sigma_ra = astrometry::Angle::from_arcsec(std::sqrt(cov(0, 0)));
    obs.sigma_dec = astrometry::Angle::from_arcsec(std::sqrt(cov(1, 1)));
    return obs;
}

} // namespace

std::vector<OpticalObservation> read_gaia_sso_csv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Gaia SSO: cannot open " + path);
    std::string line;
    if (!std::getline(in, line)) return {};
    const auto header = split_csv(line);
    std::vector<Record> records;
    std::map<std::string, int> count_per_transit;
    while (std::getline(in, line)) {
        const auto fields = split_csv(line);
        Record r;
        for (size_t i = 0; i < header.size() && i < fields.size(); ++i) r[header[i]] = fields[i];
        if (r["is_rejected"] == "true") continue;
        ++count_per_transit[r["transit_id"]];
        records.push_back(std::move(r));
    }
    std::vector<OpticalObservation> out;
    for (const auto& r : records) out.push_back(to_observation(r, count_per_transit[r.at("transit_id")]));
    return out;
}

} // namespace astdyn::observations
