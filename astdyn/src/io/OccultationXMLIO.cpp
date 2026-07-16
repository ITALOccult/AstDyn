/**
 * @file OccultationXMLIO.cpp
 * @brief Reader/writer for Occult4's "occelmnt" XML format.
 *
 * Field order follows Occult4's published specification; see OccultationEvent.hpp.
 */
#include "astdyn/io/OccultationXMLIO.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace astdyn::io {
namespace {

/// Value of <tag>...</tag> within @p src, or an empty string when absent.
std::string extract_tag(const std::string& src, const std::string& tag) {
    const std::string open = "<" + tag + ">";
    const std::string close = "</" + tag + ">";
    const auto a = src.find(open);
    if (a == std::string::npos) return {};
    const auto b = src.find(close, a);
    if (b == std::string::npos) return {};
    return src.substr(a + open.size(), b - a - open.size());
}

std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return {};
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

/// Fields are comma-separated and may legitimately be empty (e.g. the K2 flag).
std::vector<std::string> split(const std::string& csv) {
    std::vector<std::string> out;
    std::stringstream ss(csv);
    std::string item;
    while (std::getline(ss, item, ',')) out.push_back(trim(item));
    return out;
}

/// Field @p i as a double, or @p fallback when missing or unparsable.
double num(const std::vector<std::string>& f, size_t i, double fallback = 0.0) {
    if (i >= f.size() || f[i].empty()) return fallback;
    try { return std::stod(f[i]); } catch (...) { return fallback; }
}

int inum(const std::vector<std::string>& f, size_t i, int fallback = 0) {
    if (i >= f.size() || f[i].empty()) return fallback;
    try { return std::stoi(f[i]); } catch (...) { return fallback; }
}

std::string str(const std::vector<std::string>& f, size_t i) {
    return i < f.size() ? f[i] : std::string{};
}

/// Occult4 writes the JWST flag as True/False; accept 1/0 as the spec allows.
bool flag(const std::vector<std::string>& f, size_t i) {
    const std::string v = str(f, i);
    return v == "True" || v == "true" || v == "1";
}

} // namespace

OccultationEvent OccultationXMLIO::parse_event_node(const std::string& xml) {
    OccultationEvent ev;

    {   // <Elements> source, duration, y, m, d, UT, x, y, dX, dY, d2X, d2Y, d3X, d3Y
        const auto f = split(extract_tag(xml, "Elements"));
        ev.elements_source = str(f, 0);
        ev.duration_s   = num(f, 1);
        ev.year         = inum(f, 2);
        ev.month        = inum(f, 3);
        ev.day          = inum(f, 4);
        ev.ut_closest_h = num(f, 5);
        ev.x  = num(f, 6);   ev.y  = num(f, 7);
        ev.dx = num(f, 8);   ev.dy = num(f, 9);
        ev.d2x = num(f, 10); ev.d2y = num(f, 11);
        ev.d3x = num(f, 12); ev.d3y = num(f, 13);
    }
    {   // <Earth> SubstellarLong, SubstellarLat, SubsolarLong, SubsolarLat, JWST
        const auto f = split(extract_tag(xml, "Earth"));
        ev.substellar_lon_deg = num(f, 0);
        ev.substellar_lat_deg = num(f, 1);
        ev.subsolar_lon_deg   = num(f, 2);
        ev.subsolar_lat_deg   = num(f, 3);
        ev.jwst               = flag(f, 4);
    }
    {   // <Star>
        const auto f = split(extract_tag(xml, "Star"));
        ev.star_id            = str(f, 0);
        ev.star_ra_h          = num(f, 1);
        ev.star_dec_deg       = num(f, 2);
        ev.mag_b              = num(f, 3, 99.0);
        ev.mag_v              = num(f, 4, 99.0);
        ev.mag_r              = num(f, 5, 99.0);
        ev.star_diameter_mas  = num(f, 6);
        ev.double_star_code   = inum(f, 7);
        ev.k2_flag            = str(f, 8);
        ev.star_app_ra_h      = num(f, 9);
        ev.star_app_dec_deg   = num(f, 10);
        ev.mag_drop_v         = num(f, 11);
        ev.mag_drop_r         = num(f, 12);
        ev.mag_drops_adjusted = inum(f, 13);
        ev.bright_nearby_count = inum(f, 14, -1);
        ev.total_nearby_count  = inum(f, 15, -1);
    }
    {   // <Object>
        const auto f = split(extract_tag(xml, "Object"));
        ev.object_number = str(f, 0);
        ev.object_name   = str(f, 1);
        ev.object_mag    = num(f, 2);
        ev.diameter_km   = num(f, 3);
        ev.distance_au   = num(f, 4);
        ev.n_rings       = inum(f, 5);
        ev.n_moons       = inum(f, 6);
        ev.d_ra_s_hr     = num(f, 7);
        ev.d_dec_as_hr   = num(f, 8);
        ev.taxonomy      = str(f, 9);
        ev.diameter_uncertainty_km = num(f, 10);
        ev.moon_in_planet_shadow   = inum(f, 11);
        ev.mag_v_asteroid = num(f, 12);
        ev.mag_r_asteroid = num(f, 13);
    }
    {   // <Orbit> equinox, MA, yr, month, day, peri, node, i, e, a, q, H0, LogR, G
        const auto f = split(extract_tag(xml, "Orbit"));
        if (!f.empty()) {
            ev.equinox          = num(f, 0);
            ev.mean_anomaly_deg = num(f, 1);
            ev.epoch_year       = inum(f, 2);
            ev.epoch_month      = inum(f, 3);
            ev.epoch_day        = inum(f, 4);
            ev.peri_deg         = num(f, 5);
            ev.node_deg         = num(f, 6);
            ev.inclination_deg  = num(f, 7);
            ev.eccentricity     = num(f, 8);
            ev.semi_major_axis_au = num(f, 9);
            ev.perihelion_au    = num(f, 10);
            ev.h0               = num(f, 11);
            ev.coeff_log_r      = num(f, 12);
            ev.g_param          = num(f, 13, 0.15);
        }
    }
    {   // <Errors>
        const auto f = split(extract_tag(xml, "Errors"));
        ev.err_path_widths = num(f, 0);
        ev.err_major_as    = num(f, 1);
        ev.err_minor_as    = num(f, 2);
        ev.err_pa_deg      = num(f, 3);
        ev.err_1sigma_as   = num(f, 4);
        ev.error_basis     = str(f, 5);
        ev.reliability     = num(f, 6, -1.0);
        ev.duplicate_source = inum(f, 7, -1);
        ev.non_gaia_pm      = inum(f, 8, -1);
        ev.pm_added_from_ucac4 = inum(f, 9, -1);
    }
    {   // <ID> id, MJD of the prediction calculation
        const auto f = split(extract_tag(xml, "ID"));
        ev.event_id       = str(f, 0);
        ev.prediction_mjd = num(f, 1);
    }
    return ev;
}

std::string OccultationXMLIO::format_event_node(const OccultationEvent& ev) {
    std::stringstream ss;
    ss << std::fixed;
    ss << "  <Event>\n";

    ss << std::setprecision(7);
    ss << "    <Elements>" << ev.elements_source << ","
       << ev.duration_s << "," << ev.year << "," << ev.month << "," << ev.day << ","
       << ev.ut_closest_h << "," << ev.x << "," << ev.y << ","
       << ev.dx << "," << ev.dy << "," << ev.d2x << "," << ev.d2y << ","
       << ev.d3x << "," << ev.d3y << "</Elements>\n";

    ss << std::setprecision(4);
    ss << "    <Earth>" << ev.substellar_lon_deg << "," << ev.substellar_lat_deg << ","
       << ev.subsolar_lon_deg << "," << ev.subsolar_lat_deg << ","
       << (ev.jwst ? "True" : "False") << "</Earth>\n";

    ss << "    <Star>" << ev.star_id << ","
       << std::setprecision(8) << ev.star_ra_h << "," << ev.star_dec_deg << ","
       << std::setprecision(2) << ev.mag_b << "," << ev.mag_v << "," << ev.mag_r << ","
       << std::setprecision(3) << ev.star_diameter_mas << ","
       << ev.double_star_code << "," << ev.k2_flag << ","
       << std::setprecision(8) << ev.star_app_ra_h << "," << ev.star_app_dec_deg << ","
       << std::setprecision(2) << ev.mag_drop_v << "," << ev.mag_drop_r << ","
       << ev.mag_drops_adjusted << "," << ev.bright_nearby_count << ","
       << ev.total_nearby_count << "</Star>\n";

    ss << "    <Object>" << ev.object_number << "," << ev.object_name << ","
       << std::setprecision(2) << ev.object_mag << ","
       << std::setprecision(3) << ev.diameter_km << ","
       << std::setprecision(4) << ev.distance_au << ","
       << ev.n_rings << "," << ev.n_moons << ","
       << std::setprecision(3) << ev.d_ra_s_hr << "," << ev.d_dec_as_hr << ","
       << ev.taxonomy << "," << ev.diameter_uncertainty_km << ","
       << ev.moon_in_planet_shadow << ","
       << std::setprecision(2) << ev.mag_v_asteroid << "," << ev.mag_r_asteroid
       << "</Object>\n";

    ss << std::setprecision(4);
    ss << "    <Orbit>" << ev.equinox << "," << ev.mean_anomaly_deg << ","
       << ev.epoch_year << "," << ev.epoch_month << "," << ev.epoch_day << ","
       << ev.peri_deg << "," << ev.node_deg << "," << ev.inclination_deg << ","
       << std::setprecision(5) << ev.eccentricity << "," << ev.semi_major_axis_au << ","
       << ev.perihelion_au << ","
       << std::setprecision(2) << ev.h0 << "," << std::setprecision(1) << ev.coeff_log_r
       << "," << std::setprecision(2) << ev.g_param << "</Orbit>\n";

    ss << std::setprecision(4);
    ss << "    <Errors>" << std::setprecision(3) << ev.err_path_widths << ","
       << std::setprecision(4) << ev.err_major_as << "," << ev.err_minor_as << ","
       << std::setprecision(0) << ev.err_pa_deg << ","
       << std::setprecision(4) << ev.err_1sigma_as << ","
       << ev.error_basis << "," << std::setprecision(2) << ev.reliability << ","
       << ev.duplicate_source << "," << ev.non_gaia_pm << ","
       << ev.pm_added_from_ucac4 << "</Errors>\n";

    ss << "    <ID>" << ev.event_id << "," << std::setprecision(4)
       << ev.prediction_mjd << "</ID>\n";

    ss << "  </Event>\n";
    return ss.str();
}

std::vector<std::string> OccultationXMLIO::split_csv(const std::string& csv) {
    return split(csv);
}

std::vector<OccultationEvent> OccultationXMLIO::read_string(const std::string& xml) {
    std::vector<OccultationEvent> events;
    size_t pos = 0;
    const std::string open = "<Event>", close = "</Event>";
    while ((pos = xml.find(open, pos)) != std::string::npos) {
        const auto end = xml.find(close, pos);
        if (end == std::string::npos) break;
        events.push_back(parse_event_node(xml.substr(pos, end - pos + close.size())));
        pos = end + close.size();
    }
    return events;
}

std::vector<OccultationEvent> OccultationXMLIO::read_file(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) return {};
    std::stringstream buf;
    buf << in.rdbuf();
    return read_string(buf.str());
}

std::string OccultationXMLIO::write_string(const std::vector<OccultationEvent>& events) {
    std::stringstream ss;
    ss << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    ss << "<Occultations>\n";
    for (const auto& ev : events) ss << format_event_node(ev);
    ss << "</Occultations>\n";
    return ss.str();
}

bool OccultationXMLIO::write_file(const std::vector<OccultationEvent>& events,
                                  const std::string& filename) {
    std::ofstream out(filename);
    if (!out) return false;
    out << write_string(events);
    return out.good();
}

} // namespace astdyn::io
