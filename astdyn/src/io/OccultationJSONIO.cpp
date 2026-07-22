/**
 * @file OccultationJSONIO.cpp
 * @brief Implementazione della serializzazione JSON degli eventi.
 */
#include "astdyn/io/OccultationJSONIO.hpp"

#include <fstream>
#include <nlohmann/json.hpp>

namespace astdyn::io {

using json = nlohmann::json;

static json event_to_json(const OccultationEvent& ev) {
    json j;

    // ---- <ID> : identificativo dell'evento -------------------------------
    j["id"] = {
        {"event_id", ev.event_id},
        {"prediction_mjd", ev.prediction_mjd},
    };

    // ---- <Object> : l'asteroide (number esplicito per i positivi) ---------
    j["object"] = {
        {"number", ev.object_number},
        {"name", ev.object_name},
        {"mag", ev.object_mag},
        {"diameter_km", ev.diameter_km},
        {"diameter_uncertainty_km", ev.diameter_uncertainty_km},
        {"distance_au", ev.distance_au},
        {"n_rings", ev.n_rings},
        {"n_moons", ev.n_moons},
        {"d_ra_s_hr", ev.d_ra_s_hr},
        {"d_dec_as_hr", ev.d_dec_as_hr},
        {"taxonomy", ev.taxonomy},
        {"moon_in_planet_shadow", ev.moon_in_planet_shadow},
        {"mag_v", ev.mag_v_asteroid},
        {"mag_r", ev.mag_r_asteroid},
    };

    // ---- <Star> -----------------------------------------------------------
    j["star"] = {
        {"id", ev.star_id},
        {"ra_h", ev.star_ra_h},
        {"dec_deg", ev.star_dec_deg},
        {"mag_b", ev.mag_b},
        {"mag_v", ev.mag_v},
        {"mag_r", ev.mag_r},
        {"diameter_mas", ev.star_diameter_mas},
        {"double_star_code", ev.double_star_code},
        {"k2_flag", ev.k2_flag},
        {"app_ra_h", ev.star_app_ra_h},
        {"app_dec_deg", ev.star_app_dec_deg},
        {"mag_drop_v", ev.mag_drop_v},
        {"mag_drop_r", ev.mag_drop_r},
        {"mag_drops_adjusted", ev.mag_drops_adjusted},
        {"bright_nearby_count", ev.bright_nearby_count},
        {"total_nearby_count", ev.total_nearby_count},
    };

    // ---- <Event> timing + <Elements> besseliani --------------------------
    j["event"] = {
        {"year", ev.year}, {"month", ev.month}, {"day", ev.day},
        {"ut_closest_h", ev.ut_closest_h},
        {"duration_s", ev.duration_s},
    };
    j["elements"] = {
        {"source", ev.elements_source},
        {"x", ev.x}, {"y", ev.y},
        {"dx", ev.dx}, {"dy", ev.dy},
        {"d2x", ev.d2x}, {"d2y", ev.d2y},
        {"d3x", ev.d3x}, {"d3y", ev.d3y},
    };

    // ---- <Earth> : punti substellare/subsolare ---------------------------
    j["earth"] = {
        {"substellar_lon_deg", ev.substellar_lon_deg},
        {"substellar_lat_deg", ev.substellar_lat_deg},
        {"subsolar_lon_deg", ev.subsolar_lon_deg},
        {"subsolar_lat_deg", ev.subsolar_lat_deg},
        {"jwst", ev.jwst},
    };

    // ---- <Orbit> ----------------------------------------------------------
    j["orbit"] = {
        {"equinox", ev.equinox},
        {"mean_anomaly_deg", ev.mean_anomaly_deg},
        {"epoch_year", ev.epoch_year},
        {"epoch_month", ev.epoch_month},
        {"epoch_day", ev.epoch_day},
        {"peri_deg", ev.peri_deg},
        {"node_deg", ev.node_deg},
        {"inclination_deg", ev.inclination_deg},
        {"eccentricity", ev.eccentricity},
        {"semi_major_axis_au", ev.semi_major_axis_au},
        {"perihelion_au", ev.perihelion_au},
        {"h0", ev.h0},
        {"coeff_log_r", ev.coeff_log_r},
        {"g_param", ev.g_param},
    };

    // ---- <Errors> : ellisse di incertezza (covariante) -------------------
    j["errors"] = {
        {"path_widths", ev.err_path_widths},
        {"major_as", ev.err_major_as},
        {"minor_as", ev.err_minor_as},
        {"pa_deg", ev.err_pa_deg},
        {"sigma1_as", ev.err_1sigma_as},
        {"basis", ev.error_basis},
        {"reliability", ev.reliability},
        {"duplicate_source", ev.duplicate_source},
        {"non_gaia_pm", ev.non_gaia_pm},
        {"pm_added_from_ucac4", ev.pm_added_from_ucac4},
    };

    return j;
}

std::string OccultationJSONIO::write_string(const std::vector<OccultationEvent>& events) {
    json root;
    root["events"] = json::array();
    for (const auto& ev : events) {
        root["events"].push_back(event_to_json(ev));
    }
    root["count"] = static_cast<int>(events.size());
    return root.dump(2);
}

bool OccultationJSONIO::write_file(const std::vector<OccultationEvent>& events,
                                   const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) return false;
    out << write_string(events);
    return out.good();
}

} // namespace astdyn::io
