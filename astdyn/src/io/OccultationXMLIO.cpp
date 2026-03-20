/**
 * @file OccultationXMLIO.cpp
 * @brief Implementation of Occult4 XML I/O.
 */

#include "astdyn/io/OccultationXMLIO.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iostream>

namespace astdyn::io {

std::vector<OccultationEvent> OccultationXMLIO::read_string(const std::string& xml_content) {
    std::vector<OccultationEvent> events;
    size_t pos = 0;
    while ((pos = xml_content.find("<Event>", pos)) != std::string::npos) {
        size_t end_pos = xml_content.find("</Event>", pos);
        if (end_pos == std::string::npos) break;
        
        std::string event_xml = xml_content.substr(pos, end_pos - pos + 8);
        events.push_back(parse_event_node(event_xml));
        pos = end_pos + 8;
    }
    return events;
}

std::vector<OccultationEvent> OccultationXMLIO::read_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return {};
    std::stringstream buffer;
    buffer << file.rdbuf();
    return read_string(buffer.str());
}

OccultationEvent OccultationXMLIO::parse_event_node(const std::string& event_xml) {
    OccultationEvent event = {};

    auto extract_tag = [&](const std::string& tag) {
        std::string start_tag = "<" + tag + ">";
        std::string end_tag = "</" + tag + ">";
        size_t start = event_xml.find(start_tag);
        size_t end = event_xml.find(end_tag);
        if (start != std::string::npos && end != std::string::npos) {
            return event_xml.substr(start + start_tag.length(), end - (start + start_tag.length()));
        }
        return std::string("");
    };

    // 1. Elements
    std::string elements_str = extract_tag("Elements");
    auto elements_parts = split_csv(elements_str);
    if (!elements_parts.empty()) {
        event.elements_source = elements_parts[0];
        for (size_t i = 1; i < elements_parts.size(); ++i) {
            try { 
                if (!elements_parts[i].empty())
                    event.elements_data.push_back(std::stod(elements_parts[i])); 
            } catch(...) {}
        }
    }

    // 2. Earth
    std::string earth_str = extract_tag("Earth");
    auto earth_parts = split_csv(earth_str);
    if (earth_parts.size() >= 5) {
        try {
            event.longitude_deg = std::stod(earth_parts[0]);
            event.latitude_deg = std::stod(earth_parts[1]);
            event.alt_or_other = std::stod(earth_parts[2]);
            event.max_duration_sec = std::stod(earth_parts[3]);
            event.is_daylight = (earth_parts[4] == "True" || earth_parts[4] == "true" || earth_parts[4] == "1");
        } catch(...) {}
    }

    // 3. Star
    std::string star_str = extract_tag("Star");
    auto star_parts = split_csv(star_str);
    if (star_parts.size() >= 11) {
        event.star_catalog_id = star_parts[0];
        try {
            event.ra_cat_h = std::stod(star_parts[1]);
            event.dec_cat_deg = std::stod(star_parts[2]);
            event.mag_v = std::stod(star_parts[3]);
            event.mag_r = std::stod(star_parts[4]);
            event.mag_k = std::stod(star_parts[5]);
            if (star_parts.size() > 6) event.pm_ra_as_yr = std::stod(star_parts[6]);
            if (star_parts.size() > 7) event.pm_dec_as_yr = std::stod(star_parts[7]);
            if (star_parts.size() > 8 && !star_parts[8].empty()) event.parallax_as = std::stod(star_parts[8]);
            
            if (star_parts.size() >= 11) {
                event.ra_event_h = std::stod(star_parts[9]);
                event.dec_event_deg = std::stod(star_parts[10]);
                for (size_t i = 11; i < star_parts.size(); ++i) {
                    if (!star_parts[i].empty())
                        event.star_extra_data.push_back(std::stod(star_parts[i]));
                }
            }
        } catch(...) {}
    }

    // 4. Object
    std::string object_str = extract_tag("Object");
    auto object_parts = split_csv(object_str);
    if (object_parts.size() >= 10) {
        try {
            event.object_number = std::stoi(object_parts[0]);
            event.object_name = object_parts[1];
            event.h_mag = std::stod(object_parts[2]);
            event.diameter_km = std::stod(object_parts[3]);
            event.apparent_rate_arcsec_hr = std::stod(object_parts[4]);
            if (object_parts.size() > 5 && !object_parts[5].empty()) event.g_param = std::stod(object_parts[5]);
            
            if (object_parts.size() >= 10) {
                event.object_type = object_parts[9];
                for (size_t i = 10; i < object_parts.size(); ++i) {
                     if (object_parts[i].empty()) continue;
                     event.object_extra_data.push_back(std::stod(object_parts[i]));
                }
            }
        } catch(...) {}
    }

    // 5. Orbit
    std::string orbit_str = extract_tag("Orbit");
    auto orbit_parts = split_csv(orbit_str);
    if (orbit_parts.size() >= 10) {
        try {
            event.orbit_type = std::stod(orbit_parts[0]);
            event.ra_node_approx = std::stod(orbit_parts[1]);
            event.epoch_year = std::stoi(orbit_parts[2]);
            event.epoch_month = std::stoi(orbit_parts[3]);
            event.epoch_day = std::stoi(orbit_parts[4]);
            event.mean_anomaly_deg = std::stod(orbit_parts[5]);
            event.node_deg = std::stod(orbit_parts[6]);
            event.inclination_deg = std::stod(orbit_parts[7]);
            event.eccentricity = std::stod(orbit_parts[8]);
            event.semi_major_axis_au = std::stod(orbit_parts[9]);
            for (size_t i = 10; i < orbit_parts.size(); ++i) {
                 if (orbit_parts[i].empty()) continue;
                 event.orbit_extra_data.push_back(std::stod(orbit_parts[i]));
            }
        } catch(...) {}
    }

    // 6. Errors
    std::string errors_str = extract_tag("Errors");
    auto errors_parts = split_csv(errors_str);
    if (errors_parts.size() >= 6) {
        try {
            event.combined_error_arcsec = std::stod(errors_parts[0]);
            if (errors_parts.size() > 1 && !errors_parts[1].empty()) event.star_error_ra_as = std::stod(errors_parts[1]);
            if (errors_parts.size() > 2 && !errors_parts[2].empty()) event.star_error_dec_as = std::stod(errors_parts[2]);
            if (errors_parts.size() > 3 && !errors_parts[3].empty()) event.object_error_ra_as = std::stod(errors_parts[3]);
            if (errors_parts.size() > 4 && !errors_parts[4].empty()) event.object_error_dec_as = std::stod(errors_parts[4]);
            
            if (errors_parts.size() > 5) event.uncertainty_method = errors_parts[5];
            
            for (size_t i = 6; i < errors_parts.size(); ++i) {
                if (!errors_parts[i].empty())
                    event.error_extra_data.push_back(std::stod(errors_parts[i]));
            }
        } catch(...) {}
    }

    // 7. ID
    std::string id_str = extract_tag("ID");
    auto id_parts = split_csv(id_str);
    if (id_parts.size() >= 2) {
        event.event_id = id_parts[0];
        try { event.mjd = std::stod(id_parts[1]); } catch(...) {}
    }

    return event;
}

std::vector<std::string> OccultationXMLIO::split_csv(const std::string& csv) {
    std::vector<std::string> parts;
    std::stringstream ss(csv);
    std::string part;
    while (std::getline(ss, part, ',')) {
        // Trim whitespace
        part.erase(0, part.find_first_not_of(" \t\n\r"));
        auto last = part.find_last_not_of(" \t\n\r");
        if (last != std::string::npos)
            part.erase(last + 1);
        parts.push_back(part);
    }
    return parts;
}

std::string OccultationXMLIO::write_string(const std::vector<OccultationEvent>& events) {
    std::stringstream ss;
    ss << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    ss << "<Occultations>\n";
    for (const auto& event : events) {
        ss << format_event_node(event);
    }
    ss << "</Occultations>\n";
    return ss.str();
}

std::string OccultationXMLIO::format_event_node(const OccultationEvent& event) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(7);
    ss << "  <Event>\n";
    
    // Elements
    ss << "    <Elements>" << event.elements_source;
    for (double d : event.elements_data) ss << "," << d;
    ss << "</Elements>\n";

    // Earth
    ss << "    <Earth>" << event.longitude_deg << "," << event.latitude_deg << "," 
       << event.alt_or_other << "," << event.max_duration_sec << "," 
       << (event.is_daylight ? "True" : "False") << "</Earth>\n";

    // Star
    ss << "    <Star>" << event.star_catalog_id << "," << event.ra_cat_h << "," << event.dec_cat_deg << ","
       << event.mag_v << "," << event.mag_r << "," << event.mag_k << ","
       << event.pm_ra_as_yr << "," << event.pm_dec_as_yr << "," << event.parallax_as << ","
       << event.ra_event_h << "," << event.dec_event_deg;
    for (double d : event.star_extra_data) ss << "," << d;
    ss << "</Star>\n";

    // Object
    ss << "    <Object>" << event.object_number << "," << event.object_name << ","
       << event.h_mag << "," << event.diameter_km << "," << event.apparent_rate_arcsec_hr << ","
       << event.g_param << ",0,0,0," << event.object_type;
    for (double d : event.object_extra_data) ss << "," << d;
    ss << "</Object>\n";

    // Orbit
    ss << "    <Orbit>" << event.orbit_type << "," << event.ra_node_approx << ","
       << event.epoch_year << "," << event.epoch_month << "," << event.epoch_day << ","
       << event.mean_anomaly_deg << "," << event.node_deg << "," << event.inclination_deg << ","
       << event.eccentricity << "," << event.semi_major_axis_au;
    for (double d : event.orbit_extra_data) ss << "," << d;
    ss << "</Orbit>\n";

    // Errors
    ss << "    <Errors>" << event.combined_error_arcsec << ","
       << event.star_error_ra_as << "," << event.star_error_dec_as << ","
       << event.object_error_ra_as << "," << event.object_error_dec_as << ","
       << event.uncertainty_method;
    for (double d : event.error_extra_data) ss << "," << d;
    ss << "</Errors>\n";

    // ID
    ss << "    <ID>" << event.event_id << "," << event.mjd << "</ID>\n";

    ss << "  </Event>\n";
    return ss.str();
}

bool OccultationXMLIO::write_file(const std::vector<OccultationEvent>& events, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) return false;
    file << write_string(events);
    return true;
}

} // namespace astdyn::io
