#include "astdyn/AstDyn.hpp"
#include "astdyn/io/OccultationXMLIO.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/astrometry/OccultationEvent.hpp"
#include "astdyn/astrometry/ChebyshevEphemerisManager.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <sstream>
#include "astdyn/time/TimeScale.hpp"

using namespace astdyn;
using namespace astdyn::astrometry;
using namespace astdyn::io;
namespace po = boost::program_options;

/**
 * @brief Helper to convert an internal candidate to an XML-formatted event structure.
 */
OccultationEvent candidate_to_event(const OccultationCandidate& cand, const std::string& ast_id, const physics::KeplerianStateTyped<core::ECLIPJ2000>& el) {
    OccultationEvent ev;
    ev.event_id = cand.params.star_id + "_" + std::to_string((int)cand.params.t_ca.mjd());
    ev.mjd = cand.params.t_ca.mjd();
    
    ev.star_catalog_id = cand.params.star_id;
    ev.ra_event_h = cand.star.ra.to_deg() / 15.0;
    ev.dec_event_deg = cand.star.dec.to_deg();
    
    ev.ra_cat_h = ev.ra_event_h; 
    ev.dec_cat_deg = ev.dec_event_deg;
    ev.pm_ra_as_yr = cand.star.pm_ra_cosdec.to_arcsec_yr(); 
    ev.pm_dec_as_yr = cand.star.pm_dec.to_arcsec_yr();
    ev.parallax_as = cand.star.parallax.to_arcsec();
    
    ev.mag_v = cand.star.g_mag;
    ev.mag_r = cand.star.rp_mag; 
    ev.mag_k = cand.star.g_mag;  
    
    ev.object_name = ast_id;
    try { ev.object_number = std::stoi(ast_id); } catch(...) { ev.object_number = 0; }
    ev.object_type = "Asteroid";
    ev.diameter_km = cand.params.cross_track_uncertainty.to_km(); 
    ev.h_mag = cand.params.star_mag; 
    ev.apparent_rate_arcsec_hr = cand.params.total_apparent_rate;
    
    ev.longitude_deg = cand.params.center_lon.to_deg();
    ev.latitude_deg = cand.params.center_lat.to_deg();
    ev.alt_or_other = 0.0; 
    ev.max_duration_sec = cand.params.max_duration.to_seconds();
    ev.is_daylight = cand.params.is_daylight;
    
    ev.combined_error_arcsec = cand.params.impact_parameter.to_km() / 150000000.0 * 206265.0; 
    ev.star_error_ra_as = cand.star.pmra_error_mas_yr / 1000.0;
    ev.star_error_dec_as = cand.star.pmdec_error_mas_yr / 1000.0;
    ev.uncertainty_method = "AstDyn-GaiaDR3-JplDE441";
    
    ev.semi_major_axis_au = el.a.to_au();
    ev.eccentricity = el.e;
    ev.inclination_deg = el.i.to_deg();
    ev.node_deg = el.node.to_deg();
    ev.mean_anomaly_deg = el.M.to_deg(); 
    ev.ra_node_approx = el.omega.to_deg();

    auto [year, month, day, frac] = time::mjd_to_calendar(el.epoch.mjd());
    ev.epoch_year = year;
    ev.epoch_month = month;
    ev.epoch_day = day;
    
    return ev;
}

std::vector<std::string> parse_asteroid_list(const std::string& input) {
    std::vector<std::string> ids;
    if (input.empty()) return ids;
    if (input[0] == '@') {
        std::ifstream file(input.substr(1));
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty()) ids.push_back(line);
        }
    } else {
        std::stringstream ss(input);
        std::string segment;
        while (std::getline(ss, segment, ',')) {
            if (!segment.empty()) ids.push_back(segment);
        }
    }
    return ids;
}

int main(int argc, char** argv) {
    po::options_description desc("ioccultcalc: AstDyn Occultation Tool CLI\nUsage: ioccultcalc --asteroid <num1,num2...> --jd-start <jd> --duration <days>");
    desc.add_options()
        ("help,h", "produce help message")
        ("conf", po::value<std::string>(), "JSON configuration file for AstDynEngine")
        ("asteroid", po::value<std::string>(), "asteroid designation list (comma-separated or @file)")
        ("jd-start", po::value<double>(), "start Julian Date for searching (TDB)")
        ("duration", po::value<double>()->default_value(1.0), "search duration in days")
        ("mag", po::value<double>()->default_value(15.0), "magnitude limit for stars")
        ("xml-output", po::value<std::string>(), "save found occultations to XML")
        ("kml", po::value<std::string>(), "save first match to KML")
    ;

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (...) { return 1; }

    if (vm.count("help") || !vm.count("asteroid") || !vm.count("jd-start")) {
        std::cout << desc << "\n";
        return 0;
    }

    // --- 1. Engine Setup ---
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    AstDynEngine engine;
    if (vm.count("conf")) {
        engine.load_config(vm["conf"].as<std::string>());
    } else {
        AstDynConfig cfg;
        cfg.ephemeris_file = bsp_path;
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.verbose = false;
        engine.set_config(cfg);
    }
    
    try {
        ephemeris::PlanetaryEphemeris::setProvider(std::make_shared<ephemeris::DE441Provider>(engine.config().ephemeris_file));
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
    } catch (...) { return 1; }

    // --- 2. Preparing Asteroids & Polynomials ---
    std::vector<std::string> asteroid_ids = parse_asteroid_list(vm["asteroid"].as<std::string>());
    double jd_start = vm["jd-start"].as<double>();
    double duration = vm["duration"].as<double>();
    time::EpochTDB start_epoch = time::EpochTDB::from_jd(jd_start);
    time::EpochTDB end_epoch = time::EpochTDB::from_jd(jd_start + duration);

    ChebyshevEphemerisManager manager(engine.config());
    HorizonsClient horizons;
    std::map<std::string, physics::KeplerianStateTyped<core::ECLIPJ2000>> stored_elements;
    std::map<std::string, std::pair<double, double>> stored_props; // {H, D}

    std::cout << "[ioccultcalc] Pre-calculating polynomials for " << asteroid_ids.size() << " bodies over " << duration << " days...\n";
    for (const auto& id : asteroid_ids) {
        auto state_vec = horizons.query_vectors(id, start_epoch);
        if (state_vec) {
            auto elements = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(state_vec->cast_frame<core::ECLIPJ2000>());
            manager.add_asteroid(id, elements, start_epoch, end_epoch);
            stored_elements[id] = elements;
            
            auto props = horizons.query_physical_properties(id);
            if (props) {
                stored_props[id] = {props->h_mag, props->diameter_km};
            } else {
                stored_props[id] = {0.0, 100.0};
            }
        }
    }

    // --- 3. Global Search ---
    double mag_limit = vm["mag"].as<double>();
    std::cout << "[ioccultcalc] Searching occultations..." << std::endl;
    auto results = OccultationLogic::find_multi_asteroid_occultations(asteroid_ids, manager, start_epoch, end_epoch, mag_limit, engine);
    std::cout << "[ioccultcalc] Total: " << results.size() << " event(s) found.\n";

    // --- 4. Outputs ---
    if (!results.empty()) {
        std::cout << "\nLIST OF EVENTS:\n";
        for (const auto& res : results) {
            std::cout << " Body: " << std::setw(10) << res.asteroid_id 
                      << " | Star: " << std::setw(20) << res.params.star_id 
                      << " | TCA: " << std::fixed << std::setprecision(5) << res.params.t_ca.mjd() << "\n";
        }
    }

    if (vm.count("xml-output") && !results.empty()) {
        std::vector<OccultationEvent> out_events;
        for (const auto& res : results) {
            auto cand = res;
            auto prop_it = stored_props.find(res.asteroid_id);
            if (prop_it != stored_props.end()) {
                cand.params.star_mag = prop_it->second.first; // Using field as placeholder for body mag
                cand.params.cross_track_uncertainty = physics::Distance::from_km(prop_it->second.second);
            }
            
            auto it = stored_elements.find(res.asteroid_id);
            if (it != stored_elements.end()) {
                 out_events.push_back(candidate_to_event(cand, res.asteroid_id, it->second));
            }
        }
        OccultationXMLIO::write_file(out_events, vm["xml-output"].as<std::string>());
    }

    catalog::GaiaDR3Catalog::shutdown();
    return 0;
}
