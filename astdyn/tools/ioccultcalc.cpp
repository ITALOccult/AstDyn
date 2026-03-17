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
#include "astdyn/io/CovarianceIO.hpp"
#include "astdyn/math/MultivariateSampler.hpp"
#include "astdyn/core/IOCConfig.hpp"
#include "astdyn/time/TimeScale.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <sstream>
#include <filesystem>

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
        ("svg-output", po::value<std::string>(), "save world map with paths to SVG")
        ("kml", po::value<std::string>(), "save first match to KML")
        ("out-dir", po::value<std::string>(), "base directory for all output files")
        ("prefix", po::value<std::string>()->default_value("occ"), "prefix for individual output files")
        ("lat", po::value<double>(), "observer latitude (degrees) for regional search")
        ("lon", po::value<double>(), "observer longitude (degrees) for regional search")
        ("alt", po::value<double>()->default_value(0.0), "observer altitude (meters)")
        ("bsp", po::value<std::string>(), "path to satellite ephemeris file (BSP)")
        ("system-ids", po::value<std::string>(), "comma-separated NAIF IDs for system bodies (e.g. 100,201)")
        ("covariance,c", po::value<std::string>(), "path to orbital covariance file (.cor or .csv)")
        ("clones,n", po::value<int>()->default_value(0), "number of Monte Carlo clones to generate")
        ("zoom", po::value<double>()->default_value(1.0), "zoom level for SVG map")
        ("map-lat", po::value<double>(), "center latitude for SVG map")
        ("map-lon", po::value<double>(), "center longitude for SVG map")
        ("catalog", po::value<std::string>()->default_value("gaia_dr3"), "stellar catalog to use (gaia_dr3, legacy)")
    ;

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (...) { return 1; }

    // --- 1. Load Configuration ---
    core::IOCConfig adv_cfg;
    if (vm.count("conf")) {
        adv_cfg.load(vm["conf"].as<std::string>());
    }

    // --- 2. Engine Setup ---
    std::string bsp_path = adv_cfg.get<std::string>("ephemeris.file", "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp");
    
    AstDynEngine engine;
    if (vm.count("conf")) {
        engine.load_config(vm["conf"].as<std::string>());
    } else {
        AstDynConfig cfg;
        cfg.ephemeris_file = bsp_path;
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.verbose = false;
        cfg.preferred_catalog = vm["catalog"].as<std::string>();
        engine.set_config(cfg);
    }
    
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephem_ptr;
    try {
        auto provider = std::make_shared<ephemeris::DE441Provider>(engine.config().ephemeris_file);
        ephemeris::PlanetaryEphemeris::setGlobalProvider(provider);
        ephem_ptr = std::make_shared<ephemeris::PlanetaryEphemeris>(provider);
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa"})");
    } catch (...) { return 1; }

    // --- 3. Preparing Asteroids & Polynomials ---
    std::vector<std::string> asteroid_ids;
    if (vm.count("asteroid")) {
        asteroid_ids = parse_asteroid_list(vm["asteroid"].as<std::string>());
    } else if (adv_cfg.has("asteroid")) {
        asteroid_ids = parse_asteroid_list(adv_cfg.get<std::string>("asteroid", ""));
    }

    double jd_start = 0.0;
    if (vm.count("jd-start")) {
        jd_start = vm["jd-start"].as<double>();
    } else if (adv_cfg.has("jd-start")) {
        jd_start = adv_cfg.get<double>("jd-start", 0.0);
    }


    if ((asteroid_ids.empty() && !vm.count("system-ids") && !adv_cfg.has("system-ids")) || jd_start == 0.0) {
        if (!vm.count("help")) std::cerr << "Error: Missing asteroid list or jd-start (check CLI or config file).\n";
        std::cout << desc << "\n";
        return 1;
    }

    double duration = (!vm["duration"].defaulted()) ? vm["duration"].as<double>() : adv_cfg.get<double>("duration", 1.0);
    double mag_limit = (!vm["mag"].defaulted()) ? vm["mag"].as<double>() : adv_cfg.get<double>("mag", 15.0);
    
    std::string out_dir = vm.count("out-dir") ? vm["out-dir"].as<std::string>() : adv_cfg.get<std::string>("out-dir", "");
    std::string prefix = (!vm["prefix"].defaulted()) ? vm["prefix"].as<std::string>() : adv_cfg.get<std::string>("prefix", "occ");

    if (!out_dir.empty()) {
        std::filesystem::create_directories(out_dir);
    }

    time::EpochTDB start_epoch = time::EpochTDB::from_jd(jd_start);
    time::EpochTDB end_epoch = time::EpochTDB::from_jd(jd_start + duration);

    ChebyshevEphemerisManager manager(engine.config());
    HorizonsClient horizons;
    std::map<std::string, physics::CartesianStateTyped<core::GCRF>> stored_states;
    std::map<std::string, physics::KeplerianStateTyped<core::ECLIPJ2000>> stored_elements;
    std::map<std::string, std::pair<double, double>> stored_props; // {H, D}
    std::unique_ptr<io::SPKReader> system_reader;

    // --- 2a. Load Primary Asteroids ---
    if (!asteroid_ids.empty()) {
        std::cout << "[ioccultcalc] Pre-calculating polynomials for " << asteroid_ids.size() << " bodies over " << duration << " days...\n";
        for (const auto& id : asteroid_ids) {
            auto state_vec = horizons.query_vectors(id, start_epoch);
            if (state_vec) {
                auto elements = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(state_vec->cast_frame<core::ECLIPJ2000>());
                manager.add_asteroid(id, elements, start_epoch, end_epoch);
                stored_elements[id] = elements;
                stored_states[id] = state_vec->cast_frame<core::GCRF>();
                
                auto props = horizons.query_physical_properties(id);
                if (props) {
                    stored_props[id] = {props->h_mag, props->diameter_km};
                } else {
                    stored_props[id] = {0.0, 100.0};
                }
            }
        }
    }

    // --- 3b. Load System/Satellite Bodies from BSP ---
    std::string system_bsp = vm.count("bsp") ? vm["bsp"].as<std::string>() : adv_cfg.get<std::string>("bsp", "");
    std::string system_ids_str = vm.count("system-ids") ? vm["system-ids"].as<std::string>() : adv_cfg.get<std::string>("system-ids", "");

    if (!system_bsp.empty() && !system_ids_str.empty()) {
        std::vector<std::string> system_ids = parse_asteroid_list(system_ids_str);
        std::cout << "[ioccultcalc] Adding " << system_ids.size() << " bodies from BSP: " << system_bsp << "\n";
        
        system_reader = std::make_unique<io::SPKReader>(system_bsp);
        for (const auto& id_str : system_ids) {
            int naif_id = std::stoi(id_str);
            manager.add_system_body(id_str, naif_id, *system_reader, start_epoch, end_epoch);
            asteroid_ids.push_back(id_str); // Add to search list
            stored_props[id_str] = {15.0, 10.0}; // Default for satellites (H, D)
        }
    }

    // --- 3c. Global Search ---
    OccultationConfig occ_config = engine.config().occultation_settings;
    occ_config.max_mag_star = mag_limit;
    
    std::vector<OccultationCandidate> results;

    if (!system_bsp.empty() && !system_ids_str.empty()) {
        std::vector<std::string> system_ids = parse_asteroid_list(system_ids_str);
        std::cout << "[ioccultcalc] Searching system occultations (BSP: " << system_bsp << ")..." << std::endl;
        auto system_results = OccultationLogic::find_system_occultations(system_ids, system_bsp, start_epoch, end_epoch, occ_config, engine);
        
        for (const auto& res : system_results) {
            for (const auto& body : res.bodies) {
                OccultationCandidate cand;
                cand.asteroid_id = body.name;
                cand.star = res.star;
                cand.params = body.params;
                results.push_back(cand);
            }
        }
    } else {
        std::cout << "[ioccultcalc] Searching occultations with mag < " << occ_config.max_mag_star << "..." << std::endl;
        results = OccultationLogic::find_multi_asteroid_occultations(asteroid_ids, manager, start_epoch, end_epoch, occ_config, engine);
    }
    
    // --- 3d. Apply Uncertainty (if requested) ---
    std::string cov_file = vm.count("covariance") ? vm["covariance"].as<std::string>() : adv_cfg.get<std::string>("covariance", "");
    if (!cov_file.empty() && !results.empty()) {
        try {
            astdyn::Matrix6d cov_t0 = CovarianceIO::read_file(cov_file);
            for (auto& res : results) {
                if (stored_states.count(res.asteroid_id)) {
                    OccultationLogic::apply_uncertainty(res.params, res.star, cov_t0, stored_states[res.asteroid_id], engine);
                }
            }
        } catch (...) {}
    }

    // --- 3e. Regional Filter ---
    double obs_lat = 0.0, obs_lon = 0.0;
    bool has_location = false;
    if (vm.count("lat") && vm.count("lon")) {
        obs_lat = vm["lat"].as<double>();
        obs_lon = vm["lon"].as<double>();
        has_location = true;
    }

    if (has_location && !results.empty()) {
        std::vector<OccultationCandidate> filtered;
        for (const auto& res : results) {
            double d_lat = res.params.center_lat.to_deg() - obs_lat;
            double d_lon = res.params.center_lon.to_deg() - obs_lon;
            while (d_lon > 180.0) d_lon -= 360.0;
            while (d_lon < -180.0) d_lon += 360.0;
            if (std::sqrt(d_lat*d_lat + d_lon*d_lon) < 25.0) filtered.push_back(res);
        }
        results = filtered;
    }

    // --- 4. Outputs ---
    if (!results.empty()) {
        std::vector<OccultationPath> paths;
        std::vector<std::string> labels;
        std::vector<std::string> colors = {"#ef4444", "#3b82f6", "#22c55e", "#eab308", "#8b5cf6"};

        for (const auto& res : results) {
            double diam = 100.0;
            auto it = stored_props.find(res.asteroid_id);
            if (it != stored_props.end()) diam = it->second.second;
            
            time::EpochUTC tca_utc = time::to_utc(res.params.t_ca);
            auto path = OccultationMapper::compute_path(res.params, res.star.ra, res.star.dec, physics::Distance::from_km(diam), tca_utc, ephem_ptr);
            paths.push_back(path);
            labels.push_back(res.asteroid_id + " - " + std::to_string(res.star.source_id));
        }

        std::string svg_file = vm.count("svg-output") ? vm["svg-output"].as<std::string>() : adv_cfg.get<std::string>("svg-output", "");
        if (!svg_file.empty()) {
            std::filesystem::path p = std::filesystem::path(out_dir) / svg_file;
            double zoom = vm["zoom"].as<double>();
            double map_lat = vm.count("map-lat") ? vm["map-lat"].as<double>() : obs_lat;
            double map_lon = vm.count("map-lon") ? vm["map-lon"].as<double>() : obs_lon;
            
            OccultationMapper::export_global_svg(paths, labels, colors, p.string(), ephem_ptr, 
                                                 Angle::from_deg(map_lat), Angle::from_deg(map_lon), zoom);
        }

        if (!out_dir.empty()) {
            for (size_t i = 0; i < results.size(); ++i) {
                const auto& res = results[i];
                std::stringstream ss;
                ss << prefix << "_" << res.asteroid_id << "_" << res.star.source_id << "_" << (int)res.params.t_ca.mjd() << ".svg";
                std::filesystem::path ip = std::filesystem::path(out_dir) / ss.str();
                OccultationMapper::export_global_svg({paths[i]}, {labels[i]}, {colors[i % colors.size()]}, 
                                                     ip.string(), ephem_ptr, res.params.center_lat, res.params.center_lon, 4.0);
            }
        }
    }

    catalog::GaiaDR3Catalog::shutdown();
    return 0;
}
