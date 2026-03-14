/**
 * @file ioccultcalc.cpp
 * @brief Professional command-line tool for searching and comparing stellar occultations.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/io/OccultationXMLIO.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/astrometry/OccultationEvent.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <tuple>
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
    
    // Star data
    ev.star_catalog_id = cand.params.star_id;
    ev.ra_event_h = cand.star.ra.to_deg() / 15.0;
    ev.dec_event_deg = cand.star.dec.to_deg();
    ev.mag_v = cand.star.g_mag;
    ev.ra_cat_h = ev.ra_event_h; 
    ev.dec_cat_deg = ev.dec_event_deg;
    
    // Object data
    ev.object_name = ast_id;
    try { ev.object_number = std::stoi(ast_id); } catch(...) { ev.object_number = 0; }
    ev.object_type = "Asteroid";
    ev.diameter_km = 0.0; // Placeholder
    
    // Geometry
    ev.combined_error_arcsec = cand.params.impact_parameter.to_km();
    
    // Orbit elements mapping (matching library conventions)
    ev.semi_major_axis_au = el.a.to_au();
    ev.eccentricity = el.e;
    ev.inclination_deg = el.i.to_deg();
    ev.node_deg = el.node.to_deg();
    ev.mean_anomaly_deg = el.omega.to_deg(); 
    ev.ra_node_approx = el.M.to_deg();

    // Set date from epoch
    auto [year, month, day, frac] = time::mjd_to_calendar(el.epoch.mjd());
    ev.epoch_year = year;
    ev.epoch_month = month;
    ev.epoch_day = day;
    
    return ev;
}

void print_separator() {
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

int main(int argc, char** argv) {
    po::options_description desc("ioccultcalc: AstDyn Occultation Tool CLI\nUsage: ioccultcalc --asteroid <num> --jd <jd_val> [options]");
    desc.add_options()
        ("help,h", "produce help message")
        ("conf", po::value<std::string>(), "JSON configuration file for AstDynEngine")
        ("asteroid", po::value<std::string>(), "asteroid number or designation (e.g., 704 or Interamnia)")
        ("jd", po::value<double>(), "Julian Date for searching occultations (TDB)")
        ("mag", po::value<double>()->default_value(15.0), "magnitude limit for online stellar catalog search")
        ("xml-output", po::value<std::string>(), "save found occultations to the specified XML file")
        ("kml", po::value<std::string>(), "save the path of the first matching event to KML")
        ("xml-check", po::value<std::string>(), "compare search results with a reference XML file")
    ;

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const std::exception& e) {
        std::cerr << "Command line error: " << e.what() << "\n";
        return 1;
    }

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    if (!vm.count("asteroid") || !vm.count("jd")) {
        std::cerr << "Error: Mandatory parameters --asteroid and --jd are missing.\n";
        return 1;
    }

    // --- 1. System and Engine Setup ---
    // Try to find the DE441 ephemeris in standard locations
    std::string bsp_path = "de441.bsp"; // Default to current directory
    if (!std::filesystem::exists(bsp_path)) {
        bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    }
    
    AstDynEngine engine;
    
    if (vm.count("conf")) {
        engine.load_config(vm["conf"].as<std::string>());
    } else {
        AstDynConfig cfg;
        cfg.ephemeris_file = bsp_path;
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.verbose = false;
        
        // --- Full Precision Default Settings ---
        cfg.integrator_type = IntegratorType::AAS;
        cfg.aas_precision = 1e-4; 
        cfg.initial_step_size = 0.01; 
        
        cfg.propagator_settings.include_planets = true;
        cfg.propagator_settings.include_moon = true;
        cfg.propagator_settings.include_relativity = true;
        cfg.propagator_settings.include_asteroids = true;
        
        cfg.light_time_correction = true;
        cfg.aberration_correction = true;
        
        engine.set_config(cfg);
    }
    
    try {
        ephemeris::PlanetaryEphemeris::setProvider(std::make_shared<ephemeris::DE441Provider>(engine.config().ephemeris_file));
    } catch (const std::exception& e) {
        std::cerr << "Error: Ephemeris setup failed: " << e.what() << "\n";
        return 1;
    }

    // --- 2. Input Data Retrieval ---
    std::string asteroid_id = vm["asteroid"].as<std::string>();
    double jd = vm["jd"].as<double>();
    time::EpochTDB search_epoch = time::EpochTDB::from_jd(jd);

    std::cout << "[ioccultcalc] Target: " << asteroid_id << " @ JD " << std::fixed << std::setprecision(6) << jd << "\n";
    io::HorizonsClient horizons;
    auto state_vec = horizons.query_vectors(asteroid_id, search_epoch);
    if (!state_vec) {
        std::cerr << "Error: Failed to fetch Cartesian vectors from JPL Horizons.\n";
        return 1;
    }
    // Convert to Keplerian for initial state bridge (ensuring Ecliptic J2000 frame)
    auto state_vec_ecl = state_vec->cast_frame<core::ECLIPJ2000>();
    auto elements = propagation::cartesian_to_keplerian<core::ECLIPJ2000>(state_vec_ecl);
    std::cout << "[ioccultcalc] Success: Cartesian state fetched from JPL Horizons.\n";

    // --- 3. Occultation Search ---
    try {
        catalog::GaiaDR3Catalog::initialize(R"({"catalog_type":"online_esa", "timeout_seconds":60})");
    } catch (const std::exception& e) {
        std::cerr << "Error: Catalog initialization failed: " << e.what() << "\n";
        return 1;
    }

    double mag_limit = vm["mag"].as<double>();
    time::EpochTDB start_window = search_epoch - time::TimeDuration::from_days(0.5);
    time::EpochTDB end_window = search_epoch + time::TimeDuration::from_days(0.5);

    std::cout << "[ioccultcalc] Searching Gaia DR3 (Online ESA) Corridor (Mag < " << mag_limit << ")..." << std::endl;
    auto results = OccultationLogic::find_occultations(asteroid_id, elements, start_window, end_window, mag_limit, engine);
    std::cout << "[ioccultcalc] Result: " << results.size() << " event(s) found.\n";

    // --- 4. Processing and Outputs ---
    if (!results.empty()) {
        std::cout << "\nLIST OF EVENTS:\n";
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& res = results[i];
            std::cout << std::setw(2) << i + 1 << ". Star: " << std::setw(20) << res.params.star_id 
                      << " | TCA (MJD): " << std::fixed << std::setprecision(5) << res.params.t_ca.mjd()
                      << " | Impact: " << std::setw(8) << std::setprecision(1) << res.params.impact_parameter.to_km() << " km"
                      << " | Vel: " << std::setprecision(2) << res.params.shadow_velocity.to_km_s() << " km/s\n";
        }

        if (vm.count("kml")) {
             auto path_data = OccultationMapper::compute_path(results[0].params, results[0].star.ra, results[0].star.dec, 
                                                             physics::Distance::from_km(250.0), 
                                                             time::to_utc(results[0].params.t_ca));
             OccultationMapper::export_kml(path_data, vm["kml"].as<std::string>());
        }
    }

    if (vm.count("xml-output")) {
        std::vector<OccultationEvent> out_events;
        for (const auto& res : results) out_events.push_back(candidate_to_event(res, asteroid_id, elements));
        OccultationXMLIO::write_file(out_events, vm["xml-output"].as<std::string>());
    }

    // --- 5. Comparison Protocol ---
    if (vm.count("xml-check")) {
        std::string ref_path = vm["xml-check"].as<std::string>();
        auto references = OccultationXMLIO::read_file(ref_path);
        
        std::cout << "\n" << "================================================================================" << "\n";
        std::cout << " ASTDYN COMPARISON REPORT (vs " << ref_path << ")\n";
        std::cout << "================================================================================" << "\n";
        
        for (const auto& ref : references) {
            std::cout << "REFERENCE: " << ref.event_id << " [MJD " << ref.mjd << "]\n";
            
            // Re-calculate AstDyn position at the exact XML MJD time
            time::EpochTDB t_ref = time::EpochTDB::from_mjd(ref.mjd);
            auto obs = astrometry::AstrometryReducer::compute_observation(
                elements, elements.epoch, t_ref, engine.config(), astrometry::AstrometricSettings());
            
            if (obs) {
                double ra_ref = ref.ra_event_h * 15.0;
                double de_ref = ref.dec_event_deg;
                double ra_ast = (*obs).ra.value * constants::RAD_TO_DEG;
                double de_ast = (*obs).dec.value * constants::RAD_TO_DEG;
                
                double dra = (ra_ast - ra_ref) * 3600.0 * std::cos(de_ref * constants::DEG_TO_RAD);
                double dde = (de_ast - de_ref) * 3600.0;
                
                std::cout << "  - Sky Coordinates (Equatorial Apparent):\n";
                std::cout << "    - Reference:  RA " << std::setw(12) << ra_ref << " | Dec " << std::setw(12) << de_ref << "\n";
                if (std::abs(dra) > 3600.0 || std::abs(dde) > 3600.0) {
                     std::cout << "    - AstDyn:     RA " << std::setw(12) << ra_ast << " | Dec " << std::setw(12) << de_ast << " [WARNING: BIG OFFSET]\n";
                } else {
                     std::cout << "    - AstDyn:     RA " << std::setw(12) << ra_ast << " | Dec " << std::setw(12) << de_ast << "\n";
                }
                std::cout << "    - Difference: dRA*cosD: " << std::setw(10) << dra << " arcsec | dDec: " << std::setw(10) << dde << " arcsec\n";
            }
            
            // Check if any results from our search match this reference
            bool found = false;
            for (const auto& res : results) {
                if (std::abs(res.params.t_ca.mjd() - ref.mjd) < 0.04) { // within ~1 hour
                    std::cout << "  - Status: MATCH FOUND in search (MJD delta: " << (res.params.t_ca.mjd() - ref.mjd)*1440.0 << " mins)\n";
                    std::cout << "    - Impact: Ref (N/A) | AstDyn: " << res.params.impact_parameter.to_km() << " km\n";
                    found = true;
                    break;
                }
            }
            if (!found) std::cout << "  - Status: NO MATCH in search results (check Corridor/Mag limits)\n";
            print_separator();
        }
    }

    catalog::GaiaDR3Catalog::shutdown();
    return 0;
}
