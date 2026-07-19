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
#include <optional>
#include "astdyn/coordinates/EquinoctialElements.hpp"
#include "astdyn/coordinates/CelestialToTerrestrial.hpp"
#include <cstdio>
#include <chrono>

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
OccultationEvent candidate_to_event(const OccultationCandidate& cand, const std::string& ast_id,
                                    const physics::KeplerianStateTyped<core::ECLIPJ2000>& el,
                                    double diameter_km, double h_mag,
                                    const std::string& ast_name = "", double g_slope = 0.15) {
    constexpr double kEarthRadiusM = 6378137.0;
    const auto& pr = cand.params;
    OccultationEvent ev;

    const time::EpochUTC t_utc = time::to_utc(pr.t_ca);
    const time::EpochTT  t_tt  = time::to_tt(t_utc);
    auto [ey, em, ed, efrac] = time::mjd_to_calendar(t_utc.mjd());

    // ---- <Elements> : the payload Occult4 uses to draw the path -----------
    ev.elements_source = "AstDyn-AAS-GaiaDR3";
    ev.duration_s   = pr.max_duration.to_seconds();
    ev.year = ey; ev.month = em; ev.day = ed;
    ev.ut_closest_h = efrac * 24.0;
    ev.x  = pr.xi_ca.to_m()  / kEarthRadiusM;
    ev.y  = pr.eta_ca.to_m() / kEarthRadiusM;
    ev.dx = pr.dxi_dt.to_ms()  * 3600.0 / kEarthRadiusM;
    ev.dy = pr.deta_dt.to_ms() * 3600.0 / kEarthRadiusM;
    // The higher-order terms would need the Chebyshev segment differentiated a
    // second time; over a ~40 min event they are worth a few km at the path ends.
    ev.d2x = ev.d2y = ev.d3x = ev.d3y = 0.0;

    // ---- <Earth> ----------------------------------------------------------
    ev.substellar_lon_deg = pr.substar_lon.to_deg();
    ev.substellar_lat_deg = pr.substar_lat.to_deg();   // geocentric == apparent Dec
    ev.subsolar_lon_deg   = pr.subsolar_lon.to_deg();
    ev.subsolar_lat_deg   = pr.subsolar_lat.to_deg();
    ev.jwst = false;

    // ---- <Star> -----------------------------------------------------------
    ev.star_id = cand.star.source_id != 0
               ? "Gaia DR3 " + std::to_string(cand.star.source_id)
               : pr.star_id;
    // The spec asks for the BCRS position at the epoch of the event with NO
    // parallax applied, which is exactly what predict_at() returns when it is
    // called without an observer position.
    const auto s_ep = cand.star.predict_at(pr.t_ca);
    ev.star_ra_h    = s_ep.ra().to_deg() / 15.0;
    ev.star_dec_deg = s_ep.dec().to_deg();
    // Gaia BP/G/RP are the closest available proxies for B/V/R.
    ev.mag_b = cand.star.bp_mag;
    ev.mag_v = cand.star.g_mag;
    ev.mag_r = cand.star.rp_mag;
    ev.star_diameter_mas = 0.0;   // not modelled
    ev.double_star_code  = 0;
    ev.k2_flag           = "";
    {
        Angle ra_app, dec_app;
        coordinates::apparent_place(t_tt, s_ep.ra(), s_ep.dec(), ra_app, dec_app);
        ev.star_app_ra_h    = ra_app.to_deg() / 15.0;
        ev.star_app_dec_deg = dec_app.to_deg();
    }
    // Magnitude drop: during the event only the asteroid is seen, before it the
    // combined light. Occult4's 10.69 is simply 21.11 - 10.42.
    if (ev.object_mag > -5.0) {
        const double comb_v = -2.5 * std::log10(std::pow(10.0, -0.4 * ev.mag_v) +
                                                std::pow(10.0, -0.4 * ev.object_mag));
        const double comb_r = -2.5 * std::log10(std::pow(10.0, -0.4 * ev.mag_r) +
                                                std::pow(10.0, -0.4 * ev.object_mag));
        ev.mag_drop_v = ev.object_mag - comb_v;
        ev.mag_drop_r = ev.object_mag - comb_r;
    }
    ev.mag_drops_adjusted  = 0;
    // No nearby-star check is performed; Occult4/OWC expects 0 (not counted),
    // never -1, for these two fields.
    ev.bright_nearby_count = 0;
    ev.total_nearby_count  = 0;

    // ---- <Object> ---------------------------------------------------------
    ev.object_number = ast_id;
    ev.object_name   = ast_name.empty() ? ast_id : ast_name;
    // HG apparent magnitude. Occult4 writes -5.00 when it cannot compute one;
    // zero would claim an object brighter than Vega and poison the magnitude
    // drop below, so the sentinel is honoured rather than reinvented.
    ev.object_mag = hg_magnitude(h_mag, g_slope,
                                 pr.heliocentric_distance.to_au(),
                                 pr.geocentric_distance.to_au(),
                                 pr.phase_angle);
    ev.diameter_km   = diameter_km;
    ev.distance_au   = pr.geocentric_distance.to_au();
    ev.n_rings = 0;
    ev.n_moons = 0;
    // The format wants dRA in SECONDS OF TIME per hour; our rate is dRA (not
    // dRA*cos(dec)) in arcsec per hour, so only the 15 is needed.
    ev.d_ra_s_hr   = pr.d_ra_arcsec_hr / 15.0;
    ev.d_dec_as_hr = pr.d_dec_arcsec_hr;
    ev.taxonomy = "";
    ev.diameter_uncertainty_km = 0.0;
    ev.moon_in_planet_shadow   = 0;
    ev.mag_v_asteroid = 0.0;
    ev.mag_r_asteroid = 0.0;

    // ---- <Orbit> : low-precision, for plotting only -----------------------
    auto [oy, om, od, ofrac] = time::mjd_to_calendar(el.epoch.mjd());
    ev.equinox          = 0.0;      // J2000
    ev.mean_anomaly_deg = el.M.to_deg();
    ev.epoch_year = oy; ev.epoch_month = om; ev.epoch_day = od;
    ev.peri_deg           = el.omega.to_deg();
    ev.node_deg           = el.node.to_deg();
    ev.inclination_deg    = el.i.to_deg();
    ev.eccentricity       = el.e;
    ev.semi_major_axis_au = el.a.to_au();
    ev.perihelion_au      = el.a.to_au() * (1.0 - el.e);
    ev.h0          = h_mag;
    ev.coeff_log_r = 5.0;           // standard for asteroids
    ev.g_param     = 0.15;

    // ---- <Errors> : this is where SCOPE surfaces --------------------------
    // These fields are only meaningful when a covariance was supplied and
    // apply_uncertainty ran; an ellipse of zero means "not computed", not
    // "certain", and the error basis says which.
    const bool have_cov = pr.err_major.to_rad() > 0.0;

    ev.err_major_as = pr.err_major.to_arcsec();
    ev.err_minor_as = pr.err_minor.to_arcsec();
    ev.err_pa_deg   = pr.err_pa.to_deg();
    // Occult4's field 5 is the quadrature sum of the two semi-axes: its
    // 0.0790 is sqrt(0.0690^2 + 0.0387^2) = 0.0791.
    ev.err_1sigma_as = std::hypot(ev.err_major_as, ev.err_minor_as);

    // Path location uncertainty as a fraction of the path width. Only the
    // CROSS-TRACK component displaces the path, which is why apply_uncertainty
    // projects the ellipse onto the direction perpendicular to the motion.
    //
    // NOTE: reconstructing Occult4's own numbers gives 59.1 against the 59.015
    // shown on its prediction page, but its <Errors> record carries 60.015 --
    // exactly one more. The specification does not explain the offset, so the
    // physical quantity is written here and the discrepancy left visible rather
    // than papered over with a fudge.
    ev.err_path_widths = (have_cov && diameter_km > 0.0)
                       ? pr.cross_track_uncertainty.to_km() / diameter_km
                       : 0.0;

    // "Known errors" means both the star and the object have measured
    // covariances. The star half is only true for a full Gaia solution.
    // Occult4/OWC canonical error-basis strings (verified against a real
    // OccultWatcher export): only "Star+PeakEphemUncert" and "Star+Assumed"
    // are accepted. The true-vs-estimated distinction lives in our logs, not
    // in the interchange XML, which OWC parses strictly.
    ev.error_basis = have_cov ? "Star+PeakEphemUncert" : "Star+Assumed";

    // Occult4/OWC uses 0 (flag inactive), never -1: a real OWC export shows
    // these three fields as 0 in 1604/1630 events, 1 in the rest, never -1.
    ev.reliability = (cand.star.ruwe > 0.0) ? cand.star.ruwe : 0.0;
    ev.duplicate_source    = 0;
    ev.non_gaia_pm         = 0;
    ev.pm_added_from_ucac4 = 0;

    // The nonlinearity index has no home in the occelmnt format -- it is a SCOPE
    // quantity, not an Occult4 one -- so it is recorded in the source string,
    // where it travels with the prediction and stays visible.
    if (have_cov && pr.nonlinearity_index > 0.0) {
        char nbuf[64];
        std::snprintf(nbuf, sizeof(nbuf), "AstDyn-SCOPE-N%.2e", pr.nonlinearity_index);
        ev.elements_source = nbuf;
    }

    // ---- <ID> -------------------------------------------------------------
    // Format: yyyymmdd_xxxxxx -- the event date plus the tail of the star id.
    {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%04d%02d%02d", ey, em, ed);
        const std::string sid = std::to_string(cand.star.source_id);
        ev.event_id = std::string(buf) + "_" +
                      (sid.size() > 6 ? sid.substr(sid.size() - 6) : sid);
    }
    // This field is the date the prediction was COMPUTED, not the event epoch;
    // the latter lives in <Elements>.
    {
        const auto now = std::chrono::system_clock::now().time_since_epoch();
        const double unix_s = std::chrono::duration<double>(now).count();
        ev.prediction_mjd = 40587.0 + unix_s / 86400.0;   // MJD at 1970-01-01
    }

    return ev;
}

std::vector<std::string> parse_asteroid_list(const std::string& input) {
    std::vector<std::string> ids;
    if (input.empty()) return ids;
    
    std::vector<std::string> raw_segments;
    if (input[0] == '@') {
        std::ifstream file(input.substr(1));
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty()) raw_segments.push_back(line);
        }
    } else {
        std::stringstream ss(input);
        std::string segment;
        while (std::getline(ss, segment, ',')) {
            if (!segment.empty()) raw_segments.push_back(segment);
        }
    }

    for (const auto& seg : raw_segments) {
        size_t dash_pos = seg.find('-');
        if (dash_pos != std::string::npos && dash_pos > 0 && dash_pos < seg.length() - 1) {
            try {
                int start_id = std::stoi(seg.substr(0, dash_pos));
                int end_id = std::stoi(seg.substr(dash_pos + 1));
                if (start_id > end_id) std::swap(start_id, end_id);
                for (int i = start_id; i <= end_id; ++i) {
                    ids.push_back(std::to_string(i));
                }
            } catch (...) {
                ids.push_back(seg);
            }
        } else {
            ids.push_back(seg);
        }
    }
    return ids;
}

int main(int argc, char** argv) {
    po::options_description desc("ioccultcalc: AstDyn Occultation Tool CLI\nUsage: ioccultcalc --asteroid <num1,num2,start-end...> --jd-start <jd> --duration <days>");
    desc.add_options()
        ("help,h", "produce help message")
        ("conf", po::value<std::string>(), "JSON configuration file for AstDynEngine")
        ("asteroid", po::value<std::string>(), "asteroid designations (comma-separated, ranges like '1-100', or '@file')")
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
        ("max-dist-km", po::value<double>(), "maximum distance from observer to shadow centerline [km]")
        ("min-duration", po::value<double>(), "minimum event duration [seconds]")
        ("min-diameter", po::value<double>(), "minimum asteroid diameter [km]")
        ("bsp", po::value<std::string>(), "path to satellite ephemeris file (BSP)")
        ("system-ids", po::value<std::string>(), "comma-separated NAIF IDs for system bodies (e.g. 100,201)")
        ("covariance,c", po::value<std::string>(), "path to orbital covariance file (.cor or .csv)")
        ("clones,n", po::value<int>()->default_value(0), "number of Monte Carlo clones to generate")
        ("zoom", po::value<double>()->default_value(1.0), "zoom level for SVG map")
        ("map-lat", po::value<double>(), "center latitude for SVG map")
        ("map-lon", po::value<double>(), "center longitude for SVG map")
        ("max-ruwe", po::value<double>(), "maximum Gaia DR3 RUWE for stars")
        ("max-moon-phase", po::value<double>(), "maximum Moon phase [0.0 - 1.0]")
        ("min-moon-dist", po::value<double>(), "minimum angular distance from Moon [degrees]")
        ("max-shadow-dist", po::value<double>(), "maximum shadow search distance [km]")
        ("catalog", po::value<std::string>()->default_value("gaia_dr3"), "stellar catalog to use (gaia_dr3, legacy)")
        ("multibody", po::bool_switch()->default_value(false), "use high-precision multibody propagator for refinement")
        ("star-offset-mas", po::value<double>()->default_value(0.0), "apply RA offset to star in mas (e.g. for duplicity correction)")
        ("star", po::value<std::string>(), "filter by star source ID")
    ;

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const std::exception& e) { 
        std::cerr << "Error parsing options: " << e.what() << "\n";
        return 1; 
    }

    // --- 1. Load Configuration ---
    core::IOCConfig adv_cfg;
    if (vm.count("conf")) {
        adv_cfg.load(vm["conf"].as<std::string>());
    }

    // --- 2. Engine Setup ---
    // B5 fix: use $HOME env var instead of hardcoded personal path
    auto default_bsp_path = []() -> std::string {
        const char* home = std::getenv("HOME");
        return home ? std::string(home) + "/.ioccultcalc/ephemerides/de441.bsp" : "";
    };
    std::string bsp_path = adv_cfg.get<std::string>("ephemeris.file", default_bsp_path());
    
    AstDynEngine engine;
    if (vm.count("conf")) {
        engine.load_config(vm["conf"].as<std::string>());
    } else {
        AstDynConfig cfg;
        cfg.ephemeris_file = bsp_path;
        cfg.ephemeris_type = EphemerisType::DE441;
        cfg.verbose = true;
        cfg.preferred_catalog = vm["catalog"].as<std::string>();
        engine.set_config(cfg);
    }
    
    std::shared_ptr<ephemeris::PlanetaryEphemeris> ephem_ptr;
    try {
        auto provider = std::make_shared<ephemeris::DE441Provider>(engine.config().ephemeris_file);
        ephemeris::PlanetaryEphemeris::setGlobalProvider(provider);
        ephem_ptr = std::make_shared<ephemeris::PlanetaryEphemeris>(provider);
        // M2 fix: catalog config from --conf or hardcoded default
        std::string catalog_json = adv_cfg.get<std::string>("catalog_config",
            R"({"catalog_type":"sqlite_dr3","sqlite_file_path":"~/.catalog/crossreference/gaia_dr3_occult_pro.db"})");
        catalog::GaiaDR3Catalog::initialize(catalog_json);
    } catch (const std::exception& e) { 
        std::cerr << "Error initializing engine: " << e.what() << "\n";
        return 1; 
    }

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
    std::map<std::string, physics::CartesianStateTyped<core::ECLIPJ2000>> stored_states;
    std::map<std::string, physics::KeplerianStateTyped<core::ECLIPJ2000>> stored_elements;
    // M3: named constants for fallback physical properties
    constexpr double DEFAULT_ASTEROID_DIAMETER_KM = 100.0;
    constexpr double DEFAULT_SATELLITE_DIAMETER_KM = 10.0;
    constexpr double DEFAULT_SATELLITE_H_MAG       = 15.0;
    // Full properties, not just {H, D}: the designation is needed for the
    // occelmnt <Object> record, which wants "2015 BK290" and not "820987".
    std::map<std::string, io::PhysicalProperties> stored_props;
    std::unique_ptr<io::SPKReader> system_reader;

    // --- 2a. Load Primary Asteroids ---
    if (!asteroid_ids.empty()) {
        std::cout << "[ioccultcalc] Pre-calcolo polinomi per " << asteroid_ids.size() << " corpi over " << duration << " giorni...\n";
        size_t skipped_count = 0;
        for (const auto& id : asteroid_ids) {
            auto elements_opt = horizons.query_elements(id, start_epoch);
            if (!elements_opt) {
                std::cerr << "[ioccultcalc] ERRORE: JPL Horizons non ha restituito elementi per '"
                          << id << "'. L'asteroide sara' ignorato.\n";
            }
            if (elements_opt) {
                // Un asteroide problematico (integratore che diverge, effemeride
                // incoerente, ecc.) non deve abortire l'intero batch: lo si
                // salta e si prosegue. Con range grandi (es. 1-34244) e' certo
                // che alcuni falliranno.
                try {
                    auto elements = *elements_opt;
                    manager.add_asteroid(id, elements, start_epoch, end_epoch);
                    stored_elements[id] = elements;

                    // The state IS heliocentric ecliptic -- the variable is named for it.
                    auto state_eclip = propagation::keplerian_to_cartesian(elements);
                    stored_states[id] = state_eclip;

                    std::cout << "[ioccultcalc] '" << id << "' caricato da Horizons OK\n";
                    auto props = horizons.query_physical_properties(id);
                    if (props) {
                        stored_props[id] = *props;
                        manager.set_diameter(id, props->diameter_km);
                    } else {
                        stored_props[id] = io::PhysicalProperties{"", 0.0, DEFAULT_ASTEROID_DIAMETER_KM, 0.0};
                        manager.set_diameter(id, DEFAULT_ASTEROID_DIAMETER_KM);
                    }
                } catch (const std::exception& e) {
                    std::cerr << "[ioccultcalc] SKIP '" << id << "': " << e.what() << "\n";
                    ++skipped_count;
                    continue;
                }
            }
        }
        if (skipped_count > 0) {
            std::cout << "[ioccultcalc] " << skipped_count
                      << " asteroidi saltati (vedi SKIP sopra).\n";
        }
    }

    // --- 3b. Load System/Satellite Bodies from BSP ---
    std::string system_bsp = vm.count("bsp") ? vm["bsp"].as<std::string>() : adv_cfg.get<std::string>("bsp", "");
    std::string system_ids_str = vm.count("system-ids") ? vm["system-ids"].as<std::string>() : adv_cfg.get<std::string>("system-ids", "");

    if (!system_bsp.empty() && !system_ids_str.empty()) {
        std::vector<std::string> system_ids = parse_asteroid_list(system_ids_str);
        std::cout << "[ioccultcalc] Aggiunta di " << system_ids.size() << " corpi da BSP: " << system_bsp << "\n";
        
        system_reader = std::make_unique<io::SPKReader>(system_bsp);
        for (const auto& id_str : system_ids) {
            int naif_id = 0;
            try {
                naif_id = std::stoi(id_str);
            } catch (const std::exception& e) {
                std::cerr << "[ioccultcalc] Skipping invalid NAIF ID '" << id_str << "': " << e.what() << "\n";
                continue;
            }
            manager.add_system_body(id_str, naif_id, *system_reader, start_epoch, end_epoch);
            asteroid_ids.push_back(id_str);
            stored_props[id_str] = io::PhysicalProperties{"", DEFAULT_SATELLITE_H_MAG,
                                                          DEFAULT_SATELLITE_DIAMETER_KM, 0.0};
            manager.set_diameter(id_str, DEFAULT_SATELLITE_DIAMETER_KM);
        }
    }

    // --- 3c. Global Search ---
    OccultationConfig occ_config = engine.config().occultation_settings;
    occ_config.max_mag_star = mag_limit;
    
    // Apply scientific filters from CLI (overriding config if present)
    if (vm.count("min-duration")) occ_config.min_duration_s = vm["min-duration"].as<double>();
    if (vm.count("min-diameter")) occ_config.min_asteroid_diameter_km = vm["min-diameter"].as<double>();
    
    // Apply proximity filters from CLI
    if (vm.count("lat")) occ_config.obs_lat = vm["lat"].as<double>();
    if (vm.count("lon")) occ_config.obs_lon = vm["lon"].as<double>();
    if (vm.count("max-dist-km")) occ_config.max_obs_dist_km = vm["max-dist-km"].as<double>();

    // Apply scientific quality filters from CLI
    if (vm.count("max-ruwe")) occ_config.max_gaia_ruwe = vm["max-ruwe"].as<double>();
    if (vm.count("max-moon-phase")) occ_config.max_moon_phase = vm["max-moon-phase"].as<double>();
    if (vm.count("min-moon-dist")) occ_config.min_moon_dist = vm["min-moon-dist"].as<double>();
    if (vm.count("max-shadow-dist")) occ_config.max_shadow_distance = physics::Distance::from_km(vm["max-shadow-dist"].as<double>());

    std::vector<OccultationCandidate> results;

    if (!system_bsp.empty() && !system_ids_str.empty()) {
        std::vector<std::string> system_ids = parse_asteroid_list(system_ids_str);
        std::cout << "[ioccultcalc] Ricerca occultazioni di sistema (BSP: " << system_bsp << ")..." << std::endl;
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
        std::cout << "[ioccultcalc] Trovate " << results.size() << " potenziali occultazioni." << std::endl;
    } else {
        std::cout << "[ioccultcalc] Ricerca occultazioni con magnitudine < " << occ_config.max_mag_star << "..." << std::endl;
        results = OccultationLogic::find_multi_asteroid_occultations(asteroid_ids, manager, start_epoch, end_epoch, occ_config, engine);
        std::cout << "[ioccultcalc] Trovate " << results.size() << " potenziali occultazioni." << std::endl;
    }
    
    // --- 3d. Apply Uncertainty (if requested) ---
    std::string cov_file = vm.count("covariance") ? vm["covariance"].as<std::string>() : adv_cfg.get<std::string>("covariance", "");
    if (!cov_file.empty() && !results.empty()) {
        try {
            constexpr double kAuKm = 149597870.700;
            astdyn::Matrix6d cov_t0;
            std::optional<physics::CartesianStateTyped<core::ECLIPJ2000>> x0;

            if (cov_file.size() > 4 && cov_file.substr(cov_file.size() - 4) == ".eq1") {
                // AstDyS publishes the elements AND their covariance at the SAME
                // epoch, and both must be used: taking the covariance from the .eq1
                // while leaving the state at the Horizons epoch silently propagates
                // C(t0) from the wrong t0 -- 48 days out, in this case.
                auto orb = CovarianceIO::read_eq1(cov_file);
                CovarianceIO::to_si(orb);          // a -> km, lambda -> rad
                coordinates::EquinoctialElements eq(orb.elements(0), orb.elements(1),
                                                    orb.elements(2), orb.elements(3),
                                                    orb.elements(4), orb.elements(5));

                // The covariance lives in equinoctial elements; the tensor propagates
                // in Cartesian. Nothing complains if you skip this rotation, because
                // both are 6x6 -- hence the Jacobian.
                const astdyn::Matrix6d J = eq.jacobian_to_cartesian();
                const astdyn::Matrix6d C_km = J * orb.covariance * J.transpose();

                // ... and the tensor works in AU and days.
                Eigen::Matrix<double, 6, 1> d;
                d << 1.0 / kAuKm, 1.0 / kAuKm, 1.0 / kAuKm,
                     86400.0 / kAuKm, 86400.0 / kAuKm, 86400.0 / kAuKm;
                cov_t0 = d.asDiagonal() * C_km * d.asDiagonal();

                const auto st = eq.to_cartesian();
                const auto t0 = time::to_tdb(time::EpochTT::from_mjd(orb.epoch_mjd_tt));
                x0 = physics::CartesianStateTyped<core::ECLIPJ2000>::from_si(
                    t0, st.position()(0) * 1e3, st.position()(1) * 1e3, st.position()(2) * 1e3,
                        st.velocity()(0) * 1e3, st.velocity()(1) * 1e3, st.velocity()(2) * 1e3,
                    constants::GMS_SI);

                std::cout << "[ioccultcalc] AstDyS: epoca MJD " << orb.epoch_mjd_tt
                          << " TT, H=" << orb.h_mag << "\n";
            } else {
                // A bare 6x6: assumed already Cartesian, in AU and days, at the epoch
                // of the stored state.
                cov_t0 = CovarianceIO::read_file(cov_file);
            }

            for (auto& res : results) {
                if (x0) {
                    OccultationLogic::apply_uncertainty(res.params, res.star, cov_t0, *x0, engine);
                } else if (stored_states.count(res.asteroid_id)) {
                    OccultationLogic::apply_uncertainty(res.params, res.star, cov_t0,
                                                        stored_states[res.asteroid_id], engine);
                }
                if (res.params.nonlinearity_index > 0.0) {
                    std::cout << "[ioccultcalc] " << res.asteroid_id
                              << ": ellisse 1-sigma " << res.params.err_major.to_arcsec()
                              << "\" x " << res.params.err_minor.to_arcsec()
                              << "\" @ PA " << res.params.err_pa.to_deg()
                              << "   cross-track " << res.params.cross_track_uncertainty.to_km()
                              << " km   N=" << res.params.nonlinearity_index << "\n";
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "[ioccultcalc] Warning: covariance application failed (" << e.what() << "). Results will have no uncertainty.\n";
        }
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
        // B4 fix: haversine great-circle distance instead of flat Euclidean degrees
        auto haversine_deg = [](double lat1_deg, double lon1_deg, double lat2_deg, double lon2_deg) -> double {
            constexpr double DEG_TO_RAD = M_PI / 180.0;
            double dlat = (lat2_deg - lat1_deg) * DEG_TO_RAD;
            double dlon = (lon2_deg - lon1_deg) * DEG_TO_RAD;
            double a = std::sin(dlat/2)*std::sin(dlat/2)
                     + std::cos(lat1_deg*DEG_TO_RAD) * std::cos(lat2_deg*DEG_TO_RAD)
                     * std::sin(dlon/2)*std::sin(dlon/2);
            return 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0-a)) / DEG_TO_RAD;
        };
        std::vector<OccultationCandidate> filtered;
        for (const auto& res : results) {
            if (haversine_deg(obs_lat, obs_lon,
                              res.params.center_lat.to_deg(),
                              res.params.center_lon.to_deg()) < 25.0) {
                filtered.push_back(res);
            }
        }
        results = filtered;
    }

    // --- 4. Outputs ---
    if (!results.empty()) {
        // B1 fix: keep result and path paired so --star filter never creates index desync.
        // B2 fix: each event gets its own orbital elements from stored_elements.
        struct MatchedResult {
            OccultationCandidate result;
            OccultationPath      path;
            std::string          label;
        };

        const std::vector<std::string> colors = {"#ef4444", "#3b82f6", "#22c55e", "#eab308", "#8b5cf6"};
        std::vector<MatchedResult> matched;

        for (auto& res : results) {
            if (vm.count("star") && std::to_string(res.star.source_id) != vm["star"].as<std::string>()) continue;

            if (vm["star-offset-mas"].as<double>() != 0.0) {
                double offset_deg = vm["star-offset-mas"].as<double>() / 3600000.0;
                res.star.ra = RightAscension::from_deg(res.star.ra.to_deg() + offset_deg / std::cos(res.star.dec.to_rad()));
            }

            if (vm["multibody"].as<bool>()) {
                // M1: high-precision multibody refinement not yet implemented
                std::cerr << "[ioccultcalc] Warning: --multibody requested but not yet implemented; using standard propagation.\n";
            }

            double diam = DEFAULT_ASTEROID_DIAMETER_KM;
            auto it = stored_props.find(res.asteroid_id);
            if (it != stored_props.end()) diam = it->second.diameter_km;

            time::EpochUTC tca_utc = time::to_utc(res.params.t_ca);
            auto path = OccultationMapper::compute_path(res.params, res.star.ra, res.star.dec,
                                                        physics::Distance::from_km(diam), tca_utc, ephem_ptr);
            matched.push_back({res, path, res.asteroid_id + " - " + std::to_string(res.star.source_id)});
        }

        // Collect flat arrays for batch SVG export
        std::vector<OccultationPath> paths;
        std::vector<std::string> labels;
        paths.reserve(matched.size());
        labels.reserve(matched.size());
        for (const auto& m : matched) {
            paths.push_back(m.path);
            labels.push_back(m.label);
        }

        std::string xml_file = vm.count("xml-output") ? vm["xml-output"].as<std::string>() : adv_cfg.get<std::string>("xml-output", "");
        if (!xml_file.empty()) {
            std::vector<OccultationEvent> events;
            for (const auto& m : matched) {
                // B2 fix: use per-asteroid elements from stored_elements
                physics::KeplerianStateTyped<core::ECLIPJ2000> el;
                auto el_it = stored_elements.find(m.result.asteroid_id);
                if (el_it != stored_elements.end()) {
                    el = el_it->second;
                } else if (engine.has_orbit()) {
                    el = engine.orbit();
                }
                // Diameter fix: use real H mag and diameter from Horizons, not fallback
                double ev_h_mag  = 0.0;
                double ev_diam   = DEFAULT_ASTEROID_DIAMETER_KM;
                std::string ev_name;
                auto props_it = stored_props.find(m.result.asteroid_id);
                if (props_it != stored_props.end()) {
                    ev_h_mag = props_it->second.h_mag;
                    ev_diam  = props_it->second.diameter_km;
                    ev_name  = props_it->second.name;
                }
                events.push_back(candidate_to_event(m.result, m.result.asteroid_id, el,
                                                    ev_diam, ev_h_mag, ev_name));
            }
            std::filesystem::path p = std::filesystem::path(out_dir) / xml_file;
            OccultationXMLIO::write_file(events, p.string());
            std::cout << "[ioccultcalc] Saved " << events.size() << " events to " << p.string() << "\n";
        }

        std::string svg_file = vm.count("svg-output") ? vm["svg-output"].as<std::string>() : adv_cfg.get<std::string>("svg-output", "");
        if (!svg_file.empty()) {
            std::filesystem::path p = std::filesystem::path(out_dir) / svg_file;
            double zoom    = vm["zoom"].as<double>();
            double map_lat = vm.count("map-lat") ? vm["map-lat"].as<double>() : obs_lat;
            double map_lon = vm.count("map-lon") ? vm["map-lon"].as<double>() : obs_lon;
            OccultationMapper::export_global_svg(paths, labels, colors, p.string(), ephem_ptr,
                                                 Angle::from_deg(map_lat), Angle::from_deg(map_lon), zoom);
        }

        if (!out_dir.empty()) {
            for (size_t i = 0; i < matched.size(); ++i) {
                const auto& m = matched[i];
                std::stringstream ss;
                ss << prefix << "_" << m.result.asteroid_id << "_"
                   << m.result.star.source_id << "_" << (int)m.result.params.t_ca.mjd() << ".svg";
                std::filesystem::path ip = std::filesystem::path(out_dir) / ss.str();
                OccultationMapper::export_global_svg({m.path}, {m.label}, {colors[i % colors.size()]},
                                                     ip.string(), ephem_ptr,
                                                     m.result.params.center_lat, m.result.params.center_lon, 4.0);
            }
        }
    }

    catalog::GaiaDR3Catalog::shutdown();
    return 0;
}
