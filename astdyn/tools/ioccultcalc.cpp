#include "astdyn/AstDyn.hpp"
#include "astdyn/astrometry/OrbitingChebyshevEphemeris.hpp"
#include "astdyn/observations/ObservatoryDatabase.hpp"
#include "astdyn/observations/RWOReader.hpp"
#include <fstream>
#include "astdyn/io/OccultationXMLIO.hpp"
#include "astdyn/io/OccultationJSONIO.hpp"
#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/astrometry/OccultationEvent.hpp"
#include "astdyn/astrometry/ChebyshevEphemerisManager.hpp"
#include "astdyn/io/CovarianceIO.hpp"
#include "astdyn/io/AstDysOrbitFitter.hpp"
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
/// Una magnitudine e' utilizzabile? Il catalogo puo' non avere tutte le bande:
/// il nostro estratto SQLite di Gaia DR3 porta solo la G, e i campi BP/RP
/// restano a zero. Zero non e' una magnitudine plausibile — la stella piu'
/// luminosa del cielo ha G = -1.5 — quindi va trattato come "assente".
static bool magnitudine_disponibile(double m) {
    // Lo ZERO va escluso esplicitamente: e' il valore che restano i campi mai
    // popolati (GaiaStar non li inizializza e il lettore SQLite non li scrive,
    // perche' il catalogo ha la sola banda G). Escluderlo costa la possibilita'
    // di rappresentare una stella di magnitudine esattamente 0.00 — Vega ne ha
    // 0.03 — ma quel caso non si distingue comunque da un campo vuoto.
    if (m == 0.0) return false;
    return m > -2.0 && m < 30.0;
}

/// Denominazione posizionale Jhhmmss.ss+ddmmss.s a partire dalle coordinate.
///
/// Occult4 non memorizza il source_id di Gaia — il suo record e' di 42 byte e
/// contiene solo posizione, moti propri, magnitudini e incertezze — quindi
/// identifica le stelle per posizione. Un identificativo "Gaia DR3 <numero>"
/// non viene riconosciuto, e i suoi spazi rompono gli indirizzi che
/// OccultWatcher costruisce concatenandolo.
static std::string denominazione_posizionale(double ra_ore, double dec_gradi) {
    // RA in ore, minuti, secondi con due decimali
    double r = ra_ore;
    if (r < 0.0) r += 24.0;
    int rh = static_cast<int>(r);
    double rm_f = (r - rh) * 60.0;
    int rm = static_cast<int>(rm_f);
    double rs = (rm_f - rm) * 60.0;
    // il centesimo di secondo puo' arrotondare a 60: si riporta
    if (rs >= 59.995) { rs = 0.0; if (++rm == 60) { rm = 0; if (++rh == 24) rh = 0; } }

    // Dec in gradi, primi, secondi con un decimale
    const char segno = (dec_gradi < 0.0) ? '-' : '+';
    double d = std::abs(dec_gradi);
    int dd = static_cast<int>(d);
    double dm_f = (d - dd) * 60.0;
    int dm = static_cast<int>(dm_f);
    double ds = (dm_f - dm) * 60.0;
    if (ds >= 59.95) { ds = 0.0; if (++dm == 60) { dm = 0; ++dd; } }

    // Convenzione verificata sul riferimento J175336.90-214720.9, ottenuto da
    // RA 17.89358196 h (36.895056 s) e Dec -21.78915475 deg (20.957100 arcsec):
    // l'ascensione retta ARROTONDA, la declinazione TRONCA verso lo zero.
    // Troncare la declinazione evita che un raffinamento di un centesimo di
    // arcsec faccia slittare il primo e cambi il nome della stella.
    ds = std::floor(ds * 10.0) / 10.0;

    char buf[32];
    std::snprintf(buf, sizeof(buf), "J%02d%02d%05.2f%c%02d%02d%04.1f",
                  rh, rm, rs, segno, dd, dm, ds);
    return buf;
}

/// Valore che Occult4 usa per una magnitudine non nota.
constexpr double kMagAssente = 99.0;

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
    ev.elements_source = "AstDyn";
    ev.duration_s   = pr.max_duration.to_seconds();
    ev.year = ey; ev.month = em; ev.day = ed;
    ev.ut_closest_h = efrac * 24.0;
    ev.x  = pr.xi_ca.to_m()  / kEarthRadiusM;
    ev.y  = pr.eta_ca.to_m() / kEarthRadiusM;
    ev.dx = pr.dxi_dt.to_ms()  * 3600.0 / kEarthRadiusM;
    ev.dy = pr.deta_dt.to_ms() * 3600.0 / kEarthRadiusM;
    // Termini di ordine superiore, nelle unita' del formato: raggi terrestri per
    // ora elevata alla potenza corrispondente. Il polinomio di Occult4 e'
    //     x(t) = x + dx t + d2x t^2 + d3x t^3
    // quindi i coefficienti includono gia' il fattoriale.
    constexpr double ORA = 3600.0;
    ev.d2x = 0.5 * pr.d2xi_dt2  * ORA * ORA / kEarthRadiusM;
    ev.d2y = 0.5 * pr.d2eta_dt2 * ORA * ORA / kEarthRadiusM;
    ev.d3x = pr.d3xi_dt3  * ORA * ORA * ORA / (6.0 * kEarthRadiusM);
    ev.d3y = pr.d3eta_dt3 * ORA * ORA * ORA / (6.0 * kEarthRadiusM);

    // ---- <Earth> ----------------------------------------------------------
    ev.substellar_lon_deg = pr.substar_lon.to_deg();
    ev.substellar_lat_deg = pr.substar_lat.to_deg();   // geocentric == apparent Dec
    ev.subsolar_lon_deg   = pr.subsolar_lon.to_deg();
    ev.subsolar_lat_deg   = pr.subsolar_lat.to_deg();
    ev.jwst = false;

    // ---- <Star> -----------------------------------------------------------
    // The spec asks for the BCRS position at the epoch of the event with NO
    // parallax applied, which is exactly what predict_at() returns when it is
    // called without an observer position.
    const auto s_ep = cand.star.predict_at(pr.t_ca);
    ev.star_ra_h    = s_ep.ra().to_deg() / 15.0;
    ev.star_dec_deg = s_ep.dec().to_deg();
    // Identificativo posizionale: Occult4 non conosce i source_id di Gaia.
    // Il source_id resta nel JSON nativo, dove non crea problemi.
    ev.star_id = denominazione_posizionale(ev.star_ra_h, ev.star_dec_deg);
    // Gaia BP/G/RP sono i surrogati piu' vicini a B/V/R. Il catalogo puo'
    // fornirne solo alcune: quelle assenti si dichiarano, non si azzerano.
    ev.mag_b = magnitudine_disponibile(cand.star.bp_mag) ? cand.star.bp_mag : kMagAssente;
    ev.mag_v = magnitudine_disponibile(cand.star.g_mag)  ? cand.star.g_mag  : kMagAssente;
    ev.mag_r = magnitudine_disponibile(cand.star.rp_mag) ? cand.star.rp_mag : kMagAssente;
    ev.star_diameter_mas = 0.0;   // not modelled
    ev.double_star_code  = 0;
    ev.k2_flag           = "";
    {
        Angle ra_app, dec_app;
        coordinates::apparent_place(t_tt, s_ep.ra(), s_ep.dec(), ra_app, dec_app);
        ev.star_app_ra_h    = ra_app.to_deg() / 15.0;
        ev.star_app_dec_deg = dec_app.to_deg();
    }
    ev.mag_drops_adjusted  = 0;
    // Nessun controllo sulle stelle vicine viene eseguito. La specifica
    // occelmnt: "Set to -1 if a check for nearby stars has not occurred".
    // Zero direbbe che abbiamo cercato entro 15 arcsec senza trovare nulla.
    ev.bright_nearby_count = -1;
    ev.total_nearby_count  = -1;

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

    // Calo di magnitudine: durante l'evento si vede il solo asteroide, prima la
    // luce combinata. DEVE stare dopo l'assegnazione di object_mag: calcolarlo
    // prima significava usare magnitudine zero — un corpo piu' luminoso di Vega —
    // e produceva 0.00 invece di 10.69 sul caso 820987.
    // Ogni banda si calcola solo se la magnitudine stellare corrispondente c'e':
    // con mag_r assente il calcolo dava 21.11, cioe' un calo di ventun magnitudini.
    if (ev.object_mag > -5.0) {
        if (magnitudine_disponibile(ev.mag_v)) {
            const double comb_v = -2.5 * std::log10(std::pow(10.0, -0.4 * ev.mag_v) +
                                                    std::pow(10.0, -0.4 * ev.object_mag));
            ev.mag_drop_v = ev.object_mag - comb_v;
        }
        if (magnitudine_disponibile(ev.mag_r)) {
            const double comb_r = -2.5 * std::log10(std::pow(10.0, -0.4 * ev.mag_r) +
                                                    std::pow(10.0, -0.4 * ev.object_mag));
            ev.mag_drop_r = ev.object_mag - comb_r;
        }
    }
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

    // Stringhe canoniche della specifica occelmnt: Star+Assumed, Star+PEU,
    // Known errors, Assumed. "Known errors" richiede covarianze misurate per
    // ENTRAMBI, stella e asteroide; noi abbiamo quella dell'asteroide e
    // l'incertezza di picco dell'effemeride, quindi Star+PEU.
    ev.error_basis = have_cov ? "Star+PEU" : "Star+Assumed";

    // La specifica prescrive -1 per "non impostato". Scrivere 0 significa
    // affermare di aver eseguito un controllo che non facciamo: per i conteggi
    // di stelle vicine 0 vuol dire "cercate, nessuna trovata", e per il RUWE
    // vuol dire un valore nullo, impossibile (una stella singola ben misurata
    // ha RUWE ~1).
    ev.reliability = (cand.star.ruwe > 0.0) ? cand.star.ruwe : -1.0;
    ev.duplicate_source    = -1;
    ev.non_gaia_pm         = -1;
    ev.pm_added_from_ucac4 = -1;

    // L'indice di nonlinearita' e' una quantita' SCOPE, non del formato occelmnt:
    // viene riportato nel JSON nativo, in un campo proprio. Metterlo nella stringa
    // della sorgente la rendeva diversa a ogni evento, impedendo a chi legge il
    // file di riconoscere le predizioni come provenienti dalla stessa origine.
    if (have_cov && pr.nonlinearity_index > 0.0) {
        ev.nonlinearity_index = pr.nonlinearity_index;
    }

    // ---- <ID> -------------------------------------------------------------
    // Format: yyyymmdd_xxxxxx -- the event date plus the tail of the star id.
    {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%04d%02d%02d", ey, em, ed);
        // La coda va presa dall'identificativo che SCRIVIAMO, non dal source_id
        // di Gaia: e' cosi' che OccultWatcher riconosce due predizioni come
        // riferite allo stesso evento. Con J175336.90-214720.9 la coda e'
        // 4720.9, che coincide con quella delle predizioni di riferimento.
        const std::string& sid = ev.star_id;
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

// Converte una specifica temporale flessibile in JD (TDB).
// Accetta: "YYYY-MM-DD" o "YYYY MM DD.ddd" (calendario), oppure una stringa
// numerica interpretata come JD se >= 2400000, altrimenti come MJD.
// Per disambiguare esplicitamente si possono usare i prefissi "jd:" / "mjd:".
static double parse_epoch_to_jd(const std::string& spec_in) {
    std::string spec = spec_in;
    // trim
    auto a = spec.find_first_not_of(" \t");
    auto b = spec.find_last_not_of(" \t");
    if (a == std::string::npos) return 0.0;
    spec = spec.substr(a, b - a + 1);

    auto starts_with = [&](const std::string& p) {
        return spec.size() >= p.size() && spec.compare(0, p.size(), p) == 0;
    };

    if (starts_with("jd:"))  return std::stod(spec.substr(3));
    if (starts_with("mjd:")) return time::mjd_to_jd(std::stod(spec.substr(4)));

    // Calendario: contiene '-' (YYYY-MM-DD) o spazi (YYYY MM DD)
    bool looks_calendar = (spec.find('-') != std::string::npos && spec.size() >= 8)
                       || (spec.find(' ') != std::string::npos);
    if (looks_calendar) {
        std::string norm = spec;
        for (char& c : norm) if (c == '-' || c == '/') c = ' ';
        std::istringstream iss(norm);
        int year = 0, month = 0; double day = 1.0;
        iss >> year >> month >> day;
        if (year > 0 && month > 0) {
            int day_int = static_cast<int>(day);
            double frac = day - day_int;
            return time::mjd_to_jd(time::calendar_to_mjd(year, month, day_int, frac));
        }
    }

    // Numerico puro: JD se grande, altrimenti MJD.
    double v = std::stod(spec);
    return (v >= 2400000.0) ? v : time::mjd_to_jd(v);
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
        ("json-output", po::value<std::string>(), "save found occultations to native JSON")
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
    // Selezione oggetti: 'objects.asteroids' (config dinamica) ha priorita',
    // poi il flag CLI --asteroid, poi la vecchia chiave piatta 'asteroid'.
    // Tutte le forme (lista, range "100-34244", "@file") passano da
    // parse_asteroid_list.
    // "*" significa "da 1 a objects.all_max": lo traduciamo in un range e lo
    // passiamo a parse_asteroid_list, che i range gia' li gestisce. Il default
    // di all_max e' volutamente modesto (1000): "*" non deve avviare per
    // sbaglio una campagna su centinaia di migliaia di corpi (ore di calcolo).
    // Chi vuole coprire di piu' alza objects.all_max, oppure -- meglio -- usa
    // un range esplicito (es. "800000-821000") mirato alla zona di interesse.
    auto expand_star = [&](std::string spec) -> std::string {
        if (spec == "*") {
            int all_max = adv_cfg.get<int>("objects.all_max", 1000);
            return "1-" + std::to_string(all_max);
        }
        return spec;
    };

    std::vector<std::string> asteroid_ids;
    if (adv_cfg.has("objects.asteroids")) {
        asteroid_ids = parse_asteroid_list(expand_star(adv_cfg.get<std::string>("objects.asteroids", "")));
    } else if (vm.count("asteroid")) {
        asteroid_ids = parse_asteroid_list(expand_star(vm["asteroid"].as<std::string>()));
    } else if (adv_cfg.has("asteroid")) {
        asteroid_ids = parse_asteroid_list(expand_star(adv_cfg.get<std::string>("asteroid", "")));
    }

    // Finestra temporale: 'time.start' (calendario YYYY-MM-DD, "mjd:" o "jd:")
    // ha priorita', poi --jd-start CLI, poi la vecchia chiave 'jd-start'.
    double jd_start = 0.0;
    if (adv_cfg.has("time.start")) {
        jd_start = parse_epoch_to_jd(adv_cfg.get<std::string>("time.start", ""));
    } else if (vm.count("jd-start")) {
        jd_start = vm["jd-start"].as<double>();
    } else if (adv_cfg.has("jd-start")) {
        jd_start = adv_cfg.get<double>("jd-start", 0.0);
    }


    if ((asteroid_ids.empty() && !vm.count("system-ids") && !adv_cfg.has("system-ids")) || jd_start == 0.0) {
        if (!vm.count("help")) std::cerr << "Error: Missing asteroid list or jd-start (check CLI or config file).\n";
        std::cout << desc << "\n";
        return 1;
    }

    // Durata: 'time.end' (stessa sintassi di start) -> durata = end - start;
    // oppure 'time.duration_days'; poi --duration CLI; poi 'duration' piatta.
    double duration;
    if (adv_cfg.has("time.end")) {
        double jd_end = parse_epoch_to_jd(adv_cfg.get<std::string>("time.end", ""));
        duration = jd_end - jd_start;
    } else if (adv_cfg.has("time.duration_days")) {
        duration = adv_cfg.get<double>("time.duration_days", 1.0);
    } else {
        duration = (!vm["duration"].defaulted()) ? vm["duration"].as<double>() : adv_cfg.get<double>("duration", 1.0);
    }
    double mag_limit = (!vm["mag"].defaulted()) ? vm["mag"].as<double>() : adv_cfg.get<double>("mag", 15.0);
    
    std::string out_dir = vm.count("out-dir") ? vm["out-dir"].as<std::string>() : adv_cfg.get<std::string>("out-dir", "");
    std::string prefix = (!vm["prefix"].defaulted()) ? vm["prefix"].as<std::string>() : adv_cfg.get<std::string>("prefix", "occ");

    if (!out_dir.empty()) {
        std::filesystem::create_directories(out_dir);
    }

    time::EpochTDB start_epoch = time::EpochTDB::from_jd(jd_start);
    time::EpochTDB end_epoch = time::EpochTDB::from_jd(jd_start + duration);

    // Fase 0 (orchestratore): se objects.elements_dir e' impostato, gli elementi
    // di ciascun corpo vengono letti da <elements_dir>/<id>.eq1 (formato AstDyS)
    // invece che da Horizons. Il tool Python prepara gli .eq1 grezzi; qui li
    // convertiamo con la catena tipizzata (read_eq1 -> to_si ->
    // equinoctial_to_keplerian), che rispetta frame e unita'. Fallback a
    // Horizons se il file manca o e' illeggibile.
    std::string elements_dir = adv_cfg.get<std::string>("objects.elements_dir", "");
    // Fase 5: directory dei .rwo per il fit orbitale. L'orchestrator li scarica
    // per i soli positivi del secondo stadio.
    // Catalogo osservatori MPC: senza, le osservazioni vengono ridotte dal
    // geocentro e si perde la parallasse topocentrica (~4" a 2.2 AU).
    {
        std::string obs_file = adv_cfg.get<std::string>("observatories.file", "");
        if (obs_file.empty()) {
            const char* home = std::getenv("HOME");
            const std::vector<std::string> candidati = {
                std::string(home ? home : "") + "/.ioccultcalc/observatories/ObsCodes.txt",
                std::string(home ? home : "") + "/.ioccultcalc/ObsCodes.txt",
                "ObsCodes.txt"
            };
            for (const auto& c : candidati) {
                std::ifstream prova(c);
                if (prova.good()) { obs_file = c; break; }
            }
        }
        if (obs_file.empty()) {
            std::cout << "[osservatori] ATTENZIONE: catalogo MPC non trovato. Le osservazioni "
                         "saranno ridotte dal GEOCENTRO, perdendo la parallasse topocentrica. "
                         "Impostare 'observatories.file' oppure collocare ObsCodes.txt in "
                         "~/.ioccultcalc/observatories/\n";
        } else {
            try {
                size_t n = observations::ObservatoryDatabase::getInstance().loadFromMPCFile(obs_file);
                std::cout << "[osservatori] " << n << " codici da " << obs_file << "\n";
            } catch (const std::exception& e) {
                std::cout << "[osservatori] ATTENZIONE: lettura fallita (" << e.what()
                          << "): riduzione dal geocentro\n";
            }
        }
    }

    std::string observations_dir = adv_cfg.get<std::string>("objects.observations_dir", "");
    const bool fit_attivo = adv_cfg.get<bool>("diffcorr.enabled", false);
    if (fit_attivo && observations_dir.empty()) {
        std::cout << "[fit] ATTENZIONE: diffcorr.enabled=true ma objects.observations_dir "
                     "non e' impostata: nessun fit verra' eseguito.\n";
    }

    ChebyshevEphemerisManager manager(engine.config());
    HorizonsClient horizons;
    std::map<std::string, physics::CartesianStateTyped<core::ECLIPJ2000>> stored_states;
    std::map<std::string, physics::KeplerianStateTyped<core::ECLIPJ2000>> stored_elements;
    // Covarianza per-asteroide (cartesiana ECLIPJ2000, AU/day), dagli .eq1
    // che la portano (AstDyS). Assente per .eq1 senza covarianza (DB).
    std::map<std::string, astdyn::Matrix6d> stored_covariances;
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
            // Fase 0: prima l'elemento locale <elements_dir>/<id>.eq1, poi Horizons.
            // Stesso tipo di query_elements (expected), cosi' il resto del loop
            // (if(!elements_opt), *elements_opt) resta invariato.
            // Covarianza cartesiana (AU/day) calcolata dall'.eq1 se presente,
            // trasferita al blocco storage per stored_covariances[id].
            std::optional<astdyn::Matrix6d> eq1_cov_cart_opt;
                        std::expected<physics::KeplerianStateTyped<core::ECLIPJ2000>, io::HorizonsError> elements_opt =
                std::unexpected(io::HorizonsError{});
            std::string ef = elements_dir.empty() ? std::string()
                                                  : (elements_dir + "/" + id + ".eq1");
            if (!ef.empty() && std::filesystem::exists(ef)) {
                try {
                    // read_eq1 da' gli elementi equinoziali in unita' AstDyS grezze:
                    // a in AU, angoli in rad/deg. equinoctial_to_keplerian vuole
                    // proprio a in AU (converte lei internamente), quindi NON si
                    // applica to_si. Gli elementi portano la loro epoca (61200);
                    // add_asteroid/fit_chebyshev propagano da elements.epoch alla
                    // finestra usando la via gia' collaudata -- non propaghiamo qui.
                    auto orb = io::CovarianceIO::read_eq1(ef);
                    auto t0 = time::to_tdb(time::EpochTT::from_mjd(orb.epoch_mjd_tt));
                    // Se l'.eq1 porta la covarianza (AstDyS), la trasformo in
                    // cartesiana AU/day con la STESSA via del flusso --covariance:
                    // to_si (a->km, lambda->rad) -> EquinoctialElements -> jacobiana.
                    // Lavoro su una COPIA (orb_cov) per non alterare gli elementi
                    // grezzi usati subito sotto per equinoctial_to_keplerian.
                    if (!orb.covariance.isZero(0.0)) {
                        try {
                            constexpr double kAuKm_c = 149597870.700;
                            io::AstDysOrbit orb_cov = orb;
                            io::CovarianceIO::to_si(orb_cov);   // a->km, lambda->rad
                            coordinates::EquinoctialElements eqc(
                                orb_cov.elements(0), orb_cov.elements(1), orb_cov.elements(2),
                                orb_cov.elements(3), orb_cov.elements(4), orb_cov.elements(5));
                            const astdyn::Matrix6d Jc = eqc.jacobian_to_cartesian();
                            const astdyn::Matrix6d C_km_c = Jc * orb_cov.covariance * Jc.transpose();
                            Eigen::Matrix<double, 6, 1> dc;
                            dc << 1.0/kAuKm_c, 1.0/kAuKm_c, 1.0/kAuKm_c,
                                  86400.0/kAuKm_c, 86400.0/kAuKm_c, 86400.0/kAuKm_c;
                            eq1_cov_cart_opt = dc.asDiagonal() * C_km_c * dc.asDiagonal();
                        } catch (const std::exception& ce) {
                            std::cerr << "[ioccultcalc] '" << id << "': covarianza .eq1 non trasformata ("
                                      << ce.what() << ")\n";
                        }
                    }
                    // NB: read_eq1 da' a in AU e lambda in GRADI (unita' AstDyS
                    // grezze). equinoctial_to_keplerian vuole a in AU (ok) ma
                    // lambda in RADIANTI (lo usa in M = lambda - varpi con varpi
                    // in rad). to_si() faceva 'a->km, lambda->rad': serviva solo
                    // la seconda meta'. Convertiamo qui solo lambda.
                    elements_opt = io::AstDysOrbitFitter::equinoctial_to_keplerian(
                        orb.elements(0), orb.elements(1), orb.elements(2),
                        orb.elements(3), orb.elements(4),
                        orb.elements(5) * constants::DEG_TO_RAD, t0);
                    std::cout << "[ioccultcalc] '" << id << "' caricato da " << ef
                              << " (elementi @MJD " << orb.epoch_mjd_tt << " TT)\n";
                                } catch (const std::exception& e) {
                    std::cerr << "[ioccultcalc] '" << id << "': errore leggendo " << ef
                              << " (" << e.what() << "); fallback a Horizons\n";
                }
            }
            if (!elements_opt) {
                elements_opt = horizons.query_elements(id, start_epoch);
            }
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

                    // --- Fase 5: fit orbitale on-demand ---------------------
                    // Se il fit e' richiesto e per questo corpo esistono le
                    // osservazioni, si raffina l'orbita. L'orbita fittata
                    // sostituisce quella di partenza solo se il fit converge:
                    // un fit fallito non deve peggiorare un'orbita buona.
                    std::optional<astdyn::Matrix6d> cov_dal_fit;
                    if (fit_attivo && !observations_dir.empty()) {
                        const std::string rwo = observations_dir + "/" + id + ".rwo";
                        std::ifstream prova(rwo);
                        if (!prova.good()) {
                            std::cout << "[fit] " << id << ": osservazioni assenti ("
                                      << rwo << "), orbita AstDyS mantenuta\n";
                        } else {
                            prova.close();
                            try {
                                auto osservazioni = observations::RWOReader::readFile(rwo);
                                AstDynEngine motore_fit;
                                auto cfg_fit = engine.config();
                                cfg_fit.verbose = false;
                                motore_fit.set_config(cfg_fit);
                                motore_fit.set_initial_orbit(elements);
                                for (const auto& o : osservazioni) motore_fit.add_observation(o);
                                auto esito = motore_fit.fit_orbit();
                                if (esito.converged) {
                                    elements = esito.orbit;
                                    // Covarianza dal fit: sostituisce quella dell'.eq1 perche'
                                    // descrive l'orbita che stiamo effettivamente usando.
                                    // Orbita e covarianza restano accoppiate.
                                    if (esito.covariance.rows() == 6 && esito.covariance.cols() == 6) {
                                        astdyn::Matrix6d cov_fit;
                                        for (int r = 0; r < 6; ++r)
                                            for (int c = 0; c < 6; ++c)
                                                cov_fit(r, c) = esito.covariance(r, c);
                                        cov_dal_fit = cov_fit;
                                    }
                                    std::cout << "[fit] " << id << ": orbita raffinata su "
                                              << esito.num_observations << " osservazioni, RMS "
                                              << esito.rms_total_arcsec << " arcsec ("
                                              << esito.num_rejected << " scartate)"
                                              << (cov_dal_fit ? ", covarianza dal fit" : "")
                                              << "\n";
                                } else {
                                    std::cout << "[fit] " << id << ": NON convergente ("
                                              << esito.exit_reason << "), orbita AstDyS mantenuta\n";
                                }
                            } catch (const std::exception& e) {
                                std::cout << "[fit] " << id << ": errore (" << e.what()
                                          << "), orbita AstDyS mantenuta\n";
                            }
                        }
                    }
                    // --------------------------------------------------------

                    manager.add_asteroid(id, elements, start_epoch, end_epoch);
                    stored_elements[id] = elements;

                    // The state IS heliocentric ecliptic -- the variable is named for it.
                    auto state_eclip = propagation::keplerian_to_cartesian(elements);
                    stored_states[id] = state_eclip;
                    // La covarianza del fit ha la precedenza: descrive l'orbita in uso.
                    if (cov_dal_fit)              stored_covariances[id] = *cov_dal_fit;
                    else if (eq1_cov_cart_opt)    stored_covariances[id] = *eq1_cov_cart_opt;

                    // Il messaggio di provenienza (Horizons o .eq1) e' gia' stato
                    // stampato nel ramo di caricamento; qui non ripetiamo "Horizons".
                    auto props = horizons.query_physical_properties(id);
                    if (props && props->diameter_km > 0.0) {
                        stored_props[id] = *props;
                        manager.set_diameter(id, props->diameter_km);
                    } else {
                        // props assente, oppure presente ma con diametro non
                        // parsato (0): in entrambi i casi il diametro e' ignoto.
                        // Usa il default, altrimenti il filtro diametro vedrebbe
                        // 0 e non potrebbe mai scartare il corpo.
                        double d = (props && props->diameter_km > 0.0)
                                 ? props->diameter_km : DEFAULT_ASTEROID_DIAMETER_KM;
                        io::PhysicalProperties pp = props ? *props : io::PhysicalProperties{"", 0.0, d, 0.0};
                        pp.diameter_km = d;
                        stored_props[id] = pp;
                        manager.set_diameter(id, d);
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

            // Proprieta' fisiche dal blocco `system_bodies` della configurazione.
            // Il diametro determina la larghezza dell'ombra e la durata
            // dell'evento; la magnitudine determina il calo di luce. I valori
            // predefiniti (10 km, H=15) sono segnaposto: per Hi'iaka, che misura
            // ~320 km, sbaglierebbero di un fattore 32.
            const std::string chiave = "system_bodies." + id_str;
            double diam = adv_cfg.get<double>(chiave + ".diameter_km", 0.0);
            double hmag = adv_cfg.get<double>(chiave + ".H", 999.0);
            std::string nome = adv_cfg.get<std::string>(chiave + ".name", "");

            if (diam <= 0.0) {
                diam = DEFAULT_SATELLITE_DIAMETER_KM;
                std::cout << "[sistema] " << id_str << ": diametro non specificato, uso "
                          << DEFAULT_SATELLITE_DIAMETER_KM << " km (SEGNAPOSTO). "
                          << "La larghezza dell'ombra e la durata dell'evento saranno "
                          << "inattendibili: indicare il valore reale in system_bodies."
                          << id_str << ".diameter_km\n";
            }
            if (hmag > 900.0) {
                hmag = DEFAULT_SATELLITE_H_MAG;
                std::cout << "[sistema] " << id_str << ": magnitudine assoluta non "
                          << "specificata, uso H=" << DEFAULT_SATELLITE_H_MAG
                          << " (segnaposto): il calo di luce sara' indicativo\n";
            }
            if (!nome.empty() || diam != DEFAULT_SATELLITE_DIAMETER_KM) {
                std::cout << "[sistema] " << id_str
                          << (nome.empty() ? "" : " (" + nome + ")")
                          << ": diametro " << diam << " km, H=" << hmag << "\n";
            }

            stored_props[id_str] = io::PhysicalProperties{nome, hmag, diam, 0.0};
            manager.set_diameter(id_str, diam);
        }
    }

    // --- Fase 6: satelliti descritti dalla loro orbita mutua -----------------
    // Per la gran parte dei binari noti non esiste un kernel SPK: la letteratura
    // pubblica i parametri dell'orbita del satellite attorno al primario.
    {
        const auto elenco = adv_cfg.get_keys("system_bodies");
        for (const auto& sat_id : elenco) {
            const std::string base = "system_bodies." + sat_id;
            if (!adv_cfg.has(base + ".orbit.a_km")) continue;   // non e' un satellite

            const std::string primario = adv_cfg.get<std::string>(base + ".primary", "");
            if (primario.empty()) {
                std::cout << "[sistema] " << sat_id << ": manca 'primary', "
                          << "impossibile collocare l'orbita mutua. Corpo ignorato.\n";
                continue;
            }
            auto it = stored_elements.find(primario);
            if (it == stored_elements.end()) {
                std::cout << "[sistema] " << sat_id << ": il primario '" << primario
                          << "' non e' fra i corpi caricati. Corpo ignorato.\n";
                continue;
            }

            astrometry::OrbitingChebyshevEphemeris::OrbitaMutua orb;
            orb.a_km        = adv_cfg.get<double>(base + ".orbit.a_km", 0.0);
            orb.e           = adv_cfg.get<double>(base + ".orbit.e", 0.0);
            orb.i_deg       = adv_cfg.get<double>(base + ".orbit.i_deg", 0.0);
            orb.node_deg    = adv_cfg.get<double>(base + ".orbit.node_deg", 0.0);
            orb.peri_deg    = adv_cfg.get<double>(base + ".orbit.peri_deg", 0.0);
            orb.M_deg       = adv_cfg.get<double>(base + ".orbit.M_deg", 0.0);
            orb.period_days = adv_cfg.get<double>(base + ".orbit.period_days", 0.0);
            orb.epoch = time::EpochTDB::from_mjd(
                adv_cfg.get<double>(base + ".orbit.epoch_mjd", it->second.epoch.mjd()));

            const std::string piano = adv_cfg.get<std::string>(base + ".orbit.plane", "equatorial");
            orb.piano = (piano == "ecliptic")
                      ? astrometry::OrbitingChebyshevEphemeris::PianoRiferimento::Eclittico
                      : astrometry::OrbitingChebyshevEphemeris::PianoRiferimento::Equatoriale;

            if (orb.period_days <= 0.0) {
                std::cout << "[sistema] " << sat_id << ": manca 'orbit.period_days', "
                          << "necessario per il parametro gravitazionale. Corpo ignorato.\n";
                continue;
            }

            // proprieta' fisiche, con gli stessi avvisi degli altri corpi di sistema
            double diam = adv_cfg.get<double>(base + ".diameter_km", 0.0);
            double hmag = adv_cfg.get<double>(base + ".H", 999.0);
            const std::string nome = adv_cfg.get<std::string>(base + ".name", "");
            if (diam <= 0.0) {
                diam = DEFAULT_SATELLITE_DIAMETER_KM;
                std::cout << "[sistema] " << sat_id << ": diametro non specificato, uso "
                          << diam << " km (SEGNAPOSTO): larghezza dell'ombra e durata "
                          << "saranno inattendibili\n";
            }
            if (hmag > 900.0) hmag = DEFAULT_SATELLITE_H_MAG;

            try {
                manager.add_orbiting_body(sat_id, it->second, orb, start_epoch, end_epoch);
                asteroid_ids.push_back(sat_id);
                stored_props[sat_id] = io::PhysicalProperties{nome, hmag, diam, 0.0};
                manager.set_diameter(sat_id, diam);

                std::cout << "[sistema] " << sat_id
                          << (nome.empty() ? "" : " (" + nome + ")")
                          << ": satellite di " << primario
                          << ", a=" << orb.a_km << " km, P=" << orb.period_days << " d, "
                          << "diametro " << diam << " km, piano " << piano << "\n";
            } catch (const std::exception& e) {
                std::cout << "[sistema] " << sat_id << ": " << e.what()
                          << ". Corpo ignorato.\n";
            }
        }
    }

    // --- 3c. Global Search ---
    OccultationConfig occ_config = engine.config().occultation_settings;
    occ_config.max_mag_star = mag_limit;

    // --- Fase 2: filtri dal blocco 'filters' della config dinamica ---
    if (adv_cfg.has("filters.diameter_min_km"))
        occ_config.min_asteroid_diameter_km = adv_cfg.get<double>("filters.diameter_min_km", 0.0);
    if (adv_cfg.has("filters.diameter_max_km"))
        occ_config.max_asteroid_diameter_km = adv_cfg.get<double>("filters.diameter_max_km", 0.0);
    if (adv_cfg.has("filters.star_mag_max"))
        occ_config.max_mag_star = adv_cfg.get<double>("filters.star_mag_max", 14.0);
    if (adv_cfg.has("filters.event_duration_s_min"))
        occ_config.min_duration_s = adv_cfg.get<double>("filters.event_duration_s_min", 0.0);
    if (adv_cfg.has("filters.drop_mag_min"))
        occ_config.min_mag_drop = adv_cfg.get<double>("filters.drop_mag_min", 0.05);
    if (adv_cfg.has("filters.max_gaia_ruwe"))
        occ_config.max_gaia_ruwe = adv_cfg.get<double>("filters.max_gaia_ruwe", 99.0);
    
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
    // La covarianza si calcola quando richiesta: via --covariance/config (cov_file),
    // oppure quando gli .eq1 in elements_dir hanno portato covarianze per-asteroide
    // (campagna AstDyS) e il calcolo e' abilitato (config physics.covariance:true).
    bool want_covariance = !cov_file.empty()
                        || (adv_cfg.get<bool>("physics.covariance", false) && !stored_covariances.empty());
    if (want_covariance && !results.empty()) {
        try {
            constexpr double kAuKm = 149597870.700;
            astdyn::Matrix6d cov_t0;
            std::optional<physics::CartesianStateTyped<core::ECLIPJ2000>> x0;

            if (!cov_file.empty() && cov_file.size() > 4 && cov_file.substr(cov_file.size() - 4) == ".eq1") {
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
            } else if (!cov_file.empty()) {
                // A bare 6x6: assumed already Cartesian, in AU and days, at the epoch
                // of the stored state.
                cov_t0 = CovarianceIO::read_file(cov_file);
            }

            for (auto& res : results) {
                if (stored_covariances.count(res.asteroid_id) && stored_states.count(res.asteroid_id)) {
                    // Covarianza PER-ASTEROIDE (dagli .eq1 AstDyS): ellisse del corpo giusto.
                    OccultationLogic::apply_uncertainty(res.params, res.star,
                                                        stored_covariances[res.asteroid_id],
                                                        stored_states[res.asteroid_id], engine);
                } else if (x0) {
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
            std::string json_file = vm.count("json-output") ? vm["json-output"].as<std::string>() : adv_cfg.get<std::string>("json-output", "");
            if (!json_file.empty()) {
                std::filesystem::path pj = std::filesystem::path(out_dir) / json_file;
                std::filesystem::create_directories(pj.parent_path());
                OccultationJSONIO::write_file(events, pj.string());
                std::cout << "[ioccultcalc] Saved " << events.size() << " events (JSON) to " << pj.string() << "\n";
            }
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

    // output.write_empty: con 0 occultazioni il blocco sopra (if !results.empty)
    // non scrive nulla. Per un batch OWC puo' servire un file comunque, per
    // sapere che il run e' andato (vuoto != fallito). Se attivo, scrive un
    // <Occultations></Occultations> valido.
    if (results.empty() && adv_cfg.get<bool>("output.write_empty", false)) {
        std::string xml_file = vm.count("xml-output") ? vm["xml-output"].as<std::string>()
                                                      : adv_cfg.get<std::string>("xml-output", "");
        if (!xml_file.empty()) {
            std::filesystem::path pth = std::filesystem::path(out_dir) / xml_file;
            OccultationXMLIO::write_file({}, pth.string());
            std::cout << "[ioccultcalc] Saved 0 events (write_empty) to " << pth.string() << "\n";
        }
    }

    catalog::GaiaDR3Catalog::shutdown();
    return 0;
}
