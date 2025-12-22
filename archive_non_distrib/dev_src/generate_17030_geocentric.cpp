/**
 * Genera coordinate equatoriali geocentriche per asteroide 17030
 * Data: 28/11/2025 dalle 00:00 alle 24:00 UTC (ogni ora)
 * Usa: AstDyn Propagator con RKF78, tutte le perturbazioni
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/propagation/Integrator.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "astdyn/core/Constants.hpp"

using namespace astdyn;
using namespace Eigen;

// Costanti
constexpr double MJD_TO_JD = 2400000.5;
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Obliquità eclittica J2000
constexpr double EPSILON_J2000 = 23.4392911 * DEG_TO_RAD;

// Conversione UTC a MJD TDB (approssimata, ignora leap seconds ~70s)
double utc_to_mjd_tdb(int year, int month, int day, int hour, int min, double sec) {
    // Formula di Meeus per JD
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    double jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    jd += (hour - 12) / 24.0 + min / 1440.0 + sec / 86400.0;
    double mjd = jd - MJD_TO_JD;
    // TDB ≈ TT ≈ UTC + 69.184s (al 2025)
    return mjd + 69.184 / 86400.0;
}

// Rotazione ECLM J2000 → ICRF (equatoriale)
Vector3d ecliptic_to_equatorial(const Vector3d& ecl) {
    double cos_eps = std::cos(EPSILON_J2000);
    double sin_eps = std::sin(EPSILON_J2000);
    return Vector3d(
        ecl.x(),
        ecl.y() * cos_eps - ecl.z() * sin_eps,
        ecl.y() * sin_eps + ecl.z() * cos_eps
    );
}

// Conversione cartesiano equatoriale → RA/Dec
void cartesian_to_radec(const Vector3d& pos, double& ra_deg, double& dec_deg, double& dist_au) {
    dist_au = pos.norm();
    dec_deg = std::asin(pos.z() / dist_au) * RAD_TO_DEG;
    ra_deg = std::atan2(pos.y(), pos.x()) * RAD_TO_DEG;
    if (ra_deg < 0) ra_deg += 360.0;
}

// Formatta RA in HH:MM:SS.ss
std::string format_ra(double ra_deg) {
    double ra_h = ra_deg / 15.0;
    int h = static_cast<int>(ra_h);
    double rem = (ra_h - h) * 60.0;
    int m = static_cast<int>(rem);
    double s = (rem - m) * 60.0;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%02d:%02d:%05.2f", h, m, s);
    return std::string(buf);
}

// Formatta Dec in ±DD:MM:SS.s
std::string format_dec(double dec_deg) {
    char sign = dec_deg >= 0 ? '+' : '-';
    dec_deg = std::fabs(dec_deg);
    int d = static_cast<int>(dec_deg);
    double rem = (dec_deg - d) * 60.0;
    int m = static_cast<int>(rem);
    double s = (rem - m) * 60.0;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%c%02d:%02d:%04.1f", sign, d, m, s);
    return std::string(buf);
}

int main() {
    std::cout << "=== Coordinate Equatoriali Geocentriche Asteroide 17030 ===\n";
    std::cout << "Data: 28 Novembre 2025, 00:00-24:00 UTC (ogni ora)\n";
    std::cout << "Integratore: RKF78, tolleranza 1e-12 AU\n";
    std::cout << "Perturbazioni: 8 pianeti + asteroidi + relatività\n";
    std::cout << "Frame output: ICRF J2000 geocentrico\n\n";

    // 1. Carica elementi orbitali da file .eq1
    std::string eq1_path = "./astdyn/data/17030.eq1";
    io::parsers::OrbFitEQ1Parser parser;
    io::IOrbitParser::OrbitalElements elements;
    
    try {
        elements = parser.parse(eq1_path);
        std::cout << "Elementi caricati da: " << eq1_path << "\n";
        std::cout << "  Epoca elementi: MJD " << std::fixed << std::setprecision(6) 
                  << elements.epoch_mjd_tdb << "\n";
        std::cout << "  a = " << elements.semi_major_axis << " AU\n";
        std::cout << "  e = " << elements.eccentricity << "\n\n";
    } catch (const std::exception& e) {
        std::cerr << "Errore caricamento .eq1: " << e.what() << "\n";
        return 1;
    }

    // 2. Configura propagatore RKF78 con tutte le perturbazioni
    auto integrator = std::make_unique<propagation::RKF78Integrator>(0.1, 1e-12);
    
    propagation::PropagatorSettings prop_settings;
    prop_settings.include_planets = true;
    prop_settings.include_relativity = true;
    prop_settings.include_asteroids = true;
    prop_settings.perturb_mercury = true;
    prop_settings.perturb_venus = true;
    prop_settings.perturb_earth = true;
    prop_settings.perturb_mars = true;
    prop_settings.perturb_jupiter = true;
    prop_settings.perturb_saturn = true;
    prop_settings.perturb_uranus = true;
    prop_settings.perturb_neptune = true;

    std::cout << "DEBUG: Perturbazioni abilitate:\n";
    std::cout << "  include_planets = " << prop_settings.include_planets << "\n";
    std::cout << "  include_relativity = " << prop_settings.include_relativity << "\n";
    std::cout << "  include_asteroids = " << prop_settings.include_asteroids << "\n";
    std::cout << "  perturb_jupiter = " << prop_settings.perturb_jupiter << "\n\n";

    auto ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
    propagation::Propagator propagator(std::move(integrator), ephemeris, prop_settings);

    // 3. Prepara elementi kepleriani per propagazione
    propagation::KeplerianElements kep;
    kep.epoch_mjd_tdb = elements.epoch_mjd_tdb;
    kep.semi_major_axis = elements.semi_major_axis;
    kep.eccentricity = elements.eccentricity;
    kep.inclination = elements.inclination;
    kep.longitude_ascending_node = elements.longitude_asc_node;
    kep.argument_perihelion = elements.argument_perihelion;
    kep.mean_anomaly = elements.mean_anomaly;
    kep.gravitational_parameter = constants::GMS;

    // 4. Genera coordinate per ogni ora
    std::cout << std::setw(6) << "Ora" 
              << std::setw(14) << "MJD_TDB"
              << std::setw(14) << "RA (J2000)"
              << std::setw(14) << "Dec (J2000)"
              << std::setw(12) << "Dist [AU]"
              << "\n";
    std::cout << std::string(60, '-') << "\n";

    // File CSV output
    std::ofstream csv("/tmp/17030_geocentric_28nov2025.csv");
    csv << "Hour_UTC,MJD_TDB,RA_deg,Dec_deg,RA_HMS,Dec_DMS,Distance_AU,X_AU,Y_AU,Z_AU\n";
    csv << std::fixed << std::setprecision(8);

    auto t_start = std::chrono::high_resolution_clock::now();
    
    for (int hour = 0; hour <= 24; ++hour) {
        double mjd_tdb = utc_to_mjd_tdb(2025, 11, 28, hour, 0, 0.0);
        double jd_tdb = mjd_tdb + MJD_TO_JD;

        try {
            // Propaga asteroide a epoca target
            auto kep_prop = propagator.propagate_keplerian(kep, mjd_tdb);
            auto cart_ecl = propagation::keplerian_to_cartesian(kep_prop);
            
            // Posizione asteroide in ICRF (baricentrico)
            Vector3d ast_bary_icrf = ecliptic_to_equatorial(cart_ecl.position);
            
            // Posizione Terra (eliocentrica, eclittica) - converti a ICRF
            Vector3d earth_helio_ecl = ephemeris::PlanetaryEphemeris::getPosition(
                ephemeris::CelestialBody::EARTH, jd_tdb);
            Vector3d earth_helio_icrf = ecliptic_to_equatorial(earth_helio_ecl);
            
            // Vettore geocentrico = asteroide_bary - terra_helio
            // (approssimazione: ignora offset sole-baricentro ~0.01 AU)
            Vector3d ast_geo_icrf = ast_bary_icrf - earth_helio_icrf;
            
            // Converti in RA/Dec
            double ra_deg, dec_deg, dist_au;
            cartesian_to_radec(ast_geo_icrf, ra_deg, dec_deg, dist_au);
            
            std::string ra_hms = format_ra(ra_deg);
            std::string dec_dms = format_dec(dec_deg);
            
            // Output console
            std::cout << std::setw(4) << hour << ":00"
                      << std::setw(14) << std::fixed << std::setprecision(6) << mjd_tdb
                      << std::setw(14) << ra_hms
                      << std::setw(14) << dec_dms
                      << std::setw(12) << std::setprecision(6) << dist_au
                      << "\n";
            
            // Output CSV
            csv << hour << ","
                << mjd_tdb << ","
                << ra_deg << ","
                << dec_deg << ","
                << ra_hms << ","
                << dec_dms << ","
                << dist_au << ","
                << ast_geo_icrf.x() << ","
                << ast_geo_icrf.y() << ","
                << ast_geo_icrf.z() << "\n";
                
        } catch (const std::exception& e) {
            std::cerr << "Errore ora " << hour << ": " << e.what() << "\n";
        }
    }
    
    csv.close();
    
    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
    
    std::cout << "\n✓ CSV salvato: /tmp/17030_geocentric_28nov2025.csv\n";
    std::cout << "⏱️  Tempo totale propagazione (25 epoche): " << duration.count() << " ms\n";
    std::cout << "   Tempo medio per epoca: " << (duration.count() / 25.0) << " ms\n";
    
    return 0;
}
