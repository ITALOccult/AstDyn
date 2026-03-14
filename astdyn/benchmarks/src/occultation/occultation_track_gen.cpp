/**
 * @file occultation_track_gen.cpp
 * @brief Generates shadow ground tracks for comparison.
 */

#include "astdyn/AstDyn.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/io/HorizonsClient.hpp"
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/catalog/GaiaDR3Catalog.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace astdyn;
using namespace astdyn::physics;
using namespace astdyn::astrometry;
using namespace astdyn::time;

struct TrackPoint {
    double jd;
    double lat;
    double lon;
};

// Simple GCRF to WGS84 (Spherical approximation for visualization)
TrackPoint project_to_earth(const Eigen::Vector3d& p_ast, const Eigen::Vector3d& u_star, double jd) {
    // Light rays travel in direction -u_star
    // p_shadow = p_ast + k * (-u_star)  s.t. |p_shadow| = R_earth
    double R = 6371.0 * 1000.0;
    
    // Intersection of line (p_ast, -u_star) with sphere |p| = R
    // |p_ast - k*u_star|^2 = R^2
    // |p_ast|^2 - 2k(p_ast . u_star) + k^2 = R^2
    double b = -2.0 * p_ast.dot(u_star);
    double c = p_ast.squaredNorm() - R * R;
    double delta = b * b - 4.0 * c;
    
    if (delta < 0) return {jd, 0, 0}; 
    
    // We want the intersection closer to the asteroid (first hit on earth)
    double k = (-b - std::sqrt(delta)) / 2.0;
    if (k < 0) k = (-b + std::sqrt(delta)) / 2.0; // Try other side if needed
    
    Eigen::Vector3d p_geo = p_ast - k * u_star;
    
    double dist = p_geo.norm();
    double lat = std::asin(p_geo.z() / dist) * constants::RAD_TO_DEG;
    
    // Greenwich Mean Sidereal Time
    double T = (jd - 2451545.0) / 36525.0;
    double gmst_deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T * T;
    gmst_deg = std::fmod(gmst_deg, 360.0);
    if (gmst_deg < 0) gmst_deg += 360.0;
    
    double ra = std::atan2(p_geo.y(), p_geo.x()) * constants::RAD_TO_DEG;
    double lon = ra - gmst_deg;
    while (lon <= -180.0) lon += 360.0;
    while (lon > 180.0) lon -= 360.0;
    
    return {jd, lat, lon};
}

int main() {
    std::string bsp_path = "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp";
    auto ephem = std::make_shared<ephemeris::DE441Provider>(bsp_path);
    ephemeris::PlanetaryEphemeris::setProvider(ephem);

    // Star info: Gaia EDR3 3618364055431251072
    double ra_rad = (13.0 + 38.0/60.0 + 56.576/3600.0) * 15.0 * constants::DEG_TO_RAD;
    double dec_rad = -(8.0 + 48.0/60.0 + 51.12/3600.0) * constants::DEG_TO_RAD;
    Eigen::Vector3d u_star(std::cos(dec_rad)*std::cos(ra_rad), std::cos(dec_rad)*std::sin(ra_rad), std::sin(dec_rad));

    // Case 1: AstDyn Path (Current)
    EpochTDB t_tca_ast = EpochTDB::from_jd(2461112.616472); 
    
    io::HorizonsClient horizons;
    auto start_state = *horizons.query_vectors("50936", t_tca_ast - TimeDuration::from_seconds(300), "@sun");
    
    propagation::PropagatorSettings settings;
    settings.include_planets = true;
    settings.integrate_in_ecliptic = true;
    propagation::Propagator prop(std::make_shared<propagation::RKF78Integrator>(0.001), std::make_shared<ephemeris::PlanetaryEphemeris>(), settings);

    std::ofstream f_ast("track_astdyn.csv");
    f_ast << "jd,lat,lon\n";
    for (double s = -600; s <= 600; s += 2.0) {
        EpochTDB t = t_tca_ast + TimeDuration::from_seconds(s);
        auto state = prop.propagate_cartesian(start_state, t);
        auto sun_ssb = ephem->getState(ephemeris::CelestialBody::SUN, t);
        auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t);
        
        Eigen::Vector3d p_ast_geo = (state.position.to_eigen_si() + sun_ssb.position.to_eigen_si()) - earth_ssb.position.to_eigen_si();
        auto pt = project_to_earth(p_ast_geo, u_star, t.jd());
        if (std::abs(pt.lat) > 0.001)
            f_ast << std::fixed << std::setprecision(8) << pt.jd << "," << pt.lat << "," << pt.lon << "\n";
    }

    // Case 2: XML Path (Simulated to hit geocenter at 02:46:34)
    EpochTDB t_tca_xml = EpochTDB::from_jd(2461112.615671);
    auto state_tca = prop.propagate_cartesian(start_state, t_tca_xml); 
    auto sun_ssb_xml = ephem->getState(ephemeris::CelestialBody::SUN, t_tca_xml);
    auto earth_ssb_xml = ephem->getState(ephemeris::CelestialBody::EARTH, t_tca_xml);
    
    Eigen::Vector3d p_ast_geo_xml = (state_tca.position.to_eigen_si() + sun_ssb_xml.position.to_eigen_si()) - earth_ssb_xml.position.to_eigen_si();
    double dist = p_ast_geo_xml.norm();
    // Shadow center on geocenter -> Asteroid must be at k*u_star
    Eigen::Vector3d p_ast_geo_xml_target = u_star * dist; 
    
    Eigen::Vector3d offset = p_ast_geo_xml_target - p_ast_geo_xml;

    std::ofstream f_xml("track_xml.csv");
    f_xml << "jd,lat,lon\n";
    for (double s = -600; s <= 600; s += 2.0) {
        EpochTDB t = t_tca_xml + TimeDuration::from_seconds(s);
        auto state = prop.propagate_cartesian(start_state, t);
        auto sun_ssb = ephem->getState(ephemeris::CelestialBody::SUN, t);
        auto earth_ssb = ephem->getState(ephemeris::CelestialBody::EARTH, t);
        
        Eigen::Vector3d p_ast_geo = (state.position.to_eigen_si() + sun_ssb.position.to_eigen_si()) - earth_ssb.position.to_eigen_si();
        p_ast_geo += offset; 
        
        auto pt = project_to_earth(p_ast_geo, u_star, t.jd());
        if (std::abs(pt.lat) > 0.001)
            f_xml << std::fixed << std::setprecision(8) << pt.jd << "," << pt.lat << "," << pt.lon << "\n";
    }

    return 0;
}
