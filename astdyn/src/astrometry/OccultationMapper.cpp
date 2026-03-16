/**
 * @file OccultationMapper.cpp
 */

#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/astrometry/WorldMapData.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <fstream>
#include <cmath>
#include <iomanip>

namespace astdyn::astrometry {

GeoPoint OccultationMapper::project_to_earth(
    const physics::Distance& xi, 
    const physics::Distance& eta, 
    const physics::Distance& zeta, 
    const RightAscension& ra, 
    const Declination& dec, 
    const time::EpochUTC& t) 
{
    using namespace constants;

    // 1. Besselian Basis (toward star)
    double as = ra.to_rad();
    double ds = dec.to_rad();
    Eigen::Vector3d k_vec(std::cos(as) * std::cos(ds), std::sin(as) * std::cos(ds), std::sin(ds));
    Eigen::Vector3d i_vec(-std::sin(as), std::cos(as), 0.0);
    Eigen::Vector3d j_vec(-std::sin(ds) * std::cos(as), -std::sin(ds) * std::sin(as), std::cos(ds));

    // 2. Cartesian position in GCRF (m)
    Eigen::Vector3d p_eci = xi.to_m() * i_vec + eta.to_m() * j_vec + zeta.to_m() * k_vec;

    // 3. Earth Rotation Angle (ERA)
    // DU is defined in UT1 days since J2000.0
    double dut1 = time::get_dut1(t.mjd());
    double mjd_ut1 = t.mjd() + dut1 / constants::SECONDS_PER_DAY;
    double jd_ut1 = time::mjd_to_jd(mjd_ut1);
    double Du = jd_ut1 - 2451545.0;
    
    double era_frac = 0.7790572732640 + 1.00273781191135448 * Du;
    double era_rad = std::fmod(era_frac * constants::TWO_PI, constants::TWO_PI);
    if (era_rad < 0) era_rad += constants::TWO_PI;
    
    // 4. Geodetic Projection (Spherical)
    double ra_p = std::atan2(p_eci.y(), p_eci.x());
    double dec_p = std::asin(p_eci.z() / p_eci.norm());
    
    double lon_rad = ra_p - era_rad;
    while (lon_rad > PI) lon_rad -= TWO_PI;
    while (lon_rad <= -PI) lon_rad += TWO_PI;

    return {Angle::from_rad(dec_p), Angle::from_rad(lon_rad)};
}

OccultationPath OccultationMapper::compute_path(
    const OccultationParameters& params,
    const RightAscension& star_ra,
    const Declination& star_dec,
    const physics::Distance& asteroid_diameter,
    const time::EpochUTC& tca_utc,
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem,
    const time::TimeDuration& duration) 
{
    using namespace constants;
    OccultationPath path;
    path.shadow_width = asteroid_diameter;
    const int steps = 600;
    time::TimeDuration dt = duration / static_cast<double>(steps);

    double vxi = params.dxi_dt.to_ms();
    double veta = params.deta_dt.to_ms();

    auto get_pos_at_dt_from_ca = [&](const time::TimeDuration& dt_from_ca, double offset_perp_m = 0.0) -> std::optional<GeoPoint> {
        // Shadow axis position on fundamental plane
        double xi_m = params.xi_ca.to_m() + vxi * dt_from_ca.to_seconds();
        double eta_m = params.eta_ca.to_m() + veta * dt_from_ca.to_seconds();
        
        // Offset perpendicular to track (for limits and sigma)
        if (std::abs(offset_perp_m) > 1e-3) {
            double v_norm = std::sqrt(vxi*vxi + veta*veta);
            xi_m += offset_perp_m * (-veta / v_norm);
            eta_m += offset_perp_m * ( vxi  / v_norm);
        }

        double r_earth_m = constants::R_EARTH * 1000.0;
        double d2 = xi_m*xi_m + eta_m*eta_m;
        
        if (d2 >= r_earth_m * r_earth_m) return std::nullopt;
        
        double zeta_m = std::sqrt(r_earth_m * r_earth_m - d2);
        // NO redundant offset here: tca_utc is already the time of CA
        time::EpochUTC t = tca_utc + dt_from_ca;
        
        return project_to_earth(physics::Distance::from_m(xi_m), 
                                physics::Distance::from_m(eta_m), 
                                physics::Distance::from_m(zeta_m), 
                                star_ra, star_dec, t);
    };

    for (int i = -steps/2; i <= steps/2; ++i) {
        time::TimeDuration dt_from_ca = dt * static_cast<double>(i);
        
        auto pt_center = get_pos_at_dt_from_ca(dt_from_ca, 0.0);
        if (pt_center) {
            path.center_line.push_back(*pt_center);
            
            // Add time markers every 60 seconds (1 minute)
            double secs = std::round(dt_from_ca.to_seconds());
            if (std::abs(secs) < 1e-3 || static_cast<int>(std::abs(secs) + 0.5) % 60 == 0) {
                time::EpochUTC t = tca_utc + dt_from_ca;
                auto [y, mon, d, frac] = time::mjd_to_calendar(t.mjd());
                
                int hour = static_cast<size_t>(frac * 24.0) % 24;
                int minute = static_cast<size_t>(frac * 1440.0) % 60;
                
                std::stringstream ss;
                ss << std::setfill('0') << std::setw(2) << hour << ":" << std::setw(2) << minute;
                path.markers.push_back({*pt_center, ss.str()});
            }
        }
        
        double w2 = path.shadow_width.to_m() / 2.0;
        double sigma_m = params.cross_track_uncertainty.to_m();

        // Nominal Shadow Limits
        auto pt_n = get_pos_at_dt_from_ca(dt_from_ca, -w2);
        auto pt_s = get_pos_at_dt_from_ca(dt_from_ca, w2);
        if (pt_n) path.shadow_north.push_back(*pt_n);
        if (pt_s) path.shadow_south.push_back(*pt_s);

        // 1-Sigma Uncertainty Limits (Shadow + Uncertainty)
        auto pt_sn = get_pos_at_dt_from_ca(dt_from_ca, -(w2 + sigma_m));
        auto pt_ss = get_pos_at_dt_from_ca(dt_from_ca, (w2 + sigma_m));
        if (pt_sn) path.sigma1_north.push_back(*pt_sn);
        if (pt_ss) path.sigma1_south.push_back(*pt_ss);
    }

    return path;
}

void OccultationMapper::export_svg(const OccultationPath& path, const std::string& filename) {
    if (path.center_line.empty()) return;
    std::ofstream ofs(filename);
    ofs << "<svg viewBox=\"0 0 800 600\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background:#0b0e14\">\n";
    
    double center_lon = path.center_line[path.center_line.size()/2].lon.to_deg();
    double center_lat = path.center_line[path.center_line.size()/2].lat.to_deg();
    
    auto to_x = [&](double lon) { return (lon - center_lon) * 50.0 + 400.0; };
    auto to_y = [&](double lat) { return (center_lat - lat) * 50.0 + 300.0; };

    ofs << "  <polyline points=\"";
    for (const auto& pt : path.center_line) ofs << to_x(pt.lon.to_deg()) << "," << to_y(pt.lat.to_deg()) << " ";
    ofs << "\" fill=\"none\" stroke=\"#ff4444\" stroke-width=\"2\" />\n";
    ofs << "</svg>\n";
}

void OccultationMapper::export_comparison_svg(
    const OccultationPath& path1, const std::string& label1,
    const OccultationPath& path2, const std::string& label2,
    const std::string& filename) 
{
    if (path1.center_line.empty() || path2.center_line.empty()) return;
    std::ofstream ofs(filename);
    ofs << "<svg viewBox=\"0 0 800 800\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background:#0f172a\">\n";
    
    double center_lon = path2.center_line[path2.center_line.size()/2].lon.to_deg();
    double center_lat = path2.center_line[path2.center_line.size()/2].lat.to_deg();
    
    auto to_x = [&](double lon) { return (lon - center_lon) * 40.0 + 400.0; };
    auto to_y = [&](double lat) { return (center_lat - lat) * 40.0 + 400.0; };

    // Graticule
    ofs << "  <g stroke=\"#1e293b\" stroke-width=\"1\">\n";
    for (int i = -10; i <= 10; ++i) {
        ofs << "    <line x1=\"0\" y1=\"" << to_y(center_lat + i) << "\" x2=\"800\" y2=\"" << to_y(center_lat + i) << "\"/>\n";
        ofs << "    <line x1=\"" << to_x(center_lon + i) << "\" y1=\"0\" x2=\"" << to_x(center_lon + i) << "\" y2=\"800\"/>\n";
    }
    ofs << "  </g>\n";

    auto draw_path = [&](const OccultationPath& p, const std::string& color) {
        ofs << "  <path fill=\"" << color << "\" fill-opacity=\"0.15\" d=\"M ";
        for (const auto& pt : p.shadow_north) ofs << to_x(pt.lon.to_deg()) << "," << to_y(pt.lat.to_deg()) << " ";
        for (auto it = p.shadow_south.rbegin(); it != p.shadow_south.rend(); ++it) 
            ofs << "L " << to_x(it->lon.to_deg()) << "," << to_y(it->lat.to_deg()) << " ";
        ofs << "Z\" />\n";
        ofs << "  <polyline points=\"";
        for (const auto& pt : p.center_line) ofs << to_x(pt.lon.to_deg()) << "," << to_y(pt.lat.to_deg()) << " ";
        ofs << "\" fill=\"none\" stroke=\"" << color << "\" stroke-width=\"3\" />\n";
    };

    draw_path(path1, "#06b6d4");
    draw_path(path2, "#ef4444");
    ofs << "</svg>\n";
}

std::vector<GeoPoint> OccultationMapper::compute_terminator(const time::EpochTDB& t,
                                                   std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem) {
    using namespace constants;
    std::vector<GeoPoint> points;
    
    // Fallback if ephem is null to avoid crash (though it should be provided)
    auto actual_ephem = ephem ? ephem : std::make_shared<ephemeris::PlanetaryEphemeris>();

    // 1. Get Sun position in GCRF
    auto sun_ssb = actual_ephem->getState(ephemeris::CelestialBody::SUN, t);
    auto earth_ssb = actual_ephem->getState(ephemeris::CelestialBody::EARTH, t);
    Eigen::Vector3d s_vec = (sun_ssb.position.to_eigen_si() - earth_ssb.position.to_eigen_si()).normalized();
    
    // 2. Earth Rotation Angle (ERA)
    double dut1 = time::get_dut1(t.mjd());
    double mjd_ut1 = t.mjd() + dut1 / constants::SECONDS_PER_DAY;
    double Du = time::mjd_to_jd(mjd_ut1) - 2451545.0;
    double era_rad = std::fmod((0.7790572732640 + 1.00273781191135448 * Du) * constants::TWO_PI, constants::TWO_PI);
    if (era_rad < 0) era_rad += constants::TWO_PI;

    // 3. Solve for longitude at each latitude
    double sx = s_vec.x();
    double sy = s_vec.y();
    double sz = s_vec.z();
    double R = std::sqrt(sx*sx + sy*sy);
    double alpha = std::atan2(sy, sx);

    for (int ilat = -90; ilat <= 90; ++ilat) {
        double phi = ilat * constants::DEG_TO_RAD;
        double cos_val = -sz * std::tan(phi) / R;
        
        if (std::abs(cos_val) <= 1.0) {
            double d_theta = std::acos(cos_val);
            for (double sign : {1.0, -1.0}) {
                double Theta = alpha + sign * d_theta;
                double lon = Theta - era_rad;
                while (lon > PI) lon -= TWO_PI;
                while (lon <= -PI) lon += TWO_PI;
                points.push_back({Angle::from_rad(phi), Angle::from_rad(lon)});
            }
        }
    }
    // Sort points to make a coherent line (Western branch then Eastern)
    std::sort(points.begin(), points.end(), [](const GeoPoint& a, const GeoPoint& b) {
        if (std::abs(a.lon.to_deg() - b.lon.to_deg()) < 0.1) return a.lat.to_deg() < b.lat.to_deg();
        return a.lon.to_deg() < b.lon.to_deg();
    });

    return points;
}

void OccultationMapper::export_global_svg(
    const std::vector<OccultationPath>& paths,
    const std::vector<std::string>& labels,
    const std::vector<std::string>& colors,
    const std::string& filename,
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem,
    double center_lat,
    double center_lon,
    double zoom)
{
    std::ofstream ofs(filename);
    
    // Dynamic ViewBox
    double vw = 3600.0 / zoom;
    double vh = 1800.0 / zoom;
    double vx = center_lon * 10.0 - vw / 2.0;
    double vy = -center_lat * 10.0 - vh / 2.0;

    // SVG Header
    ofs << "<svg viewBox=\"" << vx << " " << vy << " " << vw << " " << vh << "\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background:#0b0e14; font-family: 'Outfit', 'Inter', sans-serif;\">\n";
    
    // Defs
    ofs << "  <defs>\n";
    ofs << "    <filter id=\"glow\" x=\"-50%\" y=\"-50%\" width=\"200%\" height=\"200%\">\n";
    ofs << "      <feGaussianBlur stdDeviation=\"" << 8.0 / zoom << "\" result=\"blur\" />\n";
    ofs << "      <feComposite in=\"SourceGraphic\" in2=\"blur\" operator=\"over\" />\n";
    ofs << "    </filter>\n";
    
    for (size_t i = 0; i < paths.size(); ++i) {
        std::string color = (i < colors.size()) ? colors[i] : "#ffffff";
        ofs << "    <linearGradient id=\"grad" << i << "\" x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">\n";
        ofs << "      <stop offset=\"0%\" stop-color=\"" << color << "\" stop-opacity=\"0.05\" />\n";
        ofs << "      <stop offset=\"50%\" stop-color=\"" << color << "\" stop-opacity=\"0.2\" />\n";
        ofs << "      <stop offset=\"100%\" stop-color=\"" << color << "\" stop-opacity=\"0.05\" />\n";
        ofs << "    </linearGradient>\n";
    }
    ofs << "  </defs>\n";

    auto lon_to_svg = [&](double lon) { return lon * 10.0; };
    auto lat_to_svg = [&](double lat) { return -lat * 10.0; };

    // 1. Sea
    ofs << "  <rect x=\"-1800\" y=\"-900\" width=\"3600\" height=\"1800\" fill=\"#0b0e14\" />\n";

    // 2. Graticule
    ofs << "  <g stroke=\"#1e293b\" stroke-width=\"" << 1.0 / zoom << "\" stroke-dasharray=\"" << 10.0/zoom << "," << 10.0/zoom << "\" opacity=\"0.5\">\n";
    for (int lat = -90; lat <= 90; lat += 15) ofs << "    <line x1=\"-1800\" y1=\"" << lat_to_svg(lat) << "\" x2=\"1800\" y2=\"" << lat_to_svg(lat) << "\" />\n";
    for (int lon = -180; lon <= 180; lon += 15) ofs << "    <line x1=\"" << lon_to_svg(lon) << "\" y1=\"-900\" x2=\"" << lon_to_svg(lon) << "\" y2=\"900\" />\n";
    ofs << "  </g>\n";

    // 3. Detailed World Map - Coastlines
    ofs << "  <g fill=\"#1e293b\" stroke=\"#334155\" stroke-width=\"" << 0.8 / zoom << "\" transform=\"scale(10, -10)\">\n";
    for (const auto& d : WorldMapData::get_coastlines()) {
        ofs << "    <path d=\"" << d << "\" />\n";
    }
    ofs << "  </g>\n";

    // 3.1 Political Boundaries
    ofs << "  <g fill=\"none\" stroke=\"#475569\" stroke-width=\"" << 0.5 / zoom << "\" stroke-dasharray=\"" << 1.0/zoom << "," << 1.0/zoom << "\" transform=\"scale(10, -10)\" opacity=\"0.7\">\n";
    for (const auto& d : WorldMapData::get_borders()) {
        ofs << "    <path d=\"" << d << "\" />\n";
    }
    ofs << "  </g>\n";

    // 3.2 Major Cities
    ofs << "  <g>\n";
    for (const auto& city : WorldMapData::get_major_cities()) {
        double cx = lon_to_svg(city.lon);
        double cy = lat_to_svg(city.lat);
        ofs << "    <circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << 5.0 / zoom << "\" fill=\"#94a3b8\" opacity=\"0.6\" />\n";
        ofs << "    <text x=\"" << cx + 10.0/zoom << "\" y=\"" << cy + 5.0/zoom << "\" fill=\"#64748b\" font-size=\"" << 20.0 / zoom << "\" opacity=\"0.5\">" << city.name << "</text>\n";
    }
    ofs << "  </g>\n";

    // 4. Day/Night Terminator
    auto term_pts = compute_terminator(time::EpochTDB::from_jd(2461121.5), ephem);
    if (!term_pts.empty()) {
        ofs << "  <polyline points=\"";
        for (const auto& pt : term_pts) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
        ofs << "\" fill=\"none\" stroke=\"#facc15\" stroke-width=\"" << 6.0 / zoom << "\" stroke-dasharray=\"" << 20.0/zoom << "," << 15.0/zoom << "\" opacity=\"0.6\" />\n";
    }

    // 5. Occultation Paths
    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& p = paths[i];
        std::string color = (i < colors.size()) ? colors[i] : "#ffffff";
        
        // --- 1-Sigma Uncertainty Region ---
        if (!p.sigma1_north.empty() && p.sigma1_north.size() == p.sigma1_south.size()) {
            ofs << "  <path fill=\"" << color << "\" fill-opacity=\"0.12\" d=\"M ";
            for (const auto& pt : p.sigma1_north) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
            for (auto it = p.sigma1_south.rbegin(); it != p.sigma1_south.rend(); ++it) 
                ofs << "L " << lon_to_svg(it->lon.to_deg()) << "," << lat_to_svg(it->lat.to_deg()) << " ";
            ofs << "Z\" />\n";
        }

        // Shadow Zone (Nominal)
        if (!p.shadow_north.empty() && p.shadow_north.size() == p.shadow_south.size()) {
            ofs << "  <path fill=\"url(#grad" << i << ")\" fill-opacity=\"0.7\" d=\"M ";
            for (const auto& pt : p.shadow_north) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
            for (auto it = p.shadow_south.rbegin(); it != p.shadow_south.rend(); ++it) 
                ofs << "L " << lon_to_svg(it->lon.to_deg()) << "," << lat_to_svg(it->lat.to_deg()) << " ";
            ofs << "Z\" />\n";
        }

        auto draw_line = [&](const std::vector<GeoPoint>& pts, const std::string& color, double width, double op = 1.0, bool dashed = false) {
            if (pts.empty()) return;
            ofs << "  <polyline points=\"";
            for (const auto& pt : pts) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
            ofs << "\" fill=\"none\" stroke=\"" << color << "\" stroke-width=\"" << width / zoom << "\" stroke-linecap=\"round\" ";
            if (dashed) ofs << "stroke-dasharray=\"" << 20.0/zoom << "," << 20.0/zoom << "\" ";
            ofs << "opacity=\"" << op << "\" />\n";
        };

        draw_line(p.sigma1_north, color, 3, 0.4, true);
        draw_line(p.sigma1_south, color, 3, 0.4, true);
        draw_line(p.shadow_north, color, 6, 0.3);
        draw_line(p.shadow_south, color, 6, 0.3);
        draw_line(p.center_line, color, 5, 1.0);
        
        for (const auto& m : p.markers) {
            double mx = lon_to_svg(m.point.lon.to_deg());
            double my = lat_to_svg(m.point.lat.to_deg());
            ofs << "  <circle cx=\"" << mx << "\" cy=\"" << my << "\" r=\"" << 8.0 / zoom << "\" fill=\"" << color << "\" />\n";
            ofs << "  <text x=\"" << mx + 15.0/zoom << "\" y=\"" << my + 5.0/zoom << "\" fill=\"#f8fafc\" font-size=\"" << 25.0 / zoom << "\" font-weight=\"bold\" stroke=\"#0b0e14\" stroke-width=\"" << 4.0/zoom << "\" paint-order=\"stroke\">" << m.label << "</text>\n";
        }
    }

    // Header and Legend (Positioned relative to ViewBox)
    double ui_scale = 1.0 / zoom;
    double margin = 50.0 * ui_scale;
    double lx = vx + margin;
    double ly = vy + margin;

    ofs << "  <rect x=\"" << lx << "\" y=\"" << ly << "\" width=\"" << 850 * ui_scale << "\" height=\"" << 200 * ui_scale << "\" rx=\"" << 25 * ui_scale << "\" fill=\"#1e293b\" fill-opacity=\"0.95\" stroke=\"#334155\" stroke-width=\"" << 2 * ui_scale << "\" />\n";
    ofs << "  <text x=\"" << lx + 40 * ui_scale << "\" y=\"" << ly + 65 * ui_scale << "\" fill=\"#f8fafc\" font-size=\"" << 65 * ui_scale << "\" font-weight=\"900\">ASTDYN PRECISION MAP</text>\n";
    ofs << "  <text x=\"" << lx + 40 * ui_scale << "\" y=\"" << ly + 130 * ui_scale << "\" fill=\"#94a3b8\" font-size=\"" << 42 * ui_scale << "\">Beta 0.5 - High-Fidelity Projection</text>\n";

    for (size_t i = 0; i < paths.size(); ++i) {
        std::string label = (i < labels.size()) ? labels[i] : "Path " + std::to_string(i);
        std::string color = (i < colors.size()) ? colors[i] : "#ffffff";
        double y_coord = ly + 220 * ui_scale + static_cast<double>(i) * 110.0 * ui_scale;
        
        ofs << "  <rect x=\"" << lx << "\" y=\"" << y_coord - 30 * ui_scale << "\" width=\"" << 45 * ui_scale << "\" height=\"" << 45 * ui_scale << "\" rx=\"" << 8 * ui_scale << "\" fill=\"" << color << "\" />\n";
        ofs << "  <text x=\"" << lx + 65 * ui_scale << "\" y=\"" << y_coord + 10 * ui_scale << "\" fill=\"#f8fafc\" font-size=\"" << 38 * ui_scale << "\">" << label << "</text>\n";
        
        ofs << "  <rect x=\"" << lx << "\" y=\"" << y_coord + 25 * ui_scale << "\" width=\"" << 45 * ui_scale << "\" height=\"" << 15 * ui_scale << "\" rx=\"" << 4 * ui_scale << "\" fill=\"" << color << "\" fill-opacity=\"0.3\" stroke=\"#f8fafc\" stroke-width=\"" << 1.0*ui_scale << "\" stroke-dasharray=\"" << 5*ui_scale << "," << 3*ui_scale << "\" />\n";
        ofs << "  <text x=\"" << lx + 65 * ui_scale << "\" y=\"" << y_coord + 40 * ui_scale << "\" fill=\"#94a3b8\" font-size=\"" << 30 * ui_scale << "\">1-sigma Uncertainty Corridor</text>\n";
    }

    ofs << "</svg>\n";
}

void OccultationMapper::export_kml(const OccultationPath& path, const std::string& filename) {
    std::ofstream ofs(filename);
    ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    ofs << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    ofs << "<Document>\n";
    ofs << "  <name>Occultation Path</name>\n";
    
    // Styles
    ofs << "  <Style id=\"centerLine\"><LineStyle><color>ff0000ff</color><width>4</width></LineStyle></Style>\n";
    ofs << "  <Style id=\"shadowLimit\"><LineStyle><color>ffff0000</color><width>2</width></LineStyle></Style>\n";
    ofs << "  <Style id=\"sigmaLimit\"><LineStyle><color>ff0000ff</color><width>1</width><dashArray>5,5</dashArray></LineStyle></Style>\n";
    ofs << "  <Style id=\"shadowArea\"><PolyStyle><color>40ffffff</color><fill>1</fill><outline>0</outline></PolyStyle></Style>\n";
    ofs << "  <Style id=\"sigmaArea\"><PolyStyle><color>200000ff</color><fill>1</fill><outline>0</outline></PolyStyle></Style>\n";

    // 1. Center Line
    ofs << "  <Placemark><name>Nominal Center</name><styleUrl>#centerLine</styleUrl>\n";
    ofs << "    <LineString><tessellate>1</tessellate><coordinates>\n";
    for (const auto& pt : path.center_line) ofs << std::fixed << std::setprecision(6) << pt.lon.to_deg() << "," << pt.lat.to_deg() << ",0 ";
    ofs << "    </coordinates></LineString></Placemark>\n";

    // 2. Shadow Limits and Area
    if (!path.shadow_north.empty()) {
        ofs << "  <Placemark><name>North Shadow Limit</name><styleUrl>#shadowLimit</styleUrl>\n";
        ofs << "    <LineString><tessellate>1</tessellate><coordinates>\n";
        for (const auto& pt : path.shadow_north) ofs << pt.lon.to_deg() << "," << pt.lat.to_deg() << ",0 ";
        ofs << "    </coordinates></LineString></Placemark>\n";
        
        ofs << "  <Placemark><name>South Shadow Limit</name><styleUrl>#shadowLimit</styleUrl>\n";
        ofs << "    <LineString><tessellate>1</tessellate><coordinates>\n";
        for (const auto& pt : path.shadow_south) ofs << pt.lon.to_deg() << "," << pt.lat.to_deg() << ",0 ";
        ofs << "    </coordinates></LineString></Placemark>\n";

        ofs << "  <Placemark><name>Shadow Path</name><styleUrl>#shadowArea</styleUrl>\n";
        ofs << "    <Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing><coordinates>\n";
        for (const auto& pt : path.shadow_north) ofs << pt.lon.to_deg() << "," << pt.lat.to_deg() << ",0 ";
        for (auto it = path.shadow_south.rbegin(); it != path.shadow_south.rend(); ++it) ofs << it->lon.to_deg() << "," << it->lat.to_deg() << ",0 ";
        ofs << path.shadow_north[0].lon.to_deg() << "," << path.shadow_north[0].lat.to_deg() << ",0 ";
        ofs << "    </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n";
    }

    // 3. Sigma Uncertainty Limits (Outer)
    if (!path.sigma1_north.empty()) {
        ofs << "  <Placemark><name>1-Sigma North (Uncertainty)</name><styleUrl>#sigmaLimit</styleUrl>\n";
        ofs << "    <LineString><tessellate>1</tessellate><coordinates>\n";
        for (const auto& pt : path.sigma1_north) ofs << pt.lon.to_deg() << "," << pt.lat.to_deg() << ",0 ";
        ofs << "    </coordinates></LineString></Placemark>\n";

        ofs << "  <Placemark><name>1-Sigma South (Uncertainty)</name><styleUrl>#sigmaLimit</styleUrl>\n";
        ofs << "    <LineString><tessellate>1</tessellate><coordinates>\n";
        for (const auto& pt : path.sigma1_south) ofs << pt.lon.to_deg() << "," << pt.lat.to_deg() << ",0 ";
        ofs << "    </coordinates></LineString></Placemark>\n";

        // Sigma Area (Outer)
        ofs << "  <Placemark><name>1-Sigma Zone</name><styleUrl>#sigmaArea</styleUrl>\n";
        ofs << "    <Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing><coordinates>\n";
        for (const auto& pt : path.sigma1_north) ofs << pt.lon.to_deg() << "," << pt.lat.to_deg() << ",0 ";
        for (auto it = path.sigma1_south.rbegin(); it != path.sigma1_south.rend(); ++it) ofs << it->lon.to_deg() << "," << it->lat.to_deg() << ",0 ";
        ofs << path.sigma1_north[0].lon.to_deg() << "," << path.sigma1_north[0].lat.to_deg() << ",0 ";
        ofs << "    </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n";
    }

    ofs << "</Document>\n</kml>\n";
}

void OccultationMapper::export_kml(
    const std::vector<OccultationPath>& paths,
    const std::vector<std::string>& labels,
    const std::string& filename) 
{
    std::ofstream ofs(filename);
    ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    ofs << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    ofs << "<Document>\n";
    ofs << "  <name>Multi-body Occultation System</name>\n";
    
    // Styles
    ofs << "  <Style id=\"centerLine\"><LineStyle><color>ff0000ff</color><width>4</width></LineStyle></Style>\n";
    ofs << "  <Style id=\"shadowLimit\"><LineStyle><color>ffff0000</color><width>2</width></LineStyle></Style>\n";
    ofs << "  <Style id=\"sigmaLimit\"><LineStyle><color>ff0000ff</color><width>1</width><dashArray>5,5</dashArray></LineStyle></Style>\n";
    ofs << "  <Style id=\"shadowArea\"><PolyStyle><color>40ffffff</color><fill>1</fill><outline>0</outline></PolyStyle></Style>\n";
    ofs << "  <Style id=\"sigmaArea\"><PolyStyle><color>200000ff</color><fill>1</fill><outline>0</outline></PolyStyle></Style>\n";

    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& path = paths[i];
        std::string label = (i < labels.size()) ? labels[i] : "Body " + std::to_string(i);
        
        ofs << "  <Folder>\n";
        ofs << "    <name>" << label << "</name>\n";
        
        // 1. Center Line
        ofs << "    <Placemark><name>Center " << label << "</name><styleUrl>#centerLine</styleUrl>\n";
        ofs << "      <LineString><tessellate>1</tessellate><coordinates>\n";
        for (const auto& pt : path.center_line) ofs << std::fixed << std::setprecision(6) << pt.lon.to_deg() << "," << pt.lat.to_deg() << ",0 ";
        ofs << "      </coordinates></LineString></Placemark>\n";

        // 2. Shadow Limits and Area
        if (!path.shadow_north.empty()) {
            ofs << "    <Placemark><name>Shadow " << label << "</name><styleUrl>#shadowArea</styleUrl>\n";
            ofs << "      <Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing><coordinates>\n";
            for (const auto& pt : path.shadow_north) ofs << pt.lon.to_deg() << "," << pt.lat.to_deg() << ",0 ";
            for (auto it = path.shadow_south.rbegin(); it != path.shadow_south.rend(); ++it) ofs << it->lon.to_deg() << "," << it->lat.to_deg() << ",0 ";
            if (!path.shadow_north.empty())
                ofs << path.shadow_north[0].lon.to_deg() << "," << path.shadow_north[0].lat.to_deg() << ",0 ";
            ofs << "      </coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n";
        }
        
        ofs << "  </Folder>\n";
    }

    ofs << "</Document>\n</kml>\n";
}

} // namespace astdyn::astrometry
