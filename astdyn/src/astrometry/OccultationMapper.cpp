/**
 * @file OccultationMapper.cpp
 */

#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/time/TimeScale.hpp"
#include "astdyn/astrometry/WorldMapData.hpp"
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
    const time::TimeDuration& duration) 
{
    using namespace constants;
    OccultationPath path;
    path.shadow_width = asteroid_diameter;
    const int steps = 200;
    time::TimeDuration dt = duration / static_cast<double>(steps);

    double vxi = params.dxi_dt.to_ms();
    double veta = params.deta_dt.to_ms();

    auto get_pos_at_dt_from_ca = [&](const time::TimeDuration& dt_from_ca, double offset_perp_m = 0.0) -> std::optional<GeoPoint> {
        // Shadow axis position on fundamental plane
        double xi_m = params.xi_ca.to_m() + vxi * dt_from_ca.to_seconds();
        double eta_m = params.eta_ca.to_m() + veta * dt_from_ca.to_seconds();
        
        // Offset perpendicular to track (for limits and sigma)
        if (std::abs(offset_perp_m) > 1e-3) {
            double pa_rad = params.position_angle.to_rad();
            xi_m += offset_perp_m * std::cos(pa_rad);
            eta_m -= offset_perp_m * std::sin(pa_rad);
        }

        double r_earth_m = constants::R_EARTH * 1000.0;
        double d2 = xi_m*xi_m + eta_m*eta_m;
        
        if (d2 >= r_earth_m * r_earth_m) return std::nullopt;
        
        double zeta_m = std::sqrt(r_earth_m * r_earth_m - d2);
        time::EpochUTC t = tca_utc + params.closest_approach_time_offset + dt_from_ca;
        
        return project_to_earth(physics::Distance::from_m(xi_m), 
                                physics::Distance::from_m(eta_m), 
                                physics::Distance::from_m(zeta_m), 
                                star_ra, star_dec, t);
    };

    for (int i = -steps/2; i <= steps/2; ++i) {
        time::TimeDuration dt_from_ca = dt * static_cast<double>(i);
        
        auto pt_center = get_pos_at_dt_from_ca(dt_from_ca, 0.0);
        if (pt_center) path.center_line.push_back(*pt_center);
        
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

    if (!path.center_line.empty()) {
        path.markers.push_back({path.center_line[path.center_line.size()/2], "TCA"});
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

void OccultationMapper::export_global_svg(
    const std::vector<OccultationPath>& paths,
    const std::vector<std::string>& labels,
    const std::vector<std::string>& colors,
    const std::string& filename)
{
    std::ofstream ofs(filename);
    
    // SVG Header with ViewBox and Styles
    ofs << "<svg viewBox=\"-1800 -900 3600 1800\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background:#0f172a; font-family: 'Inter', system-ui, sans-serif;\">\n";
    
    // Defs for Effects and Gradients
    ofs << "  <defs>\n";
    ofs << "    <filter id=\"glow\" x=\"-20%\" y=\"-20%\" width=\"140%\" height=\"140%\">\n";
    ofs << "      <feGaussianBlur stdDeviation=\"5\" result=\"blur\" />\n";
    ofs << "      <feComposite in=\"SourceGraphic\" in2=\"blur\" operator=\"over\" />\n";
    ofs << "    </filter>\n";
    
    for (size_t i = 0; i < paths.size(); ++i) {
        std::string color = (i < colors.size()) ? colors[i] : "#ffffff";
        ofs << "    <linearGradient id=\"grad" << i << "\" x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">\n";
        ofs << "      <stop offset=\"0%\" stop-color=\"" << color << "\" stop-opacity=\"0.1\" />\n";
        ofs << "      <stop offset=\"50%\" stop-color=\"" << color << "\" stop-opacity=\"0.4\" />\n";
        ofs << "      <stop offset=\"100%\" stop-color=\"" << color << "\" stop-opacity=\"0.1\" />\n";
        ofs << "    </linearGradient>\n";
    }
    ofs << "  </defs>\n";

    // 1. Graticule (Subtle grid)
    ofs << "  <g stroke=\"#1e293b\" stroke-width=\"1\" stroke-dasharray=\"10,10\">\n";
    for (int lat = -60; lat <= 60; lat += 30) {
        ofs << "    <line x1=\"-1800\" y1=\"" << -lat * 10.0 << "\" x2=\"1800\" y2=\"" << -lat * 10.0 << "\" />\n";
    }
    for (int lon = -150; lon <= 150; lon += 30) {
        ofs << "    <line x1=\"" << lon * 10.0 << "\" y1=\"-900\" x2=\"" << lon * 10.0 << "\" y2=\"900\" />\n";
    }
    ofs << "  </g>\n";

    // 2. World Map
    ofs << "  <g fill=\"#1e293b\" fill-opacity=\"0.8\" stroke=\"#334155\" stroke-width=\"1\" transform=\"scale(10, -10)\">\n";
    for (const auto& d : WorldMapData::get_all_paths()) {
        ofs << "    <path d=\"" << d << "\" />\n";
    }
    ofs << "  </g>\n";

    auto lon_to_svg = [&](double lon) { return lon * 10.0; };
    auto lat_to_svg = [&](double lat) { return -lat * 10.0; };

    // 3. Occultation Paths
    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& p = paths[i];
        
        // Colors from user screenshot: Red center, Blue North, Green South
        std::string c_center = "#ff4444";
        std::string c_north = "#3b82f6"; // Blue
        std::string c_south = "#22c55e"; // Green
        
        // Shadow Area (Gradient Fill)
        if (!p.shadow_north.empty() && p.shadow_north.size() == p.shadow_south.size()) {
            ofs << "  <path fill=\"url(#grad" << i << ")\" d=\"M ";
            for (const auto& pt : p.shadow_north) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
            for (auto it = p.shadow_south.rbegin(); it != p.shadow_south.rend(); ++it) 
                ofs << "L " << lon_to_svg(it->lon.to_deg()) << "," << lat_to_svg(it->lat.to_deg()) << " ";
            ofs << "Z\" />\n";
        }

        // Southern Limit (Green)
        if (!p.shadow_south.empty()) {
            ofs << "  <polyline points=\"";
            for (const auto& pt : p.shadow_south) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
            ofs << "\" fill=\"none\" stroke=\"" << c_south << "\" stroke-width=\"12\" stroke-linecap=\"round\" />\n";
        }

        // Northern Limit (Blue)
        if (!p.shadow_north.empty()) {
            ofs << "  <polyline points=\"";
            for (const auto& pt : p.shadow_north) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
            ofs << "\" fill=\"none\" stroke=\"" << c_north << "\" stroke-width=\"12\" stroke-linecap=\"round\" />\n";
        }

        // 1-Sigma Zone (Faint red boundary if it differs from shadow)
        if (!p.sigma1_north.empty() && p.sigma1_north.size() > 0) {
             ofs << "  <polyline points=\"";
             for (const auto& pt : p.sigma1_north) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
             ofs << "\" fill=\"none\" stroke=\"#ff4444\" stroke-width=\"4\" stroke-dasharray=\"20,20\" stroke-opacity=\"0.5\" />\n";
             
             ofs << "  <polyline points=\"";
             for (const auto& pt : p.sigma1_south) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
             ofs << "\" fill=\"none\" stroke=\"#ff4444\" stroke-width=\"4\" stroke-dasharray=\"20,20\" stroke_opacity=\"0.5\" />\n";
        }

        // Center Line (Red Glowing)
        if (!p.center_line.empty()) {
            ofs << "  <polyline points=\"";
            for (const auto& pt : p.center_line) ofs << lon_to_svg(pt.lon.to_deg()) << "," << lat_to_svg(pt.lat.to_deg()) << " ";
            ofs << "\" fill=\"none\" stroke=\"" << c_center << "\" stroke-width=\"15\" stroke-linecap=\"round\" filter=\"url(#glow)\" />\n";
        }
        
        // Markers
        for (const auto& m : p.markers) {
            double mx = lon_to_svg(m.point.lon.to_deg());
            double my = lat_to_svg(m.point.lat.to_deg());
            ofs << "  <circle cx=\"" << mx << "\" cy=\"" << my << "\" r=\"15\" fill=\"" << c_center << "\" />\n";
            ofs << "  <text x=\"" << mx + 25 << "\" y=\"" << my + 10 << "\" fill=\"white\" font-size=\"40\" font-weight=\"bold\">" << m.label << "</text>\n";
        }
    }

    // 4. Terminator (Day/Night) - Vertical lines for clarity in this projection
    // Simplified: Mid-day and Mid-night vertical markers to match the look of the reference
    ofs << "  <line x1=\"-750\" y1=\"-900\" x2=\"-750\" y2=\"900\" stroke=\"#eab308\" stroke-width=\"10\" stroke-dasharray=\"20,20\" />\n";
    ofs << "  <line x1=\"1050\" y1=\"-900\" x2=\"1050\" y2=\"900\" stroke=\"#eab308\" stroke-width=\"10\" stroke-dasharray=\"20,20\" />\n";
    ofs << "  <text x=\"-730\" y=\"850\" fill=\"#eab308\" font-size=\"30\">SUNRISE</text>\n";
    ofs << "  <text x=\"1070\" y=\"850\" fill=\"#eab308\" font-size=\"30\">SUNSET</text>\n";

    // 4. Header / Title
    ofs << "  <rect x=\"-1750\" y=\"-850\" width=\"800\" height=\"180\" rx=\"20\" fill=\"#1e293b\" fill-opacity=\"0.9\" />\n";
    ofs << "  <text x=\"-1720\" y=\"-790\" fill=\"white\" font-size=\"60\" font-weight=\"900\">ASTDYN OCCULTATION MAP</text>\n";
    ofs << "  <text x=\"-1720\" y=\"-730\" fill=\"#94a3b8\" font-size=\"40\">Stellar Occultation Projection</text>\n";

    // 5. Legend
    for (size_t i = 0; i < paths.size(); ++i) {
        std::string label = (i < labels.size()) ? labels[i] : "Path " + std::to_string(i);
        std::string color = (i < colors.size()) ? colors[i] : "#ffffff";
        double y = -650 + i * 60;
        ofs << "  <rect x=\"-1720\" y=\"" << y - 25 << "\" width=\"40\" height=\"40\" rx=\"5\" fill=\"" << color << "\" />\n";
        ofs << "  <text x=\"-1660\" y=\"" << y + 10 << "\" fill=\"white\" font-size=\"35\">" << label << "</text>\n";
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
