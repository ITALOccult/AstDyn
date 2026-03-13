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

GeoPoint OccultationMapper::project_to_earth(double xi, double eta, double /*zeta_unused*/, 
                                     const RightAscension& ra, const Declination& dec, 
                                     const time::EpochUTC& t) 
{
    using namespace constants;

    // 1. Besselian Basis (toward star)
    double as = ra.to_rad();
    double ds = dec.to_rad();
    Eigen::Vector3d k_vec(std::cos(as) * std::cos(ds), std::sin(as) * std::cos(ds), std::sin(ds));
    Eigen::Vector3d i_vec(-std::sin(as), std::cos(as), 0.0);
    Eigen::Vector3d j_vec(-std::sin(ds) * std::cos(as), -std::sin(ds) * std::sin(as), std::cos(ds));

    // 2. Zeta: Intersection with Earth surface (simplified sphere)
    double r_earth = constants::R_EARTH; 
    double diff = r_earth * r_earth - (xi * xi + eta * eta);
    double z = (diff > 0) ? std::sqrt(diff) : 0.0;
    
    // Position in GCRF/ECI [km]
    Eigen::Vector3d p_eci = xi * i_vec + eta * j_vec + z * k_vec;

    // 3. Earth Rotation Angle (ERA) - Modern IAU Standard
    // Du = days from J2000.0 (JD 2451545.0)
    double jd = time::mjd_to_jd(t.mjd());
    double Du = jd - 2451545.0;
    
    // ERA (in fractions of a rotation)
    double era_frac = 0.7790572732640 + 1.00273781191135448 * Du;
    double era_rad = std::fmod(era_frac * constants::TWO_PI, constants::TWO_PI);
    if (era_rad < 0) era_rad += constants::TWO_PI;
    
    // ECI -> ITRF (Rotation about Z by +ERA)
    double cos_era = std::cos(era_rad);
    double sin_era = std::sin(era_rad);
    double x_itrf =  p_eci.x() * cos_era + p_eci.y() * sin_era;
    double y_itrf = -p_eci.x() * sin_era + p_eci.y() * cos_era;
    double z_itrf =  p_eci.z();

    // 4. ITRF to Geodetic
    double lat = std::asin(z_itrf / r_earth) * RAD_TO_DEG;
    double lon = std::atan2(y_itrf, x_itrf) * RAD_TO_DEG;

    return {lat, lon};
}

OccultationPath OccultationMapper::compute_path(
    const OccultationParameters& params,
    const RightAscension& star_ra,
    const Declination& star_dec,
    const time::EpochUTC& tca_utc,
    double duration_sec) 
{
    using namespace constants;
    OccultationPath path;
    path.shadow_width_km = 100.0; // Typical for Nireus
    const int steps = 200;
    double dt = duration_sec / steps;

    // Movement angle on the fundamental plane
    double theta = params.position_angle_deg * DEG_TO_RAD;
    double vxi = params.shadow_velocity_kms * std::sin(theta);
    double veta = params.shadow_velocity_kms * std::cos(theta);

    auto get_pos_at = [&](double t_sec) {
        double xi = params.xi_ca_km + vxi * t_sec;
        double eta = params.eta_ca_km + veta * t_sec;
        return project_to_earth(xi, eta, 0.0, star_ra, star_dec, tca_utc + time::TimeDuration::from_seconds(t_sec));
    };

    for (int i = -steps/2; i <= steps/2; ++i) {
        double t = i * dt;
        double xi = params.xi_ca_km + vxi * t;
        double eta = params.eta_ca_km + veta * t;
        
        path.center_line.push_back(get_pos_at(t));
        
        double w = path.shadow_width_km / 2.0;
        double dxi_perp = w * std::cos(theta);
        double deta_perp = -w * std::sin(theta);

        path.sigma1_north.push_back(project_to_earth(xi - dxi_perp, eta - deta_perp, 0.0, star_ra, star_dec, tca_utc + time::TimeDuration::from_seconds(t)));
        path.sigma1_south.push_back(project_to_earth(xi + dxi_perp, eta + deta_perp, 0.0, star_ra, star_dec, tca_utc + time::TimeDuration::from_seconds(t)));
    }

    // Add Time Markers
    path.markers.push_back({get_pos_at(0.0), "TCA"});
    if (params.time_uncertainty_sec > 0.0) {
        path.markers.push_back({get_pos_at(-params.time_uncertainty_sec), "-1σ Time"});
        path.markers.push_back({get_pos_at(params.time_uncertainty_sec), "+1σ Time"});
    }
    
    return path;
}

void OccultationMapper::export_kml(const OccultationPath& path, const std::string& filename) {
    std::ofstream ofs(filename);
    ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    ofs << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    ofs << "  <Document>\n";
    ofs << "    <name>AstDyn Occultation Track</name>\n";
    
    // Style for center line
    ofs << "    <Style id=\"centerLine\">\n";
    ofs << "      <LineStyle><color>ff0000ff</color><width>3</width></LineStyle>\n";
    ofs << "    </Style>\n";
    
    // Style for shadow corridor (semi-transparent red)
    ofs << "    <Style id=\"corridorStyle\">\n";
    ofs << "      <LineStyle><color>800000ff</color><width>1</width></LineStyle>\n";
    ofs << "      <PolyStyle><color>400000ff</color></PolyStyle>\n";
    ofs << "    </Style>\n";

    // Style for markers
    ofs << "    <Style id=\"markerStyle\">\n";
    ofs << "      <IconStyle><scale>0.8</scale></IconStyle>\n";
    ofs << "      <LabelStyle><scale>1.0</scale></LabelStyle>\n";
    ofs << "    </Style>\n";

    // Draw shadow corridor
    ofs << "    <Placemark><name>Shadow Corridor</name><styleUrl>#corridorStyle</styleUrl><Polygon><outerBoundaryIs><LinearRing><coordinates>\n";
    for (const auto& p : path.sigma1_north) ofs << p.lon_deg << "," << p.lat_deg << ",0 ";
    for (auto it = path.sigma1_south.rbegin(); it != path.sigma1_south.rend(); ++it) 
        ofs << it->lon_deg << "," << it->lat_deg << ",0 ";
    ofs << path.sigma1_north[0].lon_deg << "," << path.sigma1_north[0].lat_deg << ",0 ";
    ofs << "</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n";

    // Draw center line
    ofs << "    <Placemark><name>Shadow Center</name><styleUrl>#centerLine</styleUrl><LineString><coordinates>\n";
    for (const auto& p : path.center_line) ofs << p.lon_deg << "," << p.lat_deg << ",0 ";
    ofs << "</coordinates></LineString></Placemark>\n";

    // Draw time markers
    for (const auto& m : path.markers) {
        ofs << "    <Placemark>\n";
        ofs << "      <name>" << m.label << "</name>\n";
        ofs << "      <styleUrl>#markerStyle</styleUrl>\n";
        ofs << "      <Point><coordinates>" << m.point.lon_deg << "," << m.point.lat_deg << ",0</coordinates></Point>\n";
        ofs << "    </Placemark>\n";
    }

    ofs << "  </Document>\n</kml>\n";
}

void OccultationMapper::export_svg(const OccultationPath& path, const std::string& filename) {
    if (path.center_line.empty()) return;

    // 1. Calculate bounding box for auto-scaling
    double min_lon = 180, max_lon = -180, min_lat = 90, max_lat = -90;
    auto update_bb = [&](const GeoPoint& p) {
        if (p.lon_deg < min_lon) min_lon = p.lon_deg;
        if (p.lon_deg > max_lon) max_lon = p.lon_deg;
        if (p.lat_deg < min_lat) min_lat = p.lat_deg;
        if (p.lat_deg > max_lat) max_lat = p.lat_deg;
    };
    for (const auto& p : path.center_line) update_bb(p);
    for (const auto& p : path.sigma1_north) update_bb(p);
    for (const auto& p : path.sigma1_south) update_bb(p);

    // Add 15% padding
    double dlon = max_lon - min_lon;
    double dlat = max_lat - min_lat;
    min_lon -= dlon * 0.15; max_lon += dlon * 0.15;
    min_lat -= dlat * 0.15; max_lat += dlat * 0.15;
    dlon = max_lon - min_lon;
    dlat = max_lat - min_lat;

    std::ofstream ofs(filename);
    ofs << "<svg viewBox=\"0 0 800 600\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background:#0b0e14\">\n";
    
    auto to_x = [&](double lon) { return (lon - min_lon) / dlon * 800.0; };
    auto to_y = [&](double lat) { return (max_lat - lat) / dlat * 600.0; };

    // Draw Grid
    ofs << "  <g stroke=\"#222\" stroke-width=\"1\">\n";
    for (int lo = (int)min_lon; lo <= (int)max_lon; lo++) ofs << "    <line x1=\"" << to_x(lo) << "\" y1=\"0\" x2=\"" << to_x(lo) << "\" y2=\"600\" />\n";
    for (int la = (int)min_lat; la <= (int)max_lat; la++) ofs << "    <line x1=\"0\" y1=\"" << to_y(la) << "\" x2=\"800\" y2=\"" << to_y(la) << "\" />\n";
    ofs << "  </g>\n";

    // Draw shadow corridor
    ofs << "  <path fill=\"#ff4444\" fill-opacity=\"0.15\" d=\"M ";
    for (const auto& p : path.sigma1_north) ofs << to_x(p.lon_deg) << "," << to_y(p.lat_deg) << " ";
    for (auto it = path.sigma1_south.rbegin(); it != path.sigma1_south.rend(); ++it) 
        ofs << "L " << to_x(it->lon_deg) << "," << to_y(it->lat_deg) << " ";
    ofs << "Z\" />\n";

    // Draw center line
    ofs << "  <polyline points=\"";
    for (const auto& p : path.center_line) ofs << to_x(p.lon_deg) << "," << to_y(p.lat_deg) << " ";
    ofs << "\" fill=\"none\" stroke=\"#ff4444\" stroke-width=\"2\" />\n";

    // Draw time markers
    for (const auto& m : path.markers) {
        double x = to_x(m.point.lon_deg);
        double y = to_y(m.point.lat_deg);
        ofs << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"4\" fill=\"white\" />\n";
        ofs << "  <text x=\"" << x + 8 << "\" y=\"" << y + 4 << "\" fill=\"white\" font-family=\"sans-serif\" font-size=\"10\">" << m.label << "</text>\n";
    }

    // Legend
    ofs << "  <rect x=\"20\" y=\"20\" width=\"250\" height=\"80\" rx=\"5\" fill=\"#000\" fill-opacity=\"0.7\" />\n";
    ofs << "  <text x=\"35\" y=\"45\" fill=\"white\" font-family=\"sans-serif\" font-size=\"14\" font-weight=\"bold\">AstDyn Occultation Map</text>\n";
    ofs << "  <text x=\"35\" y=\"65\" fill=\"#ccc\" font-family=\"sans-serif\" font-size=\"12\">Target: (50936) Nireus</text>\n";
    ofs << "  <text x=\"35\" y=\"83\" fill=\"#f55\" font-family=\"sans-serif\" font-size=\"11\">Uncertainty: ±1σ Corridor Shown</text>\n";

    ofs << "</svg>\n";
}
void OccultationMapper::export_comparison_svg(
    const OccultationPath& path1, const std::string& label1,
    const OccultationPath& path2, const std::string& label2,
    const std::string& filename) 
{
    // 1. Calculate bounding box for both paths
    double min_lon = 180, max_lon = -180, min_lat = 90, max_lat = -90;
    auto update_bb = [&](const GeoPoint& p) {
        if (p.lon_deg < -180 || p.lon_deg > 180) return;
        if (p.lon_deg < min_lon) min_lon = p.lon_deg;
        if (p.lon_deg > max_lon) max_lon = p.lon_deg;
        if (p.lat_deg < min_lat) min_lat = p.lat_deg;
        if (p.lat_deg > max_lat) max_lat = p.lat_deg;
    };
    
    for (const auto& p : path1.center_line) update_bb(p);
    for (const auto& p : path2.center_line) update_bb(p);
    for (const auto& p : path1.sigma1_north) update_bb(p);
    for (const auto& p : path1.sigma1_south) update_bb(p);
    for (const auto& p : path2.sigma1_north) update_bb(p);
    for (const auto& p : path2.sigma1_south) update_bb(p);

    // Expand for readability
    double dlon = max_lon - min_lon;
    double dlat = max_lat - min_lat;
    min_lon -= dlon * 0.2; max_lon += dlon * 0.2;
    min_lat -= dlat * 0.2; max_lat += dlat * 0.5; // More room at top for title
    dlon = max_lon - min_lon;
    dlat = max_lat - min_lat;

    std::ofstream ofs(filename);
    ofs << "<svg viewBox=\"0 0 1000 800\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    ofs << "  <defs>\n";
    ofs << "    <linearGradient id=\"bgGrad\" x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"100%\">\n";
    ofs << "      <stop offset=\"0%\" style=\"stop-color:#0f172a;stop-opacity:1\" />\n";
    ofs << "      <stop offset=\"100%\" style=\"stop-color:#020617;stop-opacity:1\" />\n";
    ofs << "    </linearGradient>\n";
    ofs << "    <filter id=\"glow\">\n";
    ofs << "      <feGaussianBlur stdDeviation=\"2.5\" result=\"coloredBlur\"/>\n";
    ofs << "      <feMerge><feMergeNode in=\"coloredBlur\"/><feMergeNode in=\"SourceGraphic\"/></feMerge>\n";
    ofs << "    </filter>\n";
    ofs << "  </defs>\n";

    // Background
    ofs << "  <rect width=\"1000\" height=\"800\" fill=\"url(#bgGrad)\" />\n";

    auto to_x = [&](double lon) { return (lon - min_lon) / dlon * 1000.0; };
    auto to_y = [&](double lat) { return (max_lat - lat) / dlat * 800.0; };

    // Graticule (Grid)
    ofs << "  <g stroke=\"#1e293b\" stroke-width=\"1\" stroke-dasharray=\"4\">\n";
    for (double lo = std::floor(min_lon); lo <= max_lon; lo += std::max(1.0, std::round(dlon/10.0))) {
        ofs << "    <line x1=\"" << to_x(lo) << "\" y1=\"0\" x2=\"" << to_x(lo) << "\" y2=\"800\" />\n";
        ofs << "    <text x=\"" << to_x(lo) + 2 << "\" y=\"790\" fill=\"#64748b\" font-size=\"10\" font-family=\"sans-serif\">" << (int)lo << "°</text>\n";
    }
    for (double la = std::floor(min_lat); la <= max_lat; la += std::max(1.0, std::round(dlat/10.0))) {
        ofs << "    <line x1=\"0\" y1=\"" << to_y(la) << "\" x2=\"1000\" y2=\"" << to_y(la) << "\" />\n";
        ofs << "    <text x=\"5\" y=\"" << to_y(la) - 2 << "\" fill=\"#64748b\" font-size=\"10\" font-family=\"sans-serif\">" << (int)la << "°</text>\n";
    }
    ofs << "  </g>\n";

    // --- Draw Path 2 (Occult4 Emulated) ---
    ofs << "  <path fill=\"#ef4444\" fill-opacity=\"0.1\" d=\"M ";
    for (const auto& p : path2.sigma1_north) ofs << to_x(p.lon_deg) << "," << to_y(p.lat_deg) << " ";
    for (auto it = path2.sigma1_south.rbegin(); it != path2.sigma1_south.rend(); ++it) 
        ofs << "L " << to_x(it->lon_deg) << "," << to_y(it->lat_deg) << " ";
    ofs << "Z\" />\n";
    ofs << "  <polyline points=\"";
    for (const auto& p : path2.center_line) ofs << to_x(p.lon_deg) << "," << to_y(p.lat_deg) << " ";
    ofs << "\" fill=\"none\" stroke=\"#ef4444\" stroke-width=\"3\" filter=\"url(#glow)\" />\n";

    // --- Draw Path 1 (AstDyn Nominal) ---
    ofs << "  <path fill=\"#06b6d4\" fill-opacity=\"0.15\" d=\"M ";
    for (const auto& p : path1.sigma1_north) ofs << to_x(p.lon_deg) << "," << to_y(p.lat_deg) << " ";
    for (auto it = path1.sigma1_south.rbegin(); it != path1.sigma1_south.rend(); ++it) 
        ofs << "L " << to_x(it->lon_deg) << "," << to_y(it->lat_deg) << " ";
    ofs << "Z\" />\n";
    ofs << "  <polyline points=\"";
    for (const auto& p : path1.center_line) ofs << to_x(p.lon_deg) << "," << to_y(p.lat_deg) << " ";
    ofs << "\" fill=\"none\" stroke=\"#06b6d4\" stroke-width=\"3\" filter=\"url(#glow)\" />\n";

    // Markers for both
    auto draw_markers = [&](const OccultationPath& p, const std::string& color) {
        for (const auto& m : p.markers) {
            double x = to_x(m.point.lon_deg);
            double y = to_y(m.point.lat_deg);
            ofs << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"5\" fill=\"" << color << "\" stroke=\"white\" stroke-width=\"1\" />\n";
            ofs << "  <text x=\"" << x + 8 << "\" y=\"" << y + 4 << "\" fill=\"white\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"bold\">" << m.label << "</text>\n";
        }
    };
    draw_markers(path1, "#06b6d4");
    draw_markers(path2, "#ef4444");

    // Header & Legend (Glassmorphism)
    ofs << "  <rect x=\"30\" y=\"30\" width=\"400\" height=\"180\" rx=\"15\" fill=\"#1e293b\" fill-opacity=\"0.6\" stroke=\"#334155\" stroke-width=\"1\" />\n";
    ofs << "  <text x=\"50\" y=\"70\" fill=\"white\" font-family=\"sans-serif\" font-size=\"24\" font-weight=\"bold\">Occultation Comparison</text>\n";
    ofs << "  <text x=\"50\" y=\"100\" fill=\"#94a3b8\" font-family=\"sans-serif\" font-size=\"16\">(50936) Nireus - TYC 5546-1581-1</text>\n";
    
    // Legend entry 1
    ofs << "  <rect x=\"50\" y=\"125\" width=\"15\" height=\"15\" fill=\"#06b6d4\" rx=\"3\" />\n";
    ofs << "  <text x=\"75\" y=\"138\" fill=\"white\" font-family=\"sans-serif\" font-size=\"14\">" << label1 << "</text>\n";
    
    // Legend entry 2
    ofs << "  <rect x=\"50\" y=\"155\" width=\"15\" height=\"15\" fill=\"#ef4444\" rx=\"3\" />\n";
    ofs << "  <text x=\"75\" y=\"168\" fill=\"white\" font-family=\"sans-serif\" font-size=\"14\">" << label2 << "</text>\n";

    // Metadata Footer
    ofs << "  <text x=\"50\" y=\"780\" fill=\"#475569\" font-family=\"sans-serif\" font-size=\"12\">Generated by AstDyn Engine v2.5 - ITALOccult System</text>\n";

    ofs << "</svg>\n";
}

void OccultationMapper::export_global_svg(
    const std::vector<OccultationPath>& paths,
    const std::vector<std::string>& labels,
    const std::vector<std::string>& colors,
    const std::string& filename)
{
    std::ofstream ofs(filename);
    // Standard Equirectangular Projection Space: x [-180, 180], y [-90, 90]
    // We use a larger scale (x10) for better browser compatibility and precision
    ofs << "<svg viewBox=\"-1800 -900 3600 1800\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background:#0f172a\">\n";
    
    // 1. Defs (Glow)
    ofs << "  <defs>\n";
    ofs << "    <filter id=\"glow_global\">\n";
    ofs << "      <feGaussianBlur stdDeviation=\"5\" result=\"coloredBlur\"/>\n";
    ofs << "      <feMerge><feMergeNode in=\"coloredBlur\"/><feMergeNode in=\"SourceGraphic\"/></feMerge>\n";
    ofs << "    </filter>\n";
    ofs << "  </defs>\n";

    // 2. Graticule (Grid)
    ofs << "  <g stroke=\"#1e293b\" stroke-width=\"1\">\n";
    for (int lon = -180; lon <= 180; lon += 30) {
        ofs << "    <line x1=\"" << lon*10 << "\" y1=\"-900\" x2=\"" << lon*10 << "\" y2=\"900\" />\n";
    }
    for (int lat = -90; lat <= 90; lat += 30) {
        ofs << "    <line x1=\"-1800\" y1=\"" << lat*10 << "\" x2=\"1800\" y2=\"" << lat*10 << "\" />\n";
    }
    ofs << "  </g>\n";

    // 3. Landmasses (from WorldMapData)
    ofs << "  <g fill=\"#1e293b\" stroke=\"#334155\" stroke-width=\"2\">\n";
    for (const auto& d : WorldMapData::get_all_paths()) {
        // Need to scale the d string points by 10 since our viewBox is x10
        // But the data is already (lon, -lat) which maps to SVG (x, y)
        // We'll use a transform instead of parsing the strings.
        ofs << "    <path d=\"" << d << "\" transform=\"scale(10)\" />\n";
    }
    ofs << "  </g>\n";

    auto lon_to_svg = [&](double lon) { return lon * 10.0; };
    auto lat_to_svg = [&](double lat) { return -lat * 10.0; };

    // 4. Occultation Paths
    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& p = paths[i];
        std::string color = (i < colors.size()) ? colors[i] : "#ffffff";
        
        // Corridor
        if (!p.sigma1_north.empty()) {
            ofs << "  <path fill=\"" << color << "\" fill-opacity=\"0.15\" d=\"M ";
            for (const auto& pt : p.sigma1_north) ofs << lon_to_svg(pt.lon_deg) << "," << lat_to_svg(pt.lat_deg) << " ";
            for (auto it = p.sigma1_south.rbegin(); it != p.sigma1_south.rend(); ++it) 
                ofs << "L " << lon_to_svg(it->lon_deg) << "," << lat_to_svg(it->lat_deg) << " ";
            ofs << "Z\" />\n";
        }

        // Center line
        if (!p.center_line.empty()) {
            ofs << "  <polyline points=\"";
            for (const auto& pt : p.center_line) {
                ofs << lon_to_svg(pt.lon_deg) << "," << lat_to_svg(pt.lat_deg) << " ";
            }
            ofs << "\" fill=\"none\" stroke=\"" << color << "\" stroke-width=\"15\" filter=\"url(#glow_global)\" />\n";
        }

        // Markers
        for (const auto& marker : p.markers) {
            if (marker.label == "TCA") {
                ofs << "  <circle cx=\"" << lon_to_svg(marker.point.lon_deg) << "\" cy=\"" << lat_to_svg(marker.point.lat_deg) 
                    << "\" r=\"15\" fill=\"" << color << "\" stroke=\"white\" stroke-width=\"4\" />\n";
            }
        }
    }

    // 5. Legend (Overlay)
    // Scale legend properly for the large coordinate space
    ofs << "  <g transform=\"translate(-1750, -850)\">\n";
    ofs << "    <rect width=\"1100\" height=\"" << (120 + paths.size() * 80) << "\" rx=\"40\" fill=\"#1e293b\" fill-opacity=\"0.9\" stroke=\"#334155\" stroke-width=\"5\" />\n";
    ofs << "    <text x=\"50\" y=\"80\" fill=\"white\" font-family=\"sans-serif\" font-size=\"60\" font-weight=\"bold\">Stellar Occultation Tracker</text>\n";
    for (size_t i = 0; i < paths.size(); ++i) {
        std::string color = (i < colors.size()) ? colors[i] : "#ffffff";
        ofs << "    <rect x=\"50\" y=\"" << (120 + i * 80) << "\" width=\"40\" height=\"40\" fill=\"" << color << "\" rx=\"10\" />\n";
        ofs << "    <text x=\"110\" y=\"" << (155 + i * 80) << "\" fill=\"white\" font-family=\"sans-serif\" font-size=\"45\">" << labels[i] << "</text>\n";
    }
    ofs << "  </g>\n";

    ofs << "  <text x=\"-1750\" y=\"850\" fill=\"#475569\" font-family=\"sans-serif\" font-size=\"40\">AstDyn Global Engine | Equirectangular Projection</text>\n";
    ofs << "</svg>\n";
}

} // namespace astdyn::astrometry
