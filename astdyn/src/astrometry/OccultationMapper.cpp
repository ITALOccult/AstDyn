/**
 * @file OccultationMapper.cpp
 */

#include "astdyn/astrometry/OccultationMapper.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/time/TimeScale.hpp"
#include <fstream>
#include <cmath>
#include <iomanip>

namespace astdyn::astrometry {

GeoPoint OccultationMapper::project_to_earth(double xi, double eta, double zeta, 
                                     const RightAscension& ra, const Declination& dec, 
                                     const time::EpochUTC& t) 
{
    using namespace constants;

    // 1. Rigorous Besselian to ECI transformation
    // k points towards the star
    double as = ra.to_rad();
    double ds = dec.to_rad();
    
    Eigen::Vector3d k_vec(std::cos(as) * std::cos(ds), std::sin(as) * std::cos(ds), std::sin(ds));
    Eigen::Vector3d i_vec(-std::sin(as), std::cos(as), 0.0);
    // User-suggested North vector (consistent with OccultationLogic fix)
    Eigen::Vector3d j_vec(-std::sin(ds) * std::cos(as), -std::sin(ds) * std::sin(as), std::cos(ds));

    // Calculate zeta if not provided (intersection with Earth surface)
    // zeta^2 = R_earth^2 - xi^2 - eta^2
    double r_earth = constants::R_EARTH; 
    double diff = r_earth * r_earth - (xi * xi + eta * eta);
    double z = (diff > 0) ? std::sqrt(diff) : 0.0;
    
    // Position in GCRF/ECI
    Eigen::Vector3d p_eci = xi * i_vec + eta * j_vec + z * k_vec;

    // 2. Simple GAST approximation for rotation to ITRF
    // T = centuries from J2000
    double jd = time::mjd_to_jd(t.mjd());
    double T = (jd - 2451545.0) / 36525.0;
    
    // GMST in seconds (Simplified IAU model)
    double gmst_sec = 24110.54841 + 8640184.812866 * T + 0.093104 * T*T;
    double gast_rad = std::fmod(gmst_sec * (constants::PI / 43200.0), constants::TWO_PI);
    
    // Rotate ECI -> ITRF (Rotation about Z)
    double cos_g = std::cos(gast_rad);
    double sin_g = std::sin(gast_rad);
    double x_itrf = p_eci.x() * cos_g + p_eci.y() * sin_g;
    double y_itrf = -p_eci.x() * sin_g + p_eci.y() * cos_g;
    double z_itrf = p_eci.z();

    // 3. ITRF to Geodetic (Spherical approximation)
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

    auto get_pos_at = [&](double t) {
        double s = params.shadow_velocity_kms * t;
        double xi  = params.impact_parameter_km * std::sin(theta) + s * std::cos(theta);
        double eta = params.impact_parameter_km * std::cos(theta) - s * std::sin(theta);
        return project_to_earth(xi, eta, 0.0, star_ra, star_dec, tca_utc);
    };

    for (int i = -steps/2; i <= steps/2; ++i) {
        double t = i * dt;
        double s = params.shadow_velocity_kms * t;

        // Path on fundamental plane
        double xi  = params.impact_parameter_km * std::sin(theta) + s * std::cos(theta);
        double eta = params.impact_parameter_km * std::cos(theta) - s * std::sin(theta);
        
        path.center_line.push_back(project_to_earth(xi, eta, 0.0, star_ra, star_dec, tca_utc));
        
        double w = path.shadow_width_km / 2.0;
        path.sigma1_north.push_back(project_to_earth(xi - w*std::sin(theta), eta + w*std::cos(theta), 0.0, star_ra, star_dec, tca_utc));
        path.sigma1_south.push_back(project_to_earth(xi + w*std::sin(theta), eta - w*std::sin(theta), 0.0, star_ra, star_dec, tca_utc));
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

} // namespace astdyn::astrometry
