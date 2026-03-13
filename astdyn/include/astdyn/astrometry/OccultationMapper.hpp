/**
 * @file OccultationMapper.hpp
 * @brief Geography mapping and visualization for stellar occultations.
 */

#ifndef ASTDYN_ASTROMETRY_OCCULTATION_MAPPER_HPP
#define ASTDYN_ASTROMETRY_OCCULTATION_MAPPER_HPP

#include "astdyn/astrometry/OccultationLogic.hpp"
#include <vector>
#include <string>

namespace astdyn::astrometry {

/**
 * @brief Point on Earth's surface.
 */
struct GeoPoint {
    double lat_deg; ///< Latitude [-90, 90]
    double lon_deg; ///< Longitude [-180, 180]
};

struct TimeMarker {
    GeoPoint point;
    std::string label;
};

/**
 * @brief Representation of a predicted occultation path on a map.
 */
struct OccultationPath {
    std::vector<GeoPoint> center_line;    ///< Path of the shadow center
    std::vector<GeoPoint> sigma1_north;    ///< Northern 1-sigma boundary
    std::vector<GeoPoint> sigma1_south;    ///< Southern 1-sigma boundary
    std::vector<TimeMarker> markers;       ///< Time markers (TCA, +/- sigma)
    double shadow_width_km;                ///< Width of the asteroid shadow
};

/**
 * @brief Generates geographical data and visual files for occultations.
 */
class OccultationMapper {
public:
    /**
     * @brief Computes the geographical path of the shadow.
     * 
     * @param params Calculated occultation parameters.
     * @param tca_utc Time of Closest Approach.
     * @param duration_sec Length of path to generate (default 300s).
     * @return Full path data.
     */
    static OccultationPath compute_path(
        const OccultationParameters& params,
        const RightAscension& star_ra,
        const Declination& star_dec,
        const time::EpochUTC& tca_utc,
        double duration_sec = 300.0);

    /**
     * @brief Generates an SVG map representation.
     * @param path Calculated path data.
     * @param filename Output path.
     */
    static void export_svg(const OccultationPath& path, const std::string& filename);

    /**
     * @brief Generates a high-quality comparison SVG with two paths.
     */
    static void export_comparison_svg(
        const OccultationPath& path1, const std::string& label1,
        const OccultationPath& path2, const std::string& label2,
        const std::string& filename);

    /**
     * @brief Generates a global SVG map with world background and multiple paths.
     */
    static void export_global_svg(
        const std::vector<OccultationPath>& paths,
        const std::vector<std::string>& labels,
        const std::vector<std::string>& colors,
        const std::string& filename);

    /**
     * @brief Generates a KML file for Google Earth.
     * @param path Calculated path data.
     * @param filename Output path.
     */
    static void export_kml(const OccultationPath& path, const std::string& filename);

private:
    // Internal coordinate transformation helpers
    static GeoPoint project_to_earth(double xi, double eta, double zeta, 
                                     const RightAscension& ra, const Declination& dec, 
                                     const time::EpochUTC& t);
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_OCCULTATION_MAPPER_HPP
