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
    Angle lat; ///< Latitude
    Angle lon; ///< Longitude
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
    std::vector<GeoPoint> shadow_north;   ///< Northern shadow boundary
    std::vector<GeoPoint> shadow_south;   ///< Southern shadow boundary
    std::vector<GeoPoint> sigma1_north;    ///< Northern 1-sigma boundary (shadow + uncertainty)
    std::vector<GeoPoint> sigma1_south;    ///< Southern 1-sigma boundary (shadow + uncertainty)
    std::vector<TimeMarker> markers;       ///< Time markers (TCA, +/- sigma)
    physics::Distance shadow_width;        ///< Width of the asteroid shadow
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
     * @param star_ra Star Right Ascension
     * @param star_dec Star Declination
     * @param asteroid_diameter Nominal diameter for corridor calculation
     * @param tca_utc Time of Closest Approach.
     * @param duration Length of path to generate (default 300s).
     * @return Full path data.
     */
    static OccultationPath compute_path(
        const OccultationParameters& params,
        const RightAscension& star_ra,
        const Declination& star_dec,
        const physics::Distance& asteroid_diameter,
        const time::EpochUTC& tca_utc,
        std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem,
        const time::TimeDuration& duration = time::TimeDuration::from_seconds(3600.0));

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
        const std::string& filename,
        std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem,
        Angle center_lat = Angle::from_deg(0.0),
        Angle center_lon = Angle::from_deg(0.0),
        double zoom = 1.0);

    /**
     * @brief Generates a KML file for Google Earth.
     * @param path Calculated path data.
     * @param filename Output path.
     */
    static void export_kml(const OccultationPath& path, const std::string& filename);

    /**
     * @brief Generates a KML file with multiple paths (system or comparison).
     */
    static void export_kml(
        const std::vector<OccultationPath>& paths,
        const std::vector<std::string>& labels,
        const std::string& filename);

private:
    // Internal coordinate transformation helpers
    static GeoPoint project_to_earth(const physics::Distance& xi, 
                                     const physics::Distance& eta, 
                                     const physics::Distance& zeta, 
                                     const RightAscension& ra, 
                                     const Declination& dec, 
                                     const time::EpochUTC& t);

    /**
     * @brief Computes the points of the day/night terminator.
     */
    static std::vector<GeoPoint> compute_terminator(const time::EpochTDB& t, 
                                                   std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephem);

    // CFIYH Helpers
    static std::optional<GeoPoint> calculate_geopoint_at_epoch(
        const OccultationParameters& params,
        const RightAscension& star_ra,
        const Declination& star_dec,
        const time::TimeDuration& dt_from_ca,
        const physics::Distance& offset_perp);

    static void add_time_marker(
        OccultationPath& path,
        const GeoPoint& pt,
        const time::TimeDuration& dt_from_ca,
        const time::EpochUTC& tca_utc);
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_OCCULTATION_MAPPER_HPP
