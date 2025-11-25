/**
 * @file CloseApproach.hpp
 * @brief Close approach detection and analysis
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 * 
 * Detects and characterizes close approaches between small bodies and planets.
 * Implements b-plane analysis for impact probability assessment.
 * 
 * Key concepts:
 * - Close approach: When an object passes within a specified distance of a planet
 * - b-plane: Target plane perpendicular to approach velocity at planetary distance
 * - MOID: Minimum Orbit Intersection Distance between two orbits
 * 
 * Mathematical foundation:
 * - b-plane coordinates (ξ, ζ) characterize approach geometry
 * - Miss distance: b = √(ξ² + ζ²)
 * - Impact parameter mapping for uncertainty propagation
 */

#ifndef ASTDYN_CLOSE_APPROACH_HPP
#define ASTDYN_CLOSE_APPROACH_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/propagation/Propagator.hpp"
#include "astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <vector>
#include <memory>
#include <optional>

namespace astdyn::close_approach {

/**
 * @brief Type of celestial body for close approach
 */
enum class BodyType {
    MERCURY = 1,
    VENUS = 2,
    EARTH = 3,
    MARS = 4,
    JUPITER = 5,
    SATURN = 6,
    URANUS = 7,
    NEPTUNE = 8,
    MOON = 10
};

/**
 * @brief B-plane coordinates
 * 
 * The b-plane is perpendicular to the incoming asymptote at the
 * target body's distance. Coordinates (ξ, ζ) define the miss distance.
 */
struct BPlaneCoordinates {
    double xi;                          ///< ξ coordinate [AU]
    double zeta;                        ///< ζ coordinate [AU]
    double b_magnitude;                 ///< Miss distance |b| = √(ξ²+ζ²) [AU]
    double theta;                       ///< Angle arctan2(ζ,ξ) [rad]
    
    /**
     * @brief Compute miss distance
     */
    double miss_distance() const {
        return std::sqrt(xi * xi + zeta * zeta);
    }
};

/**
 * @brief Single close approach event
 */
struct CloseApproach {
    double mjd_tdb;                     ///< Time of closest approach [MJD TDB]
    BodyType body;                      ///< Target body
    
    // Geometric quantities at closest approach
    double distance;                    ///< Distance from body center [AU]
    double relative_velocity;           ///< |v_rel| [AU/day]
    astdyn::Vector3d position_object;   ///< Object heliocentric position [AU]
    astdyn::Vector3d velocity_object;   ///< Object heliocentric velocity [AU/day]
    astdyn::Vector3d position_body;     ///< Body heliocentric position [AU]
    astdyn::Vector3d velocity_body;     ///< Body heliocentric velocity [AU/day]
    
    // Relative coordinates
    astdyn::Vector3d rel_position;      ///< Object position relative to body [AU]
    astdyn::Vector3d rel_velocity;      ///< Object velocity relative to body [AU/day]
    
    // B-plane analysis
    std::optional<BPlaneCoordinates> b_plane; ///< Target plane coordinates
    
    /**
     * @brief Check if approach is within specified threshold
     */
    bool is_close(double threshold_au) const {
        return distance < threshold_au;
    }
    
    /**
     * @brief Get distance in planetary radii
     */
    double distance_in_radii(double planet_radius_au) const {
        return distance / planet_radius_au;
    }
};

/**
 * @brief Settings for close approach detection
 */
struct CloseApproachSettings {
    double detection_distance = 0.05;   ///< Detection threshold [AU] (default ~7.5 million km)
    double min_distance = 1e-6;         ///< Minimum distance for numerical stability [AU]
    bool compute_b_plane = true;        ///< Compute b-plane coordinates
    bool refine_time = true;            ///< Refine time of closest approach
    double time_tolerance = 1e-6;       ///< Time refinement tolerance [days]
    int max_refinement_iter = 10;       ///< Maximum iterations for time refinement
    
    // Bodies to check (empty = check all)
    std::vector<BodyType> bodies_to_check;
    
    /**
     * @brief Check if should monitor body
     */
    bool should_check_body(BodyType body) const {
        if (bodies_to_check.empty()) return true;
        return std::find(bodies_to_check.begin(), bodies_to_check.end(), body) 
               != bodies_to_check.end();
    }
};

/**
 * @brief Close approach detector
 * 
 * Monitors propagation to detect close approaches to specified bodies.
 * Uses distance monitoring and root finding for accurate timing.
 */
class CloseApproachDetector {
public:
    /**
     * @brief Construct detector
     * 
     * @param propagator Orbital propagator
     * @param settings Detection settings
     */
    CloseApproachDetector(
        std::shared_ptr<propagation::Propagator> propagator,
        const CloseApproachSettings& settings = CloseApproachSettings());
    
    /**
     * @brief Find close approaches in time interval
     * 
     * Propagates orbit and detects all close approaches to monitored bodies.
     * 
     * @param initial_state Initial orbital state
     * @param t_start Start time [MJD TDB]
     * @param t_end End time [MJD TDB]
     * @return Vector of detected close approaches (sorted by time)
     */
    std::vector<CloseApproach> find_approaches(
        const propagation::CartesianElements& initial_state,
        double t_start,
        double t_end);
    
    /**
     * @brief Find close approaches (Keplerian initial state)
     */
    std::vector<CloseApproach> find_approaches(
        const propagation::KeplerianElements& initial_orbit,
        double t_start,
        double t_end);
    
    /**
     * @brief Compute b-plane coordinates for a close approach
     * 
     * @param ca Close approach (must have position/velocity filled)
     * @return B-plane coordinates
     */
    BPlaneCoordinates compute_b_plane(const CloseApproach& ca) const;
    
    /**
     * @brief Get current settings
     */
    const CloseApproachSettings& settings() const { return settings_; }
    
    /**
     * @brief Update settings
     */
    void set_settings(const CloseApproachSettings& settings) { settings_ = settings; }

private:
    std::shared_ptr<propagation::Propagator> propagator_;
    CloseApproachSettings settings_;
    
    /**
     * @brief Monitor distance during propagation step
     * 
     * @param t Time [MJD TDB]
     * @param state Current state
     * @param body Body to check
     * @return Distance to body [AU]
     */
    double compute_distance(
        double t,
        const propagation::CartesianElements& state,
        BodyType body) const;
    
    /**
     * @brief Refine time of closest approach
     * 
     * Uses root finding (bisection/Brent) to accurately determine
     * the time when distance is minimized.
     * 
     * @param state_t1 State at time t1
     * @param state_t2 State at time t2
     * @param t1 Earlier time [MJD TDB]
     * @param t2 Later time [MJD TDB]
     * @param body Target body
     * @return Time of closest approach [MJD TDB]
     */
    double refine_approach_time(
        const propagation::CartesianElements& state_t1,
        const propagation::CartesianElements& state_t2,
        double t1,
        double t2,
        BodyType body) const;
    
    /**
     * @brief Build CloseApproach structure
     * 
     * @param t Time [MJD TDB]
     * @param state Object state
     * @param body Target body
     * @return Complete CloseApproach structure
     */
    CloseApproach build_approach(
        double t,
        const propagation::CartesianElements& state,
        BodyType body) const;
    
    /**
     * @brief Get body from celestial body enum
     */
    ephemeris::CelestialBody to_celestial_body(BodyType body) const;
    
    /**
     * @brief Get planet radius in AU
     */
    double get_planet_radius(BodyType body) const;
};

/**
 * @brief Minimum Orbit Intersection Distance (MOID) calculator
 * 
 * Computes the minimum possible distance between two orbits.
 * Useful for preliminary impact risk assessment.
 */
class MOIDCalculator {
public:
    /**
     * @brief Compute MOID between object and planet orbits
     * 
     * Uses optimization to find the minimum distance between
     * the two orbital paths (not time-dependent).
     * 
     * @param object_orbit Small body orbit
     * @param planet Planet identifier
     * @return MOID [AU]
     */
    static double compute_moid(
        const propagation::KeplerianElements& object_orbit,
        BodyType planet);
    
    /**
     * @brief Compute MOID between two arbitrary orbits
     * 
     * @param orbit1 First orbit
     * @param orbit2 Second orbit
     * @return MOID [AU]
     */
    static double compute_moid(
        const propagation::KeplerianElements& orbit1,
        const propagation::KeplerianElements& orbit2);

private:
    /**
     * @brief Objective function for MOID optimization
     * 
     * Computes distance² between two points on orbits given
     * their true anomalies.
     * 
     * @param f1 True anomaly on orbit 1 [rad]
     * @param f2 True anomaly on orbit 2 [rad]
     * @param orbit1 First orbit
     * @param orbit2 Second orbit
     * @return Distance² [AU²]
     */
    static double distance_squared(
        double f1, double f2,
        const propagation::KeplerianElements& orbit1,
        const propagation::KeplerianElements& orbit2);
};

} // namespace astdyn::close_approach

#endif // ASTDYN_CLOSE_APPROACH_HPP
