/**
 * @file CloseApproach.cpp
 * @brief Implementation of close approach detection
 * @author OrbFit C++ Conversion Team
 * @date 2025-11-24
 */

#include "orbfit/close_approach/CloseApproach.hpp"
#include "orbfit/core/Constants.hpp"
#include "orbfit/propagation/OrbitalElements.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace orbfit::close_approach {

using namespace orbfit::propagation;
using namespace orbfit::ephemeris;
using namespace orbfit::constants;

// ============================================================================
// Planet Physical Data
// ============================================================================

namespace {
    // Planet radii in AU
    constexpr double MERCURY_RADIUS_AU = 2439.7 / AU_TO_KM;
    constexpr double VENUS_RADIUS_AU = 6051.8 / AU_TO_KM;
    constexpr double EARTH_RADIUS_AU = 6371.0 / AU_TO_KM;
    constexpr double MARS_RADIUS_AU = 3389.5 / AU_TO_KM;
    constexpr double JUPITER_RADIUS_AU = 69911.0 / AU_TO_KM;
    constexpr double SATURN_RADIUS_AU = 58232.0 / AU_TO_KM;
    constexpr double URANUS_RADIUS_AU = 25362.0 / AU_TO_KM;
    constexpr double NEPTUNE_RADIUS_AU = 24622.0 / AU_TO_KM;
    constexpr double MOON_RADIUS_AU = 1737.4 / AU_TO_KM;
    
    /**
     * @brief Get mean orbital elements for planets (J2000.0 epoch)
     * 
     * Returns mean elements at epoch J2000.0 from JPL approximations.
     * These are suitable for MOID calculations and general orbit analysis.
     * 
     * For osculating elements at specific epoch, use:
     *   auto mean_elem = get_planet_mean_elements(planet);
     *   auto osc_elem = propagation::mean_to_osculating(mean_elem);
     * 
     * Note: For heliocentric planetary orbits, mean ≈ osculating to high precision
     * since J2_sun ≈ 2e-7 is negligible compared to J2_Earth ≈ 1e-3.
     * 
     * Source: JPL Planetary Satellite Mean Elements
     * https://ssd.jpl.nasa.gov/planets/approx_pos.html
     * 
     * @param planet Planet identifier
     * @return Mean Keplerian elements at J2000.0
     */
    KeplerianElements get_planet_mean_elements(BodyType planet) {
        KeplerianElements elements;
        elements.epoch_mjd_tdb = 0.0; // MJD J2000.0
        elements.gravitational_parameter = 0.295912208; // Sun GM [AU³/day²]
        
        switch (planet) {
            case BodyType::MERCURY:
                elements.semi_major_axis = 0.38709927;
                elements.eccentricity = 0.20563593;
                elements.inclination = 7.00497902 * DEG_TO_RAD;
                elements.longitude_ascending_node = 48.33076593 * DEG_TO_RAD;
                elements.argument_perihelion = 77.45779628 * DEG_TO_RAD;
                elements.mean_anomaly = 252.25032350 * DEG_TO_RAD;
                break;
            case BodyType::VENUS:
                elements.semi_major_axis = 0.72333566;
                elements.eccentricity = 0.00677672;
                elements.inclination = 3.39467605 * DEG_TO_RAD;
                elements.longitude_ascending_node = 76.67984255 * DEG_TO_RAD;
                elements.argument_perihelion = 131.60246718 * DEG_TO_RAD;
                elements.mean_anomaly = 181.97909950 * DEG_TO_RAD;
                break;
            case BodyType::EARTH:
                elements.semi_major_axis = 1.00000261;
                elements.eccentricity = 0.01671123;
                elements.inclination = -0.00001531 * DEG_TO_RAD;
                elements.longitude_ascending_node = 0.0;
                elements.argument_perihelion = 102.93768193 * DEG_TO_RAD;
                elements.mean_anomaly = 100.46457166 * DEG_TO_RAD;
                break;
            case BodyType::MARS:
                elements.semi_major_axis = 1.52371034;
                elements.eccentricity = 0.09339410;
                elements.inclination = 1.84969142 * DEG_TO_RAD;
                elements.longitude_ascending_node = 49.55953891 * DEG_TO_RAD;
                elements.argument_perihelion = -23.94362959 * DEG_TO_RAD;
                elements.mean_anomaly = -4.55343205 * DEG_TO_RAD;
                break;
            case BodyType::JUPITER:
                elements.semi_major_axis = 5.20288700;
                elements.eccentricity = 0.04838624;
                elements.inclination = 1.30439695 * DEG_TO_RAD;
                elements.longitude_ascending_node = 100.47390909 * DEG_TO_RAD;
                elements.argument_perihelion = -85.78939509 * DEG_TO_RAD;
                elements.mean_anomaly = 34.39644051 * DEG_TO_RAD;
                break;
            case BodyType::SATURN:
                elements.semi_major_axis = 9.53667594;
                elements.eccentricity = 0.05386179;
                elements.inclination = 2.48599187 * DEG_TO_RAD;
                elements.longitude_ascending_node = 113.66242448 * DEG_TO_RAD;
                elements.argument_perihelion = -21.06494908 * DEG_TO_RAD;
                elements.mean_anomaly = 49.95424423 * DEG_TO_RAD;
                break;
            case BodyType::URANUS:
                elements.semi_major_axis = 19.18916464;
                elements.eccentricity = 0.04725744;
                elements.inclination = 0.77263783 * DEG_TO_RAD;
                elements.longitude_ascending_node = 74.01692503 * DEG_TO_RAD;
                elements.argument_perihelion = 96.99852891 * DEG_TO_RAD;
                elements.mean_anomaly = 313.23810451 * DEG_TO_RAD;
                break;
            case BodyType::NEPTUNE:
                elements.semi_major_axis = 30.06992276;
                elements.eccentricity = 0.00859048;
                elements.inclination = 1.77004347 * DEG_TO_RAD;
                elements.longitude_ascending_node = 131.78422574 * DEG_TO_RAD;
                elements.argument_perihelion = -55.12002969 * DEG_TO_RAD;
                elements.mean_anomaly = -55.12002969 * DEG_TO_RAD;
                break;
            default:
                throw std::invalid_argument("Unsupported planet for MOID calculation");
        }
        
        return elements;
    }
}

// ============================================================================
// CloseApproachDetector Implementation
// ============================================================================

CloseApproachDetector::CloseApproachDetector(
    std::shared_ptr<Propagator> propagator,
    const CloseApproachSettings& settings)
    : propagator_(propagator), settings_(settings)
{
    if (!propagator_) {
        throw std::invalid_argument("CloseApproachDetector: null propagator");
    }
}

std::vector<CloseApproach> CloseApproachDetector::find_approaches(
    const CartesianElements& initial_state,
    double t_start,
    double t_end)
{
    std::vector<CloseApproach> approaches;
    
    if (t_end <= t_start) {
        throw std::invalid_argument("t_end must be > t_start");
    }
    
    // Determine which bodies to check
    std::vector<BodyType> bodies;
    if (settings_.bodies_to_check.empty()) {
        // Check all major planets
        bodies = {BodyType::MERCURY, BodyType::VENUS, BodyType::EARTH, 
                  BodyType::MARS, BodyType::JUPITER, BodyType::SATURN,
                  BodyType::URANUS, BodyType::NEPTUNE};
    } else {
        bodies = settings_.bodies_to_check;
    }
    
    // Propagate with fine time steps to monitor distance
    double dt = 1.0;  // 1 day steps for monitoring
    double t = t_start;
    CartesianElements prev_state = initial_state;
    
    // Track previous distances for detecting minima
    std::map<BodyType, double> prev_distances;
    for (auto body : bodies) {
        prev_distances[body] = compute_distance(t, prev_state, body);
    }
    
    while (t < t_end) {
        double t_next = std::min(t + dt, t_end);
        
        // Propagate to next step
        CartesianElements current_state = propagator_->propagate_cartesian(prev_state, t_next);
        
        // Check each body
        for (auto body : bodies) {
            if (!settings_.should_check_body(body)) continue;
            
            double dist_current = compute_distance(t_next, current_state, body);
            double dist_prev = prev_distances[body];
            
            // Check if we passed through a local minimum (distance decreased then increased)
            if (dist_current > dist_prev && dist_prev < settings_.detection_distance) {
                // Found potential close approach - refine time
                double t_ca = t_next;
                if (settings_.refine_time) {
                    t_ca = refine_approach_time(prev_state, current_state, t, t_next, body);
                }
                
                // Get state at refined time
                CartesianElements state_ca = propagator_->propagate_cartesian(prev_state, t_ca);
                
                // Build close approach structure
                CloseApproach ca = build_approach(t_ca, state_ca, body);
                
                // Add to list if within threshold
                if (ca.distance < settings_.detection_distance) {
                    approaches.push_back(ca);
                }
            }
            
            // Update previous distance
            prev_distances[body] = dist_current;
        }
        
        // Move to next step
        t = t_next;
        prev_state = current_state;
    }
    
    // Sort by time
    std::sort(approaches.begin(), approaches.end(),
              [](const CloseApproach& a, const CloseApproach& b) {
                  return a.mjd_tdb < b.mjd_tdb;
              });
    
    return approaches;
}

std::vector<CloseApproach> CloseApproachDetector::find_approaches(
    const KeplerianElements& initial_orbit,
    double t_start,
    double t_end)
{
    // Convert to Cartesian and call main method
    CartesianElements cart = propagation::keplerian_to_cartesian(initial_orbit);
    return find_approaches(cart, t_start, t_end);
}

double CloseApproachDetector::compute_distance(
    double t,
    const CartesianElements& state,
    BodyType body) const
{
    // Get planet position
    CelestialBody cb = to_celestial_body(body);
    Vector3d planet_pos = PlanetaryEphemeris::getPosition(cb, t);
    
    // Compute relative position
    Vector3d rel_pos = state.position - planet_pos;
    
    return rel_pos.norm();
}

double CloseApproachDetector::refine_approach_time(
    const CartesianElements& state_t1,
    const CartesianElements& state_t2,
    double t1,
    double t2,
    BodyType body) const
{
    // Use golden section search to find minimum distance
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;  // Golden ratio
    const double resphi = 2.0 - phi;
    
    double a = t1;
    double b = t2;
    double tol = settings_.time_tolerance;
    
    // Initial bracketing
    double x1 = a + resphi * (b - a);
    double x2 = b - resphi * (b - a);
    
    CartesianElements state_x1 = propagator_->propagate_cartesian(state_t1, x1);
    CartesianElements state_x2 = propagator_->propagate_cartesian(state_t1, x2);
    
    double f1 = compute_distance(x1, state_x1, body);
    double f2 = compute_distance(x2, state_x2, body);
    
    int iter = 0;
    while (std::abs(b - a) > tol && iter < settings_.max_refinement_iter) {
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + resphi * (b - a);
            state_x1 = propagator_->propagate_cartesian(state_t1, x1);
            f1 = compute_distance(x1, state_x1, body);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - resphi * (b - a);
            state_x2 = propagator_->propagate_cartesian(state_t1, x2);
            f2 = compute_distance(x2, state_x2, body);
        }
        iter++;
    }
    
    return (a + b) / 2.0;
}

CloseApproach CloseApproachDetector::build_approach(
    double t,
    const CartesianElements& state,
    BodyType body) const
{
    CloseApproach ca;
    ca.mjd_tdb = t;
    ca.body = body;
    
    // Object state
    ca.position_object = state.position;
    ca.velocity_object = state.velocity;
    
    // Planet state
    CelestialBody cb = to_celestial_body(body);
    ca.position_body = PlanetaryEphemeris::getPosition(cb, t);
    ca.velocity_body = PlanetaryEphemeris::getVelocity(cb, t);
    
    // Relative quantities
    ca.rel_position = ca.position_object - ca.position_body;
    ca.rel_velocity = ca.velocity_object - ca.velocity_body;
    
    ca.distance = ca.rel_position.norm();
    ca.relative_velocity = ca.rel_velocity.norm();
    
    // Compute b-plane if requested
    if (settings_.compute_b_plane && ca.distance > settings_.min_distance) {
        ca.b_plane = compute_b_plane(ca);
    }
    
    return ca;
}

BPlaneCoordinates CloseApproachDetector::compute_b_plane(const CloseApproach& ca) const {
    BPlaneCoordinates b_plane;
    
    // B-plane coordinate system:
    // - Origin at planet center
    // - z-axis along incoming velocity (v_rel)
    // - x-axis (ξ) in plane of relative position and velocity
    // - y-axis (ζ) completes right-handed system
    
    Vector3d v_rel = ca.rel_velocity;
    Vector3d r_rel = ca.rel_position;
    
    double v_norm = v_rel.norm();
    if (v_norm < 1e-10) {
        // Degenerate case
        b_plane.xi = 0.0;
        b_plane.zeta = 0.0;
        b_plane.b_magnitude = 0.0;
        b_plane.theta = 0.0;
        return b_plane;
    }
    
    // z-axis along velocity
    Vector3d z_hat = v_rel / v_norm;
    
    // x-axis perpendicular to z in (r,v) plane
    Vector3d h = r_rel.cross(v_rel);  // Angular momentum direction
    Vector3d x_hat = z_hat.cross(h);
    double x_norm = x_hat.norm();
    if (x_norm > 1e-10) {
        x_hat /= x_norm;
    } else {
        // Use arbitrary perpendicular if degenerate
        x_hat = Vector3d(1, 0, 0);
        if (std::abs(z_hat.dot(x_hat)) > 0.9) {
            x_hat = Vector3d(0, 1, 0);
        }
        Vector3d temp = z_hat.cross(x_hat);
        x_hat = temp.cross(z_hat).normalized();
    }
    
    // y-axis completes system
    Vector3d y_hat = z_hat.cross(x_hat);
    
    // Project relative position onto b-plane (perpendicular to z_hat)
    Vector3d r_perp = r_rel - z_hat * (r_rel.dot(z_hat));
    
    // Compute ξ and ζ
    b_plane.xi = r_perp.dot(x_hat);
    b_plane.zeta = r_perp.dot(y_hat);
    b_plane.b_magnitude = r_perp.norm();
    b_plane.theta = std::atan2(b_plane.zeta, b_plane.xi);
    
    return b_plane;
}

CelestialBody CloseApproachDetector::to_celestial_body(BodyType body) const {
    switch (body) {
        case BodyType::MERCURY: return CelestialBody::MERCURY;
        case BodyType::VENUS: return CelestialBody::VENUS;
        case BodyType::EARTH: return CelestialBody::EARTH;
        case BodyType::MARS: return CelestialBody::MARS;
        case BodyType::JUPITER: return CelestialBody::JUPITER;
        case BodyType::SATURN: return CelestialBody::SATURN;
        case BodyType::URANUS: return CelestialBody::URANUS;
        case BodyType::NEPTUNE: return CelestialBody::NEPTUNE;
        case BodyType::MOON: return CelestialBody::MOON;
        default: return CelestialBody::EARTH;
    }
}

double CloseApproachDetector::get_planet_radius(BodyType body) const {
    switch (body) {
        case BodyType::MERCURY: return MERCURY_RADIUS_AU;
        case BodyType::VENUS: return VENUS_RADIUS_AU;
        case BodyType::EARTH: return EARTH_RADIUS_AU;
        case BodyType::MARS: return MARS_RADIUS_AU;
        case BodyType::JUPITER: return JUPITER_RADIUS_AU;
        case BodyType::SATURN: return SATURN_RADIUS_AU;
        case BodyType::URANUS: return URANUS_RADIUS_AU;
        case BodyType::NEPTUNE: return NEPTUNE_RADIUS_AU;
        case BodyType::MOON: return MOON_RADIUS_AU;
        default: return EARTH_RADIUS_AU;
    }
}

// ============================================================================
// MOIDCalculator Implementation
// ============================================================================

double MOIDCalculator::compute_moid(
    const KeplerianElements& object_orbit,
    BodyType planet)
{
    // Get planet mean orbital elements at J2000.0
    KeplerianElements planet_orbit = get_planet_mean_elements(planet);
    
    // Use two-orbit MOID calculation
    // Note: This uses mean elements, so accuracy is limited.
    // For precise MOID, should use osculating elements at specific epoch.
    return compute_moid(object_orbit, planet_orbit);
}

double MOIDCalculator::compute_moid(
    const KeplerianElements& orbit1,
    const KeplerianElements& orbit2)
{
    // MOID calculation via optimization
    // This is a simplified placeholder - full implementation would use
    // sophisticated optimization (e.g., Nelder-Mead, genetic algorithm)
    
    double min_dist_sq = 1e100;
    
    // Grid search over true anomalies
    const int n_points = 360;
    for (int i = 0; i < n_points; ++i) {
        double f1 = i * TWO_PI / n_points;
        for (int j = 0; j < n_points; ++j) {
            double f2 = j * TWO_PI / n_points;
            double dist_sq = distance_squared(f1, f2, orbit1, orbit2);
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
            }
        }
    }
    
    return std::sqrt(min_dist_sq);
}

double MOIDCalculator::distance_squared(
    double f1, double f2,
    const KeplerianElements& orbit1,
    const KeplerianElements& orbit2)
{
    // Compute position on orbit 1 at true anomaly f1
    // r = a(1-e²)/(1+e*cos(f))
    double r1 = orbit1.semi_major_axis * (1.0 - orbit1.eccentricity * orbit1.eccentricity)
              / (1.0 + orbit1.eccentricity * std::cos(f1));
    
    // Position in orbital plane coordinates
    double x1_orb = r1 * std::cos(f1);
    double y1_orb = r1 * std::sin(f1);
    
    // Rotate to inertial frame
    double cos_w1 = std::cos(orbit1.argument_perihelion);
    double sin_w1 = std::sin(orbit1.argument_perihelion);
    double cos_Om1 = std::cos(orbit1.longitude_ascending_node);
    double sin_Om1 = std::sin(orbit1.longitude_ascending_node);
    double cos_i1 = std::cos(orbit1.inclination);
    double sin_i1 = std::sin(orbit1.inclination);
    
    double x1 = (cos_Om1 * cos_w1 - sin_Om1 * sin_w1 * cos_i1) * x1_orb
              + (-cos_Om1 * sin_w1 - sin_Om1 * cos_w1 * cos_i1) * y1_orb;
    double y1 = (sin_Om1 * cos_w1 + cos_Om1 * sin_w1 * cos_i1) * x1_orb
              + (-sin_Om1 * sin_w1 + cos_Om1 * cos_w1 * cos_i1) * y1_orb;
    double z1 = sin_i1 * sin_w1 * x1_orb + sin_i1 * cos_w1 * y1_orb;
    
    // Same for orbit 2
    double r2 = orbit2.semi_major_axis * (1.0 - orbit2.eccentricity * orbit2.eccentricity)
              / (1.0 + orbit2.eccentricity * std::cos(f2));
    double x2_orb = r2 * std::cos(f2);
    double y2_orb = r2 * std::sin(f2);
    
    double cos_w2 = std::cos(orbit2.argument_perihelion);
    double sin_w2 = std::sin(orbit2.argument_perihelion);
    double cos_Om2 = std::cos(orbit2.longitude_ascending_node);
    double sin_Om2 = std::sin(orbit2.longitude_ascending_node);
    double cos_i2 = std::cos(orbit2.inclination);
    double sin_i2 = std::sin(orbit2.inclination);
    
    double x2 = (cos_Om2 * cos_w2 - sin_Om2 * sin_w2 * cos_i2) * x2_orb
              + (-cos_Om2 * sin_w2 - sin_Om2 * cos_w2 * cos_i2) * y2_orb;
    double y2 = (sin_Om2 * cos_w2 + cos_Om2 * sin_w2 * cos_i2) * x2_orb
              + (-sin_Om2 * sin_w2 + cos_Om2 * cos_w2 * cos_i2) * y2_orb;
    double z2 = sin_i2 * sin_w2 * x2_orb + sin_i2 * cos_w2 * y2_orb;
    
    // Distance squared
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    
    return dx*dx + dy*dy + dz*dz;
}

} // namespace orbfit::close_approach
