/**
 * @file CloseApproach.cpp
 * @brief Implementation of close approach detection
 * @author ITALOccult AstDyn Team
 * @date 2025-11-24
 */

#include "astdyn/close_approach/CloseApproach.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <map>

namespace astdyn::close_approach {

using namespace astdyn::propagation;
using namespace astdyn::ephemeris;
using namespace astdyn::constants;

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
    
    KeplerianElements get_planet_mean_elements(BodyType planet) {
        KeplerianElements elements;
        elements.epoch = time::EpochTDB::from_mjd(0.0); // MJD J2000.0
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
    std::shared_ptr<PlanetaryEphemeris> ephemeris,
    const CloseApproachSettings& settings)
    : propagator_(propagator), ephemeris_(ephemeris), settings_(settings)
{
    if (!propagator_) {
        throw std::invalid_argument("CloseApproachDetector: null propagator");
    }
}

std::vector<CloseApproach> CloseApproachDetector::detect(
    const physics::CartesianStateTyped<core::GCRF>& initial_state,
    time::EpochTDB t_start,
    time::EpochTDB t_end)
{
    std::vector<CloseApproach> approaches;
    
    if (t_end.mjd() <= t_start.mjd()) {
        throw std::invalid_argument("t_end must be > t_start");
    }
    
    // Determine which bodies to check
    std::vector<BodyType> bodies;
    if (settings_.bodies_to_check.empty()) {
        bodies = {BodyType::MERCURY, BodyType::VENUS, BodyType::EARTH, 
                  BodyType::MARS, BodyType::JUPITER, BodyType::SATURN,
                  BodyType::URANUS, BodyType::NEPTUNE};
    } else {
        bodies = settings_.bodies_to_check;
    }
    
    double dt = 1.0;
    double t = t_start.mjd();
    physics::CartesianStateTyped<core::GCRF> prev_state = initial_state;
    
    std::map<BodyType, double> prev_distances;
    for (auto body : bodies) {
        prev_distances[body] = compute_distance(t_start, prev_state, body);
    }
    
    while (t < t_end.mjd()) {
        double t_next_val = std::min(t + dt, t_end.mjd());
        time::EpochTDB t_next = time::EpochTDB::from_mjd(t_next_val);
        
        auto current_state = propagator_->propagate_cartesian(prev_state, t_next);
        
        for (auto body : bodies) {
            if (!settings_.should_check_body(body)) continue;
            
            double dist_current = compute_distance(t_next, current_state, body);
            double dist_prev = prev_distances[body];
            
            if (dist_current > dist_prev && dist_prev < settings_.detection_distance) {
                time::EpochTDB t_ca = t_next;
                if (settings_.refine_time) {
                    t_ca = detect_approach_time(prev_state, current_state, 
                                                time::EpochTDB::from_mjd(t), 
                                                t_next, body);
                }
                
                auto state_ca = propagator_->propagate_cartesian(prev_state, t_ca);
                CloseApproach ca = build_approach(t_ca, state_ca, body);
                
                if (ca.distance < settings_.detection_distance) {
                    approaches.push_back(ca);
                }
            }
            prev_distances[body] = dist_current;
        }
        
        t = t_next_val;
        prev_state = current_state;
    }
    
    std::sort(approaches.begin(), approaches.end(),
              [](const CloseApproach& a, const CloseApproach& b) {
                  return a.time.mjd() < b.time.mjd();
              });
    
    return approaches;
}

std::vector<CloseApproach> CloseApproachDetector::detect(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_orbit,
    time::EpochTDB t_start,
    time::EpochTDB t_end)
{
    // Convert orbit to Cartesian GCRF using typed converters
    auto initial_ecl = propagation::keplerian_to_cartesian<core::ECLIPJ2000>(initial_orbit);
    
    auto pos_gcrf = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(initial_ecl.position, initial_orbit.epoch);
    auto vel_gcrf = coordinates::ReferenceFrame::transform_vel<core::ECLIPJ2000, core::GCRF>(initial_ecl.position, initial_ecl.velocity, initial_orbit.epoch);
    
    auto cart_typed = physics::CartesianStateTyped<core::GCRF>(
        initial_orbit.epoch, pos_gcrf, vel_gcrf, initial_orbit.gm
    );
    
    return detect(cart_typed, t_start, t_end);
}

double CloseApproachDetector::compute_distance(
    time::EpochTDB t,
    const physics::CartesianStateTyped<core::GCRF>& state,
    BodyType body) const
{
    CelestialBody cb = to_celestial_body(body);
    auto actual_ephem = ephemeris_ ? ephemeris_ : std::make_shared<PlanetaryEphemeris>();
    auto planet_pos = actual_ephem->getPosition(cb, t);
    return (state.position.to_eigen_si() - planet_pos.to_eigen_si()).norm();
}

time::EpochTDB CloseApproachDetector::detect_approach_time(
    const physics::CartesianStateTyped<core::GCRF>& state_t1,
    const physics::CartesianStateTyped<core::GCRF>& state_t2,
    time::EpochTDB t1,
    time::EpochTDB t2,
    BodyType body) const
{
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    const double resphi = 2.0 - phi;
    
    double a = t1.mjd();
    double b = t2.mjd();
    double tol = settings_.time_tolerance;
    
    double x1_val = a + resphi * (b - a);
    double x2_val = b - resphi * (b - a);
    
    time::EpochTDB x1 = time::EpochTDB::from_mjd(x1_val);
    time::EpochTDB x2 = time::EpochTDB::from_mjd(x2_val);
    
    auto state_x1 = propagator_->propagate_cartesian(state_t1, x1);
    auto state_x2 = propagator_->propagate_cartesian(state_t1, x2);
    
    double f1 = compute_distance(x1, state_x1, body);
    double f2 = compute_distance(x2, state_x2, body);
    
    int iter = 0;
    while (std::abs(b - a) > tol && iter < settings_.max_refinement_iter) {
        if (f1 < f2) {
            b = x2_val;
            x2_val = x1_val;
            f2 = f1;
            x1_val = a + resphi * (b - a);
            x1 = time::EpochTDB::from_mjd(x1_val);
            state_x1 = propagator_->propagate_cartesian(state_t1, x1);
            f1 = compute_distance(x1, state_x1, body);
        } else {
            a = x1_val;
            x1_val = x2_val;
            f1 = f2;
            x2_val = b - resphi * (b - a);
            x2 = time::EpochTDB::from_mjd(x2_val);
            state_x2 = propagator_->propagate_cartesian(state_t1, x2);
            f2 = compute_distance(x2, state_x2, body);
        }
        iter++;
    }
    
    return time::EpochTDB::from_mjd((a + b) / 2.0);
}

CloseApproach CloseApproachDetector::build_approach(
    time::EpochTDB t,
    const physics::CartesianStateTyped<core::GCRF>& state,
    BodyType body) const
{
    CloseApproach ca;
    ca.time = t;
    ca.body = body;
    
    ca.position_object = state.position;
    ca.velocity_object = state.velocity;
    
    CelestialBody cb = to_celestial_body(body);
    auto actual_ephem = ephemeris_ ? ephemeris_ : std::make_shared<PlanetaryEphemeris>();
    ca.position_body = actual_ephem->getPosition(cb, t);
    ca.velocity_body = actual_ephem->getVelocity(cb, t);
    
    ca.rel_position = ca.position_object - ca.position_body;
    ca.rel_velocity = ca.velocity_object - ca.velocity_body;
    
    ca.distance = ca.rel_position.norm().to_m();
    ca.relative_velocity_mag = ca.rel_velocity.norm().to_ms();
    
    if (settings_.compute_b_plane && ca.distance > settings_.min_distance) {
        ca.b_plane = compute_b_plane(ca);
    }
    
    return ca;
}

BPlaneCoordinates CloseApproachDetector::compute_b_plane(const CloseApproach& ca) const {
    BPlaneCoordinates b_plane;
    Vector3d v_rel = ca.rel_velocity.to_eigen_si();
    Vector3d r_rel = ca.rel_position.to_eigen_si();
    
    double v_norm = v_rel.norm();
    if (v_norm < 1e-10) {
        b_plane.xi = 0.0; b_plane.zeta = 0.0; b_plane.b_magnitude = 0.0; b_plane.theta = 0.0;
        return b_plane;
    }
    
    Vector3d z_hat = v_rel / v_norm;
    Vector3d h = r_rel.cross(v_rel);
    Vector3d x_hat = z_hat.cross(h);
    double x_norm = x_hat.norm();
    if (x_norm > 1e-10) x_hat /= x_norm;
    else {
        x_hat = Vector3d(1, 0, 0);
        if (std::abs(z_hat.dot(x_hat)) > 0.9) x_hat = Vector3d(0, 1, 0);
        x_hat = z_hat.cross(x_hat).cross(z_hat).normalized();
    }
    
    Vector3d y_hat = z_hat.cross(x_hat);
    Vector3d r_perp = r_rel - z_hat * (r_rel.dot(z_hat));
    
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

double MOIDCalculator::compute_moid(const KeplerianElements& object_orbit, BodyType planet) {
    KeplerianElements planet_orbit = get_planet_mean_elements(planet);
    return compute_moid(object_orbit, planet_orbit);
}

double MOIDCalculator::compute_moid(const KeplerianElements& orbit1, const KeplerianElements& orbit2) {
    double min_dist_sq = 1e100;
    const int n_points = 360;
    for (int i = 0; i < n_points; ++i) {
        double f1 = i * TWO_PI / n_points;
        for (int j = 0; j < n_points; ++j) {
            double f2 = j * TWO_PI / n_points;
            double dist_sq = distance_squared(f1, f2, orbit1, orbit2);
            if (dist_sq < min_dist_sq) min_dist_sq = dist_sq;
        }
    }
    return std::sqrt(min_dist_sq);
}

double MOIDCalculator::distance_squared(double f1, double f2, const KeplerianElements& orbit1, const KeplerianElements& orbit2) {
    double r1 = orbit1.semi_major_axis * (1.0 - orbit1.eccentricity * orbit1.eccentricity) / (1.0 + orbit1.eccentricity * std::cos(f1));
    double x1_orb = r1 * std::cos(f1);
    double y1_orb = r1 * std::sin(f1);
    double cos_w1 = std::cos(orbit1.argument_perihelion); double sin_w1 = std::sin(orbit1.argument_perihelion);
    double cos_Om1 = std::cos(orbit1.longitude_ascending_node); double sin_Om1 = std::sin(orbit1.longitude_ascending_node);
    double cos_i1 = std::cos(orbit1.inclination); double sin_i1 = std::sin(orbit1.inclination);
    double x1 = (cos_Om1 * cos_w1 - sin_Om1 * sin_w1 * cos_i1) * x1_orb + (-cos_Om1 * sin_w1 - sin_Om1 * cos_w1 * cos_i1) * y1_orb;
    double y1 = (sin_Om1 * cos_w1 + cos_Om1 * sin_w1 * cos_i1) * x1_orb + (-sin_Om1 * sin_w1 + cos_Om1 * cos_w1 * cos_i1) * y1_orb;
    double z1 = sin_i1 * sin_w1 * x1_orb + sin_i1 * cos_w1 * y1_orb;
    double r2 = orbit2.semi_major_axis * (1.0 - orbit2.eccentricity * orbit2.eccentricity) / (1.0 + orbit2.eccentricity * std::cos(f2));
    double x2_orb = r2 * std::cos(f2); double y2_orb = r2 * std::sin(f2);
    double cos_w2 = std::cos(orbit2.argument_perihelion); double sin_w2 = std::sin(orbit2.argument_perihelion);
    double cos_Om2 = std::cos(orbit2.longitude_ascending_node); double sin_Om2 = std::sin(orbit2.longitude_ascending_node);
    double cos_i2 = std::cos(orbit2.inclination); double sin_i2 = std::sin(orbit2.inclination);
    double x2 = (cos_Om2 * cos_w2 - sin_Om2 * sin_w2 * cos_i2) * x2_orb + (-cos_Om2 * sin_w2 - sin_Om2 * cos_w2 * cos_i2) * y2_orb;
    double y2 = (sin_Om2 * cos_w2 + cos_Om2 * sin_w2 * cos_i2) * x2_orb + (-sin_Om2 * sin_w2 + cos_Om2 * cos_w2 * cos_i2) * y2_orb;
    double z2 = sin_i2 * sin_w2 * x2_orb + sin_i2 * cos_w2 * y2_orb;
    double dx = x2 - x1; double dy = y2 - y1; double dz = z2 - z1;
    return dx*dx + dy*dy + dz*dz;
}

} // namespace astdyn::close_approach
