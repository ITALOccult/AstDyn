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
    time::EpochTDB t_start, time::EpochTDB t_end)
{
    std::vector<CloseApproach> approaches;
    if (t_end.mjd() <= t_start.mjd()) throw std::invalid_argument("t_end must be > t_start");
    
    std::vector<BodyType> bodies = settings_.bodies_to_check.empty() ? 
        std::vector<BodyType>{BodyType::MERCURY, BodyType::VENUS, BodyType::EARTH, BodyType::MARS,
                               BodyType::JUPITER, BodyType::SATURN, BodyType::URANUS, BodyType::NEPTUNE} : 
        settings_.bodies_to_check;
    
    double dt = 1.0; // Initial search step (days)
    double t = t_start.mjd();
    physics::CartesianStateTyped<core::GCRF> prev_state = initial_state;
    
    // Initialize distances for each body at start
    std::map<BodyType, double> prev_distances;
    for (auto b : bodies) {
        prev_distances[b] = compute_distance(t_start, prev_state, b);
    }
    
    while (t < t_end.mjd()) {
        double t_next_val = std::min(t + dt, t_end.mjd());
        time::EpochTDB t_next = time::EpochTDB::from_mjd(t_next_val);
        
        auto current_state = propagator_->propagate_cartesian(prev_state, t_next);
        
        for (auto b : bodies) {
            check_body_approaches(b, t_next, current_state, prev_state, t, prev_distances, approaches);
        }
        
        t = t_next_val;
        prev_state = current_state;
    }
    
    // Final sort to guarantee chronological order
    std::sort(approaches.begin(), approaches.end(), [](const auto& a, const auto& b){ 
        return a.time.mjd() < b.time.mjd(); 
    });
    
    return approaches;
}

void CloseApproachDetector::check_body_approaches(BodyType body, time::EpochTDB t_next,
    const physics::CartesianStateTyped<core::GCRF>& current_state,
    const physics::CartesianStateTyped<core::GCRF>& prev_state,
    double t_prev_val, std::map<BodyType, double>& prev_distances,
    std::vector<CloseApproach>& approaches)
{
    if (!settings_.should_check_body(body)) return;
    
    double d_curr = compute_distance(t_next, current_state, body);
    double d_prev = prev_distances[body];

    // Local minimum detected if slope changes from negative to positive
    // OR if we are inside the detection distance and distances are decreasing
    if (d_curr > d_prev && d_prev < settings_.detection_distance) {
        // Refine the exact time of closest approach
        time::EpochTDB t_ca = settings_.refine_time ? 
            detect_approach_time(prev_state, current_state, time::EpochTDB::from_mjd(t_prev_val), t_next, body) : t_next;
        
        auto state_ca = propagator_->propagate_cartesian(prev_state, t_ca);
        auto ca = build_approach(t_ca, state_ca, body);
        
        if (ca.distance < settings_.detection_distance) {
            approaches.push_back(ca);
        }
    }
    
    // Keep track of previous distance for slope detection
    prev_distances[body] = d_curr;
}

std::vector<CloseApproach> CloseApproachDetector::detect(
    const physics::KeplerianStateTyped<core::ECLIPJ2000>& initial_orbit,
    time::EpochTDB t_start,
    time::EpochTDB t_end)
{
    // Convert orbit to Cartesian state for propagation using existing converters
    auto cart_legacy = propagation::keplerian_to_cartesian<core::GCRF>(initial_orbit);
    return detect(cart_legacy, t_start, t_end);
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
    const physics::CartesianStateTyped<core::GCRF>& s1, const physics::CartesianStateTyped<core::GCRF>& s2,
    time::EpochTDB t1, time::EpochTDB t2, BodyType body) const
{
    // Golden-section search for minimum of distance(t)
    const double PHI = (1.0 + std::sqrt(5.0)) / 2.0;
    double a = t1.mjd(), b = t2.mjd(), tol = settings_.time_tolerance;
    double x1 = b - (b - a) / PHI, x2 = a + (b - a) / PHI;
    
    auto f = [&](double mjd) {
        auto t = time::EpochTDB::from_mjd(mjd);
        // We use the propagator to interpolate/propagate between the two steps
        return compute_distance(t, propagator_->propagate_cartesian(s1, t), body);
    };

    double f1 = f(x1), f2 = f(x2);
    for (int i = 0; i < settings_.max_refinement_iter && std::abs(b - a) > tol; ++i) {
        if (f1 < f2) { 
            b = x2; x2 = x1; f2 = f1; 
            x1 = b - (b - a) / PHI; f1 = f(x1); 
        } else { 
            a = x1; x1 = x2; f1 = f2; 
            x2 = a + (b - a) / PHI; f2 = f(x2); 
        }
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
    Eigen::Vector3d v_rel = ca.rel_velocity.to_eigen_si();
    Eigen::Vector3d r_rel = ca.rel_position.to_eigen_si();
    
    double v_norm = v_rel.norm();
    if (v_norm < 1e-10) {
        b_plane.xi = 0.0; b_plane.zeta = 0.0; b_plane.b_magnitude = 0.0; b_plane.theta = 0.0;
        return b_plane;
    }
    
    // Basis vectors for b-plane
    Eigen::Vector3d k_hat = v_rel / v_norm; // velocity direction
    Eigen::Vector3d h = r_rel.cross(v_rel);
    
    if (h.norm() < 1e-10) {
        // Direct radial approach
        b_plane.xi = 0.0; b_plane.zeta = 0.0; b_plane.b_magnitude = 0.0; b_plane.theta = 0.0;
        return b_plane;
    }

    // Basis according to Öpik/Kizner theory
    Eigen::Vector3d x_hat = k_hat.cross(Eigen::Vector3d::UnitZ());
    if (x_hat.norm() < 1e-4) x_hat = k_hat.cross(Eigen::Vector3d::UnitX());
    x_hat.normalize();
    Eigen::Vector3d y_hat = k_hat.cross(x_hat);

    // Vector b points from focal point to the point of closest approach
    // In GCRF, we project onto the x_hat, y_hat basis
    b_plane.xi = r_rel.dot(x_hat);
    b_plane.zeta = r_rel.dot(y_hat);
    b_plane.b_magnitude = std::sqrt(b_plane.xi * b_plane.xi + b_plane.zeta * b_plane.zeta);
    b_plane.theta = std::atan2(b_plane.zeta, b_plane.xi);
    
    return b_plane;
}

CelestialBody CloseApproachDetector::to_celestial_body(BodyType body) const {
    return static_cast<CelestialBody>(body);
}

double CloseApproachDetector::get_planet_radius(BodyType body) const {
    return PlanetaryData::getRadius(static_cast<CelestialBody>(body)) * 1000.0;
}

// ============================================================================
// MOIDCalculator Implementation
// ============================================================================

double MOIDCalculator::compute_moid(const KeplerianElements& object_orbit, BodyType planet) {
    KeplerianElements planet_orbit = get_planet_mean_elements(planet);
    return compute_moid(object_orbit, planet_orbit);
}

double MOIDCalculator::compute_moid(const KeplerianElements& orbit1, const KeplerianElements& orbit2) {
    double min_dist_sq = 1e30;
    
    // Coarse search grid (72 x 72 -> 5 deg step)
    const int N_coarse = 72;
    double f1_min = 0, f2_min = 0;
    
    for (int i = 0; i < N_coarse; ++i) {
        double f1 = i * TWO_PI / N_coarse;
        for (int j = 0; j < N_coarse; ++j) {
            double f2 = j * TWO_PI / N_coarse;
            double d2 = distance_squared(f1, f2, orbit1, orbit2);
            if (d2 < min_dist_sq) {
                min_dist_sq = d2;
                f1_min = f1;
                f2_min = f2;
            }
        }
    }
    
    // Fine refinement around candidate (12 x 12 covering ±5 deg)
    const int N_fine = 12;
    double step_fine = (5.0 * DEG_TO_RAD) / N_fine;
    for (int i = -N_fine; i <= N_fine; ++i) {
        double f1 = f1_min + i * (step_fine / 2.0);
        for (int j = -N_fine; j <= N_fine; ++j) {
            double f2 = f2_min + j * (step_fine / 2.0);
            min_dist_sq = std::min(min_dist_sq, distance_squared(f1, f2, orbit1, orbit2));
        }
    }
    
    return std::sqrt(min_dist_sq);
}

double MOIDCalculator::distance_squared(double f1, double f2, const KeplerianElements& o1, const KeplerianElements& o2) {
    // We treat f1/f2 as mean anomalies for path distance search
    KeplerianElements tmp1 = o1; tmp1.mean_anomaly = f1;
    KeplerianElements tmp2 = o2; tmp2.mean_anomaly = f2;
    
    auto p1 = propagation::keplerian_to_cartesian(tmp1).position;
    auto p2 = propagation::keplerian_to_cartesian(tmp2).position;
    
    return (p1.to_eigen_si() - p2.to_eigen_si()).squaredNorm();
}

KeplerianElements MOIDCalculator::get_planet_mean_elements(BodyType planet) {
    KeplerianElements els;
    els.epoch = time::EpochTDB::from_mjd(MJD2000);
    els.gravitational_parameter = GMS;
    
    // Get basic radius and semi-major axis from repository
    auto data = PlanetaryData::getBodyData(static_cast<CelestialBody>(planet));
    els.semi_major_axis = data.semi_major_axis;
    
    // The following are J2000 approximations
    switch (planet) {
        case BodyType::MERCURY: els.eccentricity = 0.2056; els.inclination = 7.005 * DEG_TO_RAD; break;
        case BodyType::VENUS:   els.eccentricity = 0.0068; els.inclination = 3.394 * DEG_TO_RAD; break;
        case BodyType::EARTH:   els.eccentricity = 0.0167; els.inclination = 0.000 * DEG_TO_RAD; break;
        case BodyType::MARS:    els.eccentricity = 0.0934; els.inclination = 1.850 * DEG_TO_RAD; break;
        case BodyType::JUPITER: els.eccentricity = 0.0485; els.inclination = 1.303 * DEG_TO_RAD; break;
        case BodyType::SATURN:  els.eccentricity = 0.0555; els.inclination = 2.489 * DEG_TO_RAD; break;
        case BodyType::URANUS:  els.eccentricity = 0.0469; els.inclination = 0.773 * DEG_TO_RAD; break;
        case BodyType::NEPTUNE: els.eccentricity = 0.0090; els.inclination = 1.770 * DEG_TO_RAD; break;
        default: break;
    }
    return els;
}

} // namespace astdyn::close_approach
