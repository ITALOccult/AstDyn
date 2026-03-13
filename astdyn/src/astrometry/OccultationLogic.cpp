/**
 * @file OccultationLogic.cpp
 * @brief Implementation of occultation physical parameters computation.
 */

#include "astdyn/astrometry/OccultationLogic.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>

namespace astdyn::astrometry {

OccultationParameters OccultationLogic::compute_parameters(
    const RightAscension& star_ra, const Declination& star_dec,
    const RightAscension& ast_ra, const Declination& ast_dec,
    double ast_dist_m,
    double ast_dra_dt_rad_s, double ast_ddec_dt_rad_s,
    double ast_ddist_dt_m_s) 
{
    using namespace constants;

    double as = star_ra.to_rad();
    double ds = star_dec.to_rad();
    double a = ast_ra.to_rad();
    double d = ast_dec.to_rad();
    double rho = ast_dist_m;

    // 1. Besselian Basis (toward star)
    Eigen::Vector3d k(std::cos(as) * std::cos(ds), std::sin(as) * std::cos(ds), std::sin(ds));
    Eigen::Vector3d i(-std::sin(as), std::cos(as), 0.0);
    // BUG-3: Use the rigorous North vector formula provided by the user
    // j = (-sin(ds)cos(as), -sin(ds)sin(as), cos(ds))
    Eigen::Vector3d j(-std::sin(ds) * std::cos(as), -std::sin(ds) * std::sin(as), std::cos(ds));

    // 2. Asteroid and its velocity in 3D (GCRF)
    Eigen::Vector3d r_ast_vec(
        rho * std::cos(a) * std::cos(d),
        rho * std::sin(a) * std::cos(d),
        rho * std::sin(d)
    );

    // Derivatives of r_ast vector
    double cos_a = std::cos(a);
    double sin_a = std::sin(a);
    double cos_d = std::cos(d);
    double sin_d = std::sin(d);

    Eigen::Vector3d v_ast_vec(
        ast_ddist_dt_m_s * cos_a * cos_d - rho * sin_a * cos_d * ast_dra_dt_rad_s - rho * cos_a * sin_d * ast_ddec_dt_rad_s,
        ast_ddist_dt_m_s * sin_a * cos_d + rho * cos_a * cos_d * ast_dra_dt_rad_s - rho * sin_a * sin_d * ast_ddec_dt_rad_s,
        ast_ddist_dt_m_s * sin_d + rho * cos_d * ast_ddec_dt_rad_s
    );

    // 3. Project to Fundamental Plane (xi, eta)
    double xi = r_ast_vec.dot(i);
    double eta = r_ast_vec.dot(j);
    double dxi = v_ast_vec.dot(i);
    double deta = v_ast_vec.dot(j);

    // 4. Closest Approach Analysis (Linear approximation on the plane)
    // Distance squared: f(t) = (xi + dxi*t)^2 + (eta + deta*t)^2
    // f'(t) = 2(xi + dxi*t)dxi + 2(eta + deta*t)deta = 0
    // xi*dxi + dxi^2*t + eta*deta + deta^2*t = 0
    // t_min = -(xi*dxi + eta*deta) / (dxi^2 + deta^2)
    
    double v2 = dxi * dxi + deta * deta;
    double t_ca = 0.0;
    if (v2 > 1e-18) {
        t_ca = -(xi * dxi + eta * deta) / v2;
    }

    double xi_ca = xi + dxi * t_ca;
    double eta_ca = eta + deta * t_ca;
    double b = std::sqrt(xi_ca * xi_ca + eta_ca * eta_ca);

    // 5. Build Result
    OccultationParameters params;
    params.impact_parameter_km = b / 1000.0;
    params.shadow_velocity_kms = std::sqrt(v2) / 1000.0;
    
    // Position angle of the track (direction of velocity on plane)
    // atan2(E, N)
    params.position_angle_deg = std::atan2(dxi, deta) * RAD_TO_DEG;
    if (params.position_angle_deg < 0) params.position_angle_deg += 360.0;

    // Angular velocities in mas/s
    params.d_ra_cos_dec_mas_sec = ast_dra_dt_rad_s * std::cos(d) * (RAD_TO_DEG * 3600.0 * 1000.0);
    params.d_dec_mas_sec = ast_ddec_dt_rad_s * (RAD_TO_DEG * 3600.0 * 1000.0);
    
    params.closest_approach_time_offset_sec = t_ca;
    params.time_uncertainty_sec = 0.0;
    params.cross_track_uncertainty_km = 0.0;

    return params;
}

} // namespace astdyn::astrometry
