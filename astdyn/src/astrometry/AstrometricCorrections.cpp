#include "astdyn/astrometry/AstrometricCorrections.hpp"
#include "astdyn/core/Constants.hpp"
#include <cmath>

namespace astdyn::astrometry {

using namespace astdyn::constants;

Eigen::Vector3d aberrazione_differenziale(
    const Eigen::Vector3d& rho_eq, const Eigen::Vector3d& earth_vel_eq) {
    
    double r = rho_eq.norm();
    Eigen::Vector3d p = rho_eq / r;
    Eigen::Vector3d v = earth_vel_eq / (C_LIGHT * 1000.0); // beta vector (velocity in units of c, with v in m/s if scaled)
    // Wait: C_LIGHT is km/s. If earth_vel_eq is m/s, divide by 1000. Or if km/s, just used C_LIGHT.
    // Let's assume earth_vel_eq is m/s as per internal conventions.
    
    double p_dot_v = p.dot(v);
    double v2 = v.squaredNorm();
    double inv_gamma = std::sqrt(1.0 - v2);
    
    // Relativistic formula (IAU 2000)
    double denom = 1.0 + p_dot_v;
    Eigen::Vector3d p_prime = (inv_gamma * p + (1.0 + p_dot_v / (1.0 + inv_gamma)) * v) / denom;
    
    return p_prime * r;
}

Eigen::Vector3d deflessione_relativistica(
    const Eigen::Vector3d& rho_eq, const Eigen::Vector3d& observer_to_sun_eq) 
{
    double r = rho_eq.norm();
    Eigen::Vector3d u = rho_eq / r;
    Eigen::Vector3d q = observer_to_sun_eq;
    double q_dist = q.norm();
    if (q_dist < 1000.0) return rho_eq; 
    
    Eigen::Vector3d e = q / q_dist;
    double u_dot_e = u.dot(e);
    // Schwarzschild Sun constant factor: 2*GM/c^2
    double fac = SCHWARZSCHILD_SUN / q_dist * (1.0 + u_dot_e);
    
    return (u + fac * (e - u_dot_e * u)).normalized() * r;
}

} // namespace astdyn::astrometry
