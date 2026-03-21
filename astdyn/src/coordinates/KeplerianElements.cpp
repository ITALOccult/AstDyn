/**
 * @file KeplerianElements.cpp
 * @brief Implementation of Keplerian elements conversions
 */

#include "astdyn/coordinates/KeplerianElements.hpp"
#include "astdyn/math/MathUtils.hpp"
#include <algorithm>
#include <stdexcept>

namespace astdyn {
namespace coordinates {

CartesianState KeplerianElements::to_cartesian() const {
    // Convert Keplerian elements to Cartesian state
    // Reference: Vallado, "Fundamentals of Astrodynamics and Applications"
    
    // Get true anomaly
    double nu = true_anomaly();
    
    // Compute position and velocity in orbital plane (perifocal frame)
    double p = semi_latus_rectum();
    double r = p / (1.0 + e_ * std::cos(nu));
    
    // Position in perifocal frame
    Vector3d r_pqw;
    r_pqw << r * std::cos(nu),
             r * std::sin(nu),
             0.0;
    
    // Velocity in perifocal frame
    double sqrt_mu_p = std::sqrt(mu_ / p);
    Vector3d v_pqw;
    v_pqw << -sqrt_mu_p * std::sin(nu),
              sqrt_mu_p * (e_ + std::cos(nu)),
              0.0;
    
    // Rotation matrix from perifocal to inertial frame
    // R = R3(-Ω) * R1(-i) * R3(-ω)
    double cos_O = std::cos(Omega_);
    double sin_O = std::sin(Omega_);
    double cos_i = std::cos(i_);
    double sin_i = std::sin(i_);
    double cos_o = std::cos(omega_);
    double sin_o = std::sin(omega_);
    
    Matrix3d R;
    R(0, 0) =  cos_O * cos_o - sin_O * sin_o * cos_i;
    R(0, 1) = -cos_O * sin_o - sin_O * cos_o * cos_i;
    R(0, 2) =  sin_O * sin_i;
    
    R(1, 0) =  sin_O * cos_o + cos_O * sin_o * cos_i;
    R(1, 1) = -sin_O * sin_o + cos_O * cos_o * cos_i;
    R(1, 2) = -cos_O * sin_i;
    
    R(2, 0) =  sin_o * sin_i;
    R(2, 1) =  cos_o * sin_i;
    R(2, 2) =  cos_i;
    
    // Transform to inertial frame
    Vector3d position = R * r_pqw;
    Vector3d velocity = R * v_pqw;
    
    return CartesianState(position, velocity, mu_);
}

KeplerianElements KeplerianElements::from_cartesian(const CartesianState& state) {
    // Convert Cartesian state to Keplerian elements
    // Reference: Vallado, "Fundamentals of Astrodynamics and Applications"
    
    Vector3d r = state.position();
    Vector3d v = state.velocity();
    double mu = state.mu();
    
    double r_mag = r.norm();
    double v_mag = v.norm();
    
    // Angular momentum
    Vector3d h = r.cross(v);
    double h_mag = h.norm();
    
    // Node vector
    Vector3d k(0, 0, 1);
    Vector3d n = k.cross(h);
    double n_mag = n.norm();
    
    // Eccentricity vector
    Vector3d e_vec = ((v_mag * v_mag - mu / r_mag) * r - r.dot(v) * v) / mu;
    double e = e_vec.norm();
    
    // Specific energy
    double energy = v_mag * v_mag / 2.0 - mu / r_mag;
    
    // Semi-major axis
    double a;
    if (std::abs(energy) < 1e-10) {
        // Parabolic orbit
        a = std::numeric_limits<double>::infinity();
    } else {
        a = -mu / (2.0 * energy);
    }
    
    // Inclination
    double i = std::acos(std::clamp(h.z() / h_mag, -1.0, 1.0));
    
    // RAAN (Right Ascension of Ascending Node)
    double Omega;
    if (n_mag > 1e-10) {
        Omega = std::acos(std::clamp(n.x() / n_mag, -1.0, 1.0));
        if (n.y() < 0.0) {
            Omega = 2.0 * constants::PI - Omega;
        }
    } else {
        // Equatorial orbit
        Omega = 0.0;
    }
    
    // Argument of periapsis
    double omega;
    if (n_mag > 1e-10 && e > 1e-10) {
        omega = std::acos(std::clamp(n.dot(e_vec) / (n_mag * e), -1.0, 1.0));
        if (e_vec.z() < 0.0) {
            omega = 2.0 * constants::PI - omega;
        }
    } else if (e > 1e-10) {
        // Equatorial orbit with eccentricity
        omega = std::atan2(e_vec.y(), e_vec.x());
    } else {
        // Circular orbit
        omega = 0.0;
    }
    
    // True anomaly
    double nu;
    if (e > 1e-10) {
        nu = std::acos(std::clamp(e_vec.dot(r) / (e * r_mag), -1.0, 1.0));
        if (r.dot(v) < 0.0) {
            nu = 2.0 * constants::PI - nu;
        }
    } else {
        // Circular orbit - use argument of latitude
        if (n_mag > 1e-10) {
            nu = std::acos(std::clamp(n.dot(r) / (n_mag * r_mag), -1.0, 1.0));
            if (r.z() < 0.0) {
                nu = 2.0 * constants::PI - nu;
            }
        } else {
            nu = std::atan2(r.y(), r.x());
        }
    }
    
    // Convert true anomaly to mean anomaly
    double M = true_to_mean_anomaly(nu, e);
    
    // Normalize angles to [0, 2π)
    auto normalize = [](double angle) {
        angle = std::fmod(angle, 2.0 * constants::PI);
        if (angle < 0.0) angle += 2.0 * constants::PI;
        return angle;
    };
    
    Omega = normalize(Omega);
    omega = normalize(omega);
    M = normalize(M);
    
    return KeplerianElements(a, e, i, Omega, omega, M, mu);
}

// ============================================================================
// Jacobian Matrices
// ============================================================================

Matrix6d KeplerianElements::jacobian_to_cartesian() const {
    /**
     * Computes ∂(r,v)/∂(a,e,i,Ω,ω,M)
     * 
     * Strategy: Use finite differences for numerical Jacobian
     * This is more robust than analytical derivatives for all cases
     */
    
    Matrix6d J = Matrix6d::Zero();
    
    // Perturbation sizes (adaptive based on element magnitudes)
    double da = std::max(a_ * 1e-8, 1e-3);       // km
    double de = std::max(e_ * 1e-8, 1e-10);      // dimensionless
    double di = 1e-8;                             // rad (~0.2 arcsec)
    double dOmega = 1e-8;                         // rad
    double domega = 1e-8;                         // rad
    double dM = 1e-8;                             // rad
    
    std::array<double, 6> perturbations = {da, de, di, dOmega, domega, dM};
    
    // Nominal state
    CartesianState nominal = this->to_cartesian();
    Vector6d state_nominal;
    state_nominal << nominal.position(), nominal.velocity();
    
    // Perturb each element and compute derivative
    for (int i = 0; i < 6; ++i) {
        // Create perturbed elements
        KeplerianElements perturbed = *this;
        double* elem_ptr = nullptr;
        
        switch (i) {
            case 0: elem_ptr = &perturbed.a_; break;
            case 1: elem_ptr = &perturbed.e_; break;
            case 2: elem_ptr = &perturbed.i_; break;
            case 3: elem_ptr = &perturbed.Omega_; break;
            case 4: elem_ptr = &perturbed.omega_; break;
            case 5: elem_ptr = &perturbed.M_; break;
        }
        
        *elem_ptr += perturbations[i];
        
        // Compute perturbed state
        CartesianState state_pert = perturbed.to_cartesian();
        Vector6d state_perturbed;
        state_perturbed << state_pert.position(), state_pert.velocity();
        
        // Finite difference
        J.col(i) = (state_perturbed - state_nominal) / perturbations[i];
    }
    
    return J;
}

Matrix6d KeplerianElements::jacobian_from_cartesian(const CartesianState& state) {
    /**
     * Computes ∂(a,e,i,Ω,ω,M)/∂(r,v)
     * 
     * Strategy: Use finite differences for numerical Jacobian
     */
    
    Matrix6d J = Matrix6d::Zero();
    
    // Perturbation sizes
    double dr = std::max(state.radius() * 1e-8, 1e-3);   // km
    double dv = std::max(state.speed() * 1e-8, 1e-6);    // km/s
    
    // Nominal elements
    KeplerianElements nominal = KeplerianElements::from_cartesian(state);
    Vector6d elem_nominal;
    elem_nominal << nominal.a_, nominal.e_, nominal.i_, 
                    nominal.Omega_, nominal.omega_, nominal.M_;
    
    // Perturb each Cartesian component
    for (int i = 0; i < 6; ++i) {
        CartesianState perturbed = state;
        double perturbation = (i < 3) ? dr : dv;
        
        // Modify the appropriate component
        if (i < 3) {
            Vector3d pos = perturbed.position();
            pos[i] += perturbation;
            perturbed = CartesianState(pos, perturbed.velocity(), perturbed.mu());
        } else {
            Vector3d vel = perturbed.velocity();
            vel[i-3] += perturbation;
            perturbed = CartesianState(perturbed.position(), vel, perturbed.mu());
        }
        
        // Compute perturbed elements
        KeplerianElements elem_pert = KeplerianElements::from_cartesian(perturbed);
        Vector6d elem_perturbed;
        elem_perturbed << elem_pert.a_, elem_pert.e_, elem_pert.i_,
                         elem_pert.Omega_, elem_pert.omega_, elem_pert.M_;
        
        // Handle angle wraparound for Omega, omega, M
        for (int j = 3; j < 6; ++j) {
            double diff = elem_perturbed[j] - elem_nominal[j];
            // Unwrap angles
            if (diff > constants::PI) diff -= 2.0 * constants::PI;
            if (diff < -constants::PI) diff += 2.0 * constants::PI;
            elem_perturbed[j] = elem_nominal[j] + diff;
        }
        
        // Finite difference
        J.col(i) = (elem_perturbed - elem_nominal) / perturbation;
    }
    
    return J;
}

} // namespace coordinates
} // namespace astdyn
