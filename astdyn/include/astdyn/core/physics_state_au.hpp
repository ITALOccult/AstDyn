#ifndef ASTDYN_CORE_PHYSICS_STATE_AU_HPP
#define ASTDYN_CORE_PHYSICS_STATE_AU_HPP

#include "astdyn/core/physics_types.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include "src/core/frame_tags.hpp"

namespace astdyn::physics {

/**
 * @brief State derivative in AU units.
 * Separates velocity [AU/day] and acceleration [AU/day^2].
 */
struct DerivativeAU {
    VelocityAUD dx, dy, dz;
    AccelerationAUD2 dvx, dvy, dvz;

    // Vector space operations for RK/Adams integrators
    DerivativeAU operator+(const DerivativeAU& o) const {
        return {dx + o.dx, dy + o.dy, dz + o.dz, dvx + o.dvx, dvy + o.dvy, dvz + o.dvz};
    }
    DerivativeAU operator-(const DerivativeAU& o) const {
        return {dx - o.dx, dy - o.dy, dz - o.dz, dvx - o.dvx, dvy - o.dvy, dvz - o.dvz};
    }
    DerivativeAU operator*(double s) const {
        return {dx * s, dy * s, dz * s, dvx * s, dvy * s, dvz * s};
    }
    friend DerivativeAU operator*(double s, const DerivativeAU& d) {
        return d * s;
    }
};

/**
 * @brief High-precision Cartesian state in AU and AU/day.
 * Specifically designed for numerical integration and differential correction.
 * Forced to ECLIPJ2000 frame for mathematical consistency in the core.
 */
template <typename Frame = core::ECLIPJ2000>
struct CartesianStateAU {
    time::EpochTDB epoch;
    
    // Position [AU]
    double x, y, z;
    // Velocity [AU/day]
    double vx, vy, vz;
    
    GravitationalParameter gm;

    // Vector space operations: State + (Derivative * dt)
    CartesianStateAU operator+(const DerivativeAU& d) const {
        CartesianStateAU res = *this;
        res.x += d.dx.to_au_d(); // Logic: x is AU, dx is AU/d, but integrators use raw dt
        res.y += d.dy.to_au_d();
        res.z += d.dz.to_au_d();
        res.vx += d.dvx.to_au_d2();
        res.vy += d.dvy.to_au_d2();
        res.vz += d.dvz.to_au_d2();
        return res;
    }

    CartesianStateAU operator+(const CartesianStateAU& o) const {
        CartesianStateAU res = *this;
        res.x += o.x; res.y += o.y; res.z += o.z;
        res.vx += o.vx; res.vy += o.vy; res.vz += o.vz;
        return res;
    }
    CartesianStateAU operator-(const CartesianStateAU& o) const {
        CartesianStateAU res = *this;
        res.x -= o.x; res.y -= o.y; res.z -= o.z;
        res.vx -= o.vx; res.vy -= o.vy; res.vz -= o.vz;
        return res;
    }
    CartesianStateAU operator*(double s) const {
        CartesianStateAU res = *this;
        res.x *= s; res.y *= s; res.z *= s;
        res.vx *= s; res.vy *= s; res.vz *= s;
        return res;
    }
    friend CartesianStateAU operator*(double s, const CartesianStateAU& d) {
        return d * s;
    }

    /// Factory from standard SI state
    static CartesianStateAU from_si(const CartesianStateTyped<Frame>& si_state) {
        CartesianStateAU res;
        res.epoch = si_state.epoch;
        res.gm = si_state.gm;
        
        const double au_m = constants::AU * 1000.0;
        const double aud_ms = au_m / 86400.0;
        
        res.x = si_state.position.x_si() / au_m;
        res.y = si_state.position.y_si() / au_m;
        res.z = si_state.position.z_si() / au_m;
        
        res.vx = si_state.velocity.x_si() / aud_ms;
        res.vy = si_state.velocity.y_si() / aud_ms;
        res.vz = si_state.velocity.z_si() / aud_ms;
        
        return res;
    }

    /// Convert back to SI state
    CartesianStateTyped<Frame> to_si() const {
        const double au_m = constants::AU * 1000.0;
        const double aud_ms = au_m / 86400.0;
        return CartesianStateTyped<Frame>::from_si(
            epoch,
            x * au_m, y * au_m, z * au_m,
            vx * aud_ms, vy * aud_ms, vz * aud_ms,
            gm.to_m3_s2()
        );
    }
    
    /// Export to raw Eigen vector [x, y, z, vx, vy, vz] for integrators
    Eigen::VectorXd to_eigen() const {
        Eigen::VectorXd res(6);
        res << x, y, z, vx, vy, vz;
        return res;
    }
    
    /// Import from raw Eigen vector [AU, AU/day]
    static CartesianStateAU from_eigen(const Eigen::VectorXd& v, time::EpochTDB t, GravitationalParameter mu) {
        CartesianStateAU res;
        res.epoch = t;
        res.gm = mu;
        res.x = v[0]; res.y = v[1]; res.z = v[2];
        res.vx = v[3]; res.vy = v[4]; res.vz = v[5];
        return res;
    }
};

/**
 * @brief High-precision Covariance matrix in AU units.
 * Dimension: 6x6.
 * Blocks: Position [AU^2], Velocity [(AU/day)^2].
 */
template <typename Frame = core::ECLIPJ2000>
class CovarianceMatrixAU {
public:
    CovarianceMatrixAU() { matrix_.setZero(); }
    explicit CovarianceMatrixAU(const Matrix6d& mat) : matrix_(mat) {}
    
    const Matrix6d& matrix() const { return matrix_; }
    
    /// Combined Position uncertainty (RSS) in Kilometers
    double sigma_pos_km() const {
        return std::sqrt(matrix_.block<3,3>(0,0).trace()) * constants::AU;
    }
    
    /// Combined Velocity uncertainty (RSS) in Meters per Second
    double sigma_vel_ms() const {
        const double aud_ms = (constants::AU * 1000.0) / 86400.0;
        return std::sqrt(matrix_.block<3,3>(3,3).trace()) * aud_ms;
    }

private:
    Matrix6d matrix_;
};

} // namespace astdyn::physics

#endif // ASTDYN_CORE_PHYSICS_STATE_AU_HPP
