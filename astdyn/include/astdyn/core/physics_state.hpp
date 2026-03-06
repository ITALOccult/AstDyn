#ifndef ASTDYN_CORE_PHYSICS_STATE_HPP
#define ASTDYN_CORE_PHYSICS_STATE_HPP

#include "astdyn/math/frame_algebra.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/astrometry/sky_types.hpp"

namespace astdyn::physics {

/**
 * @brief Strongly-typed Cartesian State Vector
 * 
 * Enforces correct frame and units for position (Distance),
 * velocity (Velocity), and gravitational parameter.
 * Used for interfacing with ephemeris, propagators, etc.
 */
template <typename Frame>
struct CartesianStateTyped {
    time::EpochTDB epoch;
    math::Vector3<Frame, Distance> position;
    math::Vector3<Frame, Velocity> velocity;
    GravitationalParameter gm = GravitationalParameter::sun();

    CartesianStateTyped() = default;

    CartesianStateTyped(time::EpochTDB t, 
                        math::Vector3<Frame, Distance> p, 
                        math::Vector3<Frame, Velocity> v, 
                        GravitationalParameter g)
        : epoch(t), position(p), velocity(v), gm(g) {}

    // Factory method for creating from SI components
    [[nodiscard]] static CartesianStateTyped from_si(
        time::EpochTDB t, 
        double rx, double ry, double rz,
        double vx, double vy, double vz,
        double gm_val = constants::GM_SUN * 1e9) 
    {
        return CartesianStateTyped(
            t,
            math::Vector3<Frame, Distance>::from_si(rx, ry, rz),
            math::Vector3<Frame, Velocity>::from_si(vx, vy, vz),
            GravitationalParameter::from_si(gm_val)
        );
    }
    
    // Convert to un-typed Eigen state vector (r, v) in SI
    [[nodiscard]] Eigen::Matrix<double, 6, 1> to_eigen_si() const {
        Eigen::Matrix<double, 6, 1> state;
        state << position.x_si(), position.y_si(), position.z_si(),
                 velocity.x_si(), velocity.y_si(), velocity.z_si();
        return state;
    }
};

/**
 * @brief Strongly-typed Keplerian Elements
 */
template <typename Frame>
struct KeplerianStateTyped {
    time::EpochTDB epoch;
    Distance a;               // semi-major axis
    double e;                 // eccentricity
    astrometry::Angle i;      // inclination
    astrometry::Angle node;   // longitude of ascending node
    astrometry::Angle omega;  // argument of periapsis
    astrometry::Angle M;      // mean anomaly
    GravitationalParameter gm = GravitationalParameter::sun();

    KeplerianStateTyped() = default;

    KeplerianStateTyped(time::EpochTDB t, Distance semi, double ecc, 
                        astrometry::Angle inc, astrometry::Angle raan, 
                        astrometry::Angle argp, astrometry::Angle ma, 
                        GravitationalParameter g)
        : epoch(t), a(semi), e(ecc), i(inc), node(raan), omega(argp), M(ma), gm(g) {}

    /**
     * @brief Factory for importing from traditional astronomical databases (MPC, AstDys, etc.)
     * 
     * Since databases use raw doubles (usually AU, degrees, MJD), this explicitly forces 
     * the caller to acknowledge the units. The Frame is dictated by the template parameter.
     * 
     * Example for AstDys (Ecliptic J2000):
     *   auto k = KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
     *       EpochTDB::from_mjd(mjd), a_au, e, i_deg, node_deg, omega_deg, M_deg);
     */
    [[nodiscard]] static KeplerianStateTyped from_traditional(
        time::EpochTDB t, 
        double a_au, double ecc, 
        double i_deg, double node_deg, double omega_deg, double M_deg,
        GravitationalParameter gm_val = GravitationalParameter::sun()) 
    {
        return KeplerianStateTyped(
            t,
            Distance::from_au(a_au),
            ecc,
            astrometry::Angle::from_deg(i_deg),
            astrometry::Angle::from_deg(node_deg),
            astrometry::Angle::from_deg(omega_deg),
            astrometry::Angle::from_deg(M_deg),
            gm_val
        );
    }

    [[nodiscard]] double period_days() const {
        double a_au = a.to_au();
        double mu_au = gm.to_au3_d2();
        if (a_au <= 0 || mu_au <= 0) return 0.0;
        return 2.0 * constants::PI * std::sqrt(a_au * a_au * a_au / mu_au);
    }
    
    [[nodiscard]] double perihelion_au() const { return a.to_au() * (1.0 - e); }
    [[nodiscard]] double aphelion_au() const { return a.to_au() * (1.0 + e); }
};

} // namespace astdyn::physics

#endif // ASTDYN_CORE_PHYSICS_STATE_HPP
