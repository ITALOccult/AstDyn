/**
 * @file ReferenceFrame.hpp
 * @brief Reference frame definitions and transformations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Provides transformation matrices between various astronomical reference frames:
 * - J2000.0 (FK5): Mean equator and equinox at J2000.0
 * - ICRS: International Celestial Reference System
 * - Ecliptic: Mean ecliptic and equinox at J2000.0
 * - ITRF: International Terrestrial Reference Frame (Earth-fixed)
 * 
 * Reference: 
 * - IERS Conventions (2010)
 * - Seidelmann, "Explanatory Supplement to the Astronomical Almanac"
 * - Vallado, "Fundamentals of Astrodynamics and Applications"
 */

#ifndef ASTDYN_COORDINATES_REFERENCEFRAME_HPP
#define ASTDYN_COORDINATES_REFERENCEFRAME_HPP

#include "astdyn/core/Types.hpp"
#include "astdyn/core/Constants.hpp"
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/math/frame_algebra.hpp"
#include "astdyn/core/physics_types.hpp"
#include "astdyn/time/epoch.hpp"
#include "astdyn/core/frame_tags.hpp"
#include <cmath>
#include <string>

namespace astdyn {
namespace coordinates {

/**
 * @brief Enumeration of supported reference frames
 */
enum class FrameType {
    J2000,      ///< Mean equator and equinox at J2000.0 (FK5)
    ICRS,       ///< International Celestial Reference System
    ECLIPTIC,   ///< Mean ecliptic and equinox at J2000.0
    ITRF,       ///< International Terrestrial Reference Frame (Earth-fixed)
    MOD,        ///< Mean Of Date equator and equinox
    TOD         ///< True Of Date equator and equinox
};

/**
 * @brief Convert frame type to string
 */
inline std::string frame_type_to_string(FrameType type) {
    switch (type) {
        case FrameType::J2000: return "J2000";
        case FrameType::ICRS: return "ICRS";
        case FrameType::ECLIPTIC: return "ECLIPTIC";
        case FrameType::ITRF: return "ITRF";
        case FrameType::MOD: return "MOD";
        case FrameType::TOD: return "TOD";
        default: return "UNKNOWN";
    }
}

/**
 * @brief Reference frame transformation utilities
 * 
 * Provides static methods for transforming coordinates between
 * different astronomical reference frames.
 */
class ReferenceFrame {
public:
    // ========================================================================
    // Rotation Matrix Generators
    // ========================================================================
    
    /**
     * @brief Generate rotation matrix about X-axis
     * @param angle Rotation angle [rad]
     * @return 3x3 rotation matrix
     */
    static Matrix3d rotation_x(double angle) {
        double c = std::cos(angle);
        double s = std::sin(angle);
        
        Matrix3d R;
        R << 1.0, 0.0, 0.0,
             0.0,   c,   s,
             0.0,  -s,   c;
        return R;
    }
    
    /**
     * @brief Generate rotation matrix about Y-axis
     * @param angle Rotation angle [rad]
     * @return 3x3 rotation matrix
     */
    static Matrix3d rotation_y(double angle) {
        double c = std::cos(angle);
        double s = std::sin(angle);
        
        Matrix3d R;
        R <<   c, 0.0,  -s,
             0.0, 1.0, 0.0,
               s, 0.0,   c;
        return R;
    }
    
    /**
     * @brief Generate rotation matrix about Z-axis
     * @param angle Rotation angle [rad]
     * @return 3x3 rotation matrix
     */
    static Matrix3d rotation_z(double angle) {
        double c = std::cos(angle);
        double s = std::sin(angle);
        
        Matrix3d R;
        R <<   c,   s, 0.0,
              -s,   c, 0.0,
             0.0, 0.0, 1.0;
        return R;
    }
    
    // ========================================================================
    // J2000 ↔ ICRS Transformation
    // ========================================================================
    
    /**
     * @brief Get transformation matrix from J2000 to ICRS
     * 
     * The frame bias matrix accounts for the small offset between
     * the J2000 (FK5) and ICRS frames (~0.02 arcsec).
     * 
     * @return 3x3 rotation matrix
     */
    static Matrix3d j2000_to_icrs() {
        // Frame bias angles from IERS Conventions 2010
        // ξ₀ = -0.0166170 arcsec (x-axis offset)
        // η₀ = -0.0068192 arcsec (y-axis offset)  
        // da₀ = -0.014506 arcsec (z-axis offset)
        
        constexpr double xi0  = -0.0166170 * constants::ARCSEC_TO_RAD;
        constexpr double eta0 = -0.0068192 * constants::ARCSEC_TO_RAD;
        constexpr double da0  = -0.014506  * constants::ARCSEC_TO_RAD;
        
        // Frame bias matrix (simplified for small angles)
        Matrix3d bias;
        bias << 1.0,  da0,  -xi0,
               -da0,  1.0, -eta0,
                xi0, eta0,   1.0;
        
        return bias;
    }
    
    /**
     * @brief Get transformation matrix from ICRS to J2000
     * @return 3x3 rotation matrix (transpose of j2000_to_icrs)
     */
    static Matrix3d icrs_to_j2000() {
        return j2000_to_icrs().transpose();
    }
    
    // ========================================================================
    // J2000 ↔ Ecliptic Transformation
    // ========================================================================
    
    static double mean_obliquity(time::EpochTDB t) {
        // T is Julian centuries from J2000.0
        double T = (t.mjd() - constants::MJD2000) / constants::DAYS_PER_CENTURY;
        
        // IAU 2006 terms (in arcseconds)
        double eps_arcsec = 84381.406 
                          - 46.836769 * T 
                          - 0.0001831 * T * T 
                          + 0.00200340 * T * T * T 
                          - 0.000000576 * T * T * T * T 
                          - 0.0000000434 * T * T * T * T * T;
                          
        return eps_arcsec * constants::ARCSEC_TO_RAD;
    }

    /**
     * @brief Get transformation matrix from J2000 to Ecliptic
     * 
     * Rotation about X-axis by the mean obliquity of the ecliptic (ε)
     * 
     * @param t Epoch (defaults to J2000.0)
     * @return 3x3 rotation matrix
     */
    static Matrix3d j2000_to_ecliptic(time::EpochTDB t = time::EpochTDB::from_mjd(constants::MJD2000)) {
        return rotation_x(mean_obliquity(t));
    }
    
    /**
     * @brief Get transformation matrix from Ecliptic to J2000
     * @param t Epoch (defaults to J2000.0)
     * @return 3x3 rotation matrix
     */
    static Matrix3d ecliptic_to_j2000(time::EpochTDB t = time::EpochTDB::from_mjd(constants::MJD2000)) {
        return rotation_x(-mean_obliquity(t));
    }
    
    // ========================================================================
    // J2000 ↔ ITRF Transformation (simplified)
    // ========================================================================
    
    /**
     * @brief Compute Earth Rotation Angle (ERA)
     * @param t_utc Epoch in UTC
     * @param dut1 UT1 - UTC in seconds (from IERS data)
     * @return ERA in radians
     */
    static double computeERA(time::EpochUTC t_utc, double dut1 = 0.0) {
        double mjd_ut1 = t_utc.mjd() + dut1 / 86400.0;
        
        constexpr double era_0 = 0.7790572732640;
        constexpr double era_1 = 1.00273781191135448;
        
        double era = 2.0 * constants::PI * (era_0 + era_1 * (mjd_ut1 - 51544.5));
        
        // Normalize to [0, 2π)
        era = std::fmod(era, 2.0 * constants::PI);
        if (era < 0.0) era += 2.0 * constants::PI;
        
        return era;
    }

    /**
     * @brief Get transformation matrix from J2000 to ITRF (simplified)
     * 
     * This is a simplified transformation that only accounts for Earth rotation.
     * Full transformation requires:
     * - Precession (IAU 2006/2000A)
     * - Nutation (IAU 2000A)
     * - Earth rotation angle
     * - Polar motion
     * 
     * @param t_utc Epoch in UTC
     * @return 3x3 rotation matrix
     */
    static Matrix3d j2000_to_itrf_simple(time::EpochUTC t_utc) {
        double era = computeERA(t_utc, 0.0);
        
        // Rotation about Z-axis (negative for celestial to terrestrial)
        return rotation_z(-era);
    }
    
    /**
     * @brief Get transformation matrix from ITRF to J2000 (simplified)
     * @param mjd_ut1 Modified Julian Date in UT1 time scale
     * @return 3x3 rotation matrix
     */
    static Matrix3d itrf_to_j2000_simple(time::EpochUTC t_utc) {
        return j2000_to_itrf_simple(t_utc).transpose();
    }
    
    // ========================================================================
    // Generic Frame Transformation
    // ========================================================================
    
    /**
     * @brief Get transformation matrix between two frames
     * @param from Source frame
     * @param to Target frame
     * @param t Epoch (time::EpochTDB)
     * @return 3x3 rotation matrix
     */
    static Matrix3d get_transformation(FrameType from, FrameType to, 
                                       time::EpochTDB t = time::EpochTDB::from_mjd(constants::MJD2000)) {
        // Same frame, return identity
        if (from == to) {
            return Matrix3d::Identity();
        }
        
        // Direct transformations
        if (from == FrameType::J2000 && to == FrameType::ICRS) {
            return j2000_to_icrs();
        }
        if (from == FrameType::ICRS && to == FrameType::J2000) {
            return icrs_to_j2000();
        }
        if (from == FrameType::J2000 && to == FrameType::ECLIPTIC) {
            return j2000_to_ecliptic(t);
        }
        if (from == FrameType::ECLIPTIC && to == FrameType::J2000) {
            return ecliptic_to_j2000(t);
        }
        if (from == FrameType::J2000 && to == FrameType::ITRF) {
            return j2000_to_itrf_simple(time::EpochUTC::from_mjd(t.mjd()));
        }
        if (from == FrameType::ITRF && to == FrameType::J2000) {
            return itrf_to_j2000_simple(time::EpochUTC::from_mjd(t.mjd()));
        }
        
        // Chain transformations through J2000
        // Example: ICRS → ECLIPTIC = ICRS → J2000 → ECLIPTIC
        if (from == FrameType::ICRS && to == FrameType::ECLIPTIC) {
            return j2000_to_ecliptic(t) * icrs_to_j2000();
        }
        if (from == FrameType::ECLIPTIC && to == FrameType::ICRS) {
            return j2000_to_icrs() * ecliptic_to_j2000(t);
        }
        if (from == FrameType::ICRS && to == FrameType::ITRF) {
            return j2000_to_itrf_simple(time::EpochUTC::from_mjd(t.mjd())) * icrs_to_j2000();
        }
        if (from == FrameType::ITRF && to == FrameType::ICRS) {
            return j2000_to_icrs() * itrf_to_j2000_simple(time::EpochUTC::from_mjd(t.mjd()));
        }
        if (from == FrameType::ECLIPTIC && to == FrameType::ITRF) {
            return j2000_to_itrf_simple(time::EpochUTC::from_mjd(t.mjd())) * ecliptic_to_j2000(t);
        }
        if (from == FrameType::ITRF && to == FrameType::ECLIPTIC) {
            return j2000_to_ecliptic(t) * itrf_to_j2000_simple(time::EpochUTC::from_mjd(t.mjd()));
        }
        
        // Default: return identity (should not reach here)
        return Matrix3d::Identity();
    }
    
    // ========================================================================
    // State Vector Transformations
    // ========================================================================
    
    /**
     * @brief Transform position vector between frames
     * @param pos Position vector in source frame [km]
     * @param from Source frame
     * @param to Target frame
     * @param t Epoch (time::EpochTDB)
     * @return Position vector in target frame [km]
     */
    [[deprecated("Use typed transform_pos<FromFrame,ToFrame>() instead")]]
    static Vector3d transform_position(const Vector3d& pos, 
                                      FrameType from, FrameType to,
                                      time::EpochTDB t = time::EpochTDB::from_mjd(constants::MJD2000)) {
        Matrix3d R = get_transformation(from, to, t);
        return R * pos;
    }
    
    /**
     * @brief Transform velocity vector between frames
     * 
     * For rotating frames (ITRF), the velocity transformation includes
     * the Coriolis term: v_new = R·v_old + ω × (R·r_old)
     * 
     * @param pos Position vector in source frame [km]
     * @param vel Velocity vector in source frame [km/s]
     * @param from Source frame
     * @param to Target frame
     * @param t Epoch (time::EpochTDB)
     * @return Velocity vector in target frame [km/s]
     */
    [[deprecated("Use typed transform_vel<FromFrame,ToFrame>() instead")]]
    static Vector3d transform_velocity(const Vector3d& pos, const Vector3d& vel,
                                      FrameType from, FrameType to,
                                      time::EpochTDB t = time::EpochTDB::from_mjd(constants::MJD2000)) {
        Matrix3d R = get_transformation(from, to, t);
        Vector3d vel_rotated = R * vel;
        
        // Add Coriolis term for rotating frame transformations
        bool from_rotating = (from == FrameType::ITRF);
        bool to_rotating = (to == FrameType::ITRF);
        
        if (from_rotating != to_rotating) {
            // Earth rotation rate [rad/s]
            constexpr double omega_earth = 7.292115e-5;
            Vector3d omega(0.0, 0.0, omega_earth);
            
            Vector3d pos_rotated = R * pos;
            
            if (to_rotating) {
                // Inertial to rotating: subtract Coriolis
                vel_rotated -= omega.cross(pos_rotated);
            } else {
                // Rotating to inertial: add Coriolis
                vel_rotated += omega.cross(pos_rotated);
            }
        }
        
        return vel_rotated;
    }
    
    /**
     * @brief Transform Cartesian state between frames
     * @param state Cartesian state in source frame
     * @param from Source frame
     * @param to Target frame
     * @param t Epoch (time::EpochTDB)
     * @return Cartesian state in target frame
     */
    [[deprecated("Use typed transform_pos/transform_vel<FromFrame,ToFrame>() instead")]]
    static CartesianState transform_state(const CartesianState& state,
                                         FrameType from, FrameType to,
                                         time::EpochTDB t = time::EpochTDB::from_mjd(constants::MJD2000)) {
        Vector3d pos_new = transform_position(state.position(), from, to, t);
        Vector3d vel_new = transform_velocity(state.position(), state.velocity(),
                                              from, to, t);
        
        return CartesianState(pos_new, vel_new, state.mu());
    }
    
    // ========================================================================
    // Utility Functions
    // ========================================================================
    
    /**
     * @brief Check if a frame is inertial (non-rotating)
     */
    static bool is_inertial(FrameType frame) {
        return frame != FrameType::ITRF;
    }
    
    /**
     * @brief Check if a frame is rotating (Earth-fixed)
     */
    static bool is_rotating(FrameType frame) {
        return frame == FrameType::ITRF;
    }
    
    /**
     * @brief Get Greenwich Mean Sidereal Time (GMST)
     * @param mjd_ut1 Modified Julian Date in UT1
     * @return GMST [rad]
     */
    static double gmst(time::EpochUTC t_utc) {
        double mjd_ut1 = t_utc.mjd();
        double T = (mjd_ut1 - constants::MJD2000) / 36525.0;
        
        // GMST at 0h UT1 (IAU 2000)
        double gmst0 = 24110.54841 + 8640184.812866 * T 
                     + 0.093104 * T * T - 6.2e-6 * T * T * T;
        
        // Convert to days
        gmst0 /= constants::SECONDS_PER_DAY;
        
        // Add fraction of day
        double frac_day = mjd_ut1 - std::floor(mjd_ut1);
        double gmst_days = gmst0 + frac_day * 1.00273790935;
        
        // Convert to radians
        double gmst_rad = 2.0 * constants::PI * (gmst_days - std::floor(gmst_days));
        
        return gmst_rad;
    }
    // ========================================================================
    // Type-Safe Frame Transformations (AstDyn 3.0)
    // ========================================================================
    // These template methods encode the frame in the types, so the compiler
    // prevents mixing frames at compile-time.
    //
    // Usage:
    //   auto R = ReferenceFrame::get_rotation<core::ECLIPJ2000, core::GCRF>();
    //   auto pos_gcrf = R * pos_ecliptic;  // compile-time frame check!
    //
    //   auto pos_gcrf = ReferenceFrame::transform_pos<core::ECLIPJ2000, core::GCRF>(pos_ecl);
    // ========================================================================

    // --- Rotation matrix factories (compile-time frame pair) ---

    /**
     * @brief Get typed rotation matrix between two frame tags.
     * 
     * Specializations exist for all supported pairs.
     * ITRF transformations require a time parameter.
     */
    template <typename FromFrame, typename ToFrame>
    static math::RotationMatrix<FromFrame, ToFrame> get_rotation(
        time::EpochTDB t = time::EpochTDB::from_mjd(constants::MJD2000));

    // --- Position transformation ---

    /**
     * @brief Transform a typed position vector between frames.
     * Returns a Vector3 tagged with ToFrame.
     */
    template <typename FromFrame, typename ToFrame, typename PhysUnit>
    static math::Vector3<ToFrame, PhysUnit> transform_pos(
        const math::Vector3<FromFrame, PhysUnit>& pos,
        time::EpochTDB t = time::EpochTDB::from_mjd(constants::MJD2000))
    {
        auto R = get_rotation<FromFrame, ToFrame>(t);
        return R * pos;
    }

    // --- Velocity transformation (with Coriolis for rotating frames) ---

    /**
     * @brief Transform typed velocity between frames. 
     * Handles Coriolis term for ITRF.
     */
    template <typename FromFrame, typename ToFrame>
    static math::Vector3<ToFrame, physics::Velocity> transform_vel(
        const math::Vector3<FromFrame, physics::Distance>& pos,
        const math::Vector3<FromFrame, physics::Velocity>& vel,
        time::EpochTDB t = time::EpochTDB::from_mjd(constants::MJD2000))
    {
        auto R = get_rotation<FromFrame, ToFrame>(t);
        return R * vel; // Default: pure rotation (inertial-to-inertial)
    }
};

// ============================================================================
// Template specializations for get_rotation
// ============================================================================

// --- Identity (same frame) ---
template<> inline
math::RotationMatrix<core::GCRF, core::GCRF> ReferenceFrame::get_rotation<core::GCRF, core::GCRF>(time::EpochTDB) {
    return math::RotationMatrix<core::GCRF, core::GCRF>::from_eigen(Matrix3d::Identity());
}
template<> inline
math::RotationMatrix<core::ECLIPJ2000, core::ECLIPJ2000> ReferenceFrame::get_rotation<core::ECLIPJ2000, core::ECLIPJ2000>(time::EpochTDB) {
    return math::RotationMatrix<core::ECLIPJ2000, core::ECLIPJ2000>::from_eigen(Matrix3d::Identity());
}
template<> inline
math::RotationMatrix<core::ITRF, core::ITRF> ReferenceFrame::get_rotation<core::ITRF, core::ITRF>(time::EpochTDB) {
    return math::RotationMatrix<core::ITRF, core::ITRF>::from_eigen(Matrix3d::Identity());
}

// --- ECLIPJ2000 ↔ GCRF ---
template<> inline
math::RotationMatrix<core::ECLIPJ2000, core::GCRF> ReferenceFrame::get_rotation<core::ECLIPJ2000, core::GCRF>(time::EpochTDB) {
    // ECLIPJ2000 is the Mean Ecliptic of J2000.0. It MUST NOT depend on the current epoch t.
    return math::RotationMatrix<core::ECLIPJ2000, core::GCRF>::from_eigen(
        j2000_to_icrs() * ecliptic_to_j2000()); // Uses default J2000 epoch inside
}
template<> inline
math::RotationMatrix<core::GCRF, core::ECLIPJ2000> ReferenceFrame::get_rotation<core::GCRF, core::ECLIPJ2000>(time::EpochTDB) {
    return math::RotationMatrix<core::GCRF, core::ECLIPJ2000>::from_eigen(
        j2000_to_ecliptic() * icrs_to_j2000()); // Uses default J2000 epoch inside
}

// --- GCRF ↔ ITRF ---
template<> inline
math::RotationMatrix<core::GCRF, core::ITRF> ReferenceFrame::get_rotation<core::GCRF, core::ITRF>(time::EpochTDB t) {
    return math::RotationMatrix<core::GCRF, core::ITRF>::from_eigen(
        j2000_to_itrf_simple(time::EpochUTC::from_mjd(t.mjd())) * icrs_to_j2000());
}
template<> inline
math::RotationMatrix<core::ITRF, core::GCRF> ReferenceFrame::get_rotation<core::ITRF, core::GCRF>(time::EpochTDB t) {
    return math::RotationMatrix<core::ITRF, core::GCRF>::from_eigen(
        j2000_to_icrs() * itrf_to_j2000_simple(time::EpochUTC::from_mjd(t.mjd())));
}

// --- ECLIPJ2000 ↔ ITRF (chained via GCRF) ---
template<> inline
math::RotationMatrix<core::ECLIPJ2000, core::ITRF> ReferenceFrame::get_rotation<core::ECLIPJ2000, core::ITRF>(time::EpochTDB t) {
    return math::RotationMatrix<core::ECLIPJ2000, core::ITRF>::from_eigen(
        j2000_to_itrf_simple(time::EpochUTC::from_mjd(t.mjd())) * ecliptic_to_j2000(t));
}
template<> inline
math::RotationMatrix<core::ITRF, core::ECLIPJ2000> ReferenceFrame::get_rotation<core::ITRF, core::ECLIPJ2000>(time::EpochTDB t) {
    return math::RotationMatrix<core::ITRF, core::ECLIPJ2000>::from_eigen(
        j2000_to_ecliptic(t) * itrf_to_j2000_simple(time::EpochUTC::from_mjd(t.mjd())));
}

// ============================================================================
// Velocity specializations for ITRF (adds Coriolis term)
// ============================================================================

template<> inline
math::Vector3<core::ITRF, physics::Velocity> ReferenceFrame::transform_vel<core::GCRF, core::ITRF>(
    const math::Vector3<core::GCRF, physics::Distance>& pos,
    const math::Vector3<core::GCRF, physics::Velocity>& vel,
    time::EpochTDB t)
{
    auto R = get_rotation<core::GCRF, core::ITRF>(t);
    auto vel_rotated = R * vel;
    auto pos_rotated = R * pos;
    
    // ω_earth cross r (in ITRF)
    constexpr double omega_earth = 7.292115e-5; // rad/s
    // ω = (0, 0, ω_earth) cross pos_rotated
    double wx = -omega_earth * pos_rotated.y_si();
    double wy =  omega_earth * pos_rotated.x_si();
    
    return vel_rotated - math::Vector3<core::ITRF, physics::Velocity>::from_si(wx, wy, 0.0);
}

template<> inline
math::Vector3<core::GCRF, physics::Velocity> ReferenceFrame::transform_vel<core::ITRF, core::GCRF>(
    const math::Vector3<core::ITRF, physics::Distance>& pos,
    const math::Vector3<core::ITRF, physics::Velocity>& vel,
    time::EpochTDB t)
{
    auto R = get_rotation<core::ITRF, core::GCRF>(t);
    auto vel_rotated = R * vel;
    auto pos_rotated = R * pos;
    
    // ω_earth cross r (in GCRF)
    constexpr double omega_earth = 7.292115e-5;
    double wx = -omega_earth * pos_rotated.y_si();
    double wy =  omega_earth * pos_rotated.x_si();
    
    return vel_rotated + math::Vector3<core::GCRF, physics::Velocity>::from_si(wx, wy, 0.0);
}

} // namespace coordinates
} // namespace astdyn

#endif // ASTDYN_COORDINATES_REFERENCEFRAME_HPP
