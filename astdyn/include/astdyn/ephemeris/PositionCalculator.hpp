// PositionCalculator.hpp
// Utility class for fast computation of heliocentric position from orbital elements.
// Supports Keplerian elements expressed in either ecliptic or equatorial reference frames.
// The implementation is lightweight and header‑only for maximum performance.

#ifndef ASTDYN_EPHEMERIS_POSITIONCALCULATOR_HPP
#define ASTDYN_EPHEMERIS_POSITIONCALCULATOR_HPP

#include <Eigen/Dense>
#include <cmath>

namespace astdyn {
namespace ephemeris {

// Simple struct for Keplerian elements (mean values).
struct KeplerianElements {
    double a;      // semi‑major axis (AU)
    double e;      // eccentricity
    double i;      // inclination (rad)
    double Omega;  // longitude of ascending node (rad)
    double omega;  // argument of periapsis (rad)
    double M;      // mean anomaly at epoch (rad)
    // Optional: flag indicating whether the elements are expressed in the
    // equatorial J2000 frame. If false they are assumed to be ecliptic J2000.
    bool equatorial = false;
};

class PositionCalculator {
public:
    // Compute heliocentric position (AU) at a given MJD.
    // The function solves Kepler's equation iteratively (Newton‑Raphson) and
    // transforms the orbital plane coordinates to the requested reference frame.
    static Eigen::Vector3d computePosition(const KeplerianElements& elem,
                                           double mjd_target,
                                           double mjd_epoch,
                                           bool outputEquatorial = false) {
        // Time since epoch in days (not used for pure Keplerian propagation).
        // For a simple two‑body problem the mean motion n = sqrt(mu / a^3).
        // Using Gaussian gravitational constant k = 0.01720209895 AU^{3/2} / day.
        constexpr double k = 0.01720209895;
        double n = k * std::sqrt(1.0 / (elem.a * elem.a * elem.a)); // rad/day
        double dt = mjd_target - mjd_epoch;
        double M = elem.M + n * dt; // propagate mean anomaly
        // Normalize M to [0, 2π)
        M = std::fmod(M, 2.0 * M_PI);
        if (M < 0) M += 2.0 * M_PI;
        // Solve Kepler's equation: M = E - e sin E
        double E = solveKepler(M, elem.e);
        // True anomaly
        double nu = 2.0 * std::atan2(std::sqrt(1.0 + elem.e) * std::sin(E / 2.0),
                                      std::sqrt(1.0 - elem.e) * std::cos(E / 2.0));
        // Distance r
        double r = elem.a * (1.0 - elem.e * std::cos(E));
        // Position in orbital plane (PQW)
        double x_orb = r * std::cos(nu);
        double y_orb = r * std::sin(nu);
        // Rotate to inertial frame
        Eigen::Vector3d pos = orbitalToInertial(x_orb, y_orb, elem);
        // If the input was equatorial, we need to rotate back to ecliptic for a
        // heliocentric ecliptic coordinates (the library works internally in the
        // ecliptic frame). The obliquity of the ecliptic J2000 is ε = 23.43928°.
        if (elem.equatorial) {
            constexpr double epsilon = 23.43928 * M_PI / 180.0; // rad
            Eigen::Matrix3d R;
            R << 1, 0, 0,
                 0, std::cos(epsilon), std::sin(epsilon),
                 0, -std::sin(epsilon), std::cos(epsilon);
            pos = R * pos; // convert equatorial -> ecliptic
        }
        if (outputEquatorial) {
            // Convert back to equatorial frame
            constexpr double epsilon = 23.43928 * M_PI / 180.0; // rad
            Eigen::Matrix3d R;
            R << 1, 0, 0,
                 0, std::cos(epsilon), -std::sin(epsilon),
                 0, std::sin(epsilon),  std::cos(epsilon);
            pos = R * pos; // ecliptic -> equatorial
        }
        return pos; // AU, heliocentric
    }

    // Compute position directly from a Cartesian state vector (first three components are position).
    // The vector is assumed to be in heliocentric ecliptic coordinates.
    // If outputEquatorial is true, the result is transformed to the equatorial J2000 frame.
    static Eigen::Vector3d computePositionFromState(const Eigen::VectorXd& state,
                                                   bool outputEquatorial = false) {
        // Ensure the state has at least 3 components.
        assert(state.size() >= 3 && "State vector must contain at least 3 elements (x, y, z)");
        Eigen::Vector3d pos = state.head<3>();
        if (outputEquatorial) {
            constexpr double epsilon = 23.43928 * M_PI / 180.0; // rad
            Eigen::Matrix3d R;
            R << 1, 0, 0,
                 0, std::cos(epsilon), -std::sin(epsilon),
                 0, std::sin(epsilon),  std::cos(epsilon);
            pos = R * pos; // ecliptic -> equatorial
        }
        return pos; // AU, heliocentric (or equatorial if requested)

    }


    // Newton‑Raphson solver for Kepler's equation.
    static double solveKepler(double M, double e) {
        // Initial guess
        double E = (e < 0.8) ? M : M_PI;
        for (int iter = 0; iter < 10; ++iter) {
            double f = E - e * std::sin(E) - M;
            double fprime = 1.0 - e * std::cos(E);
            double delta = -f / fprime;
            E += delta;
            if (std::abs(delta) < 1e-12) break;
        }
        return E;
    }

    // Transform from orbital plane (PQW) to inertial (ecliptic) coordinates.
    static Eigen::Vector3d orbitalToInertial(double x_orb, double y_orb,
                                             const KeplerianElements& elem) {
        // Rotation matrices (3‑1‑3 sequence)
        double cosO = std::cos(elem.Omega);
        double sinO = std::sin(elem.Omega);
        double cosi = std::cos(elem.i);
        double sini = std::sin(elem.i);
        double cosw = std::cos(elem.omega);
        double sinw = std::sin(elem.omega);
        Eigen::Matrix3d R;
        R << cosO * cosw - sinO * sinw * cosi, -cosO * sinw - sinO * cosw * cosi, sinO * sini,
             sinO * cosw + cosO * sinw * cosi, -sinO * sinw + cosO * cosw * cosi, -cosO * sini,
             sinw * sini,                     cosw * sini,                      cosi;
        Eigen::Vector3d vec_orb(x_orb, y_orb, 0.0);
        return R * vec_orb;
    }
};

} // namespace ephemeris
} // namespace astdyn

#endif // ASTDYN_EPHEMERIS_POSITIONCALCULATOR_HPP
