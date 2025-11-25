/**
 * @file coordinates_example.cpp
 * @brief Demonstrate coordinate systems and reference frame transformations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * This example shows:
 * 1. Creating orbital states in different coordinate representations
 * 2. Converting between Cartesian, Keplerian, Equinoctial, and Cometary elements
 * 3. Transforming between reference frames (J2000, ICRS, Ecliptic, ITRF)
 * 4. Real-world example: Comet 1P/Halley orbit
 */

#include <iostream>
#include <iomanip>
#include "astdyn/coordinates/CartesianState.hpp"
#include "astdyn/coordinates/KeplerianElements.hpp"
#include "astdyn/coordinates/EquinoctialElements.hpp"
#include "astdyn/coordinates/CometaryElements.hpp"
#include "astdyn/coordinates/ReferenceFrame.hpp"
#include "astdyn/core/Constants.hpp"

using namespace astdyn;
using namespace astdyn::coordinates;
using namespace astdyn::constants;

// Helper function for formatted output
void print_separator(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(70, '=') << "\n";
}

void print_vector(const std::string& label, const Vector3d& vec, const std::string& unit = "") {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  " << std::left << std::setw(20) << label << ": ["
              << std::setw(15) << vec.x() << ", "
              << std::setw(15) << vec.y() << ", "
              << std::setw(15) << vec.z() << "]";
    if (!unit.empty()) std::cout << " " << unit;
    std::cout << "\n";
}

/**
 * Example 1: ISS Low Earth Orbit
 * Demonstrate conversions between coordinate representations
 */
void example_iss_orbit() {
    print_separator("Example 1: ISS Orbit - Coordinate Conversions");
    
    // ISS typical orbital parameters (approximate)
    double a = 6778.0;              // Semi-major axis [km]
    double e = 0.0005;              // Nearly circular
    double i = 51.6 * DEG_TO_RAD;   // Inclination [rad]
    double Omega = 45.0 * DEG_TO_RAD;
    double omega = 30.0 * DEG_TO_RAD;
    double M = 90.0 * DEG_TO_RAD;
    
    std::cout << "\n1. Starting with Keplerian Elements:\n";
    KeplerianElements kep(a, e, i, Omega, omega, M, GM_EARTH);
    std::cout << "  a = " << a << " km\n";
    std::cout << "  e = " << e << "\n";
    std::cout << "  i = " << i * RAD_TO_DEG << " deg\n";
    std::cout << "  Period = " << kep.period() / 60.0 << " minutes\n";
    
    // Convert to Cartesian
    std::cout << "\n2. Convert to Cartesian State:\n";
    CartesianState cart = kep.to_cartesian();
    print_vector("Position", cart.position(), "km");
    print_vector("Velocity", cart.velocity(), "km/s");
    std::cout << "  Altitude = " << (cart.radius() - 6371.0) << " km\n";
    std::cout << "  Speed = " << cart.speed() << " km/s\n";
    
    // Convert to Equinoctial (non-singular for low e)
    std::cout << "\n3. Convert to Equinoctial Elements (non-singular):\n";
    EquinoctialElements eq = EquinoctialElements::from_cartesian(cart);
    std::cout << "  a = " << eq.semi_major_axis() << " km\n";
    std::cout << "  h = " << eq.h() << " (e·sin(ω+Ω))\n";
    std::cout << "  k = " << eq.k() << " (e·cos(ω+Ω))\n";
    std::cout << "  λ = " << eq.mean_longitude() * RAD_TO_DEG << " deg\n";
    
    // Round-trip verification
    std::cout << "\n4. Round-trip verification (Kep → Cart → Eq → Kep):\n";
    KeplerianElements kep_final = eq.to_keplerian();
    std::cout << "  Δa = " << std::abs(kep_final.semi_major_axis() - a) << " km\n";
    std::cout << "  Δe = " << std::abs(kep_final.eccentricity() - e) << "\n";
    std::cout << "  Δi = " << std::abs(kep_final.inclination() - i) * RAD_TO_DEG << " deg\n";
}

/**
 * Example 2: Comet 1P/Halley
 * Demonstrate cometary elements for high eccentricity orbit
 */
void example_halley_comet() {
    print_separator("Example 2: Comet 1P/Halley - Cometary Elements");
    
    // Halley's orbital elements (epoch 1986 perihelion)
    double q = 0.5871 * AU;         // Perihelion distance
    double e = 0.9673;              // High eccentricity
    double i = 162.26 * DEG_TO_RAD; // Retrograde orbit
    double Omega = 58.42 * DEG_TO_RAD;
    double omega = 111.33 * DEG_TO_RAD;
    double T = 2446470.0;           // JD of perihelion passage (Feb 9, 1986)
    
    std::cout << "\n1. Cometary Elements (perihelion-based):\n";
    CometaryElements halley(q, e, i, Omega, omega, T, GM_SUN);
    std::cout << "  q = " << q / AU << " AU (perihelion)\n";
    std::cout << "  e = " << e << "\n";
    std::cout << "  i = " << i * RAD_TO_DEG << " deg (retrograde!)\n";
    std::cout << "  a = " << halley.semi_major_axis() / AU << " AU\n";
    std::cout << "  Q = " << halley.aphelion_distance() / AU << " AU (aphelion)\n";
    std::cout << "  Period = " << halley.period() / YEAR << " years\n";
    
    // State at perihelion
    std::cout << "\n2. State at perihelion (T = " << T << " JD):\n";
    CartesianState state_peri = halley.to_cartesian(T);
    std::cout << "  Distance = " << state_peri.radius() / AU << " AU\n";
    std::cout << "  Speed = " << state_peri.speed() << " km/s\n";
    
    // State 1 year after perihelion
    std::cout << "\n3. State 1 year after perihelion:\n";
    CartesianState state_1yr = halley.to_cartesian(T + 365.25);
    std::cout << "  Distance = " << state_1yr.radius() / AU << " AU\n";
    std::cout << "  Speed = " << state_1yr.speed() << " km/s\n";
    
    // Convert to standard Keplerian
    std::cout << "\n4. Convert to Keplerian at perihelion:\n";
    KeplerianElements kep_halley = halley.to_keplerian(T);
    std::cout << "  Mean anomaly = " << kep_halley.mean_anomaly() * RAD_TO_DEG 
              << " deg (should be ~0 at perihelion)\n";
}

/**
 * Example 3: Reference Frame Transformations
 * Demonstrate transformations between J2000, ICRS, Ecliptic, ITRF
 */
void example_reference_frames() {
    print_separator("Example 3: Reference Frame Transformations");
    
    // Earth satellite at specific epoch
    Vector3d pos_j2000(7000.0, 0.0, 0.0);  // km
    Vector3d vel_j2000(0.0, 7.5, 0.0);     // km/s
    CartesianState state_j2000(pos_j2000, vel_j2000, GM_EARTH);
    
    double mjd = MJD2000 + 365.25;  // One year after J2000
    
    std::cout << "\n1. Original state in J2000 frame:\n";
    print_vector("Position", state_j2000.position(), "km");
    print_vector("Velocity", state_j2000.velocity(), "km/s");
    
    // Transform to ICRS (small change)
    std::cout << "\n2. Transform to ICRS (frame bias ~0.02 arcsec):\n";
    CartesianState state_icrs = ReferenceFrame::transform_state(
        state_j2000, FrameType::J2000, FrameType::ICRS);
    print_vector("Position", state_icrs.position(), "km");
    Vector3d pos_diff = state_icrs.position() - state_j2000.position();
    std::cout << "  Difference = " << pos_diff.norm() * 1000.0 << " meters\n";
    
    // Transform to Ecliptic
    std::cout << "\n3. Transform to Ecliptic frame:\n";
    CartesianState state_ecl = ReferenceFrame::transform_state(
        state_j2000, FrameType::J2000, FrameType::ECLIPTIC);
    print_vector("Position", state_ecl.position(), "km");
    std::cout << "  (Note Z-component changed due to obliquity rotation)\n";
    
    // Transform to ITRF (Earth-fixed)
    std::cout << "\n4. Transform to ITRF (Earth-fixed) at MJD = " << mjd << ":\n";
    CartesianState state_itrf = ReferenceFrame::transform_state(
        state_j2000, FrameType::J2000, FrameType::ITRF, mjd);
    print_vector("Position", state_itrf.position(), "km");
    print_vector("Velocity", state_itrf.velocity(), "km/s");
    std::cout << "  (Velocity includes Coriolis effect)\n";
    
    // GMST at this epoch
    double gmst = ReferenceFrame::gmst(mjd);
    std::cout << "\n5. GMST = " << gmst * RAD_TO_DEG << " deg = " 
              << (gmst * RAD_TO_DEG / 15.0) << " hours\n";
    
    // Round-trip verification
    std::cout << "\n6. Round-trip verification (J2000 → ECLIPTIC → J2000):\n";
    CartesianState state_final = ReferenceFrame::transform_state(
        state_ecl, FrameType::ECLIPTIC, FrameType::J2000);
    Vector3d roundtrip_error = state_final.position() - state_j2000.position();
    std::cout << "  Position error = " << roundtrip_error.norm() * 1e6 << " microns\n";
}

/**
 * Example 4: Ground Station Visibility
 * Practical application: when is a satellite visible from a ground station?
 */
void example_ground_station_visibility() {
    print_separator("Example 4: Ground Station Visibility Calculation");
    
    // Satellite in J2000
    Vector3d sat_pos_j2000(7000.0, 0.0, 0.0);
    Vector3d sat_vel_j2000(0.0, 7.5, 0.0);
    CartesianState sat_j2000(sat_pos_j2000, sat_vel_j2000, GM_EARTH);
    
    // Ground station: Greenwich Observatory (simplified)
    // Latitude = 51.48° N, Longitude = 0° E
    double lat = 51.48 * DEG_TO_RAD;
    double lon = 0.0;
    double alt = 0.047;  // km above sea level
    double R_earth = 6371.0;  // km
    
    std::cout << "\n1. Ground Station: Royal Observatory Greenwich\n";
    std::cout << "  Latitude = " << lat * RAD_TO_DEG << " deg N\n";
    std::cout << "  Longitude = " << lon * RAD_TO_DEG << " deg E\n";
    
    double mjd = MJD2000;
    double gmst = ReferenceFrame::gmst(mjd);
    
    std::cout << "\n2. At epoch MJD = " << mjd << ":\n";
    std::cout << "  GMST = " << (gmst * RAD_TO_DEG / 15.0) << " hours\n";
    
    // Transform satellite to ITRF
    CartesianState sat_itrf = ReferenceFrame::transform_state(
        sat_j2000, FrameType::J2000, FrameType::ITRF, mjd);
    
    std::cout << "\n3. Satellite in Earth-fixed frame (ITRF):\n";
    print_vector("Position", sat_itrf.position(), "km");
    
    // Ground station position in ITRF
    double local_sidereal_time = gmst + lon;
    Vector3d station_itrf;
    station_itrf.x() = (R_earth + alt) * std::cos(lat) * std::cos(local_sidereal_time);
    station_itrf.y() = (R_earth + alt) * std::cos(lat) * std::sin(local_sidereal_time);
    station_itrf.z() = (R_earth + alt) * std::sin(lat);
    
    std::cout << "\n4. Ground station in ITRF:\n";
    print_vector("Position", station_itrf, "km");
    
    // Relative position
    Vector3d relative = sat_itrf.position() - station_itrf;
    double range = relative.norm();
    double elevation_angle = std::asin(relative.dot(station_itrf) / 
                                       (range * station_itrf.norm()));
    
    std::cout << "\n5. Visibility:\n";
    std::cout << "  Range = " << range << " km\n";
    std::cout << "  Elevation angle = " << elevation_angle * RAD_TO_DEG << " deg\n";
    if (elevation_angle > 0) {
        std::cout << "  Status: VISIBLE above horizon\n";
    } else {
        std::cout << "  Status: BELOW horizon (not visible)\n";
    }
}

/**
 * Example 5: Orbit Classification
 */
void example_orbit_classification() {
    print_separator("Example 5: Orbit Classification");
    
    std::cout << "\nClassifying different orbit types:\n\n";
    
    // 1. Circular LEO
    KeplerianElements leo(6778.0, 0.0001, 51.6*DEG_TO_RAD, 0, 0, 0, GM_EARTH);
    std::cout << "1. LEO (e=" << leo.eccentricity() << "):\n";
    std::cout << "   Circular: " << (leo.is_circular() ? "YES" : "NO") << "\n";
    std::cout << "   Elliptic: " << (leo.is_elliptic() ? "YES" : "NO") << "\n\n";
    
    // 2. Molniya orbit
    KeplerianElements molniya(26554.0, 0.74, 63.4*DEG_TO_RAD, 0, 0, 0, GM_EARTH);
    std::cout << "2. Molniya (e=" << molniya.eccentricity() << "):\n";
    std::cout << "   Circular: " << (molniya.is_circular() ? "YES" : "NO") << "\n";
    std::cout << "   Elliptic: " << (molniya.is_elliptic() ? "YES" : "NO") << "\n";
    std::cout << "   Period: " << molniya.period() / 3600.0 << " hours\n\n";
    
    // 3. Parabolic escape
    CometaryElements escape(7000.0, 1.0, 0, 0, 0, 0, GM_EARTH);
    std::cout << "3. Parabolic (e=" << escape.eccentricity() << "):\n";
    std::cout << "   Parabolic: " << (escape.is_parabolic() ? "YES" : "NO") << "\n";
    std::cout << "   Bound: " << (escape.is_bound() ? "YES" : "NO") << "\n\n";
    
    // 4. Hyperbolic flyby
    CometaryElements hyperbolic(7000.0, 1.5, 0, 0, 0, 0, GM_EARTH);
    std::cout << "4. Hyperbolic (e=" << hyperbolic.eccentricity() << "):\n";
    std::cout << "   Hyperbolic: " << (hyperbolic.is_hyperbolic() ? "YES" : "NO") << "\n";
    std::cout << "   Unbound: " << (hyperbolic.is_unbound() ? "YES" : "NO") << "\n";
}

int main() {
    std::cout << std::string(70, '=') << "\n";
    std::cout << "  AstDyn C++ Coordinates & Reference Frames Examples\n";
    std::cout << std::string(70, '=') << "\n";
    
    try {
        example_iss_orbit();
        example_halley_comet();
        example_reference_frames();
        example_ground_station_visibility();
        example_orbit_classification();
        
        print_separator("All Examples Completed Successfully!");
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
