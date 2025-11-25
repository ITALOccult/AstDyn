/**
 * @file ephemeris_example.cpp
 * @brief Demonstration of planetary ephemeris calculations
 * @author ITALOccult AstDyn Team
 * @date 2025-11-23
 * 
 * Examples:
 * 1. Query planetary positions at specific epochs
 * 2. Compute planetary distances and velocities
 * 3. Calculate barycentric corrections
 * 4. Demonstrate orbital periods validation
 */

#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/ephemeris/PlanetaryData.hpp>
#include <astdyn/core/Constants.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace astdyn;
using namespace astdyn::ephemeris;
using namespace astdyn::constants;

void printHeader(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(70, '=') << "\n";
}

void printVector(const std::string& label, const Eigen::Vector3d& vec, const std::string& unit) {
    std::cout << std::setw(20) << std::left << label
              << std::fixed << std::setprecision(6)
              << "[" << std::setw(12) << vec.x() << ", "
              << std::setw(12) << vec.y() << ", "
              << std::setw(12) << vec.z() << "] " << unit << "\n";
}

/**
 * Example 1: Solar System Snapshot at J2000.0
 * 
 * Computes and displays positions of all major planets at the J2000.0 epoch.
 */
void example1_solar_system_snapshot() {
    printHeader("Example 1: Solar System Snapshot at J2000.0");
    
    double jd = JD2000;
    std::cout << "Epoch: JD " << jd << " (January 1, 2000, 12:00 TT)\n\n";
    
    std::cout << std::setw(12) << std::left << "Body"
              << std::setw(15) << "Distance (AU)"
              << std::setw(20) << "Velocity (AU/day)"
              << std::setw(15) << "Period (days)\n";
    std::cout << std::string(70, '-') << "\n";
    
    CelestialBody planets[] = {
        CelestialBody::MERCURY, CelestialBody::VENUS, CelestialBody::EARTH,
        CelestialBody::MARS, CelestialBody::JUPITER, CelestialBody::SATURN,
        CelestialBody::URANUS, CelestialBody::NEPTUNE
    };
    
    for (auto planet : planets) {
        auto pos = PlanetaryEphemeris::getPosition(planet, jd);
        auto vel = PlanetaryEphemeris::getVelocity(planet, jd);
        auto data = PlanetaryData::getBodyData(planet);
        
        std::cout << std::setw(12) << std::left << data.name
                  << std::setw(15) << std::fixed << std::setprecision(6) << pos.norm()
                  << std::setw(20) << std::fixed << std::setprecision(6) << vel.norm()
                  << std::setw(15) << std::fixed << std::setprecision(1) << data.period << "\n";
    }
}

/**
 * Example 2: Earth-Mars Distance Over Time
 * 
 * Calculates the distance between Earth and Mars over one synodic period.
 * Identifies opposition (closest approach) and conjunction (farthest).
 */
void example2_earth_mars_distance() {
    printHeader("Example 2: Earth-Mars Distance Over Time");
    
    std::cout << "Computing Earth-Mars distance for one synodic period (~780 days)\n\n";
    
    double jd_start = JD2000;
    double min_dist = 1e10, max_dist = 0.0;
    double jd_min = jd_start, jd_max = jd_start;
    
    std::cout << std::setw(20) << std::left << "Date (JD)"
              << std::setw(20) << "Distance (AU)"
              << std::setw(20) << "Distance (Gm)\n";
    std::cout << std::string(60, '-') << "\n";
    
    for (int day = 0; day <= 780; day += 100) {
        double jd = jd_start + day;
        auto pos_earth = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd);
        auto pos_mars = PlanetaryEphemeris::getPosition(CelestialBody::MARS, jd);
        
        double dist_au = (pos_mars - pos_earth).norm();
        double dist_gm = dist_au * AU / 1e9; // Convert to Gm
        
        if (dist_au < min_dist) {
            min_dist = dist_au;
            jd_min = jd;
        }
        if (dist_au > max_dist) {
            max_dist = dist_au;
            jd_max = jd;
        }
        
        std::cout << std::setw(20) << std::fixed << std::setprecision(1) << jd
                  << std::setw(20) << std::fixed << std::setprecision(3) << dist_au
                  << std::setw(20) << std::fixed << std::setprecision(1) << dist_gm << "\n";
    }
    
    std::cout << "\nOpposition (closest): " << min_dist << " AU at JD " << jd_min << "\n";
    std::cout << "Conjunction (farthest): " << max_dist << " AU at JD " << jd_max << "\n";
}

/**
 * Example 3: Jupiter's Perturbation on Earth
 * 
 * Demonstrates the gravitational influence of Jupiter on Earth's orbit
 * by computing the acceleration due to Jupiter at different times.
 */
void example3_jupiter_perturbation() {
    printHeader("Example 3: Jupiter's Perturbation on Earth");
    
    std::cout << "Computing Jupiter's gravitational acceleration on Earth\n\n";
    
    double jd_start = JD2000;
    double gm_jupiter = PlanetaryData::GM_JUPITER; // km³/s²
    double au_km = AU; // AU to km
    
    std::cout << std::setw(15) << std::left << "Days from J2000"
              << std::setw(20) << "Jupiter Dist (AU)"
              << std::setw(25) << "Acceleration (km/s²)\n";
    std::cout << std::string(60, '-') << "\n";
    
    for (int day = 0; day <= 365 * 4; day += 365) {
        double jd = jd_start + day;
        auto pos_earth = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd);
        auto pos_jupiter = PlanetaryEphemeris::getPosition(CelestialBody::JUPITER, jd);
        
        Eigen::Vector3d delta = pos_jupiter - pos_earth;
        double dist_au = delta.norm();
        double dist_km = dist_au * au_km;
        
        // Gravitational acceleration: a = GM/r²
        double accel = gm_jupiter / (dist_km * dist_km);
        
        std::cout << std::setw(15) << day
                  << std::setw(20) << std::fixed << std::setprecision(3) << dist_au
                  << std::setw(25) << std::scientific << std::setprecision(6) << accel << "\n";
    }
    
    std::cout << "\nNote: This perturbation affects Earth's orbit determination at mm/s² level\n";
}

/**
 * Example 4: Barycentric Correction
 * 
 * Shows the difference between heliocentric and barycentric coordinates.
 * Important for high-precision astrometry and pulsar timing.
 */
void example4_barycentric_correction() {
    printHeader("Example 4: Solar System Barycentric Correction");
    
    std::cout << "Sun's position relative to Solar System Barycenter\n\n";
    
    double jd_start = JD2000;
    
    std::cout << std::setw(20) << std::left << "Date (JD)"
              << std::setw(25) << "SSB Offset (AU)"
              << std::setw(25) << "SSB Offset (km)\n";
    std::cout << std::string(70, '-') << "\n";
    
    for (int year = 0; year <= 12; year += 3) {
        double jd = jd_start + year * 365.25;
        auto r_sun_ssb = PlanetaryEphemeris::getSunBarycentricPosition(jd);
        
        double offset_au = r_sun_ssb.norm();
        double offset_km = offset_au * AU;
        
        std::cout << std::setw(20) << std::fixed << std::setprecision(1) << jd
                  << std::setw(25) << std::fixed << std::setprecision(6) << offset_au
                  << std::setw(25) << std::fixed << std::setprecision(0) << offset_km << "\n";
    }
    
    std::cout << "\nTypical offset: ~0.005 AU (~750,000 km), dominated by Jupiter\n";
    std::cout << "This is larger than the Sun's radius (695,700 km)!\n";
}

/**
 * Example 5: Planetary Configurations
 * 
 * Identifies interesting planetary alignments and configurations.
 */
void example5_planetary_configurations() {
    printHeader("Example 5: Planetary Configurations at J2000.0");
    
    double jd = JD2000;
    
    // Get all planetary positions
    auto pos_mercury = PlanetaryEphemeris::getPosition(CelestialBody::MERCURY, jd);
    auto pos_venus = PlanetaryEphemeris::getPosition(CelestialBody::VENUS, jd);
    auto pos_earth = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd);
    auto pos_mars = PlanetaryEphemeris::getPosition(CelestialBody::MARS, jd);
    auto pos_jupiter = PlanetaryEphemeris::getPosition(CelestialBody::JUPITER, jd);
    
    std::cout << "\nHeliocentric Ecliptic Coordinates (J2000.0):\n\n";
    
    printVector("Mercury:", pos_mercury, "AU");
    printVector("Venus:", pos_venus, "AU");
    printVector("Earth:", pos_earth, "AU");
    printVector("Mars:", pos_mars, "AU");
    printVector("Jupiter:", pos_jupiter, "AU");
    
    // Compute elongations as seen from Earth
    std::cout << "\nElongations (as seen from Earth):\n\n";
    
    auto compute_elongation = [&](const Eigen::Vector3d& pos_planet, const std::string& name) {
        Eigen::Vector3d earth_to_sun = -pos_earth;
        Eigen::Vector3d earth_to_planet = pos_planet - pos_earth;
        
        double cos_elong = earth_to_sun.dot(earth_to_planet) /
                          (earth_to_sun.norm() * earth_to_planet.norm());
        double elong_deg = std::acos(cos_elong) * RAD_TO_DEG;
        
        std::cout << std::setw(12) << std::left << name
                  << std::setw(20) << std::fixed << std::setprecision(2) << elong_deg << "°\n";
    };
    
    compute_elongation(pos_mercury, "Mercury");
    compute_elongation(pos_venus, "Venus");
    compute_elongation(pos_mars, "Mars");
    compute_elongation(pos_jupiter, "Jupiter");
    
    std::cout << "\nNote: Elongation is the angle between Sun and planet as seen from Earth\n";
    std::cout << "      0° = conjunction, 90° = quadrature, 180° = opposition\n";
}

int main() {
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         PLANETARY EPHEMERIS EXAMPLES - ASTDYN C++             ║\n";
    std::cout << "║                                                               ║\n";
    std::cout << "║  Demonstrates calculation of planetary positions/velocities   ║\n";
    std::cout << "║  using simplified VSOP87 analytical approximations            ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════╝\n";
    
    try {
        example1_solar_system_snapshot();
        example2_earth_mars_distance();
        example3_jupiter_perturbation();
        example4_barycentric_correction();
        example5_planetary_configurations();
        
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "All examples completed successfully!\n";
        std::cout << std::string(70, '=') << "\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
