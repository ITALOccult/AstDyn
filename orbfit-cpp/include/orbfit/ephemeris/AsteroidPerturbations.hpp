/**
 * @file AsteroidPerturbations.hpp
 * @brief Perturbations from massive asteroids
 * 
 * Implements gravitational perturbations from the 16 most massive asteroids
 * using simplified orbital elements (AST17 model or similar).
 * 
 * Data sources:
 * - AST17: Baer et al. (2011) - 16 asteroids with GM determinations
 * - JPL Small-Body Database: https://ssd.jpl.nasa.gov/sbdb.cgi
 * - Carry (2012): Density measurements
 * 
 * Included asteroids (GM > 1 km³/s²):
 * 1 Ceres, 2 Pallas, 4 Vesta, 10 Hygiea, 15 Eunomia, 16 Psyche,
 * 31 Euphrosyne, 52 Europa, 65 Cybele, 87 Sylvia, 88 Thisbe,
 * 511 Davida, 704 Interamnia, 324 Bamberga, 451 Patientia, 107 Camilla
 */

#ifndef ORBFIT_ASTEROID_PERTURBATIONS_HPP
#define ORBFIT_ASTEROID_PERTURBATIONS_HPP

#include <orbfit/coordinates/CartesianState.hpp>
#include <orbfit/core/Constants.hpp>
#include <Eigen/Dense>
#include <vector>
#include <string>

namespace orbfit {
namespace ephemeris {

/**
 * @brief Asteroid data for perturbation calculations
 */
struct AsteroidData {
    int number;              ///< Asteroid number
    std::string name;        ///< Name
    double gm;               ///< Gravitational parameter [km³/s²]
    double a;                ///< Semi-major axis [AU]
    double e;                ///< Eccentricity
    double i;                ///< Inclination [degrees]
    double omega;            ///< Argument of perihelion [degrees]
    double Omega;            ///< Longitude of ascending node [degrees]
    double M0;               ///< Mean anomaly at epoch [degrees]
    double epoch_mjd;        ///< Epoch [MJD TDB]
    double mean_motion;      ///< Mean motion [deg/day]
    
    /**
     * @brief Compute mean anomaly at given time
     */
    double meanAnomalyAt(double mjd_tdb) const {
        return M0 + mean_motion * (mjd_tdb - epoch_mjd);
    }
};

/**
 * @brief Provider for massive asteroid perturbations
 */
class AsteroidPerturbations {
public:
    /**
     * @brief Default constructor - loads 16 most massive asteroids
     */
    AsteroidPerturbations();
    
    /**
     * @brief Load asteroid data from file
     * @param filename Path to asteroid orbital elements file
     * 
     * File format (CSV):
     * number,name,GM,a,e,i,omega,Omega,M0,epoch_mjd,n
     */
    explicit AsteroidPerturbations(const std::string& filename);
    
    /**
     * @brief Compute position of asteroid at given time
     * @param asteroid Asteroid data
     * @param mjd_tdb Time [MJD TDB]
     * @return Heliocentric position [AU] in J2000 ecliptic
     */
    Eigen::Vector3d getPosition(const AsteroidData& asteroid, double mjd_tdb) const;
    
    /**
     * @brief Compute perturbation acceleration from all asteroids
     * @param position Spacecraft/asteroid position [AU]
     * @param mjd_tdb Time [MJD TDB]
     * @return Acceleration [AU/day²]
     */
    Eigen::Vector3d computePerturbation(const Eigen::Vector3d& position, 
                                        double mjd_tdb) const;
    
    /**
     * @brief Compute perturbation from single asteroid
     * @param position Spacecraft position [AU]
     * @param asteroid_pos Asteroid position [AU]
     * @param gm Asteroid GM [km³/s²]
     * @return Acceleration [AU/day²]
     */
    static Eigen::Vector3d computeSinglePerturbation(
        const Eigen::Vector3d& position,
        const Eigen::Vector3d& asteroid_pos,
        double gm);
    
    /**
     * @brief Get list of included asteroids
     */
    const std::vector<AsteroidData>& getAsteroids() const { return asteroids_; }
    
    /**
     * @brief Enable/disable specific asteroid
     */
    void setAsteroidEnabled(int number, bool enabled);
    
    /**
     * @brief Check if asteroid is enabled
     */
    bool isAsteroidEnabled(int number) const;
    
    /**
     * @brief Get total mass of all asteroids [solar masses]
     */
    double getTotalMass() const;
    
    /**
     * @brief Load default AST17 asteroid data
     */
    void loadDefaultAsteroids();
    
private:
    
    /**
     * @brief Convert Keplerian elements to Cartesian
     */
    Eigen::Vector3d keplerianToCartesian(
        double a, double e, double i, 
        double omega, double Omega, double M,
        double gm_sun) const;
    
    std::vector<AsteroidData> asteroids_;
    std::vector<bool> enabled_flags_;
};

/**
 * @brief AST17 default asteroid data (Baer et al. 2011)
 * 
 * Reference: Baer, J., et al. (2011). "Astrometric masses of 21 asteroids,
 * and an integrated asteroid ephemeris". Celestial Mechanics and Dynamical
 * Astronomy, 100(1), 27-42.
 */
namespace ast17 {
    /**
     * @brief Get default AST17 asteroid parameters
     * @return Vector of 16 most massive asteroids with orbital elements
     */
    std::vector<AsteroidData> getDefaultAsteroids();
    
    // Constants for most massive asteroids (GM in km³/s²)
    constexpr double GM_CERES = 62.6284;         // (1) Ceres (dwarf planet)
    constexpr double GM_PALLAS = 13.8;           // (2) Pallas
    constexpr double GM_VESTA = 17.8;            // (4) Vesta
    constexpr double GM_HYGIEA = 5.78;           // (10) Hygiea
    constexpr double GM_EUNOMIA = 2.1;           // (15) Eunomia
    constexpr double GM_PSYCHE = 1.8;            // (16) Psyche
    constexpr double GM_EUPHROSYNE = 1.7;        // (31) Euphrosyne
    constexpr double GM_EUROPA = 1.59;           // (52) Europa
    constexpr double GM_CYBELE = 1.58;           // (65) Cybele
    constexpr double GM_SYLVIA = 1.50;           // (87) Sylvia
    constexpr double GM_THISBE = 1.3;            // (88) Thisbe
    constexpr double GM_DAVIDA = 2.0;            // (511) Davida
    constexpr double GM_INTERAMNIA = 2.1;        // (704) Interamnia
    constexpr double GM_BAMBERGA = 0.7;          // (324) Bamberga
    constexpr double GM_PATIENTIA = 0.8;         // (451) Patientia
    constexpr double GM_CAMILLA = 1.12;          // (107) Camilla
}

} // namespace ephemeris
} // namespace orbfit

#endif // ORBFIT_ASTEROID_PERTURBATIONS_HPP
