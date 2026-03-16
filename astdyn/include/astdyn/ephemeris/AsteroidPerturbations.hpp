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

#ifndef ASTDYN_ASTEROID_PERTURBATIONS_HPP
#define ASTDYN_ASTEROID_PERTURBATIONS_HPP

#include <astdyn/coordinates/CartesianState.hpp>
#include <astdyn/core/physics_state.hpp>
#include <astdyn/core/Constants.hpp>
#include <Eigen/Dense>
#include <vector>
#include <vector>
#include <string>
#include <memory>

// Forward declaration
namespace astdyn { namespace io { class SPKReader; } }

#include <astdyn/core/physics_types.hpp>
#include <astdyn/astrometry/sky_types.hpp>
#include <astdyn/time/epoch.hpp>

namespace astdyn {
namespace ephemeris {

/**
 * @brief Asteroid data for perturbation calculations
 */
struct AsteroidData {
    int number;                           ///< Asteroid number (displayed/config)
    int naif_id;                          ///< NAIF ID for SPK lookup
    std::string name;                     ///< Name
    physics::GravitationalParameter gm;    ///< Gravitational parameter
    physics::Distance a;                  ///< Semi-major axis
    double e;                             ///< Eccentricity
    astrometry::Angle i;                  ///< Inclination
    astrometry::Angle omega;              ///< Argument of perihelion
    astrometry::Angle Omega;              ///< Longitude of ascending node
    astrometry::Angle M0;                 ///< Mean anomaly at epoch
    time::EpochTDB epoch;                 ///< Epoch
    double mean_motion_deg_day;           ///< Mean motion [deg/day]
    
    /**
     * @brief Construct from raw values (for compatibility with existing data tables)
     */
    AsteroidData(int num, std::string nam, double gm_val, double a_au, double ecc, 
                 double i_deg, double om_deg, double Om_deg, double m0_deg, 
                 double mjd, double n_deg_day, int naif = -1)
        : number(num), naif_id(naif == -1 ? (num < 1000000 ? 2000000 + num : num) : naif), 
          name(std::move(nam)), 
          gm(physics::GravitationalParameter::from_km3_s2(gm_val)),
          a(physics::Distance::from_au(a_au)), e(ecc),
          i(astrometry::Angle::from_deg(i_deg)),
          omega(astrometry::Angle::from_deg(om_deg)),
          Omega(astrometry::Angle::from_deg(Om_deg)),
          M0(astrometry::Angle::from_deg(m0_deg)),
          epoch(time::EpochTDB::from_mjd(mjd)),
          mean_motion_deg_day(n_deg_day) {}

    /**
     * @brief Default constructor
     */
    AsteroidData() = default;

    /**
     * @brief Compute mean anomaly at given time
     */
    astrometry::Angle meanAnomalyAt(time::EpochTDB t) const {
        double dt = t.mjd() - epoch.mjd();
        return astrometry::Angle::from_deg(M0.to_deg() + mean_motion_deg_day * dt).wrap_0_2pi();
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
    ~AsteroidPerturbations(); // Defined in cpp
    
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
     * @param t Time
     * @return Heliocentric position in J2000 ecliptic
     */
    math::Vector3<core::ECLIPJ2000, physics::Distance> getPosition(const AsteroidData& asteroid, time::EpochTDB t) const;
    
    /**
     * @brief Compute total perturbation acceleration from all loaded asteroids
     * 
     * @param state Heliocentric state of interest
     * @param sun_pos_bary Position of Sun relative to SSB (optional)
     * 
     * @return Acceleration vector in Frame
     */
    template <typename Frame>
    math::Vector3<Frame, physics::Acceleration> computePerturbation(
        const physics::CartesianStateTyped<Frame>& state, 
        const math::Vector3<core::GCRF, physics::Distance>& sun_pos_bary = math::Vector3<core::GCRF, physics::Distance>::from_si(0,0,0)) const;
    
    /**
     * @brief Compute perturbation from single asteroid
     * @param target_pos Position of object being perturbed
     * @param asteroid_pos Asteroid position
     * @param gm Asteroid GM
     * @return Acceleration
     */
    template <typename Frame>
    static math::Vector3<Frame, physics::Acceleration> computeSinglePerturbation(
        const math::Vector3<Frame, physics::Distance>& target_pos,
        const math::Vector3<Frame, physics::Distance>& asteroid_pos,
        const physics::GravitationalParameter& gm);
    
    /**
     * @brief Compute perturbation using raw Eigen vectors (AU for pos, MJD for time).
     * Returns acceleration in AU/day^2.
     */
    Eigen::Vector3d computePerturbationRaw(
        const Eigen::Vector3d& pos_au,
        double mjd,
        const Eigen::Vector3d& sun_pos_bary_au,
        bool in_ecliptic) const;
    
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

    /**
     * @brief Load the specific AstDyn default set of 17 asteroids (Pluto + 16 massive asteroids)
     * with custom numbering 10-26 as requested.
     */
    void loadAstDynDefaultSet();

    /**
     * @brief Load the top 30 most massive asteroids
     */
    void loadDefault30Asteroids();

    /**
     * @brief Load asteroid SPK file for precise positions
     */
    void loadSPK(const std::string& filename);
    
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
    
    // Optional SPK Reader for precise positions
    std::unique_ptr<io::SPKReader> spk_reader_;
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
} // namespace astdyn

#endif // ASTDYN_ASTEROID_PERTURBATIONS_HPP
