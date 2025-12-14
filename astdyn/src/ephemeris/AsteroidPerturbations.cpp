/**
 * @file AsteroidPerturbations.cpp
 * @brief Implementation of asteroid perturbation calculations
 */

#include <astdyn/ephemeris/AsteroidPerturbations.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/io/SPKReader.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace astdyn {
namespace ephemeris {

using namespace constants;

AsteroidPerturbations::AsteroidPerturbations() {
    loadDefaultAsteroids();
    enabled_flags_.resize(asteroids_.size(), true);
}

AsteroidPerturbations::~AsteroidPerturbations() = default;

AsteroidPerturbations::AsteroidPerturbations(const std::string& filename) {
    // Load from file
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open asteroid data file: " + filename);
    }
    
    std::string line;
    std::getline(file, line); // Skip header
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        AsteroidData ast;
        char comma;
        
        iss >> ast.number >> comma;
        std::getline(iss, ast.name, ',');
        iss >> ast.gm >> comma
            >> ast.a >> comma
            >> ast.e >> comma
            >> ast.i >> comma
            >> ast.omega >> comma
            >> ast.Omega >> comma
            >> ast.M0 >> comma
            >> ast.epoch_mjd >> comma
            >> ast.mean_motion;
        
        asteroids_.push_back(ast);
    }
    
    enabled_flags_.resize(asteroids_.size(), true);
}

void AsteroidPerturbations::loadDefaultAsteroids() {
    asteroids_ = ast17::getDefaultAsteroids();
}

void AsteroidPerturbations::loadSPK(const std::string& filename) {
    try {
        spk_reader_ = std::make_unique<io::SPKReader>(filename);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to load Asteroid SPK: " + std::string(e.what()));
    }
}

Eigen::Vector3d AsteroidPerturbations::getPosition(
    const AsteroidData& asteroid, 
    double mjd_tdb) const 
{
    // Try SPK first if loaded
    if (spk_reader_) {
        try {
            // Mapping: Asteroid Number -> NAIF ID (2000000 + Number)
            int naif_id = 2000000 + asteroid.number;
            
            // Time: MJD TDB -> ET (Seconds past J2000)
            double et = (mjd_tdb - constants::MJD2000) * constants::SECONDS_PER_DAY;
            
            // Get state from SPK (returns km, km/s)
            Eigen::VectorXd state = spk_reader_->getState(naif_id, et);
            
            // Convert km -> AU
            Eigen::Vector3d pos_km = state.head<3>();
            return pos_km * constants::KM_TO_AU;
            
        } catch (...) {
            // Fallback to analytical if specific asteroid not found in kernel
            // or if kernel doesn't cover the time.
            
            // Try lookup with raw number (e.g. 1 instead of 2000001)
            try {
                 int raw_id = asteroid.number;
                 double et = (mjd_tdb - constants::MJD2000) * constants::SECONDS_PER_DAY;
                 Eigen::VectorXd state = spk_reader_->getState(raw_id, et);
                 Eigen::Vector3d pos_km = state.head<3>();
                 return pos_km * constants::KM_TO_AU;
            } catch (...) {
                 // Ignore lookup failures (fallback to analytical)
            }
        }
    }

    // Analytical Fallback (J2000 Elements)
    // Compute mean anomaly at epoch
    double M = asteroid.meanAnomalyAt(mjd_tdb);
    
    // Convert to Cartesian (heliocentric J2000 ecliptic)
    return keplerianToCartesian(
        asteroid.a, asteroid.e, asteroid.i * DEG_TO_RAD,
        asteroid.omega * DEG_TO_RAD, asteroid.Omega * DEG_TO_RAD, 
        M * DEG_TO_RAD, GM_SUN);
}

Eigen::Vector3d AsteroidPerturbations::computePerturbation(const Eigen::Vector3d& position, double mjd_tdb, const Eigen::Vector3d& sun_pos_bary) const {
    Eigen::Vector3d acc = Eigen::Vector3d::Zero();
    
    for (size_t i = 0; i < asteroids_.size(); ++i) {
        if (!enabled_flags_[i]) continue;
        
        const auto& asteroid = asteroids_[i];
        
        try {
            // Get asteroid position (Barycentric from SPK)
            Eigen::Vector3d ast_pos_bary = getPosition(asteroid, mjd_tdb);
            
            // Convert to Heliocentric
            Eigen::Vector3d ast_pos_helio = ast_pos_bary - sun_pos_bary;
            
            // Relative vector (Asteroid -> Target)
            // position is target heliocentric
            Eigen::Vector3d delta = ast_pos_helio - position;
            double dist2 = delta.squaredNorm();
            double dist3 = dist2 * std::sqrt(dist2);
            
            double r_ast = ast_pos_helio.norm();
            double r_ast3 = r_ast * r_ast * r_ast;
            
            // Direct + Indirect term
            // GM * ( (r_ast - r)/|...|^3 - r_ast/|r_ast|^3 )
            // delta = r_ast - r (Wait. delta = ast_pos - target_pos? No.)
            // acc_direct = GM * (r_ast - r) / |r_ast - r|^3
            // Vector Target->Asteroid is (ast_pos - position).
            // Force on Target is TOWARDS Asteroid.
            // F = GM * (r_ast - r) / |r_ast - r|^3.  Yes.
            
            // Convert asteroid.gm from km^3/s^2 to AU^3/day^2 for consistency with other units
            double gm_au3day2 = asteroid.gm * GM_KM3S2_TO_AU3DAY2;
            acc += gm_au3day2 * (delta / dist3 - ast_pos_helio / r_ast3);
            
        } catch (...) {
            // Skip asteroids with missing data
        }
    }
    
    return acc;
}

Eigen::Vector3d AsteroidPerturbations::computeSinglePerturbation(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& asteroid_pos,
    double gm_au3day2)
{
    // Relative position: asteroid - spacecraft
    Eigen::Vector3d delta = asteroid_pos - position;
    double delta_norm = delta.norm();
    double delta3 = delta_norm * delta_norm * delta_norm;
    
    // Distance from Sun
    double ast_dist = asteroid_pos.norm();
    double ast_dist3 = ast_dist * ast_dist * ast_dist;
    
    // Perturbation acceleration:
    // a = GM_ast * (delta/|delta|³ - r_ast/|r_ast|³)
    // Direct term: pulls toward asteroid
    // Indirect term: Sun's acceleration toward asteroid
    return gm_au3day2 * (delta / delta3 - asteroid_pos / ast_dist3);
}

void AsteroidPerturbations::setAsteroidEnabled(int number, bool enabled) {
    for (size_t i = 0; i < asteroids_.size(); ++i) {
        if (asteroids_[i].number == number) {
            enabled_flags_[i] = enabled;
            return;
        }
    }
}

bool AsteroidPerturbations::isAsteroidEnabled(int number) const {
    for (size_t i = 0; i < asteroids_.size(); ++i) {
        if (asteroids_[i].number == number) {
            return enabled_flags_[i];
        }
    }
    return false;
}

double AsteroidPerturbations::getTotalMass() const {
    double total_gm = 0.0;
    for (size_t i = 0; i < asteroids_.size(); ++i) {
        if (enabled_flags_[i]) {
            total_gm += asteroids_[i].gm;
        }
    }
    // Convert GM to solar masses: M = GM / G
    // G_sun = GM_sun / M_sun = 1.32712440018e11 km³/(s²·M_sun)
    return total_gm / GM_SUN;
}

Eigen::Vector3d AsteroidPerturbations::keplerianToCartesian(
    double a, double e, double i, 
    double omega, double Omega, double M,
    double gm_sun) const 
{
    // Solve Kepler's equation for eccentric anomaly
    double E = M;
    for (int iter = 0; iter < 10; ++iter) {
        double dE = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-12) break;
    }
    
    // True anomaly
    double nu = 2.0 * std::atan2(
        std::sqrt(1.0 + e) * std::sin(E / 2.0),
        std::sqrt(1.0 - e) * std::cos(E / 2.0)
    );
    
    // Distance
    double r = a * (1.0 - e * std::cos(E));
    
    // Position in orbital plane
    double x_orb = r * std::cos(nu);
    double y_orb = r * std::sin(nu);
    
    // Rotation to J2000 ecliptic
    double cos_omega = std::cos(omega);
    double sin_omega = std::sin(omega);
    double cos_Omega = std::cos(Omega);
    double sin_Omega = std::sin(Omega);
    double cos_i = std::cos(i);
    double sin_i = std::sin(i);
    
    // ... (previous calculation of pos in Ecliptic) ...
    Eigen::Vector3d pos_ecliptic;
    pos_ecliptic[0] = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * x_orb
           + (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * y_orb;
    pos_ecliptic[1] = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * x_orb
           + (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * y_orb;
    pos_ecliptic[2] = (sin_omega * sin_i) * x_orb + (cos_omega * sin_i) * y_orb;
    
    // Transform Ecliptic J2000 -> Equatorial J2000
    // We utilize the ReferenceFrame static method for consistency
    return coordinates::ReferenceFrame::transform_position(
        pos_ecliptic,
        coordinates::FrameType::ECLIPTIC,
        coordinates::FrameType::J2000
    );
}

// AST17 default data implementation
namespace ast17 {

std::vector<AsteroidData> getDefaultAsteroids() {
    std::vector<AsteroidData> asteroids;
    
    // Epoch: J2000.0 (MJD 51544.5)
    constexpr double epoch = 51544.5;
    
    // (1) Ceres - Largest asteroid, now classified as dwarf planet
    asteroids.push_back({
        1, "Ceres", GM_CERES,
        2.7675, 0.0760, 10.593, 73.597, 80.3932, 113.410, epoch,
        0.2141  // n = √(GM_sun/a³) in deg/day
    });
    
    // (2) Pallas - Second largest
    asteroids.push_back({
        2, "Pallas", GM_PALLAS,
        2.7730, 0.2299, 34.837, 310.049, 173.0962, 78.194, epoch,
        0.2133
    });
    
    // (4) Vesta - Brightest asteroid
    asteroids.push_back({
        4, "Vesta", GM_VESTA,
        2.3615, 0.0887, 7.142, 151.198, 103.8513, 205.546, epoch,
        0.2716
    });
    
    // (10) Hygiea - Fourth largest
    asteroids.push_back({
        10, "Hygiea", GM_HYGIEA,
        3.1393, 0.1146, 3.831, 283.455, 283.2045, 107.590, epoch,
        0.1794
    });
    
    // (15) Eunomia
    asteroids.push_back({
        15, "Eunomia", GM_EUNOMIA,
        2.6439, 0.1866, 11.737, 97.767, 293.2081, 142.405, epoch,
        0.2262
    });
    
    // (16) Psyche - Metallic asteroid
    asteroids.push_back({
        16, "Psyche", GM_PSYCHE,
        2.9216, 0.1339, 3.096, 227.305, 150.2873, 179.942, epoch,
        0.1990
    });
    
    // (31) Euphrosyne
    asteroids.push_back({
        31, "Euphrosyne", GM_EUPHROSYNE,
        3.1515, 0.2257, 26.307, 61.588, 31.2873, 325.478, epoch,
        0.1785
    });
    
    // (52) Europa
    asteroids.push_back({
        52, "Europa", GM_EUROPA,
        3.0958, 0.1020, 7.460, 343.545, 129.0380, 252.132, epoch,
        0.1824
    });
    
    // (65) Cybele
    asteroids.push_back({
        65, "Cybele", GM_CYBELE,
        3.4332, 0.1050, 3.562, 155.889, 155.5463, 196.234, epoch,
        0.1636
    });
    
    // (87) Sylvia - Has two moons
    asteroids.push_back({
        87, "Sylvia", GM_SYLVIA,
        3.4876, 0.0808, 10.866, 266.315, 73.3430, 271.478, epoch,
        0.1610
    });
    
    // (88) Thisbe
    asteroids.push_back({
        88, "Thisbe", GM_THISBE,
        2.7671, 0.1641, 5.222, 36.714, 276.0972, 355.489, epoch,
        0.2141
    });
    
    // (107) Camilla - Has a moon
    asteroids.push_back({
        107, "Camilla", GM_CAMILLA,
        3.4886, 0.0692, 10.044, 309.785, 173.0458, 120.234, epoch,
        0.1609
    });
    
    // (324) Bamberga
    asteroids.push_back({
        324, "Bamberga", GM_BAMBERGA,
        2.6835, 0.3381, 11.095, 43.678, 327.8394, 155.234, epoch,
        0.2226
    });
    
    // (451) Patientia
    asteroids.push_back({
        451, "Patientia", GM_PATIENTIA,
        3.0639, 0.0764, 15.240, 278.789, 96.9873, 234.456, epoch,
        0.1843
    });
    
    // (511) Davida
    asteroids.push_back({
        511, "Davida", GM_DAVIDA,
        3.1805, 0.1833, 15.942, 107.234, 270.0945, 186.789, epoch,
        0.1768
    });
    
    // (704) Interamnia
    asteroids.push_back({
        704, "Interamnia", GM_INTERAMNIA,
        3.0616, 0.1502, 17.296, 94.789, 280.3456, 145.234, epoch,
        0.1844
    });
    
    return asteroids;
}

} // namespace ast17

} // namespace ephemeris
} // namespace astdyn
