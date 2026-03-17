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
#ifdef _OPENMP
#include <omp.h>
#endif

namespace astdyn {
namespace ephemeris {

using namespace constants;

AsteroidPerturbations::AsteroidPerturbations() {
    // Note: loadDefaultAsteroids is now called on demand by computePerturbation if vector is empty
    // to avoid potential hangs during object construction or shared_ptr allocation.
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
        
        double gm_val, a_au, i_deg, om_deg, Om_deg, M0_deg, mjd;
        iss >> ast.number >> comma;
        std::getline(iss, ast.name, ',');
        iss >> gm_val >> comma >> a_au >> comma >> ast.e >> comma 
            >> i_deg >> comma >> om_deg >> comma >> Om_deg >> comma 
            >> M0_deg >> comma >> mjd >> comma >> ast.mean_motion_deg_day;
        
        ast.gm = physics::GravitationalParameter::from_km3_s2(gm_val);
        ast.a = physics::Distance::from_au(a_au);
        ast.i = astrometry::Angle::from_deg(i_deg);
        ast.omega = astrometry::Angle::from_deg(om_deg);
        ast.Omega = astrometry::Angle::from_deg(Om_deg);
        ast.M0 = astrometry::Angle::from_deg(M0_deg);
        ast.epoch = time::EpochTDB::from_mjd(mjd);
        
        asteroids_.push_back(ast);
    }
    
    enabled_flags_.resize(asteroids_.size(), true);
}

void AsteroidPerturbations::loadDefaultAsteroids() {
    asteroids_ = ast17::getDefaultAsteroids();
    enabled_flags_.assign(asteroids_.size(), true);
}

void AsteroidPerturbations::loadDefault30Asteroids() {
    asteroids_.clear();
    constexpr double epoch = 51544.5;
    
    // Top 30 massive asteroids (BC405-like set)
    // number, name, GM, a, e, i, omega, Omega, M0, epoch, n, naif_id
    asteroids_.push_back({1, "Ceres", 62.6284, 2.7675, 0.0760, 10.593, 73.597, 80.3932, 113.410, epoch, 0.2141, 2000001});
    asteroids_.push_back({2, "Pallas", 13.80, 2.7730, 0.2299, 34.837, 310.049, 173.0962, 78.194, epoch, 0.2133, 2000002});
    asteroids_.push_back({3, "Juno", 1.8, 2.668, 0.258, 12.98, 248.1, 170.3, 34.2, epoch, 0.224, 2000003});
    asteroids_.push_back({4, "Vesta", 17.8, 2.3615, 0.0887, 7.142, 151.198, 103.8513, 205.546, epoch, 0.2716, 2000004});
    asteroids_.push_back({6, "Hebe", 0.8, 2.425, 0.201, 14.75, 239.1, 138.7, 120.5, epoch, 0.262, 2000006});
    asteroids_.push_back({7, "Iris", 1.2, 2.385, 0.231, 5.52, 145.23, 259.56, 18.3, epoch, 0.267, 2000007});
    asteroids_.push_back({8, "Flora", 0.5, 2.201, 0.156, 5.88, 285.4, 110.9, 15.6, epoch, 0.301, 2000008});
    asteroids_.push_back({9, "Metis", 0.1, 2.387, 0.122, 5.58, 5.7, 68.9, 140.2, epoch, 0.266, 2000009});
    asteroids_.push_back({10, "Hygiea", 5.78, 3.1393, 0.1146, 3.831, 283.455, 283.2045, 107.590, epoch, 0.1794, 2000010});
    asteroids_.push_back({13, "Egeria", 0.5, 2.576, 0.085, 16.54, 80.2, 43.3, 230.1, epoch, 0.238, 2000013});
    asteroids_.push_back({15, "Eunomia", 2.1, 2.6439, 0.1866, 11.737, 97.767, 293.2081, 142.405, epoch, 0.2262, 2000015});
    asteroids_.push_back({16, "Psyche", 1.8, 2.9216, 0.1339, 3.096, 227.305, 150.2873, 179.942, epoch, 0.1990, 2000016});
    asteroids_.push_back({19, "Fortuna", 0.4, 2.442, 0.159, 1.57, 182.1, 211.3, 15.6, epoch, 0.259, 2000019});
    asteroids_.push_back({20, "Massalia", 0.4, 2.409, 0.143, 0.71, 256.4, 206.5, 120.1, epoch, 0.263, 2000020});
    asteroids_.push_back({24, "Themis", 0.3, 3.129, 0.132, 0.76, 107.8, 35.9, 15.6, epoch, 0.180, 2000024});
    asteroids_.push_back({29, "Amphitrite", 0.4, 2.554, 0.073, 6.10, 63.4, 356.5, 15.6, epoch, 0.242, 2000029});
    asteroids_.push_back({31, "Euphrosyne", 1.7, 3.1515, 0.2257, 26.307, 61.588, 31.2873, 325.478, epoch, 0.1785, 2000031});
    asteroids_.push_back({45, "Eugenia", 0.4, 2.720, 0.082, 6.61, 204.3, 147.9, 15.6, epoch, 0.221, 2000045});
    asteroids_.push_back({48, "Doris", 0.2, 3.110, 0.075, 6.55, 250.2, 177.3, 15.6, epoch, 0.181, 2000048});
    asteroids_.push_back({52, "Europa", 1.59, 3.0958, 0.1020, 7.460, 343.545, 129.0380, 252.132, epoch, 0.1824, 2000052});
    asteroids_.push_back({65, "Cybele", 1.58, 3.4332, 0.1050, 3.562, 155.889, 155.5463, 196.234, epoch, 0.1636, 2000065});
    asteroids_.push_back({87, "Sylvia", 1.50, 3.4876, 0.0808, 10.866, 266.315, 73.3430, 271.478, epoch, 0.1610, 2000087});
    asteroids_.push_back({88, "Thisbe", 1.3, 2.7671, 0.1641, 5.222, 36.714, 276.0972, 355.489, epoch, 0.2141, 2000088});
    asteroids_.push_back({94, "Aurora", 0.2, 3.161, 0.089, 7.97, 10.2, 2.4, 15.6, epoch, 0.177, 2000094});
    asteroids_.push_back({107, "Camilla", 1.12, 3.4886, 0.0692, 10.044, 309.785, 173.0458, 120.234, epoch, 0.1609, 2000107});
    asteroids_.push_back({121, "Hermione", 0.3, 3.446, 0.135, 7.59, 240.1, 137.9, 15.6, epoch, 0.163, 2000121});
    asteroids_.push_back({324, "Bamberga", 0.7, 2.6835, 0.3381, 11.095, 43.678, 327.8394, 155.234, epoch, 0.2226, 2000324});
    asteroids_.push_back({451, "Patientia", 0.8, 3.0639, 0.0764, 15.240, 278.789, 96.9873, 234.456, epoch, 0.1843, 2000451});
    asteroids_.push_back({511, "Davida", 2.0, 3.1805, 0.1833, 15.942, 107.234, 270.0945, 186.789, epoch, 0.1768, 2000511});
    asteroids_.push_back({704, "Interamnia", 2.1, 3.0616, 0.1502, 17.296, 94.789, 280.3456, 145.234, epoch, 0.1844, 2000704});

    enabled_flags_.assign(asteroids_.size(), true);
}

void AsteroidPerturbations::loadAstDynDefaultSet() {
    asteroids_.clear();
    constexpr double epoch = 51544.5;
    
    // 10: Pluto
    asteroids_.push_back({10, "Pluto", constants::GM_PLUTO, 39.48, 0.248, 17.14, 113.83, 110.29, 14.53, epoch, 0.0039, 999});
    
    // 11: Camilla
    asteroids_.push_back({11, "Camilla", ast17::GM_CAMILLA, 3.4886, 0.0692, 10.044, 309.785, 173.0458, 120.234, epoch, 0.1609, 2000107});
    
    // 12: Ceres
    asteroids_.push_back({12, "Ceres", ast17::GM_CERES, 2.7675, 0.0760, 10.593, 73.597, 80.3932, 113.410, epoch, 0.2141, 2000001});
    
    // 13: Cybele
    asteroids_.push_back({13, "Cybele", ast17::GM_CYBELE, 3.4332, 0.1050, 3.562, 155.889, 155.5463, 196.234, epoch, 0.1636, 2000065});
    
    // 14: Davida
    asteroids_.push_back({14, "Davida", ast17::GM_DAVIDA, 3.1805, 0.1833, 15.942, 107.234, 270.0945, 186.789, epoch, 0.1768, 2000511});
    
    // 15: Eunomia
    asteroids_.push_back({15, "Eunomia", ast17::GM_EUNOMIA, 2.6439, 0.1866, 11.737, 97.767, 293.2081, 142.405, epoch, 0.2262, 2000015});
    
    // 16: Euphrosyne
    asteroids_.push_back({16, "Euphrosyne", ast17::GM_EUPHROSYNE, 3.1515, 0.2257, 26.307, 61.588, 31.2873, 325.478, epoch, 0.1785, 2000031});
    
    // 17: Europa
    asteroids_.push_back({17, "Europa", ast17::GM_EUROPA, 3.0958, 0.1020, 7.460, 343.545, 129.0380, 252.132, epoch, 0.1824, 2000052});
    
    // 18: Hygiea
    asteroids_.push_back({18, "Hygiea", ast17::GM_HYGIEA, 3.1393, 0.1146, 3.831, 283.455, 283.2045, 107.590, epoch, 0.1794, 2000010});
    
    // 19: Interamnia
    asteroids_.push_back({19, "Interamnia", ast17::GM_INTERAMNIA, 3.0616, 0.1502, 17.296, 94.789, 280.3456, 145.234, epoch, 0.1844, 2000704});
    
    // 20: Iris
    asteroids_.push_back({20, "Iris", 1.2, 2.385, 0.231, 5.52, 145.23, 259.56, 18.3, epoch, 0.267, 2000007});
    
    // 21: Juno
    asteroids_.push_back({21, "Juno", 1.8, 2.668, 0.258, 12.98, 248.1, 170.3, 34.2, epoch, 0.224, 2000003});
    
    // 22: Pallas
    asteroids_.push_back({22, "Pallas", ast17::GM_PALLAS, 2.7730, 0.2299, 34.837, 310.049, 173.0962, 78.194, epoch, 0.2133, 2000002});
    
    // 23: Psyche
    asteroids_.push_back({23, "Psyche", ast17::GM_PSYCHE, 2.9216, 0.1339, 3.096, 227.305, 150.2873, 179.942, epoch, 0.1990, 2000016});
    
    // 24: Sylvia
    asteroids_.push_back({24, "Sylvia", ast17::GM_SYLVIA, 3.4876, 0.0808, 10.866, 266.315, 73.3430, 271.478, epoch, 0.1610, 2000087});
    
    // 25: Thisbe
    asteroids_.push_back({25, "Thisbe", ast17::GM_THISBE, 2.7671, 0.1641, 5.222, 36.714, 276.0972, 355.489, epoch, 0.2141, 2000088});
    
    // 26: Vesta
    asteroids_.push_back({26, "Vesta", ast17::GM_VESTA, 2.3615, 0.0887, 7.142, 151.198, 103.8513, 205.546, epoch, 0.2716, 2000004});

    enabled_flags_.assign(asteroids_.size(), true);
}


void AsteroidPerturbations::loadSPK(const std::string& filename) {
    try {
        spk_reader_ = std::make_unique<io::SPKReader>(filename);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to load Asteroid SPK: " + std::string(e.what()));
    }
}

math::Vector3<core::ECLIPJ2000, physics::Distance> AsteroidPerturbations::getPosition(
    const AsteroidData& asteroid, 
    time::EpochTDB t) const 
{
    // Try SPK first if loaded
    if (spk_reader_) {
        try {
            int naif_id = asteroid.naif_id;
            double et = (t.mjd() - 51544.5) * 86400.0;
            
            Eigen::VectorXd state = spk_reader_->getState(naif_id, et);
            // SPK returns km -> convert to SI meters
            return math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(
                state[0] * 1000.0, state[1] * 1000.0, state[2] * 1000.0);
            
        } catch (...) {
            // Fallback for some weird cases where only number works
            try {
                 int raw_id = asteroid.number;
                 double et = (t.mjd() - 51544.5) * 86400.0;
                 Eigen::VectorXd state = spk_reader_->getState(raw_id, et);
                 return math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(
                    state[0] * 1000.0, state[1] * 1000.0, state[2] * 1000.0);
            } catch (...) { }
        }
    }

    // Analytical Fallback (J2000 Elements)
    astrometry::Angle M = asteroid.meanAnomalyAt(t);
    
    Eigen::Vector3d pos_eq = keplerianToCartesian(
        asteroid.a.to_au(), asteroid.e, asteroid.i.to_rad(),
        asteroid.omega.to_rad(), asteroid.Omega.to_rad(), 
        M.to_rad(), constants::GM_SUN);
    
    // Transform back to Ecliptic for the return type (since keplerianToCartesian returns Equatorial)
    Eigen::Matrix3d R = coordinates::ReferenceFrame::j2000_to_ecliptic();
    Eigen::Vector3d pos_ecl = R * pos_eq;
    
    return math::Vector3<core::ECLIPJ2000, physics::Distance>::from_si(
        physics::Distance::from_au(pos_ecl.x()).to_m(),
        physics::Distance::from_au(pos_ecl.y()).to_m(),
        physics::Distance::from_au(pos_ecl.z()).to_m()
    );
}

template <typename Frame>
math::Vector3<Frame, physics::Acceleration> AsteroidPerturbations::computePerturbation(
    const physics::CartesianStateTyped<Frame>& state, 
    const math::Vector3<core::GCRF, physics::Distance>& sun_pos_bary) const 
{
    math::Vector3<Frame, physics::Acceleration> total_acc = math::Vector3<Frame, physics::Acceleration>::from_si(0,0,0);
    
    // Transform Sun position to integration frame (not used currently, but kept for future use)
    // auto sun_pos_frame = coordinates::ReferenceFrame::transform_pos<core::GCRF, Frame>(sun_pos_bary, state.epoch);

    for (size_t i = 0; i < asteroids_.size(); ++i) {
        if (!enabled_flags_[i]) continue;
        
        const auto& asteroid = asteroids_[i];
        
        try {
            auto ast_pos_ecl = getPosition(asteroid, state.epoch);
            auto ast_pos_frame = coordinates::ReferenceFrame::transform_pos<core::ECLIPJ2000, Frame>(ast_pos_ecl, state.epoch);
            total_acc = total_acc + computeSinglePerturbation<Frame>(state.position, ast_pos_frame, asteroid.gm);
        } catch (...) { }
    }
    
    return total_acc;
}

template <typename Frame>
math::Vector3<Frame, physics::Acceleration> AsteroidPerturbations::computeSinglePerturbation(
    const math::Vector3<Frame, physics::Distance>& target_pos,
    const math::Vector3<Frame, physics::Distance>& asteroid_pos,
    const physics::GravitationalParameter& gm)
{
    // Relative position: asteroid - target
    auto delta = asteroid_pos - target_pos;
    double dist_m = delta.norm().to_m();
    double dist3 = dist_m * dist_m * dist_m;
    
    // Acceleration on target: a_dir = GM * delta / |delta|^3
    double gm_si = gm.to_m3_s2();
    auto acc_dir = math::Vector3<Frame, physics::Acceleration>::from_si(
        gm_si * delta.x_si() / dist3,
        gm_si * delta.y_si() / dist3,
        gm_si * delta.z_si() / dist3
    );

    // Indirect term (acceleration of Sun toward asteroid)
    double r_ast_m = asteroid_pos.norm().to_m();
    double r_ast3 = r_ast_m * r_ast_m * r_ast_m;
    auto acc_ind = math::Vector3<Frame, physics::Acceleration>::from_si(
        gm_si * asteroid_pos.x_si() / r_ast3,
        gm_si * asteroid_pos.y_si() / r_ast3,
        gm_si * asteroid_pos.z_si() / r_ast3
    );

    return acc_dir - acc_ind;
}

template math::Vector3<core::GCRF, physics::Acceleration> AsteroidPerturbations::computePerturbation(const physics::CartesianStateTyped<core::GCRF>&, const math::Vector3<core::GCRF, physics::Distance>&) const;
template math::Vector3<core::ECLIPJ2000, physics::Acceleration> AsteroidPerturbations::computePerturbation(const physics::CartesianStateTyped<core::ECLIPJ2000>&, const math::Vector3<core::GCRF, physics::Distance>&) const;

Eigen::Vector3d AsteroidPerturbations::computePerturbationRaw(
    const Eigen::Vector3d& pos_au,
    double mjd,
    const Eigen::Vector3d& sun_pos_bary_au,
    bool in_ecliptic) const
{
    Eigen::Vector3d total_acc_au_d2 = Eigen::Vector3d::Zero();
    time::EpochTDB t = time::EpochTDB::from_mjd(mjd);

    // Update cache if time changed
    if (std::abs(mjd - last_mjd_cached_) > 1e-14 || pos_cache_ecl_au_.size() != asteroids_.size()) {
        last_mjd_cached_ = mjd;
        pos_cache_ecl_au_.resize(asteroids_.size());
        for (size_t i = 0; i < asteroids_.size(); ++i) {
            if (!enabled_flags_[i]) continue;
            auto ast_pos_ecl_math = getPosition(asteroids_[i], t);
            pos_cache_ecl_au_[i] = ast_pos_ecl_math.to_eigen_si() / (constants::AU * 1000.0);
        }
    }

    for (size_t i = 0; i < asteroids_.size(); ++i) {
        if (!enabled_flags_[i]) continue;
        
        const auto& asteroid = asteroids_[i];
        
        try {
            Eigen::Vector3d ast_pos_frame_au = pos_cache_ecl_au_[i];
            if (!in_ecliptic) {
                ast_pos_frame_au = coordinates::ReferenceFrame::ecliptic_to_j2000() * ast_pos_frame_au;
            }

            Eigen::Vector3d delta = ast_pos_frame_au - pos_au;
            double d_norm = delta.norm();
            if (d_norm < 1e-4) continue; // Singularity protection (self-gravity or very close approach)
            double d_norm2 = d_norm * d_norm;
            double d3 = d_norm2 * d_norm;
            double gm_au_d2 = asteroid.gm.to_au3_d2();
            
            total_acc_au_d2 += gm_au_d2 * delta / d3;
            
            double r_ast = ast_pos_frame_au.norm();
            total_acc_au_d2 -= gm_au_d2 * ast_pos_frame_au / (r_ast * r_ast * r_ast);
            
        } catch (...) { }
    }
    
    return total_acc_au_d2;
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
            total_gm += asteroids_[i].gm.to_km3_s2();
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
    // Use rotation matrix directly (non-deprecated path)
    Eigen::Matrix3d R = coordinates::ReferenceFrame::ecliptic_to_j2000();
    return R * pos_ecliptic;
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
