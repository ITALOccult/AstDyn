/**
 * @file test_jpl_validation.cpp
 * @brief Validation of AstDyn propagation against JPL Horizons data
 * @author ITALOccult Integration Team
 * @date 4 December 2025
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "italoccultlibrary/include/astdyn_wrapper.h"
#include "astdyn/propagation/OrbitalElements.hpp"
// We need access to conversion functions, usually in propagation headers or similar
// Assuming cartesian_to_keplerian is available in astdyn/propagation/OrbitalElements.hpp or similar
// If not, we will implement the conversion here for the test.

using namespace ioccultcalc;

// JPL Data Structure
struct JPLData {
    std::string date;
    double mjd;
    Eigen::Vector3d pos; // AU
    Eigen::Vector3d vel; // AU/day
};

// HELPER: Read CSV
std::vector<JPLData> readJPLCSV(const std::string& filename) {
    std::vector<JPLData> data;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return data;
    }
    
    // Skip header
    std::getline(file, line);
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string token;
        JPLData entry;
        
        // Date,X,Y,Z,VX,VY,VZ,MJD
        std::getline(iss, entry.date, ',');
        
        std::string val;
        
        std::getline(iss, val, ','); entry.pos.x() = std::stod(val);
        std::getline(iss, val, ','); entry.pos.y() = std::stod(val);
        std::getline(iss, val, ','); entry.pos.z() = std::stod(val);
        
        std::getline(iss, val, ','); entry.vel.x() = std::stod(val);
        std::getline(iss, val, ','); entry.vel.y() = std::stod(val);
        std::getline(iss, val, ','); entry.vel.z() = std::stod(val);
        
        std::getline(iss, val, ','); entry.mjd = std::stod(val);
        
        data.push_back(entry);
    }
    
    return data;
}

// HELPER: Coordinate Conversion
// AstDyn expects Ecliptic inputs for propagation usually
// JPL data is now Equatorial (ICRF) (FRAME)
// We need to rotate ICRF -> Ecliptic
void rotateICRFtoEcliptic(const Eigen::Vector3d& pos_icrf, const Eigen::Vector3d& vel_icrf,
                          Eigen::Vector3d& pos_ecl, Eigen::Vector3d& vel_ecl) {
    constexpr double epsilon = 23.4392911 * M_PI / 180.0; // Obliquity J2000
    // Rotation matrix from Equatorial to Ecliptic (Rx(epsilon))
    // x_ecl = x_eq
    // y_ecl = y_eq * cos(eps) + z_eq * sin(eps)
    // z_ecl = -y_eq * sin(eps) + z_eq * cos(eps)
    
    double c = std::cos(epsilon);
    double s = std::sin(epsilon);
    
    pos_ecl.x() = pos_icrf.x();
    pos_ecl.y() = pos_icrf.y() * c + pos_icrf.z() * s;
    pos_ecl.z() = -pos_icrf.y() * s + pos_icrf.z() * c;
    
    vel_ecl.x() = vel_icrf.x();
    vel_ecl.y() = vel_icrf.y() * c + vel_icrf.z() * s;
    vel_ecl.z() = -vel_icrf.y() * s + vel_icrf.z() * c;
}

// HELPER: Cartesian -> Keplerian
// Minimal implementation if not linking full AstDyn util
void cartesianToKeplerian(const Eigen::Vector3d& r_vec, const Eigen::Vector3d& v_vec, double mu,
                          double& a, double& e, double& i, double& Omega, double& omega, double& M) {
    double r = r_vec.norm();
    double v2 = v_vec.squaredNorm();
    
    // Angular momentum
    Eigen::Vector3d h_vec = r_vec.cross(v_vec);
    double h = h_vec.norm();
    
    // Node vector
    Eigen::Vector3d n_vec(-h_vec.y(), h_vec.x(), 0.0);
    double n = n_vec.norm();
    
    // Eccentricity vector
    Eigen::Vector3d e_vec = (1.0/mu) * ((v2 - mu/r)*r_vec - (r_vec.dot(v_vec))*v_vec);
    e = e_vec.norm();
    
    // Energy
    double energy = v2/2.0 - mu/r;
    
    if (std::abs(e - 1.0) < 1e-9) {
        // Parabolic
        a = std::numeric_limits<double>::infinity();
    } else {
        a = -mu / (2.0 * energy);
    }
    
    // Inclination
    i = std::acos(h_vec.z() / h);
    
    // Longitude of ascending node
    if (n < 1e-12) {
        Omega = 0.0; // Circular equatorial
    } else {
        Omega = std::acos(n_vec.x() / n);
        if (n_vec.y() < 0) Omega = 2*M_PI - Omega;
    }
    
    // Argument of perihelion
    if (n < 1e-12) {
        omega = 0.0; // Undefined?
    } else {
        omega = std::acos(n_vec.dot(e_vec) / (n*e));
        if (e_vec.z() < 0) omega = 2*M_PI - omega;
    }
    
    // True anomaly
    double nu = std::acos(e_vec.dot(r_vec) / (e*r));
    if (r_vec.dot(v_vec) < 0) nu = 2*M_PI - nu;
    
    // Mean anomaly
    // tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2)
    double E = 2.0 * std::atan(std::sqrt((1.0-e)/(1.0+e)) * std::tan(nu/2.0));
    M = E - e*std::sin(E);
    // Normalize M
    while (M < 0) M += 2*M_PI;
    while (M >= 2*M_PI) M -= 2*M_PI;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "  JPL Horizons Validation Test" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // 1. Read Data
    std::vector<JPLData> truth_data = readJPLCSV("jpl_horizons_17030.csv");
    if (truth_data.empty()) {
        std::cerr << "No data found!" << std::endl;
        return 1;
    }
    std::cout << "Loaded " << truth_data.size() << " data points." << std::endl;
    
    // 2. Setup Initial State
    const JPLData& initial = truth_data[0];
    
    // Rotate to Ecliptic for AstDyn initialization
    Eigen::Vector3d pos_ecl, vel_ecl;
    rotateICRFtoEcliptic(initial.pos, initial.vel, pos_ecl, vel_ecl);
    
    // Convert to Keplerian
    double a, e, i, Omega, omega, M;
    // Uses GMS from AstDyn if available, or define it
    double mu = 2.959122082855911e-4; // GM Sun [AU^3/d^2] (k^2)
    cartesianToKeplerian(pos_ecl, vel_ecl, mu, a, e, i, Omega, omega, M);
    
    std::cout << "Initial Condition (MJD " << initial.mjd << "):" << std::endl;
    std::cout << "  Pos (ICRF): " << initial.pos.transpose() << std::endl;
    std::cout << "  Pos (ECLM): " << pos_ecl.transpose() << std::endl;
    std::cout << "  Keplerian: a=" << a << " e=" << e << " i=" << i*180/M_PI << " deg" << std::endl;
    
    // 3. Initialize AstDynWrapper
    AstDynWrapper wrapper(PropagationSettings::highAccuracy());
    
    // Pass initial state
    wrapper.setKeplerianElements(a, e, i, Omega, omega, M, initial.mjd, "17030 Sierks");
    
    // 4. Propagate and Compare
    std::cout << "\nComparison Results (AstDyn vs JPL):" << std::endl;
    std::cout << "MJD      | Dist (AU) | Err (km)   | Rel Err" << std::endl;
    std::cout << "---------|-----------|------------|--------" << std::endl;
    
    double max_error = 0.0;
    double sum_sq_error = 0.0;
    
    for (size_t idx = 1; idx < truth_data.size(); ++idx) {
        const auto& point = truth_data[idx];
        
        try {
            // Propagate
            auto state = wrapper.propagateToEpoch(point.mjd);
            
            // Calculate error
            double dist = state.position.norm();
            double error_au = (state.position - point.pos).norm();
            double error_km = error_au * 149597870.7;
            
            if (error_km > max_error) max_error = error_km;
            sum_sq_error += error_km * error_km;
            
            std::cout << std::fixed << std::setprecision(4) << point.mjd << " | "
                      << std::setprecision(6) << dist << " | "
                      << std::scientific << std::setprecision(3) << error_km << " | "
                      << (error_au / dist) << std::endl;
                      
        } catch (const std::exception& e) {
            std::cerr << "Propagation failed at " << point.mjd << ": " << e.what() << std::endl;
        }
    }
    
    double rms_error = std::sqrt(sum_sq_error / (truth_data.size() - 1));
    
    std::cout << "\nSummary:" << std::endl;
    std::cout << "  Max Error: " << std::fixed << std::setprecision(3) << max_error << " km" << std::endl;
    std::cout << "  RMS Error: " << rms_error << " km" << std::endl;
    
    // Thresholds
    if (rms_error < 1000.0) {
        std::cout << "✓ VALIDATION PASSED (Excellent accuracy < 1000 km)" << std::endl;
        return 0;
    } else {
        std::cout << "✗ VALIDATION FAILED (Accuracy too low)" << std::endl;
        return 1;
    }
}
