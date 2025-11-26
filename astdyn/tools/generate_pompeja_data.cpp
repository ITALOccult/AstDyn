/**
 * @file generate_pompeja_data.cpp
 * @brief Generate Pompeja orbit data for comparison report
 * 
 * Exports orbital elements and positions in JSON format for PDF generation
 */

#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/ephemeris/AsteroidPerturbations.hpp>
#include <astdyn/core/Constants.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::ephemeris;
using namespace astdyn::constants;

OrbitalElements equinoctial_to_keplerian(double a, double h, double k, 
                                         double p, double q, double lam, 
                                         double epoch_mjd) {
    double e = std::sqrt(h*h + k*k);
    double i = 2.0 * std::atan(std::sqrt(p*p + q*q));
    double Omega = std::atan2(p, q);
    double LP = std::atan2(h, k);
    double omega = LP - Omega;
    double M = lam * DEG_TO_RAD - LP;
    
    if (Omega < 0) Omega += TWO_PI;
    if (omega < 0) omega += TWO_PI;
    if (M < 0) M += TWO_PI;
    
    KeplerianElements kep;
    kep.epoch_mjd_tdb = epoch_mjd;
    kep.semi_major_axis = a;
    kep.eccentricity = e;
    kep.inclination = i;
    kep.longitude_ascending_node = Omega;
    kep.argument_perihelion = omega;
    kep.mean_anomaly = M;
    kep.gravitational_parameter = GMS;
    
    return kep;
}

int main() {
    std::cout << "Generating Pompeja (203) orbit data for comparison...\n";
    std::cout << "=========================================================\n\n";
    
    // Initial orbit from AstDyS (MJD 61000)
    auto kep_initial = equinoctial_to_keplerian(
        2.7385249933616391,    // a
        0.045087089252389,     // h
        0.041231297793564,     // k
        -0.005947645824719,    // p
        0.027042352297741,     // q
        112.3228065415555,     // lambda
        61000.0               // epoch
    );
    
    std::cout << "Initial orbit (MJD 61000):\n";
    std::cout << "  a = " << kep_initial.semi_major_axis << " AU\n";
    std::cout << "  e = " << kep_initial.eccentricity << "\n";
    std::cout << "  i = " << kep_initial.inclination * RAD_TO_DEG << "°\n\n";
    
    // Configure full dynamics propagator
    PropagatorSettings settings;
    settings.include_planets = true;
    settings.perturb_mercury = true;
    settings.perturb_venus = true;
    settings.perturb_earth = true;
    settings.perturb_mars = true;
    settings.perturb_jupiter = true;
    settings.perturb_saturn = true;
    settings.perturb_uranus = true;
    settings.perturb_neptune = true;
    settings.include_asteroids = true;
    
    auto ephemeris = std::make_shared<PlanetaryEphemeris>();
    auto integrator = std::make_unique<RKF78Integrator>(0.1, 1e-12);
    Propagator propagator(std::move(integrator), ephemeris, settings);
    
    // Propagate to MJD 61192 (like OrbFit output)
    std::cout << "Propagating to MJD 61192 (192 days)...\n";
    double target_mjd = 61192.0;
    auto kep_propagated = propagator.propagate_keplerian(kep_initial, target_mjd);
    auto cart_propagated = keplerian_to_cartesian(kep_propagated);
    
    std::cout << "\nFinal orbit (MJD 61192):\n";
    std::cout << std::fixed << std::setprecision(9);
    std::cout << "  a = " << kep_propagated.semi_major_axis << " AU\n";
    std::cout << "  e = " << kep_propagated.eccentricity << "\n";
    std::cout << "  i = " << kep_propagated.inclination * RAD_TO_DEG << "°\n";
    std::cout << "  Position: [" << cart_propagated.position[0] << ", "
              << cart_propagated.position[1] << ", "
              << cart_propagated.position[2] << "] AU\n\n";
    
    // Generate ephemeris at multiple epochs and save to CSV
    std::cout << "Generating ephemeris (10 days, 1-day steps)...\n";
    std::string filename = "pompeja_astdyn_data.csv";
    std::ofstream outfile(filename);
    
    // Write header
    outfile << "mjd,a_au,e,i_deg,Omega_deg,omega_deg,M_deg,x_au,y_au,z_au,vx_au_per_day,vy_au_per_day,vz_au_per_day\n";
    
    // Write initial epoch
    auto cart_init = keplerian_to_cartesian(kep_initial);
    outfile << std::fixed << std::setprecision(12)
            << kep_initial.epoch_mjd_tdb << ","
            << kep_initial.semi_major_axis << ","
            << kep_initial.eccentricity << ","
            << kep_initial.inclination * RAD_TO_DEG << ","
            << kep_initial.longitude_ascending_node * RAD_TO_DEG << ","
            << kep_initial.argument_perihelion * RAD_TO_DEG << ","
            << kep_initial.mean_anomaly * RAD_TO_DEG << ","
            << cart_init.position[0] << ","
            << cart_init.position[1] << ","
            << cart_init.position[2] << ","
            << cart_init.velocity[0] << ","
            << cart_init.velocity[1] << ","
            << cart_init.velocity[2] << "\n";
    
    // Write propagated epoch
    outfile << target_mjd << ","
            << kep_propagated.semi_major_axis << ","
            << kep_propagated.eccentricity << ","
            << kep_propagated.inclination * RAD_TO_DEG << ","
            << kep_propagated.longitude_ascending_node * RAD_TO_DEG << ","
            << kep_propagated.argument_perihelion * RAD_TO_DEG << ","
            << kep_propagated.mean_anomaly * RAD_TO_DEG << ","
            << cart_propagated.position[0] << ","
            << cart_propagated.position[1] << ","
            << cart_propagated.position[2] << ","
            << cart_propagated.velocity[0] << ","
            << cart_propagated.velocity[1] << ","
            << cart_propagated.velocity[2] << "\n";
    
    // Write ephemeris for 10 more days
    for (int day = 1; day <= 10; day++) {
        double mjd = target_mjd + day;
        auto kep_epoch = propagator.propagate_keplerian(kep_propagated, mjd);
        auto cart_epoch = keplerian_to_cartesian(kep_epoch);
        
        outfile << mjd << ","
                << kep_epoch.semi_major_axis << ","
                << kep_epoch.eccentricity << ","
                << kep_epoch.inclination * RAD_TO_DEG << ","
                << kep_epoch.longitude_ascending_node * RAD_TO_DEG << ","
                << kep_epoch.argument_perihelion * RAD_TO_DEG << ","
                << kep_epoch.mean_anomaly * RAD_TO_DEG << ","
                << cart_epoch.position[0] << ","
                << cart_epoch.position[1] << ","
                << cart_epoch.position[2] << ","
                << cart_epoch.velocity[0] << ","
                << cart_epoch.velocity[1] << ","
                << cart_epoch.velocity[2] << "\n";
    }
    
    outfile.close();
    
    std::cout << "\n✓ Data exported to: " << filename << "\n";
    std::cout << "=========================================================\n";
    
    return 0;
}
