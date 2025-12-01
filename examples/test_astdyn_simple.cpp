/**
 * @file test_astdyn_simple.cpp
 * @brief Test diretto di AstDyn con file .eq1
 * @author Michele Bigi - ITALOccultLibrary
 * @date 2025-12-01
 * 
 * Test semplificato che usa direttamente le classi AstDyn
 * per leggere elementi .eq1 e propagare un'orbita.
 */

#include <iostream>
#include <iomanip>
#include <chrono>
#include <astdyn/AstDyn.hpp>
#include <astdyn/io/parsers/OrbFitEQ1Parser.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/coordinates/CartesianState.hpp>
#include <astdyn/coordinates/KeplerianElements.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>

using namespace astdyn;
using namespace astdyn::propagation;
using namespace astdyn::io::parsers;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <eq1_file> <target_mjd_tdb>" << std::endl;
        std::cerr << "Example: " << argv[0] << " ../astdyn/data/17030.eq1 2460643.77083" << std::endl;
        return 1;
    }

    std::string eq1_file = argv[1];
    double target_mjd_tdb = std::stod(argv[2]);

    std::cout << "========================================" << std::endl;
    std::cout << "  AstDyn Orbit Propagation Test" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;

    try {
        // ===== FASE 1: Lettura elementi orbitali .eq1 =====
        std::cout << "[1/4] Lettura file .eq1..." << std::endl;
        OrbFitEQ1Parser parser;
        auto orbital_elements = parser.parse(eq1_file);
        
        std::cout << "  ✓ File letto: " << eq1_file << std::endl;
        std::cout << "  ✓ Oggetto: " << orbital_elements.object_name << std::endl;
        std::cout << "  ✓ Epoca (MJD TDB): " << std::fixed << std::setprecision(5) 
                  << orbital_elements.epoch_mjd_tdb << std::endl;
        std::cout << std::endl;

        // ===== FASE 2: Configurazione propagatore =====
        std::cout << "[2/4] Configurazione propagatore..." << std::endl;
        
        // Crea effemeridi planetarie
        auto ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
        
        // Crea integratore RKF78 (Runge-Kutta-Fehlberg 7/8)
        auto integrator = std::make_unique<RKF78Integrator>(
            0.1,     // step iniziale [giorni]
            1e-12    // tolleranza
        );
        
        // Configura perturbazioni (tutte attive)
        PropagatorSettings settings;
        settings.include_planets = true;
        settings.include_relativity = true;
        settings.include_asteroids = true;
        settings.perturb_mercury = true;
        settings.perturb_venus = true;
        settings.perturb_earth = true;
        settings.perturb_mars = true;
        settings.perturb_jupiter = true;
        settings.perturb_saturn = true;
        settings.perturb_uranus = true;
        settings.perturb_neptune = true;
        
        // Crea propagatore
        Propagator propagator(std::move(integrator), ephemeris, settings);
        
        std::cout << "  ✓ Integratore: RKF78" << std::endl;
        std::cout << "  ✓ Tolleranza: 1e-12" << std::endl;
        std::cout << "  ✓ Perturbazioni: 8 pianeti + relativistic + asteroids" << std::endl;
        std::cout << std::endl;

        // ===== FASE 3: Conversione e propagazione =====
        std::cout << "[3/4] Propagazione orbita..." << std::endl;
        std::cout << "  Epoca iniziale: " << orbital_elements.epoch_mjd_tdb << " MJD TDB" << std::endl;
        std::cout << "  Epoca target:   " << target_mjd_tdb << " MJD TDB" << std::endl;
        std::cout << "  Delta tempo:    " << (target_mjd_tdb - orbital_elements.epoch_mjd_tdb) 
                  << " giorni" << std::endl;
        std::cout << std::endl;

        // Converti OrbitalElements in KeplerianElements
        KeplerianElements kep_initial;
        kep_initial.semi_major_axis = orbital_elements.semi_major_axis;
        kep_initial.eccentricity = orbital_elements.eccentricity;
        kep_initial.inclination = orbital_elements.inclination;
        kep_initial.longitude_ascending_node = orbital_elements.longitude_asc_node;
        kep_initial.argument_perihelion = orbital_elements.argument_perihelion;
        kep_initial.mean_anomaly = orbital_elements.mean_anomaly;
        kep_initial.epoch_mjd_tdb = orbital_elements.epoch_mjd_tdb;
        kep_initial.gravitational_parameter = constants::GMS; // GM del Sole
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Propagazione
        KeplerianElements kep_final = propagator.propagate_keplerian(kep_initial, target_mjd_tdb);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << "  ✓ Propagazione completata in " << duration.count() << " ms" << std::endl;
        std::cout << std::endl;

        // ===== FASE 4: Risultati =====
        std::cout << "[4/4] Risultati alla epoch " << target_mjd_tdb << " MJD TDB:" << std::endl;
        std::cout << std::endl;
        
        // Converti in cartesiano per visualizzare posizione/velocità
        CartesianElements cart_final = keplerian_to_cartesian(kep_final);
        
        // IMPORTANTE: Conversione ECLM J2000 → ICRF (equatoriale J2000)
        // Gli elementi .eq1 sono in frame eclittico, dobbiamo ruotare
        double epsilon = 23.4393 * M_PI / 180.0;  // Obliquità eclittica J2000
        double cos_eps = std::cos(epsilon);
        double sin_eps = std::sin(epsilon);
        
        // Ruota da eclittico a equatoriale
        double x_ecl = cart_final.position.x();
        double y_ecl = cart_final.position.y();
        double z_ecl = cart_final.position.z();
        
        double x_icrf = x_ecl;
        double y_icrf = y_ecl * cos_eps - z_ecl * sin_eps;
        double z_icrf = y_ecl * sin_eps + z_ecl * cos_eps;
        
        double vx_ecl = cart_final.velocity.x();
        double vy_ecl = cart_final.velocity.y();
        double vz_ecl = cart_final.velocity.z();
        
        double vx_icrf = vx_ecl;
        double vy_icrf = vy_ecl * cos_eps - vz_ecl * sin_eps;
        double vz_icrf = vy_ecl * sin_eps + vz_ecl * cos_eps;
        
        std::cout << "Posizione ICRF (AU):" << std::endl;
        std::cout << "  X = " << std::fixed << std::setprecision(12) << x_icrf << std::endl;
        std::cout << "  Y = " << std::fixed << std::setprecision(12) << y_icrf << std::endl;
        std::cout << "  Z = " << std::fixed << std::setprecision(12) << z_icrf << std::endl;
        std::cout << std::endl;
        
        std::cout << "Velocità ICRF (AU/day):" << std::endl;
        std::cout << "  VX = " << std::fixed << std::setprecision(12) << vx_icrf << std::endl;
        std::cout << "  VY = " << std::fixed << std::setprecision(12) << vy_icrf << std::endl;
        std::cout << "  VZ = " << std::fixed << std::setprecision(12) << vz_icrf << std::endl;
        std::cout << std::endl;
        
        std::cout << "Elementi Kepleriani:" << std::endl;
        std::cout << "  a = " << std::fixed << std::setprecision(9) << kep_final.semi_major_axis << " AU" << std::endl;
        std::cout << "  e = " << std::fixed << std::setprecision(9) << kep_final.eccentricity << std::endl;
        std::cout << "  i = " << std::fixed << std::setprecision(6) << kep_final.inclination * 180.0 / M_PI << "°" << std::endl;
        std::cout << "  Ω = " << std::fixed << std::setprecision(6) << kep_final.longitude_ascending_node * 180.0 / M_PI << "°" << std::endl;
        std::cout << "  ω = " << std::fixed << std::setprecision(6) << kep_final.argument_perihelion * 180.0 / M_PI << "°" << std::endl;
        std::cout << "  M = " << std::fixed << std::setprecision(6) << kep_final.mean_anomaly * 180.0 / M_PI << "°" << std::endl;
        std::cout << std::endl;

        // Statistiche integrazione
        const auto& stats = propagator.statistics();
        std::cout << "Statistiche integrazione:" << std::endl;
        std::cout << "  Numero di step: " << stats.num_steps << std::endl;
        std::cout << "  Valutazioni funzione: " << stats.num_function_evals << std::endl;
        std::cout << "  Step rifiutati: " << stats.num_rejected_steps << std::endl;
        std::cout << std::endl;

        std::cout << "========================================" << std::endl;
        std::cout << "  Test completato con successo!" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << std::endl;
        
        std::cout << "Confronta questi risultati con JPL Horizons:" << std::endl;
        std::cout << "  https://ssd.jpl.nasa.gov/horizons/app.html" << std::endl;
        std::cout << "  Target: Asteroid " << orbital_elements.object_name << std::endl;
        std::cout << "  Epoch: " << target_mjd_tdb << " MJD TDB" << std::endl;
        std::cout << "  Reference frame: ICRF" << std::endl;
        std::cout << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << std::endl;
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
}
