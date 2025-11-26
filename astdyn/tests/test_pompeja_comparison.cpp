/**
 * @file test_pompeja_comparison.cpp
 * @brief Confronto completo OrbFit vs AstDyn vs Horizons per 203 Pompeja
 * 
 * Confronta:
 * - Elementi equinoziali (EQU)
 * - Elementi eclittici (ECL) 
 * - Vettori cartesiani (posizione e velocità)
 * - Coordinate equatoriali (RA, Dec)
 * 
 * All'epoca MJD 61192.0 TDT
 */

#include <gtest/gtest.h>
#include <astdyn/orbit_determination/DifferentialCorrector.hpp>
#include <astdyn/orbit_determination/Residuals.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <astdyn/coordinates/KeplerianElements.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace astdyn;
using namespace astdyn::observations;
using namespace astdyn::orbit_determination;
using namespace astdyn::coordinates;

class PompejaComparisonTest : public ::testing::Test {
protected:
    struct EquinoctialElements {
        double a, h, k, p, q, lambda;
        double mjd_tdt;
    };
    
    // Parse OrbFit .eq1 file
    EquinoctialElements parse_eq1_file(const std::string& filename) {
        EquinoctialElements elem{};
        std::ifstream file(filename);
        std::string line;
        bool in_elements = false;
        
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '!') continue;
            if (line.find("END_OF_HEADER") != std::string::npos) {
                in_elements = true;
                continue;
            }
            if (!in_elements) continue;
            if (line.find("203") != std::string::npos && line.size() < 10) continue;
            
            if (line.find("EQU") != std::string::npos) {
                std::istringstream iss(line);
                std::string equ_tag;
                iss >> equ_tag >> elem.a >> elem.h >> elem.k >> elem.p >> elem.q >> elem.lambda;
                elem.lambda *= constants::DEG_TO_RAD;
            }
            
            if (line.find("MJD") != std::string::npos && line.find("TDT") != std::string::npos) {
                std::istringstream iss(line);
                std::string mjd_tag;
                double mjd_value;
                std::string tdt_tag;
                iss >> mjd_tag >> mjd_value >> tdt_tag;
                elem.mjd_tdt = mjd_value;
            }
        }
        return elem;
    }
    
    // Parse OrbFit .oel file
    EquinoctialElements parse_oel_file(const std::string& filename) {
        EquinoctialElements elem{};
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        std::string line;
        
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '!') continue;
            
            std::istringstream iss(line);
            std::string key;
            if (iss >> key) {
                // Skip lines with '=' (format lines)
                if (key.find('=') != std::string::npos || key == "format" || key == "recordtype") continue;
                
                double value;
                if (iss >> value) {
                    if (key == "MJD") elem.mjd_tdt = value;
                    else if (key == "a") elem.a = value;
                    else if (key == "h") elem.h = value;
                    else if (key == "k") elem.k = value;
                    else if (key == "p") elem.p = value;
                    else if (key == "q") elem.q = value;
                    else if (key == "lambda") elem.lambda = value;
                }
            }
        }
        return elem;
    }
    
    // Convert equinoctial to Keplerian
    propagation::KeplerianElements equinoctial_to_keplerian(const EquinoctialElements& eq) {
        propagation::KeplerianElements kep;
        kep.semi_major_axis = eq.a;
        kep.eccentricity = std::sqrt(eq.h * eq.h + eq.k * eq.k);
        kep.inclination = 2.0 * std::atan(std::sqrt(eq.p * eq.p + eq.q * eq.q));
        kep.longitude_ascending_node = std::atan2(eq.p, eq.q);
        if (kep.longitude_ascending_node < 0) kep.longitude_ascending_node += 2.0 * constants::PI;
        double omega_plus_Omega = std::atan2(eq.h, eq.k);
        if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * constants::PI;
        kep.argument_perihelion = omega_plus_Omega - kep.longitude_ascending_node;
        if (kep.argument_perihelion < 0) kep.argument_perihelion += 2.0 * constants::PI;
        kep.mean_anomaly = eq.lambda - omega_plus_Omega;
        while (kep.mean_anomaly < 0) kep.mean_anomaly += 2.0 * constants::PI;
        while (kep.mean_anomaly >= 2.0 * constants::PI) kep.mean_anomaly -= 2.0 * constants::PI;
        kep.epoch_mjd_tdb = eq.mjd_tdt;
        kep.gravitational_parameter = constants::GMS;
        return kep;
    }
    
    void print_header(const std::string& title) {
        std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
        std::cout << "║ " << std::left << std::setw(59) << title << "║\n";
        std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    }
};

TEST_F(PompejaComparisonTest, FullComparison) {
    print_header("POMPEJA 203 - CONFRONTO ORBFIT vs ASTDYN vs HORIZONS");
    
    std::cout << "Epoca di confronto: MJD 61192.0 TDT (2027-07-11)\n";
    std::cout << "Osservazioni: 100 recenti (2025-01 a 2025-10)\n\n";
    
    // ========================================================================
    // 1. ORBFIT - Carica elementi fittati
    // ========================================================================
    print_header("1. ORBFIT - Elementi fittati");
    
    auto orbfit_eq = parse_oel_file("astdyn/tools/203.oel");
    auto orbfit_kep = equinoctial_to_keplerian(orbfit_eq);
    auto orbfit_cart = keplerian_to_cartesian(orbfit_kep);
    
    std::cout << "Elementi Equinoziali (EQU):\n";
    std::cout << "  a      = " << std::fixed << std::setprecision(10) << orbfit_eq.a << " AU\n";
    std::cout << "  h      = " << orbfit_eq.h << "\n";
    std::cout << "  k      = " << orbfit_eq.k << "\n";
    std::cout << "  p      = " << orbfit_eq.p << "\n";
    std::cout << "  q      = " << orbfit_eq.q << "\n";
    std::cout << "  λ      = " << orbfit_eq.lambda << " rad\n\n";
    
    std::cout << "Elementi Kepleriani (ECL J2000):\n";
    std::cout << "  a = " << orbfit_kep.semi_major_axis << " AU\n";
    std::cout << "  e = " << orbfit_kep.eccentricity << "\n";
    std::cout << "  i = " << std::setprecision(6) << orbfit_kep.inclination * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  Ω = " << orbfit_kep.longitude_ascending_node * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  ω = " << orbfit_kep.argument_perihelion * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  M = " << orbfit_kep.mean_anomaly * constants::RAD_TO_DEG << " deg\n\n";
    
    std::cout << "Vettore di Stato (ECL J2000):\n";
    std::cout << std::setprecision(12);
    std::cout << "  r = [" << orbfit_cart.position[0] << ", " 
              << orbfit_cart.position[1] << ", " << orbfit_cart.position[2] << "] AU\n";
    std::cout << "  v = [" << orbfit_cart.velocity[0] << ", " 
              << orbfit_cart.velocity[1] << ", " << orbfit_cart.velocity[2] << "] AU/day\n";
    
    // Convert to equatorial
    auto R_ecl_to_eq = coordinates::ReferenceFrame::ecliptic_to_j2000().transpose();
    Eigen::Vector3d orbfit_pos_eq = R_ecl_to_eq * orbfit_cart.position;
    double orbfit_ra = std::atan2(orbfit_pos_eq[1], orbfit_pos_eq[0]);
    if (orbfit_ra < 0) orbfit_ra += 2 * constants::PI;
    double orbfit_dec = std::asin(orbfit_pos_eq[2] / orbfit_pos_eq.norm());
    
    std::cout << "\nCoordinate Equatoriali (EQU J2000):\n";
    std::cout << "  RA  = " << std::setprecision(8) << orbfit_ra * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  Dec = " << orbfit_dec * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  r   = " << orbfit_pos_eq.norm() << " AU\n";
    
    // ========================================================================
    // 2. ASTDYN - Esegui differential correction
    // ========================================================================
    print_header("2. ASTDYN - Differential Correction");
    
    // Load initial elements
    auto initial_eq = parse_eq1_file("astdyn/tools/203_astdys.eq1");
    auto initial_kep = equinoctial_to_keplerian(initial_eq);
    auto initial_cart = keplerian_to_cartesian(initial_kep);
    
    std::cout << "Elementi iniziali: MJD " << initial_eq.mjd_tdt << "\n";
    std::cout << "Target epoca:      MJD " << orbfit_eq.mjd_tdt << "\n\n";
    
    // Load observations
    auto observations = RWOReader::readFile("astdyn/tools/203_recent100.rwo");
    std::cout << "Osservazioni caricate: " << observations.size() << "\n\n";
    
    // Setup differential corrector
    auto ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
    auto integrator = std::make_unique<propagation::RK4Integrator>(0.1);
    propagation::PropagatorSettings prop_settings;
    prop_settings.include_planets = true;
    prop_settings.perturb_venus = true;
    prop_settings.perturb_earth = true;
    prop_settings.perturb_mars = true;
    prop_settings.perturb_jupiter = true;
    prop_settings.perturb_saturn = true;
    
    auto propagator = std::make_shared<propagation::Propagator>(std::move(integrator), ephemeris, prop_settings);
    auto residual_calc = std::make_shared<ResidualCalculator>(ephemeris, propagator);
    auto stm_computer = std::make_shared<StateTransitionMatrix>(propagator);
    
    DifferentialCorrector corrector(residual_calc, stm_computer);
    
    DifferentialCorrectorSettings settings;
    settings.max_iterations = 20;
    settings.convergence_tolerance = 1.0e-6;
    settings.outlier_sigma = 3.0;
    settings.reject_outliers = true;
    settings.verbose = false;
    
    std::cout << "Running differential correction...\n";
    auto result = corrector.fit(observations, initial_cart, settings);
    
    std::cout << "  Converged: " << (result.converged ? "YES" : "NO") << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  Final RMS: " << std::setprecision(3) << result.statistics.rms_total << " arcsec\n";
    std::cout << "  Observations: " << result.statistics.num_observations << "\n";
    std::cout << "  Outliers: " << result.statistics.num_outliers << "\n\n";
    
    // Convert AstDyn result to Keplerian and Equinoctial
    auto astdyn_kep = cartesian_to_keplerian(result.final_state);
    
    // Convert to equinoctial
    double e = astdyn_kep.eccentricity;
    double tan_half_i = std::tan(astdyn_kep.inclination / 2.0);
    double omega_plus_Omega = astdyn_kep.argument_perihelion + astdyn_kep.longitude_ascending_node;
    
    EquinoctialElements astdyn_eq;
    astdyn_eq.a = astdyn_kep.semi_major_axis;
    astdyn_eq.h = e * std::sin(omega_plus_Omega);
    astdyn_eq.k = e * std::cos(omega_plus_Omega);
    astdyn_eq.p = tan_half_i * std::sin(astdyn_kep.longitude_ascending_node);
    astdyn_eq.q = tan_half_i * std::cos(astdyn_kep.longitude_ascending_node);
    astdyn_eq.lambda = astdyn_kep.mean_anomaly + omega_plus_Omega;
    astdyn_eq.mjd_tdt = result.final_state.epoch_mjd_tdb;
    
    std::cout << "Elementi Equinoziali (EQU):\n";
    std::cout << "  a      = " << std::fixed << std::setprecision(10) << astdyn_eq.a << " AU\n";
    std::cout << "  h      = " << astdyn_eq.h << "\n";
    std::cout << "  k      = " << astdyn_eq.k << "\n";
    std::cout << "  p      = " << astdyn_eq.p << "\n";
    std::cout << "  q      = " << astdyn_eq.q << "\n";
    std::cout << "  λ      = " << astdyn_eq.lambda << " rad\n\n";
    
    std::cout << "Elementi Kepleriani (ECL J2000):\n";
    std::cout << "  a = " << astdyn_kep.semi_major_axis << " AU\n";
    std::cout << "  e = " << astdyn_kep.eccentricity << "\n";
    std::cout << "  i = " << std::setprecision(6) << astdyn_kep.inclination * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  Ω = " << astdyn_kep.longitude_ascending_node * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  ω = " << astdyn_kep.argument_perihelion * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  M = " << astdyn_kep.mean_anomaly * constants::RAD_TO_DEG << " deg\n\n";
    
    std::cout << "Vettore di Stato (ECL J2000):\n";
    std::cout << std::setprecision(12);
    std::cout << "  r = [" << result.final_state.position[0] << ", " 
              << result.final_state.position[1] << ", " << result.final_state.position[2] << "] AU\n";
    std::cout << "  v = [" << result.final_state.velocity[0] << ", " 
              << result.final_state.velocity[1] << ", " << result.final_state.velocity[2] << "] AU/day\n";
    
    // Convert to equatorial
    Eigen::Vector3d astdyn_pos_eq = R_ecl_to_eq * result.final_state.position;
    double astdyn_ra = std::atan2(astdyn_pos_eq[1], astdyn_pos_eq[0]);
    if (astdyn_ra < 0) astdyn_ra += 2 * constants::PI;
    double astdyn_dec = std::asin(astdyn_pos_eq[2] / astdyn_pos_eq.norm());
    
    std::cout << "\nCoordinate Equatoriali (EQU J2000):\n";
    std::cout << "  RA  = " << std::setprecision(8) << astdyn_ra * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  Dec = " << astdyn_dec * constants::RAD_TO_DEG << " deg\n";
    std::cout << "  r   = " << astdyn_pos_eq.norm() << " AU\n";
    
    // ========================================================================
    // 3. HORIZONS - Dati di riferimento (se disponibili)
    // ========================================================================
    print_header("3. JPL HORIZONS - Elementi di riferimento");
    std::cout << "Nota: Horizons usa elementi osculanti, non fittati su osservazioni.\n";
    std::cout << "      I valori potrebbero differire leggermente.\n\n";
    std::cout << "(Dati da inserire manualmente da JPL Horizons System)\n\n";
    
    // ========================================================================
    // 4. CONFRONTO - Differenze tra OrbFit e AstDyn
    // ========================================================================
    print_header("4. CONFRONTO ORBFIT vs ASTDYN");
    
    std::cout << "ELEMENTI EQUINOZIALI:\n";
    std::cout << "  Parametro   OrbFit          AstDyn          Differenza      Diff [m]\n";
    std::cout << "  ─────────────────────────────────────────────────────────────────────\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  a      " << orbfit_eq.a << "  " << astdyn_eq.a << "  " 
              << std::setprecision(6) << std::scientific << (astdyn_eq.a - orbfit_eq.a) 
              << "  " << (astdyn_eq.a - orbfit_eq.a) * constants::AU * 1000 << " m\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  h      " << orbfit_eq.h << "  " << astdyn_eq.h << "  " 
              << std::scientific << (astdyn_eq.h - orbfit_eq.h) << "\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  k      " << orbfit_eq.k << "  " << astdyn_eq.k << "  " 
              << std::scientific << (astdyn_eq.k - orbfit_eq.k) << "\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  p      " << orbfit_eq.p << "  " << astdyn_eq.p << "  " 
              << std::scientific << (astdyn_eq.p - orbfit_eq.p) << "\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  q      " << orbfit_eq.q << "  " << astdyn_eq.q << "  " 
              << std::scientific << (astdyn_eq.q - orbfit_eq.q) << "\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  λ      " << orbfit_eq.lambda << "  " << astdyn_eq.lambda << "  " 
              << std::scientific << (astdyn_eq.lambda - orbfit_eq.lambda) << " rad\n\n";
    
    std::cout << "VETTORE DI STATO (ECL J2000):\n";
    std::cout << "  Componente  OrbFit          AstDyn          Differenza [AU]  Diff [m]\n";
    std::cout << "  ─────────────────────────────────────────────────────────────────────\n";
    std::cout << std::fixed << std::setprecision(12);
    double dr_x = result.final_state.position[0] - orbfit_cart.position[0];
    double dr_y = result.final_state.position[1] - orbfit_cart.position[1];
    double dr_z = result.final_state.position[2] - orbfit_cart.position[2];
    std::cout << "  r_x    " << orbfit_cart.position[0] << "  " << result.final_state.position[0] 
              << "  " << std::scientific << std::setprecision(6) << dr_x << "  " << dr_x * constants::AU * 1000 << " m\n";
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "  r_y    " << orbfit_cart.position[1] << "  " << result.final_state.position[1] 
              << "  " << std::scientific << std::setprecision(6) << dr_y << "  " << dr_y * constants::AU * 1000 << " m\n";
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "  r_z    " << orbfit_cart.position[2] << "  " << result.final_state.position[2] 
              << "  " << std::scientific << std::setprecision(6) << dr_z << "  " << dr_z * constants::AU * 1000 << " m\n";
    
    double dr_norm = std::sqrt(dr_x*dr_x + dr_y*dr_y + dr_z*dr_z);
    std::cout << "\n  |Δr| totale: " << dr_norm << " AU = " << dr_norm * constants::AU * 1000 << " m\n\n";
    
    double dv_x = result.final_state.velocity[0] - orbfit_cart.velocity[0];
    double dv_y = result.final_state.velocity[1] - orbfit_cart.velocity[1];
    double dv_z = result.final_state.velocity[2] - orbfit_cart.velocity[2];
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "  v_x    " << orbfit_cart.velocity[0] << "  " << result.final_state.velocity[0] 
              << "  " << std::scientific << std::setprecision(6) << dv_x << " AU/day\n";
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "  v_y    " << orbfit_cart.velocity[1] << "  " << result.final_state.velocity[1] 
              << "  " << std::scientific << std::setprecision(6) << dv_y << " AU/day\n";
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "  v_z    " << orbfit_cart.velocity[2] << "  " << result.final_state.velocity[2] 
              << "  " << std::scientific << std::setprecision(6) << dv_z << " AU/day\n";
    
    double dv_norm = std::sqrt(dv_x*dv_x + dv_y*dv_y + dv_z*dv_z);
    std::cout << "\n  |Δv| totale: " << dv_norm << " AU/day = " << dv_norm * constants::AU * 1000 / 86400 << " m/s\n\n";
    
    std::cout << "COORDINATE EQUATORIALI (EQU J2000):\n";
    std::cout << "  Coord   OrbFit          AstDyn          Differenza       Diff [arcsec]\n";
    std::cout << "  ─────────────────────────────────────────────────────────────────────\n";
    std::cout << std::fixed << std::setprecision(8);
    double dra = astdyn_ra - orbfit_ra;
    double ddec = astdyn_dec - orbfit_dec;
    std::cout << "  RA     " << orbfit_ra * constants::RAD_TO_DEG << " deg  " 
              << astdyn_ra * constants::RAD_TO_DEG << " deg  " 
              << std::setprecision(6) << std::scientific << dra << " rad  " 
              << dra * constants::RAD_TO_ARCSEC << " arcsec\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "  Dec    " << orbfit_dec * constants::RAD_TO_DEG << " deg  " 
              << astdyn_dec * constants::RAD_TO_DEG << " deg  " 
              << std::setprecision(6) << std::scientific << ddec << " rad  " 
              << ddec * constants::RAD_TO_ARCSEC << " arcsec\n\n";
    
    // ========================================================================
    // 5. CONCLUSIONI
    // ========================================================================
    print_header("5. CONCLUSIONI");
    
    std::cout << "Qualità del fit:\n";
    std::cout << "  ✓ Convergenza:     " << (result.converged ? "OK" : "FALLITA") << "\n";
    std::cout << "  ✓ RMS residui:     " << std::fixed << std::setprecision(3) 
              << result.statistics.rms_total << " arcsec (target: <1 arcsec)\n";
    std::cout << "  ✓ Iterazioni:      " << result.iterations << " (max: 20)\n";
    std::cout << "  ✓ Outliers:        " << result.statistics.num_outliers << "/" 
              << observations.size() << " (" << std::setprecision(1) 
              << (100.0 * result.statistics.num_outliers / observations.size()) << "%)\n\n";
    
    std::cout << "Accordo OrbFit vs AstDyn:\n";
    std::cout << "  • Posizione: " << std::scientific << std::setprecision(3) 
              << dr_norm * constants::AU * 1000 << " m  ";
    if (dr_norm * constants::AU * 1000 < 1000) std::cout << "✓ ECCELLENTE";
    else if (dr_norm * constants::AU * 1000 < 10000) std::cout << "✓ BUONO";
    else std::cout << "⚠ VERIFICARE";
    std::cout << "\n";
    
    std::cout << "  • Velocità:  " << std::setprecision(3) 
              << dv_norm * constants::AU * 1000 / 86400 << " m/s  ";
    if (dv_norm * constants::AU * 1000 / 86400 < 0.001) std::cout << "✓ ECCELLENTE";
    else if (dv_norm * constants::AU * 1000 / 86400 < 0.01) std::cout << "✓ BUONO";
    else std::cout << "⚠ VERIFICARE";
    std::cout << "\n";
    
    std::cout << "  • RA/Dec:    " << std::setprecision(3) 
              << std::sqrt(dra*dra + ddec*ddec) * constants::RAD_TO_ARCSEC << " arcsec  ";
    if (std::sqrt(dra*dra + ddec*ddec) * constants::RAD_TO_ARCSEC < 1.0) std::cout << "✓ ECCELLENTE";
    else if (std::sqrt(dra*dra + ddec*ddec) * constants::RAD_TO_ARCSEC < 10.0) std::cout << "✓ BUONO";
    else std::cout << "⚠ VERIFICARE";
    std::cout << "\n\n";
    
    // Assertions
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.statistics.rms_total, 1.0);
    EXPECT_LT(dr_norm * constants::AU * 1000, 10000);  // < 10 km
    EXPECT_LT(std::sqrt(dra*dra + ddec*ddec) * constants::RAD_TO_ARCSEC, 10.0);  // < 10 arcsec
}
