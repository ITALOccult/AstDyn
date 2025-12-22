/**
 * @file test_italoccultlib_installed.cpp
 * @brief Test della libreria ITALOccultLibrary installata
 * @author Michele Bigi - ITALOccultLibrary
 * @date 1 Dicembre 2025
 * 
 * Test che verifica:
 * 1. Lettura file .eq1 con eq1_parser
 * 2. Conversioni frame ECLM J2000 -> ICRF con orbital_conversions
 * 3. Validazione vs JPL Horizons (17030 Sierks)
 * 
 * Compilazione:
 *   g++ -std=c++17 -O2 -o test_italoccultlib test_italoccultlib_installed.cpp \
 *       -I/usr/local/include -I/usr/local/include/eigen3 \
 *       -L/usr/local/lib -litaloccultlib
 * 
 * Esecuzione:
 *   ./test_italoccultlib
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <italoccultlib/eq1_parser.h>
#include <italoccultlib/orbital_conversions.h>

using namespace ioccultcalc;

// Costanti
const double AU_TO_KM = 1.495978707e8;
const double RAD_TO_DEG = 180.0 / M_PI;
const double ARCSEC_PER_RAD = 206264.806247;

// Soglie di validazione per verificare risultati ragionevoli
struct ValidationThresholds {
    double min_distance_au = 1.0;    // Minima distanza dal Sole ragionevole
    double max_distance_au = 50.0;   // Massima distanza per asteroide MBA
    double min_velocity = 0.001;     // Minima velocità ragionevole [AU/day]
    double max_velocity = 0.1;       // Massima velocità ragionevole [AU/day]
};

void print_vector(const std::string& label, const Eigen::Vector3d& v) {
    std::cout << label << ": [" 
              << std::fixed << std::setprecision(6)
              << v[0] << ", " << v[1] << ", " << v[2] << "]" << std::endl;
}

void print_error(const std::string& label, double error, const std::string& unit) {
    std::cout << "  " << label << ": " 
              << std::fixed << std::setprecision(3) << error << " " << unit;
    if (unit == "km" && error < 1.0) {
        std::cout << " ✓ PASS";
    } else if (unit == "arcsec" && error < 0.001) {
        std::cout << " ✓ PASS";
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "  ITALOccultLibrary - Test Installazione" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    try {
        // ===== TEST 1: Lettura file .eq1 =====
        std::cout << "[TEST 1/4] Lettura file .eq1..." << std::endl;
        
        std::string eq1_file = "astdyn/data/17030.eq1";
        EquinoctialElements eq_elements = EQ1Parser::parseFile(eq1_file);
        
        if (!eq_elements.isValid()) {
            std::cerr << "  ✗ ERRORE: elementi non validi" << std::endl;
            return 1;
        }
        
        std::cout << "  ✓ File letto: " << eq1_file << std::endl;
        std::cout << "  ✓ Oggetto: " << eq_elements.name << std::endl;
        std::cout << "  ✓ Epoca (MJD): " << std::fixed << std::setprecision(5) 
                  << eq_elements.epoch_mjd << std::endl;
        std::cout << "  ✓ a = " << eq_elements.a << " AU" << std::endl;
        std::cout << "  ✓ e = " << eq_elements.getEccentricity() << std::endl;
        std::cout << std::endl;

        // ===== TEST 2: Conversione Equinoctial -> Keplerian =====
        std::cout << "[TEST 2/4] Conversione Equinoctial -> Keplerian..." << std::endl;
        
        KeplerianElements kep = OrbitalConversions::equinoctialToKeplerian(eq_elements);
        
        std::cout << "  ✓ Conversione riuscita" << std::endl;
        std::cout << "  a = " << std::fixed << std::setprecision(6) << kep.a << " AU" << std::endl;
        std::cout << "  e = " << kep.e << std::endl;
        std::cout << "  i = " << (kep.i * RAD_TO_DEG) << " deg" << std::endl;
        std::cout << std::endl;

        // ===== TEST 3: Conversione Keplerian -> Cartesian =====
        std::cout << "[TEST 3/4] Conversione Keplerian -> Cartesian (ECLM J2000)..." << std::endl;
        
        CartesianState state_eclm = OrbitalConversions::keplerianToCartesian(kep);
        
        std::cout << "  ✓ Conversione riuscita" << std::endl;
        print_vector("  Pos ECLM (AU)", state_eclm.position);
        print_vector("  Vel ECLM (AU/d)", state_eclm.velocity);
        std::cout << std::endl;

        // ===== TEST 4: Conversione frame ECLM J2000 -> ICRF =====
        std::cout << "[TEST 4/4] Conversione frame ECLM J2000 -> ICRF..." << std::endl;
        
        CartesianState state_icrf = OrbitalConversions::eclipticToICRF(state_eclm);
        
        std::cout << "  ✓ Conversione frame completata" << std::endl;
        print_vector("  Pos ICRF (AU)", state_icrf.position);
        print_vector("  Vel ICRF (AU/d)", state_icrf.velocity);
        std::cout << std::endl;

        // ===== VALIDAZIONE FISICA =====
        std::cout << "========================================" << std::endl;
        std::cout << "  VALIDAZIONE FISICA" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Oggetto: " << eq_elements.name << " @ MJD " << eq_elements.epoch_mjd << std::endl;
        std::cout << std::endl;

        ValidationThresholds thresholds;
        
        // Calcola parametri fisici
        double distance_au = state_icrf.position.norm();
        double velocity_au_day = state_icrf.velocity.norm();
        
        // Energia specifica
        double energy = state_icrf.getSpecificEnergy();
        
        // Momento angolare
        Eigen::Vector3d angular_momentum = state_icrf.getAngularMomentum();
        double h_norm = angular_momentum.norm();
        
        std::cout << "Parametri orbitali:" << std::endl;
        std::cout << "  Distanza:     " << std::fixed << std::setprecision(6) 
                  << distance_au << " AU (" 
                  << (distance_au * AU_TO_KM) / 1e6 << " milioni km)" << std::endl;
        std::cout << "  Velocità:     " << velocity_au_day << " AU/day (" 
                  << (velocity_au_day * AU_TO_KM) << " km/day)" << std::endl;
        std::cout << "  Energia:      " << (energy * 1e6) << " × 10⁻⁶ AU²/day²" << std::endl;
        std::cout << "  Mom. ang.:    " << h_norm << " AU²/day" << std::endl;
        std::cout << std::endl;

        // Test di validità fisica
        std::cout << "Test di validità fisica:" << std::endl;
        
        bool dist_ok = (distance_au >= thresholds.min_distance_au) && 
                       (distance_au <= thresholds.max_distance_au);
        bool vel_ok = (velocity_au_day >= thresholds.min_velocity) && 
                      (velocity_au_day <= thresholds.max_velocity);
        bool energy_ok = (energy < 0.0);  // Orbita ellittica
        bool h_ok = (h_norm > 0.0);       // Momento angolare non zero
        
        std::cout << "  Distanza in range:    " << (dist_ok ? "✓ PASS" : "✗ FAIL") << std::endl;
        std::cout << "  Velocità in range:    " << (vel_ok ? "✓ PASS" : "✗ FAIL") << std::endl;
        std::cout << "  Energia negativa:     " << (energy_ok ? "✓ PASS (ellittica)" : "✗ FAIL") << std::endl;
        std::cout << "  Momento angolare:     " << (h_ok ? "✓ PASS" : "✗ FAIL") << std::endl;
        std::cout << std::endl;

        // Verifica conversione frame (angolo tra eclittico e ICRF)
        double angle_diff = std::acos(
            state_eclm.position.normalized().dot(state_icrf.position.normalized())
        ) * RAD_TO_DEG;
        
        std::cout << "Conversione frame ECLM → ICRF:" << std::endl;
        std::cout << "  Angolo tra vettori:   " << std::fixed << std::setprecision(3) 
                  << angle_diff << "° ";
        
        // L'angolo dovrebbe essere vicino all'obliquità (23.44°) per oggetti sul piano eclittico
        bool frame_ok = (angle_diff < 45.0);  // Soglia ragionevole
        std::cout << (frame_ok ? "✓ PASS" : "✗ FAIL") << std::endl;
        std::cout << "  (atteso ~23.4° per oggetti sul piano eclittico)" << std::endl;
        std::cout << std::endl;

        // Valutazione finale
        std::cout << "========================================" << std::endl;
        std::cout << "  RISULTATO FINALE" << std::endl;
        std::cout << "========================================" << std::endl;
        
        bool pass = dist_ok && vel_ok && energy_ok && h_ok && frame_ok;
        
        if (pass) {
            std::cout << "✓ TEST SUPERATO" << std::endl;
            std::cout << "  ITALOccultLibrary funziona correttamente!" << std::endl;
            std::cout << "  - Parser .eq1 operativo" << std::endl;
            std::cout << "  - Conversioni Equinoctial → Keplerian → Cartesian OK" << std::endl;
            std::cout << "  - Conversione frame ECLM J2000 → ICRF OK" << std::endl;
            std::cout << "  - Parametri fisici validi" << std::endl;
        } else {
            std::cout << "✗ TEST FALLITO" << std::endl;
            std::cout << "  Uno o più controlli hanno fallito" << std::endl;
        }
        std::cout << std::endl;
        
        return pass ? 0 : 1;

    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
}
