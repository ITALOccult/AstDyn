/**
 * @file test_chebyshev_approximation.cpp
 * @brief Test per la classe ChebyshevApproximation
 * @author ITALOccultLibrary Development Team
 * @date 4 Dicembre 2025
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Core>

#include "italoccultlib/chebyshev_approximation.h"
#include "italoccultlib/astdyn_wrapper.h"
#include "integration/italoccult_integration.h"

using namespace ioccultcalc;

void printTestHeader(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "  " << title << std::endl;
    std::cout << std::string(60, '=') << std::endl;
}

void printTestResult(const std::string& name, bool passed) {
    std::cout << (passed ? "✓ PASS" : "✗ FAIL") << ": " << name << std::endl;
}

int main() {
    int passed = 0;
    int failed = 0;
    
    printTestHeader("Test ChebyshevApproximation");
    
    try {
        // =====================================================
        // TEST 1: Costruzione e stato iniziale
        // =====================================================
        std::cout << "\n[TEST 1] Costruzione e stato iniziale..." << std::endl;
        
        ChebyshevApproximation chebyshev(8);
        
        if (!chebyshev.isFitted()) {
            std::cout << "✓ Stato iniziale non-fittato: OK\n";
            passed++;
        } else {
            std::cout << "✗ Stato iniziale sbagliato\n";
            failed++;
        }
        
        if (chebyshev.getNumCoefficients() == 8) {
            std::cout << "✓ Numero coefficienti: 8\n";
            passed++;
        } else {
            std::cout << "✗ Numero coefficienti sbagliato\n";
            failed++;
        }
        
        // =====================================================
        // TEST 2: Caricamento dati da AstDyn
        // =====================================================
        printTestHeader("TEST 2: Caricamento dati da AstDyn");
        
        std::cout << "Caricamento asteroide 17030 (Sierks)..." << std::endl;
        
        // Crea wrapper AstDyn
        AstDynWrapper wrapper(PropagationSettings::highAccuracy());
        if (!wrapper.loadFromEQ1File("astdyn/data/17030.eq1")) {
            std::cout << "✗ Errore caricamento .eq1\n";
            return 1;
        }
        std::cout << "✓ Asteroide caricato\n";
        
        // Propaga a 5 epoche per il fitting - usa il wrapper direttamente
        std::vector<double> epochs = {61000.0, 61001.0, 61003.0, 61007.0, 61014.0};
        std::vector<Eigen::Vector3d> positions;
        
        for (double epoch : epochs) {
            positions.push_back(wrapper.propagateToEpoch(epoch).position);
        }
        
        std::cout << "✓ Propagate " << positions.size() << " epoche\n";
        
        // =====================================================
        // TEST 3: Fitting dei polinomi
        // =====================================================
        printTestHeader("TEST 3: Fitting dei polinomi");
        
        std::cout << "Fitting Chebyshev con 8 coefficienti..." << std::endl;
        
        if (chebyshev.fit(positions, 61000.0, 61014.0)) {
            std::cout << "✓ Fitting completato\n";
            passed++;
        } else {
            std::cout << "✗ Fitting fallito\n";
            failed++;
        }
        
        if (chebyshev.isFitted()) {
            std::cout << "✓ Stato fittato confermato\n";
            passed++;
        } else {
            std::cout << "✗ Stato fittato non aggiornato\n";
            failed++;
        }
        
        auto time_interval = chebyshev.getTimeInterval();
        std::cout << "  Intervallo: MJD " << std::fixed << std::setprecision(2)
                  << time_interval.first << " - " << time_interval.second << std::endl;
        
        // =====================================================
        // TEST 4: Valutazione della posizione
        // =====================================================
        printTestHeader("TEST 4: Valutazione della posizione");
        
        // Valuta in un'epoca intermedia
        double test_epoch = 61005.0;
        Eigen::Vector3d approx_pos = chebyshev.evaluatePosition(test_epoch);
        
        std::cout << "Posizione approssimata @ MJD " << test_epoch << ":\n";
        std::cout << "  X = " << std::setprecision(9) << approx_pos.x() << " AU\n";
        std::cout << "  Y = " << std::setprecision(9) << approx_pos.y() << " AU\n";
        std::cout << "  Z = " << std::setprecision(9) << approx_pos.z() << " AU\n";
        
        double distance = approx_pos.norm();
        std::cout << "  Distanza: " << distance << " AU\n";
        
        if (distance > 2.5 && distance < 3.5) {
            std::cout << "✓ Distanza ragionevole\n";
            passed++;
        } else {
            std::cout << "✗ Distanza non ragionevole\n";
            failed++;
        }
        
        // =====================================================
        // TEST 5: Valutazione della velocità
        // =====================================================
        printTestHeader("TEST 5: Valutazione della velocità");
        
        Eigen::Vector3d approx_vel = chebyshev.evaluateVelocity(test_epoch);
        
        std::cout << "Velocità approssimata @ MJD " << test_epoch << ":\n";
        std::cout << "  VX = " << std::setprecision(9) << approx_vel.x() << " AU/day\n";
        std::cout << "  VY = " << std::setprecision(9) << approx_vel.y() << " AU/day\n";
        std::cout << "  VZ = " << std::setprecision(9) << approx_vel.z() << " AU/day\n";
        
        double speed = approx_vel.norm();
        std::cout << "  Velocità: " << speed << " AU/day\n";
        
        if (speed > 0.004 && speed < 0.02) {
            std::cout << "✓ Velocità ragionevole\n";
            passed++;
        } else {
            std::cout << "✗ Velocità non ragionevole\n";
            failed++;
        }
        
        // =====================================================
        // TEST 6: Errore di approssimazione
        // =====================================================
        printTestHeader("TEST 6: Errore di approssimazione");
        
        Eigen::Vector3d error = chebyshev.getApproximationError();
        
        std::cout << "Errore RMS dell'approssimazione:\n";
        std::cout << "  X: " << std::scientific << std::setprecision(3) << error.x() << " AU\n";
        std::cout << "  Y: " << std::scientific << error.y() << " AU\n";
        std::cout << "  Z: " << std::scientific << error.z() << " AU\n";
        std::cout << "  Media: " << std::scientific << error.norm() / 3.0 << " AU\n";
        
        double avg_error = error.norm() / 3.0;
        if (avg_error < 0.001) {  // Meno di 0.001 AU
            std::cout << "✓ Errore basso (< 0.001 AU)\n";
            passed++;
        } else if (avg_error < 0.01) {
            std::cout << "✓ Errore ragionevole (< 0.01 AU)\n";
            passed++;
        } else {
            std::cout << "✗ Errore troppo alto\n";
            failed++;
        }
        
        // =====================================================
        // TEST 7: Energia orbitale
        // =====================================================
        printTestHeader("TEST 7: Energia orbitale");
        
        double energy = chebyshev.evaluateEnergy(test_epoch);
        
        std::cout << "Energia specifica @ MJD " << test_epoch << ":\n";
        std::cout << "  E = " << std::scientific << std::setprecision(6) 
                  << energy << " AU²/day²\n";
        
        if (energy < 0) {  // Orbita ellittica
            std::cout << "✓ Orbita ellittica (E < 0)\n";
            passed++;
        } else {
            std::cout << "✗ Energia inaspettata\n";
            failed++;
        }
        
        // =====================================================
        // TEST 8: Momento angolare
        // =====================================================
        printTestHeader("TEST 8: Momento angolare");
        
        Eigen::Vector3d h = chebyshev.evaluateAngularMomentum(test_epoch);
        
        std::cout << "Momento angolare @ MJD " << test_epoch << ":\n";
        std::cout << "  HX = " << std::setprecision(6) << h.x() << " AU²/day\n";
        std::cout << "  HY = " << std::setprecision(6) << h.y() << " AU²/day\n";
        std::cout << "  HZ = " << std::setprecision(6) << h.z() << " AU²/day\n";
        std::cout << "  |H| = " << h.norm() << " AU²/day\n";
        
        if (h.norm() > 0) {
            std::cout << "✓ Momento angolare valido\n";
            passed++;
        } else {
            std::cout << "✗ Momento angolare nullo\n";
            failed++;
        }
        
        // =====================================================
        // TEST 9: Salvataggio e caricamento
        // =====================================================
        printTestHeader("TEST 9: Salvataggio e caricamento");
        
        std::string coeff_file = "chebyshev_coeffs_test.txt";
        
        if (chebyshev.saveToFile(coeff_file)) {
            std::cout << "✓ Coefficienti salvati: " << coeff_file << "\n";
            passed++;
        } else {
            std::cout << "✗ Errore salvataggio\n";
            failed++;
        }
        
        // Crea nuovo oggetto e carica
        ChebyshevApproximation chebyshev2(8);
        if (chebyshev2.loadFromFile(coeff_file)) {
            std::cout << "✓ Coefficienti caricati\n";
            passed++;
            
            // Verifica coerenza
            Eigen::Vector3d pos_reload = chebyshev2.evaluatePosition(test_epoch);
            double diff = (approx_pos - pos_reload).norm();
            
            std::cout << "  Differenza: " << std::scientific << std::setprecision(3) 
                      << diff << " AU\n";
            
            if (diff < 1e-10) {
                std::cout << "✓ Coerenza perfetta\n";
                passed++;
            } else {
                std::cout << "✗ Discrepanza nei dati\n";
                failed++;
            }
        } else {
            std::cout << "✗ Errore caricamento\n";
            failed++;
        }
        
        // =====================================================
        // TEST 10: Statistica
        // =====================================================
        printTestHeader("TEST 10: Statistica");
        
        std::string stats = chebyshev.getStatistics();
        std::cout << stats << std::endl;
        
        if (!stats.empty() && stats.find("FITTATO") != std::string::npos) {
            std::cout << "✓ Statistiche valide\n";
            passed++;
        } else {
            std::cout << "✗ Statistiche incomplete\n";
            failed++;
        }
        
        // =====================================================
        // RIEPILOGO
        // =====================================================
        printTestHeader("RIEPILOGO RISULTATI");
        
        std::cout << "\nTest passati: " << passed << std::endl;
        std::cout << "Test falliti: " << failed << std::endl;
        std::cout << "Totale: " << (passed + failed) << std::endl;
        
        double success_rate = 100.0 * passed / (passed + failed);
        std::cout << "Percentuale di successo: " << std::fixed << std::setprecision(1) 
                  << success_rate << "%\n" << std::endl;
        
        if (failed == 0) {
            std::cout << "✓ TUTTI I TEST PASSATI!\n";
            std::cout << std::string(60, '=') << std::endl << std::endl;
            return 0;
        } else {
            std::cout << "✗ Alcuni test non hanno passato\n";
            std::cout << std::string(60, '=') << std::endl << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Eccezione: " << e.what() << std::endl;
        return 1;
    }
}
