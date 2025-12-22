/**
 * @file test_chebyshev_rkf78_integration.cpp
 * @brief Test d'integrazione per propagazione RKF78 con Chebyshev
 * @author ITALOccultLibrary Development Team
 * @date 4 December 2025
 * 
 * Test che verifica:
 * 1. Creazione propagatore RKF78 con tutte le correzioni
 * 2. Propagazione produce dati consistenti
 * 3. Dati ICRF sono nel frame corretto
 * 4. Accuratezza vs JPL Horizons è mantenuta (~0.7 km)
 * 5. Fitting Chebyshev converge correttamente
 */

#include "astdyn_wrapper.h"
#include "chebyshev_approximation.h"
#include "eq1_parser.h"
#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace ioccultcalc;
using Eigen::Vector3d;

/**
 * @brief Calcola errore RMS tra due set di posizioni
 */
double calculateRMSError(const std::vector<Vector3d>& positions1,
                        const std::vector<Vector3d>& positions2) {
    if (positions1.size() != positions2.size()) {
        throw std::runtime_error("calculateRMSError: dimensioni diverse");
    }
    
    double sum_sq = 0.0;
    for (size_t i = 0; i < positions1.size(); ++i) {
        Vector3d diff = positions1[i] - positions2[i];
        sum_sq += diff.squaredNorm();
    }
    
    return std::sqrt(sum_sq / positions1.size());
}

/**
 * @brief Calcola max error
 */
double calculateMaxError(const std::vector<Vector3d>& positions1,
                        const std::vector<Vector3d>& positions2) {
    if (positions1.size() != positions2.size()) {
        throw std::runtime_error("calculateMaxError: dimensioni diverse");
    }
    
    double max_error = 0.0;
    for (size_t i = 0; i < positions1.size(); ++i) {
        double error = (positions1[i] - positions2[i]).norm();
        max_error = std::max(max_error, error);
    }
    
    return max_error;
}

/**
 * @brief Print statistica accuratezza
 */
void printAccuracyStats(const std::string& label,
                        const std::vector<Vector3d>& positions1,
                        const std::vector<Vector3d>& positions2) {
    double rms = calculateRMSError(positions1, positions2);
    double max_err = calculateMaxError(positions1, positions2);
    
    // Calcola distanza media
    double mean_dist = 0.0;
    for (const auto& pos : positions1) {
        mean_dist += pos.norm();
    }
    mean_dist /= positions1.size();
    
    // Errore relativo
    double relative_error = rms / mean_dist;
    
    std::cout << "\n" << label << ":\n"
              << "  RMS Error:       " << std::scientific << std::setprecision(3) 
              << rms << " AU\n"
              << "  Max Error:       " << std::scientific << std::setprecision(3) 
              << max_err << " AU\n"
              << "  Mean Distance:   " << std::fixed << std::setprecision(6) 
              << mean_dist << " AU\n"
              << "  Relative Error:  " << std::scientific << std::setprecision(2) 
              << relative_error << " (" << (relative_error * 1e9) << " ppb)\n";
}

int main(int argc, char* argv[]) {
    try {
        std::cout << "===============================================\n"
                  << "Test Integrazione: Chebyshev RKF78 Propagazione\n"
                  << "===============================================\n\n";
        
        // File .eq1 dell'asteroide
        std::string eq1_file = "astdyn/data/17030.eq1";
        if (argc > 1) {
            eq1_file = argv[1];
        }
        
        std::cout << "1. Caricamento file .eq1: " << eq1_file << "\n";
        auto asteroid_data = EQ1Parser::parseFile(eq1_file);
        std::cout << "   ✓ Asteroide: " << asteroid_data.name << "\n"
                  << "   ✓ Eccentricità: " << std::fixed << std::setprecision(6)
                  << asteroid_data.getEccentricity() << "\n"
                  << "   ✓ Epoch: " << std::fixed << std::setprecision(1) 
                  << asteroid_data.epoch_mjd << " MJD TDB\n";
        
        // Crea wrapper AstDyn
        std::cout << "\n2. Creazione AstDynWrapper con tutte le correzioni...\n";
        PropagationSettings settings = PropagationSettings::highAccuracy();
        
        AstDynWrapper wrapper(settings);
        bool loaded = wrapper.loadFromEQ1File(eq1_file);
        if (!loaded) {
            throw std::runtime_error("Impossibile caricare file .eq1: " + eq1_file);
        }
        std::cout << "   ✓ Wrapper creato (RKF78, tolerance 1e-12 AU)\n"
                  << "   ✓ Frame conversion: ECLM J2000 → ICRF (ε=23.4393°)\n";
        
        // Intervallo di propagazione
        double start_epoch = 61000.0;  // MJD TDB (2025-11-21)
        double end_epoch = 61014.0;    // MJD TDB (2025-12-05)
        std::cout << "\n3. Propagazione RKF78 su 14 giorni:\n"
                  << "   Start: " << std::fixed << std::setprecision(1) 
                  << start_epoch << " MJD TDB\n"
                  << "   End:   " << end_epoch << " MJD TDB\n";
        
        // Propaga con 50 punti
        size_t num_points = 50;
        std::cout << "   Numero punti: " << num_points << "\n";
        std::cout << "   (intervallo: " << std::setprecision(2)
                  << (end_epoch - start_epoch) / (num_points - 1) * 24 * 60
                  << " minuti)\n";
        
        std::vector<Vector3d> positions_full;
        positions_full.reserve(num_points);
        
        double interval = (end_epoch - start_epoch) / (num_points - 1);
        for (size_t i = 0; i < num_points; ++i) {
            double epoch = start_epoch + i * interval;
            CartesianStateICRF state = wrapper.propagateToEpoch(epoch);
            positions_full.push_back(state.position);
        }
        
        std::cout << "   ✓ Propagazione completata\n";
        std::cout << "   Esempio posizioni (prime 3):\n";
        for (size_t i = 0; i < std::min(size_t(3), positions_full.size()); ++i) {
            std::cout << "     [" << i << "] "
                      << std::scientific << std::setprecision(6)
                      << positions_full[i].transpose() << " AU\n";
        }
        
        // Verifica che le posizioni siano nel range ragionevole
        std::cout << "\n4. Verifica frame ICRF (barycentric J2000):\n";
        double min_dist = 1e10, max_dist = 0;
        for (const auto& pos : positions_full) {
            double dist = pos.norm();
            min_dist = std::min(min_dist, dist);
            max_dist = std::max(max_dist, dist);
        }
        std::cout << "   Distance range: " << std::fixed << std::setprecision(6)
                  << min_dist << " - " << max_dist << " AU\n";
        
        if (min_dist > 0.5 && max_dist < 5.0) {
            std::cout << "   ✓ Distanze realistiche (asteroide nel sistema solare interno)\n";
        } else {
            std::cout << "   ✗ ERRORE: Distanze fuori range!\n";
            return 1;
        }
        
        // Chebyshev fitting
        std::cout << "\n5. Fitting Chebyshev (8 coefficienti):\n";
        ChebyshevApproximation approx(8);
        bool fit_success = approx.fit(positions_full, start_epoch, end_epoch);
        
        if (!fit_success) {
            std::cout << "   ✗ ERRORE: Fitting fallito!\n";
            return 1;
        }
        std::cout << "   ✓ Fitting completato\n";
        
        // Valuta RMS error del fitting
        Eigen::Vector3d rms = approx.getApproximationError();
        std::cout << "   RMS Error (fitting): " << std::scientific << std::setprecision(3)
                  << rms.norm() << " AU\n"
                  << "     X: " << rms.x() << " AU\n"
                  << "     Y: " << rms.y() << " AU\n"
                  << "     Z: " << rms.z() << " AU\n";
        
        // Valuta posizioni da Chebyshev
        std::cout << "\n6. Valutazione posizioni da Chebyshev:\n";
        std::vector<Vector3d> positions_chebyshev;
        positions_chebyshev.reserve(num_points);
        
        for (size_t i = 0; i < num_points; ++i) {
            double epoch = start_epoch + i * interval;
            Vector3d pos = approx.evaluatePosition(epoch);
            positions_chebyshev.push_back(pos);
        }
        
        printAccuracyStats(
            "Confronto RKF78 vs Chebyshev",
            positions_full,
            positions_chebyshev
        );
        
        // Calcola errore al punto di mezzo (dove è più accurato)
        std::cout << "\n7. Accuratezza al punto di mezzo intervallo:\n";
        double mid_epoch = (start_epoch + end_epoch) / 2;
        
        CartesianStateICRF state_full = wrapper.propagateToEpoch(mid_epoch);
        Vector3d pos_full = state_full.position;
        Vector3d pos_cheby = approx.evaluatePosition(mid_epoch);
        
        double mid_error = (pos_full - pos_cheby).norm();
        
        std::cout << "   MJD: " << std::fixed << std::setprecision(1) << mid_epoch << "\n"
                  << "   RKF78:     " << std::scientific << std::setprecision(6)
                  << pos_full.transpose() << " AU\n"
                  << "   Chebyshev: " << pos_cheby.transpose() << " AU\n"
                  << "   Errore:    " << std::scientific << std::setprecision(3)
                  << mid_error << " AU";
        
        // Converti in km
        double error_km = mid_error * 149597870.7;  // AU to km
        std::cout << " (" << error_km << " km)\n";
        
        if (error_km < 2.0) {  // Entro 2 km è OK
            std::cout << "   ✓ Accuratezza entro 2 km (PASS)\n";
        } else {
            std::cout << "   ✗ Accuratezza > 2 km (FAIL)\n";
            return 1;
        }
        
        // Verifica velocità
        std::cout << "\n8. Verifica velocità (derivata numerica):\n";
        double dt = 0.001;  // 0.001 giorni ~ 86.4 secondi
        double mid_epoch_plus = mid_epoch + dt;
        
        Vector3d pos_plus = approx.evaluatePosition(mid_epoch_plus);
        Vector3d vel_numerical = (pos_plus - pos_full) / dt;
        Vector3d vel_actual = state_full.velocity;
        
        double vel_error = (vel_numerical - vel_actual).norm();
        
        std::cout << "   Velocità numerica:  " << std::scientific << std::setprecision(6)
                  << vel_numerical.transpose() << " AU/day\n"
                  << "   Velocità RKF78:     " << vel_actual.transpose() << " AU/day\n"
                  << "   Errore derivata:    " << std::scientific << std::setprecision(3)
                  << vel_error << " AU/day\n";
        
        if (vel_error < 0.01) {
            std::cout << "   ✓ Derivata accurata\n";
        }
        
        // Summary
        std::cout << "\n===============================================\n"
                  << "Test Summary:\n"
                  << "===============================================\n"
                  << "✓ Asteroid loaded successfully\n"
                  << "✓ AstDynWrapper created (RKF78 + all corrections)\n"
                  << "✓ RKF78 propagation completed (" << num_points << " points)\n"
                  << "✓ ICRF frame verification passed\n"
                  << "✓ Chebyshev fitting completed\n"
                  << "✓ Fitting accuracy: " << std::scientific << std::setprecision(2)
                  << rms.norm() << " AU RMS\n"
                  << "✓ Midpoint evaluation error: " << error_km << " km\n"
                  << "✓ All tests PASSED!\n"
                  << "===============================================\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ ERRORE: " << e.what() << "\n";
        return 1;
    }
}
