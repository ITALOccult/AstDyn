/**
 * @file propagation_demo.cpp
 * @brief Esempio di propagazione orbitale veloce e precisa
 * 
 * Questo esempio mostra come propagare elementi orbitali per 100-200 giorni
 * utilizzando AstDynWrapper e (opzionalmente) fitting di Chebyshev per
 * massimizzare la velocità di accesso ai dati.
 *
 * Compilazione:
 * g++ -std=c++17 -O3 -I../italoccultlibrary/include -I../astdyn/include \
 *     -o propagation_demo propagation_demo.cpp \
 *     -L../build -lioccult -lastdyn
 */

#include "astdyn_wrapper.h"
#include "chebyshev_approximation.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>

using namespace ioccultcalc;

int main() {
    // -------------------------------------------------------------------------
    // 1. Configurazione Iniziale
    // -------------------------------------------------------------------------
    std::cout << "=== Demo Propagazione Orbitale AstDyn ===\n\n";

    // Imposta elementi orbitali iniziali (es. un asteroide fittizio o reale)
    // Esempio: Cerere
    double a = 2.766;           // Semiasse maggiore [AU]
    double e = 0.076;           // Eccentricità
    double i = 10.59 * M_PI / 180.0; // Inclinazione [rad]
    double Omega = 80.3 * M_PI / 180.0; // Long. nodo ascendente [rad]
    double omega = 73.0 * M_PI / 180.0; // Argomento perielio [rad]
    double M = 0.0;             // Anomalia media [rad]
    double epoch = 60000.0;     // Epoca iniziale (MJD TDB)

    // Configura il propagatore per alta precisione (RKF78 + Perturbazioni)
    // Per 100-200 giorni, RKF78 è molto veloce e garantisce precisione ~km
    auto settings = PropagationSettings::highAccuracy();
    AstDynWrapper propagator(settings);

    std::cout << "Inizializzazione propagatore..." << std::endl;
    propagator.setKeplerianElements(a, e, i, Omega, omega, M, epoch, "Ceres-Test");

    // -------------------------------------------------------------------------
    // 2. Propagazione Numerica (La parte "Precisa")
    // -------------------------------------------------------------------------
    double end_epoch = epoch + 200.0; // Propaga per 200 giorni
    int steps = 100;                  // Numero di punti per il fitting (1 ogni 2 giorni)
    std::vector<Eigen::Vector3d> positions;
    std::vector<double> times;

    std::cout << "Propagazione numerica per " << (end_epoch - epoch) << " giorni..." << std::endl;
    auto t_start = std::chrono::high_resolution_clock::now();

    double dt = (end_epoch - epoch) / (steps - 1);
    for (int k = 0; k < steps; ++k) {
        double t = epoch + k * dt;
        auto state = propagator.propagateToEpoch(t);
        positions.push_back(state.position);
        times.push_back(t);
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double duration_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    
    std::cout << "Propagazione completata in " << duration_ms << " ms." << std::endl;

    // -------------------------------------------------------------------------
    // 3. Compressione Chebyshev (La parte "Veloce" per query successive)
    // -------------------------------------------------------------------------
    // Se devi interrogare questa orbita migliaia di volte, convertila in polinomi di Chebyshev.
    // Questo permette query in < 1 microsecondo con precisione quasi identica alla propagazione numerica.
    
    std::cout << "\nFitting Polinomi di Chebyshev..." << std::endl;
    ChebyshevApproximation chebyshev(16); // 16 coefficienti sono solitamente sufficienti per 200gg
    bool fit_ok = chebyshev.fit(positions, epoch, end_epoch);
    
    if (fit_ok) {
        std::cout << "Fitting riuscito! Errore RMS: " << chebyshev.getApproximationError().norm() << " AU" << std::endl;
        
        // Esempio di query veloce
        double query_time = epoch + 123.456; // Un tempo a caso nel range
        auto pos_cheb = chebyshev.evaluatePosition(query_time);
        
        // Verifica accuratezza contro propagatore completo
        auto state_ref = propagator.propagateToEpoch(query_time);
        double err = (pos_cheb - state_ref.position).norm() * 149597870.7; // AU -> km
        
        std::cout << "Test a T+" << 123.456 << " giorni:" << std::endl;
        std::cout << "  Errore approssimazione: " << err << " km" << std::endl;
        std::cout << "  (L'errore è trascurabile per occultazioni, ma la query è 100x più veloce)" << std::endl;
    }

    return 0;
}
