/**
 * @file test_ioccultcalc_integration.cpp
 * @brief Test integrazione ITALOccultLibrary → IOccultCalc
 * @author IOccultCalc Team
 * @date 1 Dicembre 2025
 * 
 * Test completo del modulo di integrazione:
 * 1. Caricamento asteroide da .eq1
 * 2. Propagazione con ITALOccultLibrary
 * 3. Conversione in formato IOccultCalc
 * 4. Verifica compatibilità dati
 * 
 * Compilazione:
 *   g++ -std=c++17 -O2 -o test_ioccultcalc_integration \
 *       test_ioccultcalc_integration.cpp \
 *       ../integration/italoccult_integration.cpp \
 *       -I/usr/local/include -I/usr/local/include/eigen3 \
 *       -I../integration \
 *       -L/usr/local/lib -litaloccultlib -lastdyn \
 *       -L/opt/homebrew/lib -lboost_filesystem -lboost_program_options -lboost_date_time
 */

#include "italoccult_integration.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

const double AU_TO_KM = 1.495978707e8;

void print_asteroid_state(const AsteroidState& state) {
    std::cout << "Asteroide: " << state.name << std::endl;
    std::cout << "Epoca: MJD " << std::fixed << std::setprecision(2) 
              << state.epoch_mjd_tdb << " TDB" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Posizione ICRF (AU):" << std::endl;
    std::cout << "  X = " << std::setprecision(9) << state.x_au << std::endl;
    std::cout << "  Y = " << state.y_au << std::endl;
    std::cout << "  Z = " << state.z_au << std::endl;
    
    double r_au = std::sqrt(state.x_au * state.x_au + 
                           state.y_au * state.y_au + 
                           state.z_au * state.z_au);
    std::cout << "  R = " << std::setprecision(6) << r_au << " AU ("
              << (r_au * AU_TO_KM / 1e6) << " milioni km)" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Velocità ICRF (AU/day):" << std::endl;
    std::cout << "  VX = " << std::setprecision(9) << state.vx_au_day << std::endl;
    std::cout << "  VY = " << state.vy_au_day << std::endl;
    std::cout << "  VZ = " << state.vz_au_day << std::endl;
    
    double v = std::sqrt(state.vx_au_day * state.vx_au_day + 
                        state.vy_au_day * state.vy_au_day + 
                        state.vz_au_day * state.vz_au_day);
    std::cout << "  V = " << v << " AU/day" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Parametri orbitali (approssimati):" << std::endl;
    std::cout << "  a = " << state.a_au << " AU" << std::endl;
    std::cout << "  e = " << state.e << std::endl;
    std::cout << "  i = " << state.i_deg << " deg" << std::endl;
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "  Test Integrazione IOccultCalc" << std::endl;
    std::cout << "  ITALOccultLibrary → IOccultCalc" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    try {
        // ===== TEST 1: Integrazione base =====
        std::cout << "[TEST 1/5] Creazione integrator..." << std::endl;
        
        ITALOccultIntegration integrator(true);  // High accuracy
        std::cout << "  ✓ Integrator creato (high accuracy mode)" << std::endl;
        std::cout << std::endl;

        // ===== TEST 2: Caricamento asteroide =====
        std::cout << "[TEST 2/5] Caricamento asteroide da .eq1..." << std::endl;
        
        std::string eq1_file = "astdyn/data/17030.eq1";
        bool loaded = integrator.loadAsteroidFromEQ1(eq1_file);
        
        if (!loaded) {
            std::cerr << "  ✗ ERRORE: impossibile caricare " << eq1_file << std::endl;
            return 1;
        }
        
        std::cout << "  ✓ Asteroide caricato: " << integrator.getAsteroidName() << std::endl;
        std::cout << "  ✓ Epoca: MJD " << integrator.getCurrentEpoch() << std::endl;
        std::cout << std::endl;

        // ===== TEST 3: Propagazione singola =====
        std::cout << "[TEST 3/5] Propagazione singola epoca..." << std::endl;
        
        double target_mjd = integrator.getCurrentEpoch() + 7.0;
        auto state = integrator.propagateToEpoch(target_mjd);
        
        std::cout << "  ✓ Propagazione completata" << std::endl;
        std::cout << std::endl;
        
        print_asteroid_state(state);
        std::cout << std::endl;
        
        std::cout << "Info propagazione:" << std::endl;
        std::cout << integrator.getLastPropagationInfo() << std::endl;
        std::cout << std::endl;

        // ===== TEST 4: Propagazione multipla =====
        std::cout << "[TEST 4/5] Propagazione multipla epoche..." << std::endl;
        
        double epoch_0 = integrator.getCurrentEpoch();
        std::vector<double> epochs = {
            epoch_0,
            epoch_0 + 1.0,
            epoch_0 + 3.0,
            epoch_0 + 7.0,
            epoch_0 + 14.0
        };
        
        auto states = integrator.propagateToEpochs(epochs);
        
        std::cout << "  ✓ Propagate " << states.size() << " epoche" << std::endl;
        std::cout << std::endl;
        
        std::cout << "Traiettoria:" << std::endl;
        for (size_t i = 0; i < states.size(); ++i) {
            double r = std::sqrt(states[i].x_au * states[i].x_au + 
                               states[i].y_au * states[i].y_au + 
                               states[i].z_au * states[i].z_au);
            std::cout << "  MJD " << std::fixed << std::setprecision(1) 
                      << states[i].epoch_mjd_tdb 
                      << ": R = " << std::setprecision(6) << r << " AU" << std::endl;
        }
        std::cout << std::endl;

        // ===== TEST 5: Quick propagate helper =====
        std::cout << "[TEST 5/5] Test funzione helper quickPropagateFromEQ1..." << std::endl;
        
        auto quick_state = quickPropagateFromEQ1(eq1_file, target_mjd, true);
        
        std::cout << "  ✓ Quick propagate completato" << std::endl;
        
        // Verifica consistenza con metodo normale
        double diff = std::abs(quick_state.x_au - state.x_au) +
                     std::abs(quick_state.y_au - state.y_au) +
                     std::abs(quick_state.z_au - state.z_au);
        
        bool consistent = (diff < 1e-12);
        std::cout << "  Consistenza con metodo normale: " 
                  << (consistent ? "✓ PASS" : "✗ FAIL") << std::endl;
        std::cout << "  (differenza: " << std::scientific << diff << " AU)" << std::endl;
        std::cout << std::endl;

        // ===== VALIDAZIONE FINALE =====
        std::cout << "========================================" << std::endl;
        std::cout << "  VALIDAZIONE INTEGRAZIONE" << std::endl;
        std::cout << "========================================" << std::endl;
        
        // Verifica formato dati IOccultCalc
        bool has_name = !state.name.empty();
        bool has_epoch = (state.epoch_mjd_tdb > 0.0);
        bool has_position = (state.x_au != 0.0 || state.y_au != 0.0 || state.z_au != 0.0);
        bool has_velocity = (state.vx_au_day != 0.0 || state.vy_au_day != 0.0 || state.vz_au_day != 0.0);
        bool has_orbital = (state.a_au > 0.0);
        
        std::cout << "Formato AsteroidState IOccultCalc:" << std::endl;
        std::cout << "  Nome:         " << (has_name ? "✓" : "✗") << std::endl;
        std::cout << "  Epoca:        " << (has_epoch ? "✓" : "✗") << std::endl;
        std::cout << "  Posizione:    " << (has_position ? "✓" : "✗") << std::endl;
        std::cout << "  Velocità:     " << (has_velocity ? "✓" : "✗") << std::endl;
        std::cout << "  Orbitali:     " << (has_orbital ? "✓" : "✗") << std::endl;
        std::cout << std::endl;

        bool all_ok = has_name && has_epoch && has_position && 
                     has_velocity && has_orbital && consistent;
        
        std::cout << "========================================" << std::endl;
        std::cout << "  RISULTATO FINALE" << std::endl;
        std::cout << "========================================" << std::endl;
        
        if (all_ok) {
            std::cout << "✓ INTEGRAZIONE COMPLETATA CON SUCCESSO" << std::endl;
            std::cout << std::endl;
            std::cout << "ITALOccultLibrary è pronta per IOccultCalc:" << std::endl;
            std::cout << "  ✓ Caricamento .eq1 funzionante" << std::endl;
            std::cout << "  ✓ Propagazione AstDyn operativa" << std::endl;
            std::cout << "  ✓ Conversione formato IOccultCalc OK" << std::endl;
            std::cout << "  ✓ API high-level pronta" << std::endl;
            std::cout << "  ✓ Helper functions disponibili" << std::endl;
        } else {
            std::cout << "✗ INTEGRAZIONE FALLITA" << std::endl;
            std::cout << "  Verificare log sopra per dettagli" << std::endl;
        }
        std::cout << std::endl;
        
        return all_ok ? 0 : 1;

    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
}
