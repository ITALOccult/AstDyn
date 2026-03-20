/**
 * @file aas_propagation_demo.cpp
 * @brief Esempio di propagazione orbitale di (1) Ceres tramite integratore AAS
 * 
 * Questo codice dimostra come configurare l'AstDynEngine per utilizzare l'integratore
 * Adaptive Symplectic (AAS) e propagare uno stato iniziale cartesiano.
 */

#include "astdyn/AstDynEngine.hpp"
#include "astdyn/core/Constants.hpp"
#include <iostream>
#include <iomanip>

using namespace astdyn;

int main() {
    // 1. Inizializzazione motore AstDyn
    AstDynEngine engine;
    
    // 2. Configurazione per AAS Integrator (Adaptive Symplectic)
    AstDynConfig config;
    config.integrator_type = IntegratorType::AAS;
    config.aas_precision = 1e-8;   // Precisione elevata
    
    // Attivazione delle perturbazioni scientifiche (Modello AST17)
    // NOTA: include_planets richiede un file di effemeridi (es. DE441.bsp) configurato tramite PlanetaryEphemeris::setProvider.
    // In questo esempio usiamo il modello AST17 analitico per gli asteroidi che non richiede file esterni.
    config.propagator_settings.include_planets = false;    // Richiede provider esterno (DE441)
    config.propagator_settings.include_asteroids = true;   // AST17 (16 Asteroidi - Analitico)
    config.propagator_settings.include_relativity = true;  // Correzioni GR (Analitico)
    
    // Selezioniamo quali pianeti includere esplicitamente
    config.propagator_settings.perturb_venus = true;
    config.propagator_settings.perturb_earth = true;
    config.propagator_settings.perturb_mars = true;
    config.propagator_settings.perturb_jupiter = true;
    config.propagator_settings.perturb_saturn = true;
    
    engine.set_config(config);
    
    // 3. Dati orbitali di (1) Ceres (Epoch: 2024-01-01 / MJD 60310.0)
    auto epoch = time::EpochTDB::from_mjd(60310.0);
    
    // Versione Tipizzata 3.0 (Best Practice):
    // Definiamo i vettori posizione e velocità usando i factory per unità astronomiche (AU, AU/day)
    // che gestiscono automaticamente le conversioni in SI (metri, m/s).
    auto pos = math::Vector3<core::GCRF, physics::Distance>::from_si(
        physics::Distance::from_au(-1.777322920405232).to_m(),
        physics::Distance::from_au( 2.108027731776856).to_m(),
        physics::Distance::from_au( 1.611110014769046).to_m()
    );

    auto vel = math::Vector3<core::GCRF, physics::Velocity>::from_si(
        physics::Velocity::from_au_d(-8.841551069485747e-3).to_ms(),
        physics::Velocity::from_au_d(-6.012176882200350e-3).to_ms(),
        physics::Velocity::from_au_d(-3.882415174092120e-3).to_ms()
    );

    // Stato Cartesiano tipizzato completo (Epoca, Posizione, Velocità, GM)
    auto ceres_state = physics::CartesianStateTyped<core::GCRF>(
        epoch, 
        pos, 
        vel, 
        physics::GravitationalParameter::sun()
    );
    
    // Carichiamo l'orbita nell'engine (con conversione automatica interna)
    engine.set_initial_orbit(propagation::cartesian_to_keplerian(ceres_state));
    
    std::cout << "--- Stato Iniziale Ceres (MJD 60310.0) ---" << std::endl;
    std::cout << "Posizione [AU]: " << std::fixed << std::setprecision(8)
              << ceres_state.position.x_si() / (constants::AU * 1000.0) << " " 
              << ceres_state.position.y_si() / (constants::AU * 1000.0) << " " 
              << ceres_state.position.z_si() / (constants::AU * 1000.0) << std::endl;

    // 4. Propagazione di 1000 giorni (MJD 61310.0)
    auto target_epoch = time::EpochTDB::from_mjd(61310.0);
    
    std::cout << "\nPropagazione in corso con AAS..." << std::endl;
    auto final_orbit = engine.propagate_to(target_epoch);
    
    // 5. Risultato
    // Convertiamo il risultato in Cartesiano per visualizzare le coordinate finali
    auto final_cart = propagation::keplerian_to_cartesian(final_orbit);
    
    std::cout << "--- Stato Finale (MJD 61310.0) ---" << std::endl;
    std::cout << "Posizione [AU]: " << std::setprecision(10)
              << final_cart.position.x_si() / (constants::AU * 1000.0) << " " 
              << final_cart.position.y_si() / (constants::AU * 1000.0) << " " 
              << final_cart.position.z_si() / (constants::AU * 1000.0) << std::endl;

    return 0;
}
