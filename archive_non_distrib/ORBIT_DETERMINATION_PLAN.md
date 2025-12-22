# Piano Implementazione: Orbit Determination (OrbFit-compatible)

## Obiettivo
Implementare un sistema completo di determinazione orbitale compatibile con OrbFit, che:
1. Legge osservazioni da file `.rwo` (formato OrbFit)
2. Propaga orbita con calcolo State Transition Matrix (STM)
3. Calcola residui O-C (Observed - Computed)
4. Esegue fit least-squares per correggere elementi orbitali
5. Produce output compatibile con OrbFit

## Componenti Necessari

### 1. Parser RWO (Già parzialmente presente)
**File:** `astdyn/include/astdyn/parsers/OrbFitRWOParser.hpp`

**Funzionalità:**
- Legge file `.rwo` (osservazioni ottiche + radar)
- Estrae: RA, Dec, tempo, codice osservatorio, pesi
- Gestisce errori e selezione osservazioni

### 2. State Transition Matrix (STM)
**Nuovo:** `astdyn/include/astdyn/propagation/STMPropagator.hpp`

**Teoria:**
La STM $\Phi(t, t_0)$ descrive come piccole variazioni negli elementi orbitali si propagano:

$$\delta \mathbf{x}(t) = \Phi(t, t_0) \delta \mathbf{x}(t_0)$$

dove $\mathbf{x} = [\mathbf{r}, \mathbf{v}]^T$ è lo stato.

**Equazione:**
$$\dot{\Phi} = \frac{\partial \mathbf{f}}{\partial \mathbf{x}} \Phi, \quad \Phi(t_0, t_0) = \mathbf{I}_{6\times6}$$

**Implementazione:**
- Propaga simultaneamente stato (6D) e STM (36D) → sistema 42D
- Usa **Radau15** (stabile per sistemi grandi)
- Calcola Jacobiano $\partial \mathbf{f}/\partial \mathbf{x}$ numericamente o analiticamente

### 3. Residual Calculator
**Nuovo:** `astdyn/include/astdyn/orbit_determination/ResidualCalculator.hpp`

**Funzionalità:**
- Propaga orbita all'epoca di ogni osservazione
- Calcola posizione topocent rica (osservatore)
- Converte a RA/Dec
- Calcola residui: $O - C$

**Correzioni necessarie:**
- Aberrazione stellare
- Light-time
- Parallasse topocent rica
- Correzione atmosferica (opzionale)

### 4. Least Squares Fitter
**Nuovo:** `astdyn/include/astdyn/orbit_determination/LeastSquaresFitter.hpp`

**Algoritmo:**
```
1. Elementi iniziali: x₀
2. LOOP fino a convergenza:
   a. Propaga con STM: x(tᵢ), Φ(tᵢ, t₀)
   b. Calcola residui: Δρᵢ = O - C
   c. Calcola matrice design: A = ∂ρ/∂x₀ = (∂ρ/∂x) × Φ
   d. Risolvi: (AᵀWA) δx₀ = AᵀW Δρ
   e. Aggiorna: x₀ ← x₀ + δx₀
   f. Check convergenza: |δx₀| < ε
3. RETURN: x₀ ottimizzato, covarianza, RMS residui
```

**Peso osservazioni:**
$$W = \text{diag}(1/\sigma_i^2)$$

### 5. Integratore Ottimale: **NON Radau15 standard**

**Problema Radau15 attuale:**
- Troppo lento per Kepler semplice
- Newton solver non ottimizzato per OD

**Soluzione: Usa RKF78 con STM!**

**Perché RKF78 è meglio per OD:**
- ✅ Veloce (esplicito)
- ✅ Preciso (ordine 7)
- ✅ Adaptive step (gestisce close approaches)
- ✅ Già ottimizzato in AstDyn
- ✅ Usato da JPL per effemeridi

**Radau15 serve solo se:**
- Problema è stiff (rare in OD)
- Serve precisione estrema (> 1e-13)
- Close approaches molto stretti

## Architettura Proposta

```cpp
// 1. Carica osservazioni
OrbFitRWOParser parser;
auto observations = parser.parse("asteroid.rwo");

// 2. Elementi iniziali (da .eq1 o preliminary orbit)
KeplerianElements initial_elements = ...;

// 3. Setup orbit determination
OrbitDetermination od;
od.setObservations(observations);
od.setInitialElements(initial_elements);

// 4. Configura propagatore (RKF78 con STM)
RKF78Integrator integrator(0.1, 1e-12);
STMPropagator propagator(integrator);
od.setPropagator(propagator);

// 5. Esegui fit
auto result = od.fit();

// 6. Risultati
std::cout << "RMS residui: " << result.rms << " arcsec\n";
std::cout << "Elementi corretti:\n" << result.elements << "\n";
std::cout << "Covarianza:\n" << result.covariance << "\n";

// 7. Salva in formato OrbFit
result.saveToEQ1("asteroid_fitted.eq1");
```

## Timeline Implementazione

### Fase 1: STM Propagator (2-3 ore)
- [ ] Creare `STMPropagator` class
- [ ] Implementare calcolo Jacobiano
- [ ] Propagare stato + STM simultaneamente
- [ ] Test: verificare $\Phi \Phi^{-1} = I$

### Fase 2: Residual Calculator (2 ore)
- [ ] Implementare conversione Cartesian → RA/Dec
- [ ] Aggiungere correzioni (aberrazione, light-time)
- [ ] Gestire osservatori (MPC codes)
- [ ] Test: confronto con OrbFit

### Fase 3: Least Squares (3 ore)
- [ ] Implementare algoritmo LS
- [ ] Calcolare matrice design
- [ ] Risolvere sistema normale
- [ ] Calcolare covarianza
- [ ] Test: fit sintetico

### Fase 4: Integrazione (1 ora)
- [ ] Creare classe `OrbitDetermination`
- [ ] Integrare tutti i componenti
- [ ] Test end-to-end con dati reali

### Fase 5: Validazione (2 ore)
- [ ] Confronto con OrbFit su stesso dataset
- [ ] Verificare residui identici
- [ ] Verificare elementi corretti identici
- [ ] Documentazione

**Totale: ~10-12 ore**

## File da Creare

1. `astdyn/include/astdyn/propagation/STMPropagator.hpp`
2. `astdyn/src/propagation/STMPropagator.cpp`
3. `astdyn/include/astdyn/orbit_determination/ResidualCalculator.hpp`
4. `astdyn/src/orbit_determination/ResidualCalculator.cpp`
5. `astdyn/include/astdyn/orbit_determination/LeastSquaresFitter.hpp`
6. `astdyn/src/orbit_determination/LeastSquaresFitter.cpp`
7. `astdyn/include/astdyn/orbit_determination/OrbitDetermination.hpp`
8. `astdyn/src/orbit_determination/OrbitDetermination.cpp`
9. `examples/orbit_determination_example.cpp`
10. `test_orbit_determination.cpp`

## Esempio d'Uso

```cpp
// Esempio completo: fit orbita asteroide 17030 Sierks

#include "astdyn/orbit_determination/OrbitDetermination.hpp"
#include "astdyn/parsers/OrbFitRWOParser.hpp"
#include "astdyn/parsers/OrbFitEQ1Parser.hpp"

int main() {
    // 1. Carica osservazioni
    OrbFitRWOParser rwo_parser;
    auto obs = rwo_parser.parse("17030.rwo");
    std::cout << "Loaded " << obs.size() << " observations\n";
    
    // 2. Carica elementi iniziali
    OrbFitEQ1Parser eq1_parser;
    auto initial = eq1_parser.parse("17030_initial.eq1");
    
    // 3. Setup OD
    OrbitDetermination od;
    od.setObservations(obs);
    od.setInitialElements(initial.elements);
    od.setEpoch(initial.epoch_mjd);
    
    // 4. Configura opzioni
    od.setMaxIterations(10);
    od.setConvergenceTolerance(1e-6);
    od.setOutlierRejection(3.0);  // 3-sigma
    
    // 5. Esegui fit
    auto result = od.fit();
    
    // 6. Report
    std::cout << "\n=== ORBIT DETERMINATION RESULTS ===\n";
    std::cout << "Iterations: " << result.num_iterations << "\n";
    std::cout << "RMS (RA):  " << result.rms_ra << " arcsec\n";
    std::cout << "RMS (Dec): " << result.rms_dec << " arcsec\n";
    std::cout << "RMS (tot): " << result.rms_total << " arcsec\n";
    std::cout << "Observations used: " << result.num_used << "/" << obs.size() << "\n";
    
    // 7. Elementi corretti
    std::cout << "\nCorrected elements:\n";
    std::cout << "  a = " << result.elements.a << " AU\n";
    std::cout << "  e = " << result.elements.e << "\n";
    std::cout << "  i = " << result.elements.i * 180/M_PI << " deg\n";
    
    // 8. Salva risultati
    result.saveToEQ1("17030_fitted.eq1");
    result.saveResidualsToRWO("17030_residuals.rwo");
    
    return 0;
}
```

## Prossimi Passi

Vuoi che proceda con l'implementazione di Orbit Determination completo?

Posso iniziare con:
1. **STMPropagator** (fondamentale)
2. **ResidualCalculator** (calcolo O-C)
3. **LeastSquaresFitter** (fit elementi)

Oppure preferisci prima finire il benchmark e poi fare OD?
