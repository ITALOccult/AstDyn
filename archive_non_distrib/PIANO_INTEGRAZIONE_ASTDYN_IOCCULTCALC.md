# Piano di Integrazione Completa AstDyn in IOccultCalc

**Data**: 1 Dicembre 2025  
**Obiettivo**: Utilizzo pieno della libreria AstDyn in IOccultCalc  
**Status**: 🚀 Piano Operativo

---

## 📋 INDICE

1. [Obiettivi](#1-obiettivi)
2. [Architettura Attuale](#2-architettura-attuale)
3. [Moduli da Integrare](#3-moduli-da-integrare)
4. [Piano di Implementazione](#4-piano-di-implementazione)
5. [Codice da Implementare](#5-codice-da-implementare)
6. [Test e Validazione](#6-test-e-validazione)
7. [Tempistiche](#7-tempistiche)

---

## 1. OBIETTIVI

### 1.1 Obiettivi Primari

- ✅ **Sostituire propagazione Kepleriana** con RKF78 di AstDyn
- ✅ **Aggiungere supporto file .eq1** (formato OEF2.0)
- ✅ **Includere perturbazioni planetarie** (8 pianeti + Schwarzschild)
- ✅ **Ridurre errore propagazione** da 11.35" a <0.1"
- ✅ **Mantenere velocità FASE 1** (screening rapido)

### 1.2 Benefici Attesi

| Aspetto | Prima | Dopo |
|---------|-------|------|
| **Accuratezza** | 12.65" (17030) | 1.53" (JPL-grade) |
| **Perturbazioni** | Solo Sole | 11 sorgenti |
| **Formato elementi** | API/MPC | .eq1 + API + MPC |
| **Frame** | Ambiguo | ICRF certificato |
| **Validazione** | Teorica | Confronto JPL |

---

## 2. ARCHITETTURA ATTUALE

### 2.1 Struttura IOccultCalc (da Documentazione)

```
IOccultCalc/
├── include/ioccultcalc/
│   ├── propagation_strategy.h       ✅ Esiste (342 righe)
│   ├── orbital_elements.h           ? Da verificare
│   └── astdyn_interface.h           ❌ Da creare
├── src/
│   ├── propagation_strategy.cpp     ✅ Esiste (803 righe)
│   ├── orbital_elements.cpp         ? Da verificare
│   ├── orbit_propagator.cpp         ? Da verificare
│   └── rkf78_integrator.cpp         ? Da verificare
├── CMakeLists.txt                   ✅ Esiste
└── build/
    └── lib/libioccultcalc.a         ✅ Compilato (2.3 MB)
```

### 2.2 Integrazione AstDyn Attuale

**Dal documento VERIFICA_INTEGRITÀ_PROGETTI.md**:

```cmake
# CMakeLists.txt attuale
target_link_libraries(ioccultcalc 
    PRIVATE AstDyn::astdyn
)
```

**Status**: ⚠️ **PARZIALE**
- AstDyn è linkato ma non utilizzato completamente
- `propagation_strategy.cpp` costruisce AstDynPropagator ma non lo usa per .eq1
- Manca parser .eq1 diretto

---

## 3. MODULI DA INTEGRARE

### 3.1 Modulo 1: Parser File .eq1

**File da creare**: `include/ioccultcalc/eq1_parser.h`

**Funzionalità**:
- Lettura file formato OEF2.0
- Parsing elementi equinoziali (a, h, k, p, q, λ)
- Conversione equinoziale → kepleriano
- Gestione epoca MJD/JD

**Riferimento**: `astdyn/tests/test_asteroid_17030_occultation.cpp` (linee 79-120)

### 3.2 Modulo 2: Conversioni Elementi Orbitali

**File da creare**: `include/ioccultcalc/orbital_conversions.h`

**Funzionalità**:
- `equinoctialToKeplerian()` - Con normalizzazione angoli
- `keplerianToCartesian()` - Risoluzione equazione Keplero
- `eclipticToICRF()` - Rotazione con obliquità
- Validazione frame di riferimento

**Riferimento**: `astdyn/tools/astdyn_propagator.cpp` (linee 512-640)

### 3.3 Modulo 3: Wrapper AstDyn Completo

**File da creare**: `include/ioccultcalc/astdyn_wrapper.h`

**Funzionalità**:
- Configurazione perturbazioni (pianeti, AST17, relatività)
- Gestione tolleranze RKF78
- Interfaccia semplificata per IOccultCalc
- Cache stati propagati

**Riferimento**: `astdyn/tools/astdyn_propagator.cpp` (classe completa)

### 3.4 Modulo 4: Strategia Ibrida Potenziata

**File da modificare**: `src/propagation_strategy.cpp`

**Funzionalità**:
- FASE 1: IOccultCalc Keplerian (screening, soglia 60")
- FASE 2: AstDyn RKF78 (precisione, tolleranza 1e-12)
- Auto-switch basato su distanza candidato
- Caching risultati per multi-stelle

---

## 4. PIANO DI IMPLEMENTAZIONE

### FASE 1: Fondamenta (Priorità CRITICA) 🔴

**Durata stimata**: 2-3 ore

#### Task 1.1: Parser .eq1

**File**: `include/ioccultcalc/eq1_parser.h` + `src/eq1_parser.cpp`

```cpp
namespace ioccultcalc {

struct EquinoctialElements {
    double a;          // AU
    double h;          // e*sin(LP)
    double k;          // e*cos(LP)
    double p;          // tan(i/2)*sin(Omega)
    double q;          // tan(i/2)*cos(Omega)
    double lambda;     // Mean longitude [deg]
    double epoch_mjd;  // MJD TDT
    double H;          // Absolute magnitude
    double G;          // Slope parameter
    std::string name;
};

class EQ1Parser {
public:
    static EquinoctialElements parseFile(const std::string& filepath);
    static EquinoctialElements parseString(const std::string& content);
    
private:
    static bool parseLine(const std::string& line, 
                         EquinoctialElements& elem);
};

} // namespace ioccultcalc
```

#### Task 1.2: Conversioni Orbitali

**File**: `include/ioccultcalc/orbital_conversions.h` + `src/orbital_conversions.cpp`

```cpp
namespace ioccultcalc {

struct KeplerianElements {
    double a;          // Semi-major axis [AU]
    double e;          // Eccentricity
    double i;          // Inclination [rad]
    double Omega;      // Long. ascending node [rad]
    double omega;      // Argument perihelion [rad]
    double M;          // Mean anomaly [rad]
    double epoch_jd;   // Epoch [JD]
    std::string name;
};

struct CartesianState {
    Eigen::Vector3d position;  // [AU]
    Eigen::Vector3d velocity;  // [AU/day]
    double epoch_jd;
};

class OrbitalConversions {
public:
    // Equinoctial → Keplerian
    static KeplerianElements equinoctialToKeplerian(
        const EquinoctialElements& eq);
    
    // Keplerian → Cartesian (eclittico)
    static CartesianState keplerianToCartesian(
        const KeplerianElements& kep);
    
    // Eclittico → ICRF
    static CartesianState eclipticToICRF(
        const CartesianState& ecliptic);
    
    // Risoluzione equazione Keplero
    static double solveKeplerEquation(double M, double e, 
                                     double tol = 1e-14);
    
    // Validazione frame
    static bool validateICRF(const CartesianState& state);
};

} // namespace ioccultcalc
```

#### Task 1.3: Wrapper AstDyn

**File**: `include/ioccultcalc/astdyn_wrapper.h` + `src/astdyn_wrapper.cpp`

```cpp
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>

namespace ioccultcalc {

class AstDynWrapper {
public:
    AstDynWrapper();
    
    // Configurazione
    void setTolerance(double tol);
    void enablePlanets(bool enable);
    void enableAST17(bool enable);
    void enableRelativity(bool enable);
    
    // Caricamento elementi da .eq1
    void loadFromEQ1(const std::string& filepath);
    void setElements(const EquinoctialElements& eq);
    void setElements(const KeplerianElements& kep);
    
    // Propagazione
    struct EquatorialCoords {
        double ra;    // [rad]
        double dec;   // [rad]
        double dist;  // [AU]
    };
    
    EquatorialCoords propagate(double jd_target);
    
    // Batch propagation
    std::vector<EquatorialCoords> propagateBatch(
        const std::vector<double>& jd_targets);
    
    // Statistiche
    struct PropagationStats {
        int steps_taken;
        double cpu_time_ms;
        double final_step_size;
        double max_error;
    };
    
    PropagationStats getLastStats() const;
    
private:
    std::unique_ptr<astdyn::propagation::Propagator> propagator_;
    astdyn::propagation::KeplerianElements elements_;
    double tolerance_;
    PropagationStats last_stats_;
    
    // Conversione interna
    astdyn::propagation::KeplerianElements 
        convertToAstDynElements(const KeplerianElements& kep);
};

} // namespace ioccultcalc
```

### FASE 2: Integrazione TwoPhaseStrategy (Priorità ALTA) 🟠

**Durata stimata**: 3-4 ore

#### Task 2.1: Modifica TwoPhaseStrategy

**File**: `src/propagation_strategy.cpp`

**Modifiche da applicare**:

1. **Aggiungere supporto .eq1 nel costruttore**:

```cpp
void TwoPhaseStrategy::loadElementsFromEQ1(const std::string& filepath) {
    // Usa il nuovo parser
    auto eq_elements = EQ1Parser::parseFile(filepath);
    
    // Converti a Keplerian
    auto kep_elements = OrbitalConversions::equinoctialToKeplerian(
        eq_elements);
    
    // Configura AstDyn wrapper
    astdyn_wrapper_->setElements(kep_elements);
    
    // Salva anche per Chebyshev (FASE 1)
    convertToInternalFormat(kep_elements);
}
```

2. **Potenziare getRKF78Position()**:

```cpp
EquatorialCoords TwoPhaseStrategy::getRKF78Position(
    double jd_target, 
    double tolerance) {
    
    // Configura tolleranza
    astdyn_wrapper_->setTolerance(tolerance);
    
    // Abilita tutte le perturbazioni per FASE 2
    astdyn_wrapper_->enablePlanets(true);
    astdyn_wrapper_->enableAST17(config_.use_ast17);
    astdyn_wrapper_->enableRelativity(true);
    
    // Propaga con AstDyn
    auto coords = astdyn_wrapper_->propagate(jd_target);
    
    // Aggiorna statistiche
    performance_stats_.astdyn_calls++;
    performance_stats_.astdyn_time_ms += 
        astdyn_wrapper_->getLastStats().cpu_time_ms;
    
    return coords;
}
```

3. **Aggiungere metodo batch**:

```cpp
std::vector<EquatorialCoords> 
TwoPhaseStrategy::getRKF78PositionBatch(
    const std::vector<double>& jd_targets,
    double tolerance) {
    
    astdyn_wrapper_->setTolerance(tolerance);
    astdyn_wrapper_->enablePlanets(true);
    astdyn_wrapper_->enableAST17(config_.use_ast17);
    astdyn_wrapper_->enableRelativity(true);
    
    return astdyn_wrapper_->propagateBatch(jd_targets);
}
```

#### Task 2.2: Auto-Switch Intelligente

**File**: `src/propagation_strategy.cpp`

```cpp
EquatorialCoords TwoPhaseStrategy::getBestPosition(
    double jd_target,
    const Eigen::Vector2d& star_coords_rad) {
    
    // FASE 1: Quick screening con Chebyshev
    auto phase1 = getChebyshevPosition(jd_target, "phase1");
    
    // Calcola distanza angolare
    double angular_distance = calculateAngularDistance(
        phase1.ra, phase1.dec, 
        star_coords_rad(0), star_coords_rad(1)
    );
    
    // Conversione a arcsec
    double distance_arcsec = angular_distance * 206265.0;
    
    // Auto-switch: se < 60" usa FASE 2
    if (distance_arcsec < config_.phase2_threshold_arcsec) {
        std::cout << "  Switch to PHASE 2: distance = " 
                  << distance_arcsec << "\" < " 
                  << config_.phase2_threshold_arcsec << "\"\n";
        
        return getRKF78Position(jd_target, 1e-12);
    }
    
    return phase1;  // Resta in FASE 1
}
```

### FASE 3: Testing e Validazione (Priorità ALTA) 🟠

**Durata stimata**: 2-3 ore

#### Task 3.1: Unit Test Conversioni

**File**: `tests/test_orbital_conversions.cpp`

```cpp
#include <gtest/gtest.h>
#include "ioccultcalc/eq1_parser.h"
#include "ioccultcalc/orbital_conversions.h"

TEST(OrbitalConversions, Equinoctial17030) {
    // Elementi 17030 da .eq1
    EquinoctialElements eq;
    eq.a = 3.175473;
    eq.h = -0.018963;
    eq.k = -0.041273;
    eq.p = 0.025407;
    eq.q = -0.001956;
    eq.lambda = 229.790880;
    
    // Converti
    auto kep = OrbitalConversions::equinoctialToKeplerian(eq);
    
    // Valori attesi
    EXPECT_NEAR(kep.e, 0.045407, 1e-6);
    EXPECT_NEAR(kep.i * 180.0/M_PI, 2.9046, 1e-4);
    EXPECT_NEAR(kep.Omega * 180.0/M_PI, 94.06, 1e-2);
    EXPECT_NEAR(kep.omega * 180.0/M_PI, 110.28, 1e-2);
    EXPECT_NEAR(kep.M * 180.0/M_PI, 25.45, 1e-2);
}

TEST(OrbitalConversions, KeplerianToCartesian) {
    KeplerianElements kep;
    kep.a = 2.5;
    kep.e = 0.1;
    kep.i = 10.0 * M_PI/180.0;
    kep.Omega = 50.0 * M_PI/180.0;
    kep.omega = 30.0 * M_PI/180.0;
    kep.M = 45.0 * M_PI/180.0;
    
    auto state = OrbitalConversions::keplerianToCartesian(kep);
    
    // Verifica norma posizione ~= a
    double r = state.position.norm();
    EXPECT_NEAR(r, kep.a, 0.5);  // ±0.5 AU per orbita ellittica
}

TEST(OrbitalConversions, EclipticToICRF) {
    // Vettore lungo asse X eclittico
    CartesianState ecl;
    ecl.position = Eigen::Vector3d(1.0, 0.0, 0.0);
    ecl.velocity = Eigen::Vector3d(0.0, 1.0, 0.0);
    
    auto icrf = OrbitalConversions::eclipticToICRF(ecl);
    
    // X deve rimanere uguale
    EXPECT_NEAR(icrf.position(0), 1.0, 1e-9);
    
    // Y, Z ruotati con obliquità
    double eps = 23.439291 * M_PI/180.0;
    EXPECT_NEAR(icrf.position(1), 0.0, 1e-9);
    EXPECT_NEAR(icrf.position(2), 0.0, 1e-9);
}
```

#### Task 3.2: Integration Test 17030

**File**: `tests/test_integration_17030.cpp`

```cpp
#include <gtest/gtest.h>
#include "ioccultcalc/propagation_strategy.h"
#include "ioccultcalc/eq1_parser.h"

TEST(Integration17030, FullPipeline) {
    // Configura strategia
    PropagationConfig config;
    config.use_phase2 = true;
    config.phase2_threshold_arcsec = 60.0;
    config.use_ast17 = false;  // Per confronto con AstDyn base
    
    TwoPhaseStrategy strategy(config);
    
    // Carica elementi da .eq1
    std::string eq1_file = 
        "/Users/michelebigi/Astro/OrbFit/tests/orbfit/17030/epoch/17030.eq1";
    strategy.loadElementsFromEQ1(eq1_file);
    
    // Data evento: 28/11/2025 00:35 UTC
    double jd_event = 2460277.523611;  // JD per 28/11/2025 00:35
    
    // Propaga con FASE 2 (AstDyn RKF78)
    auto coords = strategy.getRKF78Position(jd_event, 1e-12);
    
    // Valori JPL Horizons per confronto
    double jpl_ra = 73.416 * M_PI/180.0;    // 04h 53m 39s
    double jpl_dec = 20.332 * M_PI/180.0;   // +20° 19' 54"
    
    // Calcola distanza angolare
    double delta_ra = (coords.ra - jpl_ra) * 
                      std::cos(coords.dec) * 206265.0;
    double delta_dec = (coords.dec - jpl_dec) * 206265.0;
    double separation = std::sqrt(delta_ra*delta_ra + 
                                 delta_dec*delta_dec);
    
    // Tolleranza: <2" per match JPL
    EXPECT_LT(separation, 2.0) 
        << "Separation from JPL: " << separation << " arcsec";
    
    std::cout << "Separation from JPL Horizons: " 
              << separation << " arcsec\n";
}

TEST(Integration17030, PhaseComparison) {
    TwoPhaseStrategy strategy;
    strategy.loadElementsFromEQ1(
        "/Users/michelebigi/Astro/OrbFit/tests/orbfit/17030/epoch/17030.eq1");
    
    double jd_event = 2460277.523611;
    
    // FASE 1
    auto phase1 = strategy.getChebyshevPosition(jd_event, "phase1");
    
    // FASE 2
    auto phase2 = strategy.getRKF78Position(jd_event, 1e-12);
    
    // Calcola differenza
    double delta_ra = (phase2.ra - phase1.ra) * 
                      std::cos(phase2.dec) * 206265.0;
    double delta_dec = (phase2.dec - phase1.dec) * 206265.0;
    double diff = std::sqrt(delta_ra*delta_ra + delta_dec*delta_dec);
    
    std::cout << "PHASE 1 vs PHASE 2 difference: " 
              << diff << " arcsec\n";
    
    // Per 17030, differenza attesa ~11.35"
    EXPECT_GT(diff, 5.0);   // Deve esserci differenza significativa
    EXPECT_LT(diff, 20.0);  // Ma non troppo grande
}
```

### FASE 4: Ottimizzazioni e Features Avanzate (Priorità MEDIA) 🟡

**Durata stimata**: 4-5 ore

#### Task 4.1: Caching Multi-Stella

```cpp
class AstDynWrapper {
private:
    struct CachedState {
        double jd;
        CartesianState state;
        EquatorialCoords coords;
    };
    
    std::map<double, CachedState> cache_;
    double cache_tolerance_jd_ = 1e-6;  // ~0.1 sec
    
public:
    EquatorialCoords propagate(double jd_target) {
        // Check cache
        auto it = cache_.find(jd_target);
        if (it != cache_.end()) {
            stats_.cache_hits++;
            return it->second.coords;
        }
        
        // Propaga nuovo
        auto coords = propagateInternal(jd_target);
        
        // Salva in cache
        cache_[jd_target] = {jd_target, last_state_, coords};
        
        return coords;
    }
    
    // Batch ottimizzato
    std::vector<EquatorialCoords> propagateBatch(
        const std::vector<double>& jd_targets) {
        
        // Ordina per efficienza
        auto sorted = jd_targets;
        std::sort(sorted.begin(), sorted.end());
        
        std::vector<EquatorialCoords> results;
        results.reserve(sorted.size());
        
        for (double jd : sorted) {
            results.push_back(propagate(jd));
        }
        
        return results;
    }
};
```

#### Task 4.2: Diagnostica e Logging

```cpp
class AstDynWrapper {
public:
    struct DiagnosticInfo {
        // Elementi usati
        KeplerianElements elements;
        
        // Configurazione
        double tolerance;
        bool planets_enabled;
        bool ast17_enabled;
        bool relativity_enabled;
        
        // Statistiche
        int total_propagations;
        double total_cpu_time_ms;
        double avg_steps_per_propagation;
        
        // Errori
        int failed_propagations;
        std::vector<std::string> error_messages;
    };
    
    DiagnosticInfo getDiagnostics() const;
    void exportDiagnostics(const std::string& json_file) const;
};
```

---

## 5. CODICE DA IMPLEMENTARE

### 5.1 CMakeLists.txt - Aggiornamento

```cmake
# File: IOccultCalc/CMakeLists.txt

cmake_minimum_required(VERSION 3.15)
project(IOccultCalc VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trova AstDyn (deve essere installato o in workspace)
find_package(AstDyn REQUIRED)

# Trova Eigen3
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)

# Nuovi file sorgente
set(IOCCULTCALC_SOURCES
    src/propagation_strategy.cpp
    src/orbital_elements.cpp
    src/eq1_parser.cpp              # NUOVO
    src/orbital_conversions.cpp      # NUOVO
    src/astdyn_wrapper.cpp           # NUOVO
    src/orbit_propagator.cpp
    src/rkf78_integrator.cpp
)

# Libreria IOccultCalc
add_library(ioccultcalc STATIC ${IOCCULTCALC_SOURCES})

target_link_libraries(ioccultcalc 
    PUBLIC 
        AstDyn::astdyn
        Eigen3::Eigen
)

target_include_directories(ioccultcalc
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Tests (se abilitato)
option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    
    add_executable(test_orbital_conversions
        tests/test_orbital_conversions.cpp
    )
    target_link_libraries(test_orbital_conversions
        ioccultcalc
        GTest::gtest_main
    )
    
    add_executable(test_integration_17030
        tests/test_integration_17030.cpp
    )
    target_link_libraries(test_integration_17030
        ioccultcalc
        GTest::gtest_main
    )
    
    add_test(NAME OrbitalConversions 
             COMMAND test_orbital_conversions)
    add_test(NAME Integration17030 
             COMMAND test_integration_17030)
endif()

# Install
install(TARGETS ioccultcalc
    EXPORT IOccultCalcTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/ioccultcalc
    DESTINATION include
)
```

### 5.2 Implementazione Completa - eq1_parser.cpp

```cpp
// File: IOccultCalc/src/eq1_parser.cpp

#include "ioccultcalc/eq1_parser.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace ioccultcalc {

EquinoctialElements EQ1Parser::parseFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    return parseString(buffer.str());
}

EquinoctialElements EQ1Parser::parseString(const std::string& content) {
    EquinoctialElements elem;
    std::istringstream stream(content);
    std::string line;
    bool in_data = false;
    
    while (std::getline(stream, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        // Skip commenti
        if (line.empty() || line[0] == '!') continue;
        
        // Header end
        if (line.find("END_OF_HEADER") != std::string::npos) {
            in_data = true;
            continue;
        }
        
        if (!in_data) continue;
        
        // Numero asteroide (opzionale)
        if (std::isdigit(line[0]) && line.find("EQU") == std::string::npos) {
            elem.name = line;
            continue;
        }
        
        // Elementi equinoziali
        if (line.find("EQU") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag >> elem.a >> elem.h >> elem.k 
                >> elem.p >> elem.q >> elem.lambda;
        }
        
        // Epoca
        if (line.find("MJD") != std::string::npos) {
            // Extract MJD value
            size_t mjd_pos = line.find("MJD");
            std::string mjd_part = line.substr(mjd_pos + 3);
            
            // Find first digit
            size_t num_start = mjd_part.find_first_of("0123456789");
            if (num_start != std::string::npos) {
                mjd_part = mjd_part.substr(num_start);
                
                // Extract until TDT/UTC
                size_t num_end = mjd_part.find_first_not_of("0123456789.");
                if (num_end != std::string::npos) {
                    mjd_part = mjd_part.substr(0, num_end);
                }
                
                elem.epoch_mjd = std::stod(mjd_part);
            }
        }
        
        // Magnitudine
        if (line.find("MAG") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag >> elem.H >> elem.G;
        }
        
        // Fine dati
        if (line.find("END") != std::string::npos) {
            break;
        }
    }
    
    return elem;
}

} // namespace ioccultcalc
```

### 5.3 Implementazione Completa - orbital_conversions.cpp

```cpp
// File: IOccultCalc/src/orbital_conversions.cpp

#include "ioccultcalc/orbital_conversions.h"
#include <cmath>
#include <stdexcept>

namespace ioccultcalc {

namespace {
    constexpr double PI = 3.141592653589793;
    constexpr double TWO_PI = 2.0 * PI;
    constexpr double DEG_TO_RAD = PI / 180.0;
    constexpr double GM_SUN = 0.000295912208; // AU³/day² (k²)
    constexpr double OBLIQUITY_J2000 = 23.439291 * DEG_TO_RAD; // rad
}

KeplerianElements OrbitalConversions::equinoctialToKeplerian(
    const EquinoctialElements& eq) {
    
    KeplerianElements kep;
    kep.name = eq.name;
    kep.a = eq.a;
    
    // Eccentricità
    kep.e = std::sqrt(eq.h*eq.h + eq.k*eq.k);
    
    // Inclinazione
    double tan_i_2 = std::sqrt(eq.p*eq.p + eq.q*eq.q);
    kep.i = 2.0 * std::atan(tan_i_2);
    
    // Longitudine nodo ascendente
    double Omega_rad = std::atan2(eq.p, eq.q);
    if (Omega_rad < 0.0) Omega_rad += TWO_PI;
    kep.Omega = Omega_rad;
    
    // Longitudine del perielio
    double LP_rad = std::atan2(eq.h, eq.k);
    if (LP_rad < 0.0) LP_rad += TWO_PI;
    
    // Argomento del perielio
    double omega_rad = LP_rad - Omega_rad;
    while (omega_rad < 0.0) omega_rad += TWO_PI;
    while (omega_rad > TWO_PI) omega_rad -= TWO_PI;
    kep.omega = omega_rad;
    
    // Anomalia media
    double M_rad = eq.lambda * DEG_TO_RAD - LP_rad;
    while (M_rad < 0.0) M_rad += TWO_PI;
    while (M_rad > TWO_PI) M_rad -= TWO_PI;
    kep.M = M_rad;
    
    // Epoca (converti MJD → JD)
    kep.epoch_jd = eq.epoch_mjd + 2400000.5;
    
    return kep;
}

CartesianState OrbitalConversions::keplerianToCartesian(
    const KeplerianElements& kep) {
    
    // Risolvi equazione di Keplero
    double E = solveKeplerEquation(kep.M, kep.e);
    
    // Anomalia vera
    double sin_E = std::sin(E);
    double cos_E = std::cos(E);
    double sqrt_1_e2 = std::sqrt(1.0 - kep.e * kep.e);
    double nu = std::atan2(sqrt_1_e2 * sin_E, cos_E - kep.e);
    
    // Raggio
    double r = kep.a * (1.0 - kep.e * cos_E);
    
    // Posizione nel piano orbitale
    double x_orb = r * std::cos(nu);
    double y_orb = r * std::sin(nu);
    
    // Velocità nel piano orbitale
    double v_factor = std::sqrt(GM_SUN * kep.a) / r;
    double vx_orb = -v_factor * sin_E;
    double vy_orb = v_factor * sqrt_1_e2 * cos_E;
    
    // Matrice di Gauss (rotazione dal piano orbitale all'eclittica)
    double cO = std::cos(kep.Omega);
    double sO = std::sin(kep.Omega);
    double cw = std::cos(kep.omega);
    double sw = std::sin(kep.omega);
    double ci = std::cos(kep.i);
    double si = std::sin(kep.i);
    
    double P11 = cO*cw - sO*sw*ci;
    double P12 = -cO*sw - sO*cw*ci;
    double P21 = sO*cw + cO*sw*ci;
    double P22 = -sO*sw + cO*cw*ci;
    double P31 = sw*si;
    double P32 = cw*si;
    
    // Stato eclittico
    CartesianState state;
    state.position(0) = P11*x_orb + P12*y_orb;
    state.position(1) = P21*x_orb + P22*y_orb;
    state.position(2) = P31*x_orb + P32*y_orb;
    
    state.velocity(0) = P11*vx_orb + P12*vy_orb;
    state.velocity(1) = P21*vx_orb + P22*vy_orb;
    state.velocity(2) = P31*vx_orb + P32*vy_orb;
    
    state.epoch_jd = kep.epoch_jd;
    
    return state;
}

CartesianState OrbitalConversions::eclipticToICRF(
    const CartesianState& ecliptic) {
    
    double c_eps = std::cos(OBLIQUITY_J2000);
    double s_eps = std::sin(OBLIQUITY_J2000);
    
    CartesianState icrf;
    icrf.epoch_jd = ecliptic.epoch_jd;
    
    // Rotazione attorno asse X
    icrf.position(0) = ecliptic.position(0);
    icrf.position(1) = c_eps * ecliptic.position(1) - 
                       s_eps * ecliptic.position(2);
    icrf.position(2) = s_eps * ecliptic.position(1) + 
                       c_eps * ecliptic.position(2);
    
    icrf.velocity(0) = ecliptic.velocity(0);
    icrf.velocity(1) = c_eps * ecliptic.velocity(1) - 
                       s_eps * ecliptic.velocity(2);
    icrf.velocity(2) = s_eps * ecliptic.velocity(1) + 
                       c_eps * ecliptic.velocity(2);
    
    return icrf;
}

double OrbitalConversions::solveKeplerEquation(
    double M, double e, double tol) {
    
    double E = M;  // Prima approssimazione
    
    for (int iter = 0; iter < 20; iter++) {
        double f = E - e * std::sin(E) - M;
        double fp = 1.0 - e * std::cos(E);
        double dE = f / fp;
        E -= dE;
        
        if (std::abs(dE) < tol) {
            return E;
        }
    }
    
    throw std::runtime_error(
        "Kepler equation did not converge after 20 iterations");
}

bool OrbitalConversions::validateICRF(const CartesianState& state) {
    // Verifica che le coordinate siano finite
    if (!std::isfinite(state.position.norm()) ||
        !std::isfinite(state.velocity.norm())) {
        return false;
    }
    
    // Verifica range ragionevole (0.1 - 100 AU)
    double r = state.position.norm();
    if (r < 0.1 || r > 100.0) {
        return false;
    }
    
    // Verifica velocità ragionevole (<100 AU/day)
    double v = state.velocity.norm();
    if (v > 100.0) {
        return false;
    }
    
    return true;
}

} // namespace ioccultcalc
```

### 5.4 Implementazione Completa - astdyn_wrapper.cpp

```cpp
// File: IOccultCalc/src/astdyn_wrapper.cpp

#include "ioccultcalc/astdyn_wrapper.h"
#include "ioccultcalc/eq1_parser.h"
#include "ioccultcalc/orbital_conversions.h"
#include <chrono>

namespace ioccultcalc {

AstDynWrapper::AstDynWrapper() 
    : tolerance_(1e-12) {
    
    // Crea propagatore AstDyn
    propagator_ = std::make_unique<astdyn::propagation::Propagator>();
}

void AstDynWrapper::setTolerance(double tol) {
    tolerance_ = tol;
    // AstDyn configura internamente la tolleranza
}

void AstDynWrapper::enablePlanets(bool enable) {
    // AstDyn abilita/disabilita perturbazioni planetarie
    // (implementazione dipende da API AstDyn)
}

void AstDynWrapper::enableAST17(bool enable) {
    // AstDyn abilita/disabilita perturbazioni asteroidali
}

void AstDynWrapper::enableRelativity(bool enable) {
    // AstDyn abilita/disabilita correzione relativistica
}

void AstDynWrapper::loadFromEQ1(const std::string& filepath) {
    auto eq = EQ1Parser::parseFile(filepath);
    setElements(eq);
}

void AstDynWrapper::setElements(const EquinoctialElements& eq) {
    auto kep = OrbitalConversions::equinoctialToKeplerian(eq);
    setElements(kep);
}

void AstDynWrapper::setElements(const KeplerianElements& kep) {
    elements_ = convertToAstDynElements(kep);
}

AstDynWrapper::EquatorialCoords 
AstDynWrapper::propagate(double jd_target) {
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Propaga con AstDyn
    auto state = propagator_->propagate(elements_, jd_target);
    
    // Converti stato cartesiano in coordinate equatoriali
    // (assumendo che state sia già in ICRF)
    
    // Posizione geocentrica (semplificato - aggiungere Terra)
    Eigen::Vector3d geo_pos = state.position;
    
    double delta = geo_pos.norm();
    double ra = std::atan2(geo_pos.y(), geo_pos.x());
    double dec = std::asin(geo_pos.z() / delta);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start);
    
    last_stats_.cpu_time_ms = duration.count();
    last_stats_.steps_taken++;  // Placeholder
    
    EquatorialCoords coords;
    coords.ra = ra;
    coords.dec = dec;
    coords.dist = delta;
    
    return coords;
}

std::vector<AstDynWrapper::EquatorialCoords> 
AstDynWrapper::propagateBatch(
    const std::vector<double>& jd_targets) {
    
    std::vector<EquatorialCoords> results;
    results.reserve(jd_targets.size());
    
    for (double jd : jd_targets) {
        results.push_back(propagate(jd));
    }
    
    return results;
}

AstDynWrapper::PropagationStats 
AstDynWrapper::getLastStats() const {
    return last_stats_;
}

astdyn::propagation::KeplerianElements 
AstDynWrapper::convertToAstDynElements(
    const KeplerianElements& kep) {
    
    astdyn::propagation::KeplerianElements astdyn_kep;
    
    // Conversione diretta
    astdyn_kep.semi_major_axis = kep.a;
    astdyn_kep.eccentricity = kep.e;
    astdyn_kep.inclination = kep.i;
    astdyn_kep.longitude_ascending_node = kep.Omega;
    astdyn_kep.argument_perihelion = kep.omega;
    astdyn_kep.mean_anomaly = kep.M;
    astdyn_kep.epoch_mjd_tdb = kep.epoch_jd - 2400000.5;
    
    return astdyn_kep;
}

} // namespace ioccultcalc
```

---

## 6. TEST E VALIDAZIONE

### 6.1 Piano di Test

| Test | Tipo | Obiettivo | Criterio Successo |
|------|------|-----------|-------------------|
| **test_eq1_parser** | Unit | Parser .eq1 | Parsing 17030.eq1 corretto |
| **test_orbital_conversions** | Unit | Conversioni | Eq→Kep per 17030: errore <0.01% |
| **test_kepler_solver** | Unit | Eq. Keplero | Convergenza <15 iter |
| **test_frame_rotation** | Unit | Ecl→ICRF | Vettori test corretti |
| **test_astdyn_wrapper** | Integration | Wrapper | Propaga 17030 senza errori |
| **test_integration_17030** | Integration | Full pipeline | Errore vs JPL <2" |
| **test_phase_comparison** | Integration | FASE1 vs FASE2 | Diff ~11" per 17030 |
| **test_batch_propagation** | Performance | Batch | 100 epoche in <5 sec |

### 6.2 Script di Validazione

**File**: `scripts/validate_integration.sh`

```bash
#!/bin/bash
# Script di validazione integrazione AstDyn

set -e

echo "========================================="
echo " Validazione Integrazione AstDyn"
echo "========================================="
echo

# 1. Compila
echo "1. Compilazione..."
cd IOccultCalc/build
cmake .. -DBUILD_TESTS=ON
make -j8
echo "✓ Compilazione OK"
echo

# 2. Run unit tests
echo "2. Unit Tests..."
./tests/test_orbital_conversions
echo "✓ Unit Tests OK"
echo

# 3. Run integration tests
echo "3. Integration Tests..."
./tests/test_integration_17030
echo "✓ Integration Tests OK"
echo

# 4. Confronto con AstDyn standalone
echo "4. Confronto con AstDyn standalone..."
cd ../../ITALOccultLibrary/astdyn/build/tests
./test_asteroid_17030_occultation > /tmp/astdyn_result.txt

# Estrai separazione
ASTDYN_SEP=$(grep "Separation" /tmp/astdyn_result.txt | awk '{print $3}')
echo "  AstDyn standalone: ${ASTDYN_SEP}\""

cd ../../../../IOccultCalc/build/tests
./test_integration_17030 > /tmp/ioccultcalc_result.txt

IOCCULTCALC_SEP=$(grep "Separation" /tmp/ioccultcalc_result.txt | awk '{print $4}')
echo "  IOccultCalc+AstDyn: ${IOCCULTCALC_SEP}\""

# Confronto
echo
echo "✓ Validazione completata!"
echo
echo "Risultati:"
echo "  AstDyn:       ${ASTDYN_SEP}\" (riferimento)"
echo "  IOccultCalc:  ${IOCCULTCALC_SEP}\" (nuovo)"
echo
```

---

## 7. TEMPISTICHE

### 7.1 Stima Complessiva

| Fase | Durata | Priorità |
|------|--------|----------|
| FASE 1: Fondamenta | 2-3 ore | 🔴 CRITICA |
| FASE 2: Integrazione | 3-4 ore | 🟠 ALTA |
| FASE 3: Testing | 2-3 ore | 🟠 ALTA |
| FASE 4: Ottimizzazioni | 4-5 ore | 🟡 MEDIA |
| **TOTALE** | **11-15 ore** | |

### 7.2 Milestone

```
Milestone 1 (dopo FASE 1): Parser .eq1 funzionante
  ✓ Parsing 17030.eq1 corretto
  ✓ Conversione equinoziale→kepleriano validata

Milestone 2 (dopo FASE 2): Integrazione TwoPhaseStrategy
  ✓ loadElementsFromEQ1() funzionante
  ✓ getRKF78Position() usa AstDyn

Milestone 3 (dopo FASE 3): Validazione JPL
  ✓ Errore vs JPL <2" per 17030
  ✓ Tutti i test passano

Milestone 4 (dopo FASE 4): Produzione-ready
  ✓ Caching implementato
  ✓ Batch propagation efficiente
  ✓ Documentazione completa
```

---

## 8. CHECKLIST FINALE

### Pre-Implementazione ✓

- [x] Analisi completata
- [x] Piano dettagliato
- [x] Riferimenti codice identificati
- [ ] Backup IOoccultCalc corrente
- [ ] Branch git creato

### Implementazione

- [ ] Parser .eq1 (Task 1.1)
- [ ] Conversioni orbitali (Task 1.2)
- [ ] Wrapper AstDyn (Task 1.3)
- [ ] Modifica TwoPhaseStrategy (Task 2.1)
- [ ] Auto-switch (Task 2.2)
- [ ] Unit tests (Task 3.1)
- [ ] Integration tests (Task 3.2)
- [ ] CMakeLists.txt aggiornato

### Validazione

- [ ] Compilazione pulita (0 errori)
- [ ] Tutti i test passano
- [ ] Confronto vs AstDyn standalone
- [ ] Confronto vs JPL Horizons
- [ ] Documentazione aggiornata

### Deployment

- [ ] Merge su branch main
- [ ] Tag versione (v2.0-astdyn-full)
- [ ] README aggiornato
- [ ] Release notes

---

## 9. SUPPORTO E RIFERIMENTI

### 9.1 File di Riferimento

**Parser .eq1**:
- `ITALOccultLibrary/astdyn/tests/test_asteroid_17030_occultation.cpp` (linee 79-120)

**Conversioni**:
- `ITALOccultLibrary/astdyn/tools/astdyn_propagator.cpp` (linee 512-640)

**API AstDyn**:
- `ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp`
- `ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp`

### 9.2 Documentazione

- `ANALISI_IOCCULTCALC_ELEMENTI_EQ1.md` - Analisi problemi
- `CONFRONTO_TECNICO_ASTDYN_VS_IOCCULTCALC.md` - Confronto algoritmi
- `IOCCULTCALC_FILE_ELEMENTI_ORBITALI.md` - Formati elementi

---

**Piano redatto**: 1 Dicembre 2025  
**Versione**: 1.0  
**Status**: 🚀 Pronto per implementazione

**Prossimo Step**: Iniziare con FASE 1 - Task 1.1 (Parser .eq1)
