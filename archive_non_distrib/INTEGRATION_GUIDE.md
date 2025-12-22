# Guida Integrazione AstDyn in IOccultCalc

**Versione**: 1.0  
**Data**: 1 Dicembre 2025  
**Autore**: IOccultCalc Integration Team  

---

## 📋 Indice

1. [Prerequisiti](#prerequisiti)
2. [Panoramica Integrazione](#panoramica-integrazione)
3. [STEP 1: Copia Template Files](#step-1-copia-template-files)
4. [STEP 2: Crea Esempio Integrazione](#step-2-crea-esempio-integrazione)
5. [STEP 3: Modifica propagation_strategy](#step-3-modifica-propagation_strategy)
6. [STEP 4: Aggiorna CMakeLists.txt](#step-4-aggiorna-cmakeliststxt)
7. [STEP 5: Build e Test](#step-5-build-e-test)
8. [STEP 6: Validazione con 17030](#step-6-validazione-con-17030)
9. [Troubleshooting](#troubleshooting)

---

## Prerequisiti

### Software Richiesto
- ✅ **CMake** 3.15+
- ✅ **C++17 compiler** (GCC 8+, Clang 7+, MSVC 2019+)
- ✅ **Eigen3** 3.3+ (`brew install eigen` su macOS)
- ✅ **AstDyn library** compilata e installata
- ✅ **Google Test** (per unit tests)

### Repository
```bash
# IOccultCalc repository
cd ~/VisualStudioCode/GitHub/IOccultCalc

# ITALOccultLibrary repository (templates)
cd ~/VisualStudioCode/GitHub/ITALOccultLibrary
```

### Verifica Installazioni
```bash
# Verifica Eigen3
pkg-config --modversion eigen3
# Expected: 3.3.x o superiore

# Verifica AstDyn
ls /usr/local/include/astdyn/
# Expected: AstDynPropagator.hpp, AstDynTypes.hpp, ...

# Verifica CMake
cmake --version
# Expected: 3.15.x o superiore
```

---

## Panoramica Integrazione

### Architettura Attuale IOccultCalc

```
IOccultCalc/
├── include/ioccultcalc/
│   ├── propagation_strategy.h    ← DA MODIFICARE
│   └── ...
├── src/
│   ├── propagation_strategy.cpp  ← DA MODIFICARE
│   └── ...
└── CMakeLists.txt                ← DA MODIFICARE
```

### Architettura Dopo Integrazione

```
IOccultCalc/
├── include/ioccultcalc/
│   ├── propagation_strategy.h    ← MODIFICATO
│   ├── eq1_parser.h              ← NUOVO
│   ├── orbital_conversions.h     ← NUOVO
│   └── astdyn_wrapper.h          ← NUOVO
├── src/
│   ├── propagation_strategy.cpp  ← MODIFICATO
│   ├── eq1_parser.cpp            ← NUOVO
│   ├── orbital_conversions.cpp   ← NUOVO
│   └── astdyn_wrapper.cpp        ← NUOVO
└── CMakeLists.txt                ← MODIFICATO
```

### Flusso di Propagazione

**PRIMA** (errore ~12"):
```
API → Keplerian propagation → RKF78 (senza perturbazioni complete)
```

**DOPO** (errore <2"):
```
.eq1 file → eq1_parser → orbital_conversions → astdyn_wrapper → RKF78 + 11 perturbations
```

---

## STEP 1: Copia Template Files

### 1.1 Copia Headers

```bash
cd ~/VisualStudioCode/GitHub/ITALOccultLibrary

# Copia headers
cp templates_ioccultcalc/include/eq1_parser.h \
   ../IOccultCalc/include/ioccultcalc/

cp templates_ioccultcalc/include/orbital_conversions.h \
   ../IOccultCalc/include/ioccultcalc/

cp templates_ioccultcalc/include/astdyn_wrapper.h \
   ../IOccultCalc/include/ioccultcalc/
```

### 1.2 Copia Source Files

```bash
# Copia implementations
cp templates_ioccultcalc/src/eq1_parser.cpp \
   ../IOccultCalc/src/

cp templates_ioccultcalc/src/orbital_conversions.cpp \
   ../IOccultCalc/src/

cp templates_ioccultcalc/src/astdyn_wrapper.cpp \
   ../IOccultCalc/src/
```

### 1.3 Verifica Copia

```bash
cd ../IOccultCalc

# Verifica headers
ls -lh include/ioccultcalc/*.h | grep -E "(eq1|orbital|astdyn)"
# Expected: eq1_parser.h, orbital_conversions.h, astdyn_wrapper.h

# Verifica sources
ls -lh src/*.cpp | grep -E "(eq1|orbital|astdyn)"
# Expected: eq1_parser.cpp, orbital_conversions.cpp, astdyn_wrapper.cpp
```

---

## STEP 2: Crea Esempio Integrazione

Prima di modificare IOccultCalc, creiamo un esempio standalone per testare l'integrazione.

### 2.1 File di Test: `test_astdyn_integration_standalone.cpp`

Questo file sarà creato nel workspace ITALOccultLibrary come esempio.

**Path**: `ITALOccultLibrary/examples/test_astdyn_integration_standalone.cpp`

Contenuto (verrà creato nel prossimo step).

---

## STEP 3: Modifica propagation_strategy

### 3.1 Modifica Header: `include/ioccultcalc/propagation_strategy.h`

**Aggiungi includes** all'inizio del file (dopo gli includes esistenti):

```cpp
// Includes esistenti...
#include <memory>
#include <vector>

// NUOVO: Aggiungi questi includes
#include "eq1_parser.h"
#include "orbital_conversions.h"
#include "astdyn_wrapper.h"
```

**Aggiungi membri privati** alla classe `TwoPhaseStrategy`:

```cpp
class TwoPhaseStrategy {
public:
    // Metodi esistenti...
    
private:
    // Membri esistenti...
    
    // NUOVO: Aggiungi questi membri
    std::unique_ptr<AstDynWrapper> astdyn_wrapper_;
    
    /**
     * @brief Converti elementi eq1 in stato cartesiano ICRF
     * @param eq Elementi equinoziali
     * @return Stato cartesiano in frame ICRF
     */
    CartesianState convertEq1ToICRF(const EquinoctialElements& eq);
    
    /**
     * @brief Propaga usando AstDyn completo
     * @param initial_state Stato iniziale ICRF
     * @param target_epoch_jd Epoca target (JD TDB)
     * @return Risultato propagazione
     */
    PropagationResult propagateWithAstDyn(
        const CartesianState& initial_state,
        double target_epoch_jd
    );
    
    /**
     * @brief Propaga da file .eq1
     * @param eq1_file Path al file .eq1
     * @param target_epoch_jd Epoca target (JD TDB)
     * @return Risultato propagazione
     */
    PropagationResult propagateFromEq1File(
        const std::string& eq1_file,
        double target_epoch_jd
    );
};
```

### 3.2 Modifica Implementation: `src/propagation_strategy.cpp`

**Nel costruttore**, inizializza il wrapper AstDyn:

```cpp
TwoPhaseStrategy::TwoPhaseStrategy(/* parametri esistenti */)
    : /* inizializzazioni esistenti */
{
    // Inizializzazioni esistenti...
    
    // NUOVO: Inizializza wrapper AstDyn con configurazione JPL-compliant
    astdyn_wrapper_ = std::make_unique<AstDynWrapper>(
        AstDynConfig::jplCompliant()
    );
}
```

**Aggiungi implementazione metodi**, alla fine del file:

```cpp
// ============================================================================
// NUOVO: Conversione eq1 → ICRF
// ============================================================================

CartesianState TwoPhaseStrategy::convertEq1ToICRF(
    const EquinoctialElements& eq)
{
    // 1. Converti equinoziale → kepleriano
    auto kep = OrbitalConversions::equinoctialToKeplerian(eq);
    
    // 2. Converti kepleriano → cartesiano (eclittico)
    auto cart_ecliptic = OrbitalConversions::keplerianToCartesian(kep);
    
    // 3. Ruota da eclittico a ICRF
    auto cart_icrf = OrbitalConversions::eclipticToICRF(cart_ecliptic);
    
    // 4. Valida risultato
    if (!OrbitalConversions::validateICRF(cart_icrf)) {
        throw std::runtime_error("Invalid ICRF state after conversion");
    }
    
    return cart_icrf;
}

// ============================================================================
// NUOVO: Propagazione con AstDyn
// ============================================================================

PropagationResult TwoPhaseStrategy::propagateWithAstDyn(
    const CartesianState& initial_state,
    double target_epoch_jd)
{
    if (!astdyn_wrapper_) {
        throw std::runtime_error("AstDyn wrapper not initialized");
    }
    
    return astdyn_wrapper_->propagate(
        initial_state.position,
        initial_state.velocity,
        initial_state.epoch_jd,
        target_epoch_jd
    );
}

// ============================================================================
// NUOVO: Propagazione da file .eq1
// ============================================================================

PropagationResult TwoPhaseStrategy::propagateFromEq1File(
    const std::string& eq1_file,
    double target_epoch_jd)
{
    try {
        // 1. Parse file .eq1
        auto eq_elements = EQ1Parser::parseFile(eq1_file);
        
        // 2. Converti in ICRF
        auto icrf_state = convertEq1ToICRF(eq_elements);
        
        // 3. Propaga con AstDyn
        auto result = propagateWithAstDyn(icrf_state, target_epoch_jd);
        
        return result;
        
    } catch (const EQ1ParseException& e) {
        PropagationResult error_result;
        error_result.success = false;
        error_result.error_message = std::string("EQ1 parse error: ") + e.what();
        return error_result;
        
    } catch (const std::exception& e) {
        PropagationResult error_result;
        error_result.success = false;
        error_result.error_message = std::string("Propagation error: ") + e.what();
        return error_result;
    }
}
```

---

## STEP 4: Aggiorna CMakeLists.txt

### 4.1 Modifica `IOccultCalc/CMakeLists.txt`

**Trova la sezione** `set(IOCCULTCALC_SOURCES ...)` e aggiungi:

```cmake
set(IOCCULTCALC_SOURCES
    # Sources esistenti...
    src/propagation_strategy.cpp
    # ... altri ...
    
    # NUOVO: Aggiungi questi source files
    src/eq1_parser.cpp
    src/orbital_conversions.cpp
    src/astdyn_wrapper.cpp
)
```

**Aggiungi dipendenze** Eigen3 e AstDyn:

```cmake
# NUOVO: Find Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# NUOVO: Find AstDyn
find_package(AstDyn REQUIRED)

# Target library
add_library(IOccultCalc ${IOCCULTCALC_SOURCES})

# NUOVO: Aggiungi include directories
target_include_directories(IOccultCalc
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# NUOVO: Link contro Eigen3 e AstDyn
target_link_libraries(IOccultCalc
    PUBLIC
        AstDyn::AstDyn
        Eigen3::Eigen
    PRIVATE
        # Link esistenti...
)

# NUOVO: Richiedi C++17
target_compile_features(IOccultCalc PUBLIC cxx_std_17)
```

### 4.2 Verifica CMakeLists.txt

```bash
cd IOccultCalc
cat CMakeLists.txt | grep -A5 "eq1_parser\|orbital_conversions\|astdyn_wrapper"
# Dovresti vedere i nuovi source files
```

---

## STEP 5: Build e Test

### 5.1 Clean Build

```bash
cd IOccultCalc

# Rimuovi build precedente
rm -rf build
mkdir build
cd build
```

### 5.2 Configure con CMake

```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr/local
```

**Output atteso**:
```
-- Found Eigen3: 3.3.x
-- Found AstDyn: /usr/local/lib/cmake/AstDyn
-- Configuring done
-- Generating done
```

### 5.3 Compile

```bash
make -j8
```

**Verifica errori**:
- ❌ Se errori di linking AstDyn → verifica installazione
- ❌ Se errori Eigen3 → `brew install eigen`
- ❌ Se errori C++17 → aggiorna compiler

### 5.4 Run Tests

```bash
ctest --output-on-failure
```

---

## STEP 6: Validazione con 17030

### 6.1 Prepara File di Test

```bash
# Copia file .eq1 di test
cp ~/VisualStudioCode/GitHub/ITALOccultLibrary/astdyn/data/17030.eq1 \
   ./test_data/
```

### 6.2 Crea Test Standalone

Crea file `test_17030_validation.cpp`:

```cpp
#include <iostream>
#include <iomanip>
#include "ioccultcalc/propagation_strategy.h"

int main() {
    try {
        // Setup
        ioccultcalc::TwoPhaseStrategy strategy;
        
        // Epoca occultazione: 26 Nov 2025, 06:30 UTC
        double target_jd = 2460643.77083;
        
        // Propaga da file .eq1
        auto result = strategy.propagateFromEq1File(
            "test_data/17030.eq1",
            target_jd
        );
        
        if (!result.success) {
            std::cerr << "Error: " << result.error_message << std::endl;
            return 1;
        }
        
        // Output risultato
        std::cout << std::fixed << std::setprecision(12);
        std::cout << "=== Asteroid 17030 Propagation ===" << std::endl;
        std::cout << "Target epoch: " << target_jd << " JD" << std::endl;
        std::cout << "Position (ICRF, AU):" << std::endl;
        std::cout << "  x = " << result.position(0) << std::endl;
        std::cout << "  y = " << result.position(1) << std::endl;
        std::cout << "  z = " << result.position(2) << std::endl;
        std::cout << "Computation time: " << result.computation_time_ms << " ms" << std::endl;
        std::cout << "Steps: " << result.num_steps << std::endl;
        
        // TODO: Confronta con JPL Horizons
        // Expected accuracy: < 2 arcsec
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
}
```

### 6.3 Compila e Run

```bash
# Aggiungi al CMakeLists.txt
add_executable(test_17030_validation test_17030_validation.cpp)
target_link_libraries(test_17030_validation IOccultCalc)

# Rebuild
make test_17030_validation

# Run
./test_17030_validation
```

**Output atteso**:
```
=== Asteroid 17030 Propagation ===
Target epoch: 2460643.770833333 JD
Position (ICRF, AU):
  x = 1.234567890123
  y = -0.987654321098
  z = 0.123456789012
Computation time: 12.345 ms
Steps: 156
```

---

## Troubleshooting

### Errore: "AstDyn not found"

```bash
# Verifica installazione AstDyn
ls /usr/local/lib/cmake/AstDyn/
# Se vuoto, reinstalla AstDyn:
cd ~/VisualStudioCode/GitHub/ITALOccultLibrary/astdyn/build
sudo make install
```

### Errore: "Eigen3 not found"

```bash
# macOS
brew install eigen

# Linux
sudo apt-get install libeigen3-dev

# Verifica
pkg-config --modversion eigen3
```

### Errore: "undefined reference to AstDynPropagator"

Verifica `CMakeLists.txt`:
```cmake
target_link_libraries(IOccultCalc
    PUBLIC AstDyn::AstDyn  # ← Deve essere presente
)
```

### Errore: "eq1_parser.h: No such file"

Verifica copia files:
```bash
ls -l include/ioccultcalc/eq1_parser.h
# Se mancante, ripeti STEP 1
```

### Errore Compilazione C++17

```cmake
# In CMakeLists.txt
target_compile_features(IOccultCalc PUBLIC cxx_std_17)
```

### Performance Basse

Se propagazione > 50ms:
```cpp
// Usa configurazione fast per screening
auto wrapper = std::make_unique<AstDynWrapper>(
    AstDynConfig::fast()
);
```

---

## Checklist Finale

Prima di committare le modifiche:

- [ ] ✅ Files copiati (3 headers + 3 sources)
- [ ] ✅ `propagation_strategy.h` modificato
- [ ] ✅ `propagation_strategy.cpp` modificato
- [ ] ✅ `CMakeLists.txt` aggiornato
- [ ] ✅ Build completa senza errori
- [ ] ✅ Test 17030 eseguito con successo
- [ ] ✅ Errore < 2" rispetto JPL Horizons
- [ ] ✅ Performance accettabile (< 20ms/propagation)

---

## Prossimi Passi

Dopo integrazione base:

1. **FASE 3**: Creare unit tests completi
2. **FASE 4**: Implementare ottimizzazioni (cache, batch, parallel)
3. **Deployment**: Integrare in pipeline produzione IOccultCalc

---

## Supporto

Per problemi o domande:
- Consulta `ANALISI_IOCCULTCALC_ELEMENTI_EQ1.md`
- Consulta `PIANO_INTEGRAZIONE_ASTDYN_IOCCULTCALC.md`
- Verifica esempi in `ITALOccultLibrary/examples/`

---

**Fine Guida Integrazione**  
*Versione 1.0 - 1 Dicembre 2025*
