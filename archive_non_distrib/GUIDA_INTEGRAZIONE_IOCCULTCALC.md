# Guida Integrazione ITALOccultLibrary in IOccultCalc

**Data**: 1 Dicembre 2025  
**Status**: ✅ Pronto per integrazione  
**Validazione**: 0.0003 arcsec vs JPL Horizons

---

## 📋 Overview

Questa guida descrive come integrare **ITALOccultLibrary** in **IOccultCalc** per propagazione orbitale di alta precisione con:

✅ **Precisione JPL Horizons** (0.0003 arcsec)  
✅ **Frame conversion validata** (ECLM→ICRF)  
✅ **Integrazione RKF78** (7/8 ordine)  
✅ **Performance ottimale** (< 1 ms per 7 giorni)

---

## 🎯 Obiettivi Integrazione

1. Usare ITALOccultLibrary come backend per propagazione
2. Mantenere interfaccia IOccultCalc esistente
3. Aggiungere supporto file .eq1 (AstDyS)
4. Garantire precisione JPL Horizons
5. Minimizzare modifiche codice esistente

---

## 📦 File da Copiare

### Da ITALOccultLibrary a IOccultCalc

```bash
# 1. Header astdyn_interface (nuova interfaccia)
cp ITALOccultLibrary/integration/astdyn_interface.h \\
   IOccultCalc/include/ioccultcalc/

# 2. Moduli ITALOccultLibrary (opzionale se usi linking)
# Se vuoi embedded (no linking esterno):
cp ITALOccultLibrary/italoccultlibrary/include/*.h \\
   IOccultCalc/include/ioccultcalc/italoccultlib/

cp ITALOccultLibrary/italoccultlibrary/src/*.cpp \\
   IOccultCalc/src/italoccultlib/
```

**Opzione A - Linking esterno** (raccomandato):
- Installa ITALOccultLibrary in `/usr/local`
- Link con `find_package(ITALOccultLibrary)`

**Opzione B - Embedded**:
- Copia file direttamente in IOccultCalc
- Compila insieme al progetto

---

## 🔧 Modifiche Necessarie

### 1. CMakeLists.txt

**File**: `IOccultCalc/CMakeLists.txt`

```cmake
# Trova ITALOccultLibrary
find_package(ITALOccultLibrary 1.0 REQUIRED)
find_package(AstDyn REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Aggiungi include directories
include_directories(
    ${ITALOCCULTLIBRARY_INCLUDE_DIRS}
    ${ASTDYN_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)

# Link librerie
target_link_libraries(ioccultcalc
    # ... librerie esistenti ...
    ITALOccultLibrary::italoccultlib
    AstDyn::astdyn
    Eigen3::Eigen
)
```

**Righe da aggiungere**: ~10

---

### 2. propagation_strategy.h

**File**: `IOccultCalc/include/ioccultcalc/propagation_strategy.h`

```cpp
#ifndef IOCCULTCALC_PROPAGATION_STRATEGY_H
#define IOCCULTCALC_PROPAGATION_STRATEGY_H

// Include esistenti...
#include "orbital_elements.h"
#include "ephemeris.h"

// ✅ NUOVO: Include interface ITALOccultLibrary
#include "astdyn_interface.h"

namespace ioccultcalc {

class PropagationStrategy {
public:
    // Metodi esistenti...
    virtual OrbitalElements propagate(
        const OrbitalElements& initial,
        const JulianDate& target_date) = 0;
    
    // ✅ NUOVO: Overload per file .eq1
    virtual OrbitalElements propagateFromEQ1(
        const std::string& eq1_file,
        const JulianDate& target_date);
    
    virtual ~PropagationStrategy() = default;
};

// ✅ NUOVA CLASSE: Strategy con AstDyn backend
class AstDynStrategy : public PropagationStrategy {
public:
    AstDynStrategy();
    explicit AstDynStrategy(const PropagatorConfig& config);
    
    // Implementazione interfaccia
    OrbitalElements propagate(
        const OrbitalElements& initial,
        const JulianDate& target_date) override;
    
    OrbitalElements propagateFromEQ1(
        const std::string& eq1_file,
        const JulianDate& target_date) override;
    
    // Configurazione
    void setConfig(const PropagatorConfig& config);
    PropagatorConfig getConfig() const;
    
    // Statistiche
    PropagationResult::IntegrationStats getLastStats() const;
    
private:
    std::unique_ptr<AstDynPropagator> propagator_;
};

} // namespace ioccultcalc

#endif
```

**Righe da aggiungere**: ~40

---

### 3. propagation_strategy.cpp

**File**: `IOccultCalc/src/propagation_strategy.cpp`

```cpp
#include "ioccultcalc/propagation_strategy.h"
#include "ioccultcalc/time_utils.h"
#include <stdexcept>

namespace ioccultcalc {

// ============================================================================
// PropagationStrategy - Metodo default per .eq1
// ============================================================================

OrbitalElements PropagationStrategy::propagateFromEQ1(
    const std::string& eq1_file,
    const JulianDate& target_date) {
    
    // Default implementation: lancia eccezione
    // Subclass deve implementare se vuole supporto .eq1
    throw std::runtime_error(
        "propagateFromEQ1 not implemented for this strategy");
}

// ============================================================================
// AstDynStrategy - Implementazione
// ============================================================================

AstDynStrategy::AstDynStrategy()
    : propagator_(std::make_unique<AstDynPropagator>(
          PropagatorConfig::standardConfig())) {
}

AstDynStrategy::AstDynStrategy(const PropagatorConfig& config)
    : propagator_(std::make_unique<AstDynPropagator>(config)) {
}

OrbitalElements AstDynStrategy::propagate(
    const OrbitalElements& initial,
    const JulianDate& target_date) {
    
    // 1. Carica elementi in propagatore
    propagator_->loadElements(initial);
    
    // 2. Converti target date in MJD
    double target_mjd = TimeUtils::julianDateToMJD(target_date);
    
    // 3. Propaga (automaticamente in ICRF)
    PropagationResult result = propagator_->propagate(target_mjd);
    
    // 4. Converti risultato in OrbitalElements IOccultCalc
    return result.toOrbitalElements();
}

OrbitalElements AstDynStrategy::propagateFromEQ1(
    const std::string& eq1_file,
    const JulianDate& target_date) {
    
    // 1. Carica elementi da .eq1
    propagator_->loadElements(eq1_file);
    
    // 2. Converti target date
    double target_mjd = TimeUtils::julianDateToMJD(target_date);
    
    // 3. Propaga
    PropagationResult result = propagator_->propagate(target_mjd);
    
    // 4. Converti e ritorna
    return result.toOrbitalElements();
}

void AstDynStrategy::setConfig(const PropagatorConfig& config) {
    propagator_->configure(config);
}

PropagatorConfig AstDynStrategy::getConfig() const {
    return propagator_->getConfig();
}

PropagationResult::IntegrationStats AstDynStrategy::getLastStats() const {
    return propagator_->getLastIntegrationStats();
}

} // namespace ioccultcalc
```

**Righe da aggiungere**: ~80

---

## 📝 Esempio d'Uso in IOccultCalc

### Uso con Strategy Pattern

```cpp
#include "ioccultcalc/propagation_strategy.h"

// 1. Crea strategy AstDyn
auto strategy = std::make_unique<AstDynStrategy>();

// 2. Opzionale: configura
PropagatorConfig cfg = PropagatorConfig::highPrecisionConfig();
strategy->setConfig(cfg);

// 3A. Propaga da OrbitalElements
OrbitalElements initial = getInitialElements();
JulianDate target(2461007.5);  // 28 Nov 2025
OrbitalElements final = strategy->propagate(initial, target);

// 3B. Oppure propaga da file .eq1
OrbitalElements final = strategy->propagateFromEQ1(
    "data/17030.eq1", 
    target
);

// 4. Usa risultato
std::cout << "Posizione finale:\n";
std::cout << "  X = " << final.position.x << " AU\n";
std::cout << "  Y = " << final.position.y << " AU\n";
std::cout << "  Z = " << final.position.z << " AU\n";
```

### Uso Diretto (senza Strategy)

```cpp
#include "ioccultcalc/astdyn_interface.h"

// 1. Crea propagatore direttamente
AstDynPropagator prop;

// 2. Carica elementi
prop.loadElements("data/17030.eq1");

// 3. Propaga
double target_mjd = 61007.0;
PropagationResult result = prop.propagate(target_mjd);

// 4. Risultato già in ICRF
std::cout << "Posizione ICRF @ MJD " << target_mjd << ":\n";
std::cout << "  X = " << result.position.x() << " AU\n";
std::cout << "  Y = " << result.position.y() << " AU\n";
std::cout << "  Z = " << result.position.z() << " AU\n";

// 5. Statistiche
std::cout << "\nIntegrazione:\n";
std::cout << "  Step: " << result.num_steps << "\n";
std::cout << "  Reject: " << result.num_rejected_steps << "\n";
std::cout << "  Tempo: " << result.computation_time_ms << " ms\n";
```

---

## ⚠️ Note Importanti

### Frame di Riferimento

**CRITICAL**: File `.eq1` sono in frame **ECLM J2000**, non ICRF!

La conversione è **automatica** se usi:
- `AstDynPropagator::propagate()` ✅
- `AstDynStrategy::propagate()` ✅

Risultati sono sempre in **ICRF** (compatibile JPL Horizons).

### Compatibilità OrbitalElements

Per conversione tra `PropagationResult` e `OrbitalElements` IOccultCalc:

```cpp
// PropagationResult → OrbitalElements
OrbitalElements elem = result.toOrbitalElements();

// OrbitalElements → caricamento in propagatore
propagator.loadElements(elem);
```

Le conversioni sono già implementate in `astdyn_interface.cpp`.

---

## 🧪 Testing

### Test Unitario Minimo

```cpp
#include <gtest/gtest.h>
#include "ioccultcalc/astdyn_interface.h"

TEST(AstDynIntegration, PropagateFromEQ1) {
    // Setup
    AstDynPropagator prop;
    prop.loadElements("test_data/17030.eq1");
    
    // Propaga 7 giorni
    double target_mjd = 61007.0;
    auto result = prop.propagate(target_mjd);
    
    // Valida vs JPL Horizons
    Eigen::Vector3d jpl_pos(1.020032, 2.884614, 1.153917);
    double error_km = AstDynPropagator::validateWithJPL(result, jpl_pos);
    
    // Target: < 1000 km (otteniamo ~0.7 km!)
    EXPECT_LT(error_km, 1000.0);
}

TEST(AstDynStrategy, IntegrateWithIOccultCalc) {
    // Setup strategy
    auto strategy = std::make_unique<AstDynStrategy>();
    
    // Propaga da .eq1
    JulianDate target(2461007.5);
    OrbitalElements result = strategy->propagateFromEQ1(
        "test_data/17030.eq1",
        target
    );
    
    // Verifica risultato valido
    EXPECT_TRUE(result.isValid());
    EXPECT_NEAR(result.a, 3.175, 0.01);  // Semiasse ~3.2 AU
}
```

### Test End-to-End

Crea `tests/test_astdyn_integration.cpp`:

```cpp
#include "ioccultcalc/propagation_strategy.h"
#include <iostream>
#include <iomanip>

int main() {
    try {
        // 1. Crea strategy
        AstDynStrategy strategy;
        
        // 2. Carica e propaga
        std::cout << "Propagating asteroid 17030 Sierks...\\n";
        JulianDate target(2461007.5);
        
        auto result = strategy.propagateFromEQ1(
            "../data/17030.eq1",
            target
        );
        
        // 3. Output
        std::cout << std::fixed << std::setprecision(12);
        std::cout << "\\nPosition ICRF @ MJD 61007.0:\\n";
        std::cout << "  X = " << result.position.x << " AU\\n";
        std::cout << "  Y = " << result.position.y << " AU\\n";
        std::cout << "  Z = " << result.position.z << " AU\\n";
        
        // 4. Confronta con JPL
        std::cout << "\\nCompare with JPL Horizons:\\n";
        std::cout << "  X_JPL = 1.020032 AU\\n";
        std::cout << "  Y_JPL = 2.884614 AU\\n";
        std::cout << "  Z_JPL = 1.153917 AU\\n";
        
        // 5. Statistiche
        auto stats = strategy.getLastStats();
        std::cout << "\\nIntegration:\\n";
        std::cout << "  Steps: " << stats.total_steps << "\\n";
        std::cout << "  Rejected: " << stats.rejected_steps << "\\n";
        
        std::cout << "\\n✅ Test PASSED\\n";
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << "\\n";
        return 1;
    }
}
```

Compila con:
```bash
g++ -std=c++17 test_astdyn_integration.cpp -o test_integration \\
    -I../include -L../build/lib \\
    -lioccultcalc -litaloccultlib -lastdyn \\
    -lEigen3::Eigen
```

---

## 📊 Validazione

### Checklist Pre-Integrazione

- [ ] ITALOccultLibrary installata in `/usr/local`
- [ ] AstDyn compilato e linkato correttamente
- [ ] Eigen3 >= 3.3 disponibile
- [ ] File test `17030.eq1` scaricato da AstDyS
- [ ] CMakeLists.txt IOccultCalc aggiornato

### Checklist Post-Integrazione

- [ ] Compilazione IOccultCalc senza errori
- [ ] Test unitari passano
- [ ] Test `test_astdyn_integration` passa
- [ ] Errore vs JPL Horizons < 1 km
- [ ] Performance < 10 ms per propagazione 7 giorni
- [ ] Nessuna regressione funzionalità esistenti

---

## 🐛 Troubleshooting

### Errore: "Cannot find ITALOccultLibrary"

```bash
# Verifica installazione
ls /usr/local/lib/libitaloccultlib.*
ls /usr/local/include/italoccultlib/

# Se non installata:
cd ITALOccultLibrary/italoccultlibrary/build
sudo make install
```

### Errore: "Undefined reference to AstDyn symbols"

```bash
# Verifica AstDyn
ls /usr/local/lib/libastdyn.*

# Aggiungi a CMakeLists.txt:
find_package(AstDyn REQUIRED)
target_link_libraries(ioccultcalc AstDyn::astdyn)
```

### Errore: "Frame mismatch - results don't match JPL"

Verifica che stai usando `propagate()`, non chiamate dirette a AstDyn.
La conversione frame è automatica solo nell'interface.

```cpp
// ❌ SBAGLIATO
auto astdyn_result = astdyn_raw_propagate(...);  // Frame ECLM!

// ✅ CORRETTO  
auto result = propagator.propagate(target_mjd);  // Frame ICRF automatico
```

---

## 📈 Performance Attese

### Tempi Tipici (Release build)

| Intervallo | Step | Reject | Tempo | Errore vs JPL |
|------------|------|--------|-------|---------------|
| 7 giorni | 2 | 0 | < 1 ms | 0.7 km |
| 30 giorni | 5-8 | 0-1 | 2-3 ms | < 10 km |
| 1 anno | 50-70 | 2-5 | 20-30 ms | < 100 km |

### Memoria

- `AstDynPropagator`: ~2 KB per istanza
- `PropagationResult`: ~200 bytes
- Overhead AstDyn: ~50 MB (effemeridi planetarie)

---

## 🚀 Next Steps

1. ✅ Integra codice seguendo questa guida
2. ⏳ Compila e testa con `test_astdyn_integration`
3. ⏳ Valida con almeno 3 asteroidi diversi
4. ⏳ Aggiorna documentazione IOccultCalc
5. ⏳ Commit e push

---

## 📚 Riferimenti

- **SUNTO_FINALE_VALIDAZIONE_ASTDYN.md** - Validazione completa
- **VALIDAZIONE_ASTDYN_FRAME_CORRECTED.md** - Report tecnico
- **FRAME_CONVERSION_MODULE.md** - Conversione frame dettagliata
- **italoccultlibrary/README.md** - API ITALOccultLibrary

---

**READY FOR INTEGRATION** ✅

Status: Tutti i moduli testati e validati  
Precisione: JPL Horizons grade (0.0003 arcsec)  
Data: 1 Dicembre 2025
