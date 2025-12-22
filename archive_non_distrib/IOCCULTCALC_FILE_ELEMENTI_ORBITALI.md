# File di Elementi Orbitali in IOoccultCalc vs AstDynPropagator

## 📁 Localizzazione File

### AstDynPropagator (ITALOccultLibrary) ✅

**File elemento per asteroide 17030:**
```
/Users/michelebigi/Astro/OrbFit/tests/orbfit/17030/epoch/17030.eq1
```

**Formato**: OrbFit OEF2.0 (Equinoziale eclittico)

**Struttura file**:
```
! Asteroid 17030 Sierks
! OEF 2.0
MJD  61000.000000 TDT
EQU   3.175473   -0.018963   -0.041273   0.025407   -0.001956   229.790880
MAG   13.290000   0.130000
END
```

**Campi**:
- `a` = 3.175473 AU (semiasse)
- `h = e·sin(ϖ)` = -0.018963
- `k = e·cos(ϖ)` = -0.041273
- `p = tan(i/2)·sin(Ω)` = 0.025407
- `q = tan(i/2)·cos(Ω)` = -0.001956
- `λ = M + ϖ` = 229.790880° (longitudine media)
- `H` = 13.29 (magnitudine assoluta)
- `G` = 0.13 (parametro di pendenza)

**Frame**: ECLM J2000 (eclittica media J2000)

---

### IOoccultCalc 🔍

IOoccultCalc **non usa direttamente file .eq1**. Invece:

1. **Carica elementi da API AstDys online** via `astdys_client.cpp`
2. **Supporta OrbFit** tramite wrapper Fortran (`orbfit_wrapper.f90`)
3. **Accetta elementi da MPC** (Minor Planet Center)
4. **Genera elementi internamente** per calcolo

**Fonte principali in IOoccultCalc**:

```cpp
// CMakeLists.txt (riga 123-130)
// Wrapper Fortran per OrbFit:
set(FORTRAN_SOURCES
    src/orbfit_c_wrapper.f90
)

// Se abilitato con: cmake -DUSE_ORBFIT=ON
if(USE_ORBFIT)
    set(ORBFIT_PATH "/Users/michelebigi/Astro/OrbFit" ...)
    # Link a librerie OrbFit
    set(ORBFIT_LIBS
        ${ORBFIT_PATH}/src/lib/libgauss.a      # Elemento iniziale
        ${ORBFIT_PATH}/src/lib/libprop.a       # Propagazione
        ${ORBFIT_PATH}/src/lib/libsuit.a       # Utility
        ${ORBFIT_PATH}/src/lib/libmoid.a       # MOID
    )
endif()
```

---

## 🔀 Differenze Principali

| Aspetto | AstDynPropagator | IOoccultCalc |
|---------|------------------|--------------|
| **File input** | `.eq1` diretto (OrbFit OEF2.0) | API AstDys + OrbFit wrapper |
| **Formato elementi** | Equinoziale (eclittico) | Kepleriano (da OrbFit) |
| **Frame** | ECLM J2000 | Dipende da sorgente |
| **Caricamento** | File statico | Online + cache |
| **Propagatore** | RKF78 C++ puro | RKF78 + OrbFit Fortran |
| **Perturbazioni** | 8 pianeti + Schwarzschild | OrbFit (dipende da config) |

---

## 📄 File di Test che Usano Elementi

### Nel workspace ITALOccultLibrary:

1. **`test_asteroid_17030_occultation.cpp`** (LEGGI QUESTO!)
   - Riga 387: Carica `/Users/michelebigi/Astro/OrbFit/tests/orbfit/17030/epoch/17030.eq1`
   - Riga 79-120: Funzione `read_orbital_elements()` che parsifica il file
   - **Formato**: Legge elementi equinoziali direttamente da .eq1

2. **`test_astdyn_17030_nov26_30.cpp`** 
   - Riga 50-70: `load_17030_elements()` - Carica da file .eq1
   - Mostra come parsificare OEF2.0

3. **`astdyn/tools/astdyn_propagator.cpp`**
   - Hardcoded nell'esempio (riga 952)
   - Non legge da file, usa costanti

### Fuori workspace (IOoccultCalc):

4. **`test_astdyn_vs_ioccultcalc_vs_jpl.cpp`** (IN IOoccultCalc)
   - Riga 37: Riferimento a `OCCULTATION_17030_RESULTS.md` (output AstDyn)
   - Confronta IOoccultCalc vs AstDynPropagator vs JPL
   - **Tipo di file elementi**: Probabilmente da API AstDys

---

## 🎯 Come IOoccultCalc Carica Elementi

Basato su CMakeLists.txt e struttura codice:

```cpp
// Pseudocodice flusso IOoccultCalc:

// Opzione 1: Da API AstDys
AstDySClient client;
OrbitalElements elem = client.getElementsFor("17030");
// Ritorna elementi Kepleriani

// Opzione 2: Da OrbFit (Fortran)
if (USE_ORBFIT) {
    OrbFitWrapper wrapper;
    OrbitalElements elem = wrapper.readFromEQ1(filename);
    // Ritorna elementi Kepleriani
}

// Opzione 3: Da MPC
MPCClient mpc;
OrbitalElements elem = mpc.getElementsFor("17030");
// Ritorna elementi MPC (Kepleriani)
```

---

## 📊 Struttura File .eq1 (OEF2.0) - Dettagli

```
Linea 1: ! Asteroid description (commento)
Linea 2: ! Formato OEF 2.0
Linea 3: MJD  61000.000000 TDT    ← Epoca MJD (Modified Julian Date)
         └─ 61000.0 = 2018-March-16.0

Linea 4: EQU   3.175473   -0.018963   -0.041273   0.025407   -0.001956   229.790880
         └─ Tag EQU (Equinoctial)
         └─ a, h, k, p, q, λ (6 elementi equinoziali)

Linea 5: MAG   13.290000   0.130000   ← Magnitudine assoluta H, Slope G

Linea 6: END   ← Fine sezione elementi
```

### Cosa significano i numeri (per 17030):

```
a        = 3.175473 AU        → Distanza media dal Sole
h        = -0.018963          → e × sin(ϖ) dove ϖ = Ω + ω
k        = -0.041273          → e × cos(ϖ)
p        = 0.025407           → tan(i/2) × sin(Ω)
q        = -0.001956          → tan(i/2) × cos(Ω)
λ        = 229.790880°        → Longitudine media (Ω + ω + M)
H        = 13.29 mag          → Magnitudine assoluta
G        = 0.13               → Parametro di pendenza (Bowell)

Da questi 6 elementi deriva:
e = √(h² + k²) = 0.045407         (eccentricità)
i = 2·atan(√(p²+q²)) = 2.9046°    (inclinazione)
Ω = atan2(p, q) = 94.06°          (nodo ascendente)
ϖ = atan2(h, k) = 204.34°         (longitudine perielio)
ω = ϖ - Ω = 110.28°               (argomento perielio)
M = λ - ϖ = 25.45°                (anomalia media)
```

---

## 🔗 File Correlati

### AstDynPropagator legge:
```
17030.eq1 (file statico OrbFit)
    ↓
Converter (Fase 1): Equinoziali → Kepleriani
    ↓
Cartesian converter (Fase 2): Kepleriani → Stato eclittico
    ↓
Rotator (Fase 3): Eclittico → ICRF
    ↓
RKF78 propagator (Fase 4): Propaga stato
    ↓
Coordinate extractor (Fase 5): RA, Dec, distanza
```

### IOoccultCalc legge:
```
API AstDys / OrbFit / MPC
    ↓
(Elementi in formato sorgente)
    ↓
Normalizzazione a Kepleriani (possibile conversione)
    ↓
Propagatore interno (RKF78 o OrbFit Fortran)
    ↓
Output coordinate
```

---

## ⚠️ Problemi da Risolvere in IOoccultCalc

Basato su CMakeLists.txt (righe 50-60):

```cmake
# File sorgente:
src/orbital_elements.cpp
src/orbit_propagator.cpp
src/rkf78_integrator.cpp

# Ma NON tutti sono abilitati/completati:
# src/orbfit_force_model.cpp  ← COMMENTATO
# src/gaia_catalog_preloader.cpp ← COMMENTATO

# Header files mancanti:
# include/ioccultcalc/ra15_integrator.hpp  # File mancante
# include/ioccultcalc/gaia_catalog_preloader.h  # File mancante
```

**Conclusione**: IOoccultCalc ha struttura simile ma **NON COMPLETA** per leggere `.eq1` direttamente.

---

## ✅ Raccomandazione per IOoccultCalc

Per usare lo stesso flusso di AstDynPropagator:

1. **Copia il parser** da `test_asteroid_17030_occultation.cpp` (riga 79-120)
   ```cpp
   OrbitalElements read_orbital_elements(const std::string& filename);
   ```

2. **Integra il flusso conversione**:
   - Equinoziali → Kepleriani
   - Kepleriani → Stato cartesiano
   - Stato cartesiano → ICRF

3. **Usa RKF78** invece di OrbFit per coerenza

4. **Carica file .eq1** direttamente:
   ```cpp
   std::string elem_file = "/Users/michelebigi/Astro/OrbFit/tests/orbfit/17030/epoch/17030.eq1";
   OrbitalElements elements = read_orbital_elements(elem_file);
   ```

---

## 📚 Riferimenti nei File

**Dove trovare il parser per .eq1**:
- `ITALOccultLibrary/astdyn/tests/test_asteroid_17030_occultation.cpp` **Riga 79-120**

**Dove trovare le conversioni**:
- `ITALOccultLibrary/astdyn/tools/astdyn_propagator.cpp` **Riga 511-630**

**Dove trovare il test completo**:
- `ITALOccultLibrary/ALGORITMO_TEST_OCCULTAZIONE.md` **Tutta la sezione FASE 1-6**
