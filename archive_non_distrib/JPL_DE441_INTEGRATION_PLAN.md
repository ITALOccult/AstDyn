# Piano di Integrazione JPL DE441

## Obiettivo
Integrare le effemeridi planetarie **JPL DE441** nella libreria AstDyn per fornire la massima precisione nelle effemeridi planetarie.

## Background

### JPL Development Ephemerides

| Versione | Epoca | Accuratezza | Dimensione | Note |
|:---------|:------|:------------|:-----------|:-----|
| **DE405** | 1600-2200 | ~1 km | ~60 MB | Standard attuale |
| **DE430** | 1550-2650 | ~10 cm | ~120 MB | Alta precisione |
| **DE441** | 1550-2650 | **~cm** | ~3.3 GB | **Massima precisione** |

### Stato Attuale AstDyn

**Implementazione corrente:**
- ✅ VSOP87 analitico (built-in, nessuna dipendenza)
- ✅ Accuratezza: 1-20 arcsec (1800-2050)
- ✅ Veloce, auto-contenuto
- ❌ Limitato per propagazioni ultra-precise

**Supporto opzionale:**
- ⚠️ JPL DE via CSPICE (menzionato ma non implementato)
- ⚠️ Richiede installazione SPICE Toolkit

## Approcci Possibili

### Opzione 1: CSPICE (NASA SPICE Toolkit)

**Pro:**
- ✅ Standard NASA/JPL
- ✅ Supporta DE405, DE430, DE441
- ✅ Libreria matura e testata
- ✅ Molte funzionalità aggiuntive

**Contro:**
- ❌ Dipendenza esterna pesante (~100 MB)
- ❌ Compilazione complessa
- ❌ Licenza NASA (open ma specifica)
- ❌ API C (richiede wrapper C++)

**Implementazione:**
```cpp
#include <cspice/SpiceUsr.h>

// Carica kernel
furnsh_c("de441.bsp");

// Ottieni posizione
SpiceDouble state[6];
SpiceDouble lt;
spkez_c(399,  // Earth
        et,   // Ephemeris time
        "J2000",
        "NONE",
        10,   // Sun
        state,
        &lt);
```

---

### Opzione 2: Implementazione Custom (Lettura BSP)

**Pro:**
- ✅ Nessuna dipendenza esterna
- ✅ Controllo completo
- ✅ Più leggero
- ✅ Integrazione nativa C++

**Contro:**
- ❌ Richiede implementare parser BSP
- ❌ Interpolazione Chebyshev complessa
- ❌ Tempo di sviluppo: ~1-2 settimane

**Formato BSP:**
- File binario con coefficienti Chebyshev
- Interpolazione polinomiale per posizioni/velocità
- Struttura DAF (Double Precision Array File)

---

### Opzione 3: Ibrido (VSOP87 + DE441 opzionale)

**Pro:**
- ✅ Mantiene VSOP87 come default
- ✅ DE441 come opzione per massima precisione
- ✅ Nessuna breaking change
- ✅ Utente sceglie il trade-off

**Contro:**
- ⚠️ Complessità architetturale
- ⚠️ Due code path da mantenere

**Architettura:**
```cpp
class PlanetaryEphemeris {
public:
    enum class Provider {
        VSOP87,    // Built-in, veloce
        JPL_DE441  // Opzionale, massima precisione
    };
    
    void setProvider(Provider p);
    Vector3d getPosition(CelestialBody body, double jd_tdb);
};
```

---

## Raccomandazione: **Opzione 3 (Ibrido)**

### Implementazione Proposta

#### Fase 1: Interfaccia Astratta
```cpp
// EphemerisProvider.hpp
class EphemerisProvider {
public:
    virtual ~EphemerisProvider() = default;
    virtual Vector3d getPosition(CelestialBody body, double jd_tdb) = 0;
    virtual Vector3d getVelocity(CelestialBody body, double jd_tdb) = 0;
};

class VSOP87Provider : public EphemerisProvider {
    // Implementazione esistente
};

class DE441Provider : public EphemerisProvider {
    // Nuova implementazione
};
```

#### Fase 2: Factory Pattern
```cpp
class PlanetaryEphemeris {
private:
    std::unique_ptr<EphemerisProvider> provider_;
    
public:
    void useVSOP87() {
        provider_ = std::make_unique<VSOP87Provider>();
    }
    
    void useDE441(const std::string& bsp_file) {
        provider_ = std::make_unique<DE441Provider>(bsp_file);
    }
};
```

#### Fase 3: Uso
```cpp
// Default: VSOP87 (nessuna dipendenza)
PlanetaryEphemeris eph;
auto pos = eph.getPosition(CelestialBody::EARTH, jd);

// Opzionale: DE441 (massima precisione)
eph.useDE441("/path/to/de441.bsp");
auto pos_precise = eph.getPosition(CelestialBody::EARTH, jd);
```

---

## Implementazione DE441Provider

### Opzione A: Via CSPICE

**File necessari:**
- `de441.bsp` (~3.3 GB, scaricabile da NAIF)
- `libcspice.a` (compilato)

**Codice:**
```cpp
class DE441Provider : public EphemerisProvider {
private:
    bool spice_loaded_ = false;
    
public:
    DE441Provider(const std::string& bsp_file) {
        furnsh_c(bsp_file.c_str());
        spice_loaded_ = true;
    }
    
    Vector3d getPosition(CelestialBody body, double jd_tdb) override {
        double et = (jd_tdb - 2451545.0) * 86400.0;  // JD to ET
        SpiceDouble state[6];
        SpiceDouble lt;
        
        int naif_id = bodyToNAIFId(body);
        spkez_c(naif_id, et, "J2000", "NONE", 10, state, &lt);
        
        // Convert km to AU
        return Vector3d(state[0], state[1], state[2]) / 149597870.691;
    }
};
```

### Opzione B: Parser BSP Custom

**Vantaggi:**
- Nessuna dipendenza CSPICE
- Più leggero

**Svantaggi:**
- Complesso da implementare
- Richiede ~1 settimana di lavoro

---

## Timeline

### Fase 1: Setup (1 giorno)
- [x] Documentare piano
- [ ] Scaricare DE441.bsp (3.3 GB)
- [ ] Decidere: CSPICE vs Custom
- [ ] Setup build system

### Fase 2: Implementazione (2-3 giorni)
- [ ] Creare interfaccia `EphemerisProvider`
- [ ] Implementare `DE441Provider`
- [ ] Testare vs VSOP87
- [ ] Validare vs JPL Horizons

### Fase 3: Integrazione (1 giorno)
- [ ] Aggiornare `PlanetaryEphemeris`
- [ ] Documentazione
- [ ] Esempi d'uso

### Fase 4: Testing (1 giorno)
- [ ] Unit tests
- [ ] Performance benchmarks
- [ ] Validazione accuratezza

**Totale:** ~5-7 giorni

---

## Dipendenze

### Se usiamo CSPICE:

**Download:**
```bash
# CSPICE Toolkit
wget https://naif.jpl.nasa.gov/pub/naif/toolkit/C/MacIntel_OSX_AppleC_64bit/packages/cspice.tar.Z

# DE441 BSP
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441.bsp
```

**Build:**
```bash
# Estrai CSPICE
tar xzf cspice.tar.Z
cd cspice
# Compila
./makeall.csh
```

**Link in CMake:**
```cmake
find_library(CSPICE_LIB cspice PATHS /path/to/cspice/lib)
target_link_libraries(astdyn ${CSPICE_LIB})
```

---

## Validazione

### Test di Accuratezza

**Confronto VSOP87 vs DE441:**
```cpp
// Test: Earth position at J2000.0
double jd = 2451545.0;

// VSOP87
eph.useVSOP87();
auto pos_vsop = eph.getPosition(CelestialBody::EARTH, jd);

// DE441
eph.useDE441("de441.bsp");
auto pos_de441 = eph.getPosition(CelestialBody::EARTH, jd);

// Differenza
double diff_km = (pos_de441 - pos_vsop).norm() * 149597870.691;
std::cout << "Difference: " << diff_km << " km\n";
```

**Atteso:**
- Differenza: ~10-100 km (VSOP87 ha errore ~20 arcsec)
- DE441 dovrebbe matchare JPL Horizons a livello di cm

---

## Prossimi Passi

1. **Decidere approccio:**
   - CSPICE (più veloce, dipendenza)
   - Custom (più lavoro, nessuna dipendenza)

2. **Download file:**
   - DE441.bsp (3.3 GB)
   - CSPICE toolkit (se necessario)

3. **Implementare:**
   - Interfaccia provider
   - DE441Provider
   - Tests

Vuoi che proceda con l'implementazione? Quale approccio preferisci?
