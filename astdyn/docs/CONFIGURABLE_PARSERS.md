# Sistema di Parser Configurabili - AstDyn

## Panoramica

Il sistema di parser configurabili permette di leggere dati orbitali e osservazioni da diversi formati in modo flessibile e estendibile. Invece di avere codice rigido per un solo formato, il sistema usa **interfacce** e **factory pattern** per supportare facilmente nuovi formati.

## Architettura

### 1. Interfacce Base (`IParser.hpp`)

Due interfacce principali definiscono il contratto per tutti i parser:

#### `IOrbitParser` - Parser per Elementi Orbitali
```cpp
class IOrbitParser {
public:
    struct OrbitalElements {
        std::string object_name;
        double epoch_mjd_tdb;
        double semi_major_axis;      // AU
        double eccentricity;
        double inclination;          // radianti
        double longitude_asc_node;   // Ω (radianti)
        double argument_perihelion;  // ω (radianti)
        double mean_anomaly;         // M (radianti)
        double magnitude, mag_slope; // H, G
    };
    
    virtual OrbitalElements parse(const std::string& filepath) = 0;
    virtual std::string name() const = 0;
    virtual bool can_handle(const std::string& filepath) const = 0;
};
```

#### `IObservationParser` - Parser per Osservazioni
```cpp
class IObservationParser {
public:
    struct OpticalObservation {
        std::string object_name;
        double mjd_utc;
        double ra, dec;          // radianti
        double mag;
        std::string obs_code;
        double sigma_ra, sigma_dec;  // arcsec
    };
    
    virtual std::vector<OpticalObservation> parse(
        const std::string& filepath, 
        size_t max_count = 0) = 0;
    virtual std::string name() const = 0;
    virtual bool can_handle(const std::string& filepath) const = 0;
};
```

### 2. Factory Pattern (`ParserFactory`)

La factory crea automaticamente il parser corretto in base al file:

```cpp
// Auto-detection basata sull'estensione
auto parser = ParserFactory::create_orbit_parser("data.eq1");

// Formato esplicito (se l'estensione non è standard)
auto parser = ParserFactory::create_orbit_parser("data.txt", "eq1");
```

### 3. Implementazioni Concrete

#### `OrbFitEQ1Parser` - Formato OrbFit .eq1 (OEF2.0)

Parser per file di elementi orbitali equinoziali da AstDyS:

**Formato file:**
```
! Object 203
 EQU   2.7385249934   0.0450870893   0.0412312978  -0.0059476458   0.0270423523  112.3228065416
 MJD 61000.000000 TDB
 MAG 8.958000 0.150000
```

**Caratteristiche:**
- Converte automaticamente da elementi equinoziali a Kepleriani
- Gestisce linee con spazi iniziali (formato OEF2.0)
- Supporta magnitude opzionale

#### `OrbFitRWOParser` - Formato OrbFit .rwo

Parser per osservazioni ottiche (formato AstDyS/OrbFit).

**Note:** Il formato RWO reale è molto complesso e richiede integrazione con il parser esistente `RWOReader`.

## Utilizzo

### Esempio 1: Auto-detection (Consigliato)

```cpp
#include "astdyn/io/IParser.hpp"

using namespace astdyn::io;

// Il sistema rileva automaticamente il formato da .eq1
auto orbit_parser = ParserFactory::create_orbit_parser("203.eq1");
auto elements = orbit_parser->parse("203.eq1");

std::cout << "a = " << elements.semi_major_axis << " AU\n";
std::cout << "e = " << elements.eccentricity << "\n";
```

### Esempio 2: Formato Esplicito

```cpp
// Se il file non ha estensione standard
auto parser = ParserFactory::create_orbit_parser("data.txt", "eq1");
auto elements = parser->parse("data.txt");
```

### Esempio 3: Controllo Capacità

```cpp
auto parser = ParserFactory::create_orbit_parser("203.eq1");
std::cout << "Parser: " << parser->name() << "\n";
std::cout << "Can handle .eq1: " << parser->can_handle("test.eq1") << "\n";
std::cout << "Can handle .mpc: " << parser->can_handle("test.mpc") << "\n";
```

### Esempio 4: Workflow Completo

```cpp
// 1. Parse elementi orbitali
auto orbit_parser = ParserFactory::create_orbit_parser("203.eq1");
auto elements = orbit_parser->parse("203.eq1");

// 2. Parse osservazioni
auto obs_parser = ParserFactory::create_observation_parser("203.rwo");
auto observations = obs_parser->parse("203.rwo", 100);  // max 100 obs

// 3. Usa i dati per orbit fitting
// ... (integrazione con AstDynEngine)
```

## Estendibilità

### Aggiungere un Nuovo Formato

1. **Crea il parser** che implementa `IOrbitParser` o `IObservationParser`:

```cpp
// include/astdyn/io/parsers/MPCParser.hpp
class MPCParser : public IOrbitParser {
public:
    OrbitalElements parse(const std::string& filepath) override {
        // Implementazione specifica per formato MPC
        ...
    }
    
    std::string name() const override {
        return "MPC Parser";
    }
    
    bool can_handle(const std::string& filepath) const override {
        return filepath.find(".mpc") != std::string::npos;
    }
};
```

2. **Aggiungi alla factory** in `ParserFactory.cpp`:

```cpp
if (lower_format == "mpc") {
    return std::make_unique<parsers::MPCParser>();
}
```

3. **Usa il nuovo parser** automaticamente:

```cpp
auto parser = ParserFactory::create_orbit_parser("data.mpc");
// Factory seleziona automaticamente MPCParser!
```

## Formati Supportati

### Attualmente Implementati

| Formato | Estensione | Parser | Note |
|---------|-----------|--------|------|
| OrbFit Equinoctial | `.eq1` | `OrbFitEQ1Parser` | Elementi equinoziali → Kepleriani |

### In Sviluppo

| Formato | Estensione | Status |
|---------|-----------|---------|
| OrbFit RWO | `.rwo` | Da integrare con `RWOReader` esistente |
| MPC | `.mpc` | Pianificato |
| OEL | `.oel` | Pianificato |

## Vantaggi

1. **Flessibilità**: Aggiungi nuovi formati senza modificare il codice esistente
2. **Auto-detection**: L'utente non deve specificare il formato
3. **Consistenza**: Tutti i parser restituiscono la stessa struttura dati
4. **Testabilità**: Ogni parser è isolato e testabile indipendentemente
5. **Estendibilità**: Factory pattern facilita l'aggiunta di nuovi formati

## Struttura dei File

```
astdyn/
├── include/astdyn/io/
│   ├── IParser.hpp                 # Interfacce base + Factory
│   ├── EQ1Parser.hpp               # Parser originale standalone
│   └── parsers/
│       ├── OrbFitEQ1Parser.hpp     # Implementazione configurabile .eq1
│       └── OrbFitRWOParser.hpp     # Implementazione .rwo (stub)
├── src/io/
│   └── ParserFactory.cpp           # Implementazione factory
└── examples/
    └── example_configurable_parsers.cpp  # Esempi d'uso
```

## Note Tecniche

### Conversione Equinoziale → Kepleriana

Gli elementi equinoziali (a, h, k, p, q, λ) vengono convertiti in Kepleriani (a, e, i, Ω, ω, M):

```
e = sqrt(h² + k²)
i = 2·atan(sqrt(p² + q²))
Ω = atan2(p, q)
ω̄ = atan2(h, k)
ω = ω̄ - Ω
M = λ - ω̄
```

Tutti gli angoli sono normalizzati in [0, 2π).

### Gestione Formato OEF2.0

Il formato OEF2.0 ha una particolarità: le linee di dati iniziano con uno **spazio**:

```
 EQU   2.7385...    (notare lo spazio iniziale)
 MJD 61000.0 TDB
```

I parser gestiscono correttamente questa caratteristica usando:
```cpp
if (line.find("EQU") == 0 || line.find(" EQU") == 0)
```

## Esempi Compilati

```bash
# Compila gli esempi
cd ITALOccultLibrary
cmake --build build --target example_configurable_parsers

# Esegui
./build/examples/example_configurable_parsers astdyn/tools/203_astdys_latest.eq1
```

Output:
```
Using parser: OrbFit EQ1 Parser (OEF2.0)

=== Orbital Elements ===
Object: 203
Epoch: MJD 61000.000000 TDB

Keplerian Elements:
  a = 2.7385249934 AU
  e = 0.06109718
  i = 3.172079°
  Ω = 347.595960°
  ω = 59.961709°
  M = 64.765138°

Magnitude: H = 8.958000, G = 0.150000
```

## Integrazione con AstDyn Engine

Il sistema di parser fornisce una **interfaccia standardizzata** che può essere facilmente integrata nel workflow di orbit fitting:

```cpp
// 1. Parser configurabile
auto orbit_parser = ParserFactory::create_orbit_parser(orbit_file);
auto initial_orbit = orbit_parser->parse(orbit_file);

// 2. Convertire a formato AstDyn Engine
OrbitState state = convert_to_orbit_state(initial_orbit);

// 3. Eseguire fit
AstDynEngine engine;
auto result = engine.run_differential_correction(state, observations);
```

Questo separa la **logica di parsing** dalla **logica di calcolo**, rendendo il codice più modulare e manutenibile.
