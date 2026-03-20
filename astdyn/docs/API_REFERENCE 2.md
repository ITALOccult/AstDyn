# AstDyn - API Guide & Reference

Benvenuto nella Guida completa all'uso delle API della libreria **AstDyn**.
Questa libreria è specializzata nei calcoli astrodinamici, validazione type-safe, time-mapping preciso, trasformazioni di frame geocentrati/eliocentrici e I/O standardizzato (MPC, AstDys, JPL Horizons).

Tutto il codice è racchiuso nel namespace base `astdyn::`.

---

## 1. Tipi Base e Sicurezza Temporale (`astdyn::types` & `astdyn::utils`)

Al cuore di AstDyn c'è una fortissima infrastruttura di "domain-driven type safety", la cui priorità è rendere impossibili a *compile-time* gli errori fisici e geometrici (es. confondere ITRF con GCRF o TT con UTC).

### Time and Epochs
- `time::EpochTDB`: Standard epoch for dynamics (Barycentric Dynamical Time).
- `time::EpochUTC`: Standard epoch for observations (Coordinated Universal Time).
- `time::EpochTT`: Terrestrial Time.
- `.mjd()`: Returns the Modified Julian Date as a double.
- `.jd()`: Returns the Julian Date as a double.

```cpp
#include <astdyn/time/epoch.hpp>

using astdyn::time::EpochTDB;
using astdyn::time::EpochUTC;

// Costruzione
auto utc_time = EpochUTC::from_mjd(60310.0);
auto tdb_time = EpochTDB::from_mjd(60310.0); // MJD TDB diretto

// Conversione Sicura
auto to_tt = utc_time.to_tt();
double tdb_mjd = utc_time.to_tdb().mjd(); // Restituisce TDB come double nativo
```

### `astdyn::types::Vector3<Frame, Unit>`
Vettori tridimensionali per posizioni e velocità, dotati di tag semantico e unità di misura rigida.

```cpp
#include <astdyn/types/vectors.hpp>
#include <astdyn/core/frame_tags.hpp>

using namespace astdyn::types;
using namespace astdyn::core;

// x, y, z in metri (GCRF - Inerziale Terrestre)
Vector3<GCRF, Meter> position(6.371e6, 0.0, 0.0);

// I calcoli matematici conservano l'Integrità (Impedisce somme frame sbagiati)
Vector3<GCRF, Meter> r2(1.0, 1.0, 1.0);
auto res = position + r2; // OK: Entrambi GCRF Meter

// Vector3<ITRF, Meter> bad = position; // ERRORE in compilazione! Frame diversi
```

---

## 2. Coordinate e Trasformazioni (`astdyn::coordinates`)

AstDyn mette a disposizione trasformazioni rigorose attraverso matrici di rotazione del frame, pilotate dai tipi `time::Epoch` (dato l'angolo siderale richiesto).

```cpp
#include <astdyn/coordinates/frame_transforms.hpp>

// Vettore statico generato da una stazione a Terra al tempo `t`
Vector3<ITRF, Meter> ground_station_pos(1.0, 2.0, 3.0); 

// Eseguiamo la "Rotazione Frame"
// Richiede Vector3 esplicito e il tempo base
Vector3<GCRF, Meter> inertial_pos = astdyn::coordinates::itrf_to_gcrf(ground_station_pos, utc_time);
```

### Conversioni di Stato (`keplerian_to_cartesian`)

Il core dinamico contiene free-functions capaci di viaggiare tra la geometria pura cartesiana e i classici parametri orbitali.

```cpp
#include <astdyn/coordinates/state_conversions.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>

// Propagare da Kepleriano a Cartesiano O viceversa
auto cartesian_state = astdyn::coordinates::keplerian_to_cartesian(kep_elements);

// Da Cartesiano (posizione/velocità in AU-Day) ad AstDyn Elements
KeplerianElements extracted = astdyn::coordinates::cartesian_to_keplerian(cartesian_state);
```

---

## 3. Gestione dell'Input/Output (`astdyn::io`)

AstDyn supporta molteplici database fisici e input di stream per leggere elementi orbitali o osservazioni (rwo, matrici mpc, file jpl).
Tutti i parser concreti espongono le opzioni:
1. `parse(filepath)` (Legge iterativamente dal FS)
2. `parse_stream(stream)` (Legge direttamente la sorgente C++ std::istream)

### Elementi Orbitali (`IOrbitParser` / `EQ1Parser` / `JPLHorizonsParser`)

Estrae `OrbitalElements` generici agnostici (MJD in TDB, Angoli e Distanze in Radianti ed AU).

```cpp
#include <astdyn/io/IParser.hpp>
#include <astdyn/io/parsers/JPLHorizonsParser.hpp>
#include <astdyn/io/EQ1Parser.hpp>

// Utilizzando i Factory in auto-detection file (rileva tipo)
auto parser = astdyn::io::ParserFactory::create_orbit_parser("pompeja.eq1");
auto kep_elements = parser->parse("pompeja.eq1"); 

// Oppure costruendo in RAM esplicitamente:
std::istringstream stream_txt("$$SOE\n2460642.5 = ... EC=0.04 ..");

astdyn::io::parsers::JPLHorizonsParser horizons;
IOrbitParser::OrbitalElements from_net = horizons.parse_stream(stream_txt);
```

### Lettura Osservazioni Ottiche (`IObservationParser` / `MPCReader` / `RWOReader`)
Decodifica la logica angolare delle stazioni. Esiste una forma raggruppata per facilitare il mapping object/osservazioni (`readFileGrouped`).

```cpp
#include <astdyn/observations/MPCReader.hpp>
#include <astdyn/observations/RWOReader.hpp>

using namespace astdyn::observations;

// MPC Classico (80 colonne) da stringa in rete
std::istringstream raw_data("00433         C2023 01 15.41667 10 34 23.45 +19 40 25.8          17.5 V      568\n");
auto  opt_observations = MPCReader::readStream(raw_data);

// Loop metriche (RA/Dec pre-calcolate in radianti)
for(const OpticalObservation& obs : opt_observations) {
    std::cout << "Osservatorio: " << obs.obs_code << "\n";
    std::cout << "Data MJD UTC: " << obs.mjd_utc << "\n";
    std::cout << "Sigma di incertezza: " << obs.sigma_ra << "\n";
}

// Analisi di residui RWO OrbFit (File System)
auto orbfit_residuals = RWOReader::readFile("/percorso/residui_ceres.rwo");
```

---

## 4. Propagazione Dinamica (`astdyn::propagation` ed `astdyn::AstDynEngine`)

Il cuore astronomico è delegabile attraverso un motore unitario ad altissimo tenore computazionale. `AstDynEngine` incapsula e nasconde RKF78, perturbazioni e gestione della pipeline orbitale.

```cpp
#include <astdyn/AstDynEngine.hpp>

using astdyn::AstDynEngine;
using astdyn::propagation::KeplerianElements;

// 1. Inizializzare Motore
AstDynEngine engine;
engine.set_verbose(true);

// 2. Elementi di Partenza
KeplerianElements initial_state;
initial_state.epoch = time::EpochTDB::from_mjd(60000.0);
initial_state.semi_major_axis = 2.76;
// ... (riempi attributi fisici) ...
engine.set_initial_orbit(initial_state);

// 3. Propagazione al Punto B (EpochTDB Forward/Backward)
time::EpochTDB target_epoch = time::EpochTDB::from_mjd(60642.0); 
KeplerianElements final_state = engine.propagate_to(target_epoch);

// 4. Oppure: Generazione Ephemeris in batch (Per grafico/plot lunghi)
// E.g Da epoca 60000 a 60300 a passi di 5 giorni
auto ephemeris = engine.compute_ephemeris(time::EpochTDB::from_mjd(60000.0), 
                                         time::EpochTDB::from_mjd(60300.0), 5.0);
std::cout << "Trovati " << ephemeris.size() << " campioni\n";
```

### Configurazione del Propagatore (Avanzato)
Il motore nasconde un factory dedicato a produrre le istanze RKF78 in alta affidabilità. Se ti serve iniettare logica extra su perturbazioni, contatta l'estrazione propagatori:

```cpp
#include <astdyn/propagation/RKF78Integrator.hpp>

// Instanzio direttamente senza proxy.
astdyn::propagation::RKF78Integrator numerical_integrator;
numerical_integrator.set_tolerance(1e-12); // override rigoroso
```

---

## 5. Integratori Numerici Sottomodulo (`astdyn::propagation::Integrator`)

AstDyn garantisce il controllo *low-level* fornendo una flotta di integratori specializzati per problemi fisici o stiff differenti, tutti estendenti la classe base `Integrator`.
Tutti supportano l'intercettamento delle metriche prestazionali attraverso `integrator.statistics()`.

### Integratori Disponibili
*   **`RK4Integrator`**: Classico Runge-Kutta a passo fisso. O(h^4). Estremamente veloce ma non adattivo. Ideale per step geometrici rigidi minimi.
*   **`RKF78Integrator`**: *(Il raccomandato e usato da `AstDynEngine`)* Runge-Kutta-Fehlberg 7(8). Adattivo, eccellente compromesso tra performance e tolleranza ($1e^{-12}$) per l'Astrodinamica standard su scale TDB di decenni.
*   **`GaussIntegrator`**: Integratore implicito di Gauss-Legendre (Simpletico, Ordine 8). Strutturato per salvaguardare la conservazione dell'energia hamiltoniana nel sistema. Obbligatorio per long-term runs (millenni, zero secular drift dell'energia orbitale).
*   **`RadauIntegrator`**: Radau IIA implicito (Ordine 15, A-stable). Sviluppato per condizioni *stiff* estreme che manderebbero in timeout l'RKF78 e per precisioni limite > $1e^{-13}$. Newton-solver based implicito.

### Esempio di utilizzo Low-Level
Invece di usare `AstDynEngine`, il backend ODE permette:

```cpp
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/GaussIntegrator.hpp>

// Sistema Conservativo, usiamo Gauss 8 implicito per zero energy drift
astdyn::propagation::GaussIntegrator gauss(
    0.1,    // step iniziale
    1e-13,  // tolleranza rigida
    1e-8,   // step sub-milli minimo
    10.0    // step max constraint
);

// Funzione differenziale f(t, y) = dy/dt definita dal problema fisico
astdyn::propagation::DerivativeFunction deriv_f = [](double t, const Eigen::VectorXd& y) {
    Eigen::VectorXd dy = Eigen::VectorXd::Zero(6);
    // ... Implementazione newtoniana N-Body accelerazione ...
    return dy;
};

// Stato Cartesiano iniziale e propagazione da 0 a 365 giorni rigida
Eigen::VectorXd initial_state(6); // riempire x,y,z,vx,vy,vz
Eigen::VectorXd final_y = gauss.integrate(deriv_f, initial_state, 0.0, 365.0);

// Metriche per diagnosticare eventuale stifness problematica
auto stats = gauss.statistics();
std::cout << "Step totali O(8): " << stats.num_steps << "\n";
```

---

## Best Practices
1. **Passa sempre un tipo `Epoch`**. Tutte le nuove interfacce preferiscono le classi specializzate in `time` per impedire ambiguità su `TDB`, `UTC`, o `TT`.
2. **Usa `Vector3` espliciti** ove possible per la geometria (ad. es geodetica e ground radar tracking), avviliti solo al fallback standard (array doppi) quando l'IO non offre il porting immediato.
3. Lo **Stream I/O** è fortemente preferito per la programmazione distribuita: genera sempre uno `std::istringstream` invece di buttare dati transitori in una directory cache.
4. L'Engine copre già *l'Energy Conservation* in propagazioni `J2` su run lunghi (l'errore di drift sull'RKF78 si assale in `< 0.001` AU/D over 5 anni testati), ma in framework di milioni di particelle a N-Body, preferire l'integrazione di `GaussIntegrator` bypassando il wrapper.
