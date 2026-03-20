# File di Configurazione .OOP - Guida Rapida

## Cos'è il File .OOP?

Il file `.oop` (**O**rbFit **O**ptions **P**arser) è un file di configurazione testuale che permette di controllare il comportamento di AstDyn senza modificare il codice sorgente.

## Formato del File

Il formato è semplice e leggibile:

```ini
! Commenti iniziano con punto esclamativo
chiave = valore

! Sezioni logiche per organizzare le opzioni
[integrator]
type = RKF78
tolerance = 1.0e-12

[propagator]
max_step = 10.0
```

## Struttura Completa

### 1. Configurazione Integratore

```ini
[integrator]
type = AAS                       ! RK4, RKF78, AAS, DOPRI, RADAU
step_size = 0.1                 ! Passo iniziale [giorni]
tolerance = 1.0e-12             ! Tolleranza errore (adattivo)
min_step = 1.0e-6              ! Passo minimo [giorni]
max_step = 100.0               ! Passo massimo [giorni]
```

**Integratori Disponibili:**
- `RK4` - Runge-Kutta 4° ordine (passo fisso)
- `RKF78` - Runge-Kutta-Fehlberg 7(8) (adattivo, alta precisione)
- `AAS` - AstDyn-Adaptive Symplectic (veloce, conserva energia, default)
- `DOPRI` - Dormand-Prince (adattivo)
- `RADAU` - Radau IIA (implicito, problemi stiff)

### 2. Perturbazioni

```ini
[perturbations]
planets = true                  ! Perturbazioni planetarie
planet_list = SUN,JUPITER,SATURN,EARTH,VENUS,MARS,URANUS,NEPTUNE
asteroids = true                ! Perturbazioni asteroidali (Default: true)
asteroid_set = 17               ! 17 (Default) o 30
relativity = true               ! Correzioni relativistiche (Default: true)
sun_j2 = true                   ! Armoniche solari J2
earth_j2 = true                 ! Armoniche terrestri J2
solar_radiation = false         ! Pressione radiazione solare
```

### 3. Determinazione Orbitale

```ini
[orbit_determination]
method = GAUSS                  ! GAUSS, LAPLACE, HERGET
diffcorr_max_iter = 20         ! Iterazioni max correzione differenziale
diffcorr_tolerance = 1.0e-6    ! Convergenza [AU]
outlier_threshold = 3.0        ! Soglia sigma per outlier
outlier_rejection = true       ! Abilita rimozione automatica
weight_model = GAUSSIAN        ! GAUSSIAN, UNIFORM, HUBER
```

### 4. Effemeridi

```ini
[ephemeris]
source = SPICE                 ! SPICE, JPL, ANALYTIC
kernel_path = /path/to/de440.bsp
preload = true                 ! Precarica kernel in memoria
```

### 5. Osservazioni

```ini
[observations]
format = MPC                   ! MPC, ADES, RWO, ORBFIT
sigma_default = 0.5            ! Incertezza default [arcsec]
catalog_bias_correction = true ! Correggi bias catalogo
time_bias_correction = false   ! Correggi bias temporale
```

### 6. Output

```ini
[output]
verbose = true                 ! Output verboso
log_file = astdyn.log         ! File di log
residuals_file = residuals.dat
covariance_matrix = true      ! Salva matrice covarianza
ephemeris_format = CSV        ! CSV, JSON, ORBFIT
```

### 7. Numeriche e Prestazioni

```ini
[numerical]
convergence_check = RELATIVE  ! RELATIVE, ABSOLUTE, BOTH
matrix_solver = QR            ! QR, SVD, CHOLESKY
condition_number_warn = 1e10  ! Avviso matrice mal condizionata

[performance]
parallel_residuals = true     ! Parallelizza calcolo residui (OpenMP)
batch_size = 100              ! Dimensione batch per elaborazione
cache_ephemeris = true        ! Cache posizioni planetarie
```

## Esempio Completo: Pompeja

```ini
! ============================================================
! Configurazione AstDyn per asteroide (203) Pompeja
! ============================================================

! --- INTEGRATORE ---
[integrator]
type = RKF78
tolerance = 1.0e-12
step_size = 0.1
min_step = 1.0e-6
max_step = 10.0

! --- FORZE ---
[perturbations]
planets = true
planet_list = SUN,JUPITER,SATURN
asteroids = false
relativity = false
solar_radiation = false

! --- DETERMINAZIONE ORBITALE ---
[orbit_determination]
method = GAUSS
diffcorr_max_iter = 20
diffcorr_tolerance = 1.0e-6
outlier_threshold = 3.0
outlier_rejection = true
weight_model = GAUSSIAN

! --- EFFEMERIDI ---
[ephemeris]
source = SPICE
kernel_path = data/de440.bsp
preload = true

! --- OSSERVAZIONI ---
[observations]
format = MPC
sigma_default = 0.5
catalog_bias_correction = true

! --- OUTPUT ---
[output]
verbose = true
log_file = pompeja.log
residuals_file = pompeja_residuals.dat
covariance_matrix = true
ephemeris_format = CSV
```

## Uso in C++

### Caricamento Configurazione

```cpp
#include <astdyn/io/AstDynConfig.hpp>

// Carica da file
auto config = astdyn::config::AstDynConfigManager::load_from_file("config.oop");

// Accedi ai valori
std::string integrator = config.get<std::string>("integrator.type");
double tolerance = config.get<double>("integrator.tolerance");
bool use_planets = config.get<bool>("perturbations.planets");

// Con valori di default
int max_iter = config.get<int>("orbit_determination.diffcorr_max_iter", 10);
```

### Uso con AstDynEngine

```cpp
#include <astdyn/AstDynEngine.hpp>

// Crea engine con configurazione
astdyn::AstDynEngine engine("pompeja.oop");

// Oppure carica dopo
astdyn::AstDynEngine engine;
engine.load_config("pompeja.oop");

// La configurazione è applicata automaticamente
auto result = engine.fit_orbit(observations);
```

### Configurazione Programmatica

```cpp
// Crea configurazione nel codice
astdyn::config::ConfigBuilder builder;

builder.set("integrator.type", "RKF78")
       .set("integrator.tolerance", 1e-12)
       .set("perturbations.planets", true)
       .set("outlier_threshold", 3.0);

auto config = builder.build();

// Salva su file
config.save_to_file("generated.oop");
```

## Validazione

Il parser valida automaticamente:
- **Tipi**: Verifica che i valori siano del tipo corretto
- **Range**: Controlla che i valori siano in intervalli validi
- **Dipendenze**: Assicura che opzioni dipendenti siano coerenti
- **Chiavi sconosciute**: Avvisa se una chiave non è riconosciuta

Esempio errori:

```
ERROR: Invalid value for 'integrator.tolerance': must be > 0
WARNING: Unknown key 'integrator.typoo' (did you mean 'integrator.type'?)
ERROR: 'ephemeris.kernel_path' is required when 'ephemeris.source = SPICE'
```

## Gerarchie e Override

Le configurazioni possono essere sovrascritte:

```cpp
// 1. Configurazione default (built-in)
auto config = astdyn::config::default_config();

// 2. Override con file globale
config.load_file("/etc/astdyn/global.oop");

// 3. Override con file utente
config.load_file("~/.astdyn/user.oop");

// 4. Override con file progetto
config.load_file("pompeja.oop");

// 5. Override programmatico
config.set("verbose", true);
```

Priorità: programmatico > progetto > utente > globale > default

## Best Practices

### 1. Commenta le Tue Configurazioni

```ini
! Tolleranza molto stretta per orbita ad alta precisione
! Necessaria per propagazioni > 1 anno
tolerance = 1.0e-14  ! [AU] - aumenta tempo calcolo ~2x
```

### 2. Usa Sezioni Logiche

```ini
[integrator]
# ...

[perturbations]
# ...

[output]
# ...
```

### 3. Documenta Scelte Non Ovvie

```ini
! Disabilitiamo Urano e Nettuno per Pompeja (main belt)
! Impatto < 1 metro su 100 anni, risparmio 15% tempo calcolo
planet_list = SUN,JUPITER,SATURN,EARTH,VENUS,MARS
```

### 4. Versionamento

```ini
! ==========================================================
! AstDyn Configuration v1.2
! Object: (203) Pompeja
! Created: 2025-11-26
! Author: Michele Bigi
! Purpose: High-precision orbit determination
! ==========================================================
config_version = 1.2
```

### 5. Profili Multipli

Crea file diversi per usi diversi:

- `pompeja_quick.oop` - Tolleranze larghe, calcolo veloce
- `pompeja_precise.oop` - Tolleranze strette, massima precisione
- `pompeja_production.oop` - Bilanciamento ottimale

## Formati Compatibili

AstDyn supporta anche:

### File .OEF (Orbital Element File)

```
203 Pompeja
KEP MEAN 60000.0 TDB ECLM J2000
2.743617  0.062415  11.7402  339.8591  258.0345  45.3215
! a[AU]  e  i[deg]  Omega[deg]  omega[deg]  M[deg]
```

### File .RWO (Residuals & Weights)

Estensione formato MPC con pesi e residui:

```
     203 C2022 01 15.12345 12 34 56.789 +12 34 56.78  17.2    F51  1.0 1.0  -0.12  +0.08
! Standard MPC + weight_RA weight_Dec residual_RA residual_Dec
```

## Riferimenti

- **Documentazione completa**: `docs/manual/it/18_parser.tex`
- **Esempi**: `examples/example_config.cpp`
- **Test**: `tests/test_config_parser.cpp`
- **Header**: `include/astdyn/io/AstDynConfig.hpp`
- **Implementazione**: `src/io/AstDynConfig.cpp`

## Quick Reference

| Categoria | Chiave | Valori | Default |
|-----------|--------|--------|---------|
| Integratore | `integrator.type` | RK4, RKF78, DOPRI | RKF78 |
| | `integrator.tolerance` | double > 0 | 1e-12 |
| Perturbazioni | `perturbations.planets` | bool | true |
| | `planet_list` | lista corpi | SUN,JUPITER,SATURN |
| OD | `diffcorr_max_iter` | int > 0 | 20 |
| | `diffcorr_tolerance` | double > 0 | 1e-6 |
| | `outlier_threshold` | double > 0 | 3.0 |
| Output | `verbose` | bool | false |
| | `log_file` | string | "" |

---

**Versione**: 1.0  
**Data**: 26 Novembre 2025  
**Progetto**: AstDyn - Modern C++ Orbit Determination Library
