# Guida: Esportare AstDynPropagator in Altri Progetti

## Metodo 1: Header-Only (Consigliato)

### Passo 1: Crea header singolo
```bash
# Copia il file principale
cp astdyn/tools/astdyn_propagator.cpp include/astdyn_propagator.hpp

# Modifica il file per renderlo header-only
```

### Passo 2: Modifica per header-only
Aggiungi all'inizio del file:
```cpp
#ifndef ASTDYN_PROPAGATOR_HPP
#define ASTDYN_PROPAGATOR_HPP

// ... tutto il codice ...

#endif // ASTDYN_PROPAGATOR_HPP
```

Commenta o rimuovi il `main()` alla fine del file:
```cpp
#ifndef ASTDYN_NO_MAIN
int main() {
    // ... esempio ...
}
#endif
```

### Passo 3: Uso nel tuo progetto
```cpp
#define ASTDYN_NO_MAIN  // Disabilita main
#include "astdyn_propagator.hpp"

using namespace astdyn;

int main() {
    // Elementi orbitali
    OrbitalElements elem;
    elem.epoch = 2458193.5;
    elem.a = 3.173489964321051;
    elem.e = 0.04796607451625862;
    elem.i = 2.904309538190326;
    elem.Omega = 104.1845838362649;
    elem.omega = 102.1497438064497;
    elem.M = 99.03517819281583;
    elem.name = "17030";
    
    // Crea propagatore
    AstDynPropagator prop(1e-12);
    prop.usePlanets(true);
    prop.useRelativity(true);
    
    // Propaga a data target
    double target_jd = 2460277.5;  // 28 Nov 2025
    EquatorialCoords coords = prop.propagateElements(elem, target_jd);
    
    std::cout << "RA:  " << coords.ra << "°\n";
    std::cout << "Dec: " << coords.dec << "°\n";
    
    return 0;
}
```

---

## Metodo 2: Libreria Statica

### Passo 1: Compila come libreria
```bash
cd astdyn/tools

# Compila oggetto
g++ -std=c++17 -O2 -c -DASTDYN_NO_MAIN \
    astdyn_propagator.cpp -o astdyn_propagator.o

# Crea libreria statica
ar rcs libastdyn.a astdyn_propagator.o

# Opzionale: installa
sudo cp libastdyn.a /usr/local/lib/
sudo cp astdyn_propagator.cpp /usr/local/include/astdyn_propagator.hpp
```

### Passo 2: Usa nel progetto
```bash
# Compila tuo progetto linkando la libreria
g++ -std=c++17 -O2 my_program.cpp -lastdyn -o my_program
```

---

## Metodo 3: Libreria Condivisa (.so/.dylib)

### Passo 1: Compila shared library
```bash
# Linux
g++ -std=c++17 -O2 -fPIC -shared -DASTDYN_NO_MAIN \
    astdyn_propagator.cpp -o libastdyn.so

# macOS
g++ -std=c++17 -O2 -fPIC -dynamiclib -DASTDYN_NO_MAIN \
    astdyn_propagator.cpp -o libastdyn.dylib
```

### Passo 2: Installa (opzionale)
```bash
# Linux
sudo cp libastdyn.so /usr/local/lib/
sudo ldconfig

# macOS
sudo cp libastdyn.dylib /usr/local/lib/
```

### Passo 3: Compila con la libreria
```bash
g++ -std=c++17 my_program.cpp -L/usr/local/lib -lastdyn -o my_program
```

---

## Metodo 4: CMake (Professionale)

### Passo 1: Crea struttura
```
my_project/
├── CMakeLists.txt
├── src/
│   └── main.cpp
└── external/
    └── astdyn/
        ├── astdyn_propagator.cpp
        └── astdyn_propagator.hpp
```

### Passo 2: CMakeLists.txt principale
```cmake
cmake_minimum_required(VERSION 3.15)
project(MyProject CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Aggiungi AstDyn come libreria
add_library(astdyn 
    external/astdyn/astdyn_propagator.cpp
)

target_compile_definitions(astdyn PRIVATE ASTDYN_NO_MAIN)
target_include_directories(astdyn PUBLIC external/astdyn)

# Tuo eseguibile
add_executable(my_program src/main.cpp)
target_link_libraries(my_program astdyn)
```

### Passo 3: Build
```bash
mkdir build && cd build
cmake ..
make
```

---

## Metodo 5: Submodule Git (Migliore per Sviluppo)

### Passo 1: Aggiungi come submodule
```bash
cd my_project
git submodule add https://github.com/manvalan/ITALOccultLibrary.git external/italoccult
git submodule update --init --recursive
```

### Passo 2: CMakeLists.txt
```cmake
# Aggiungi path
add_subdirectory(external/italoccult/astdyn)

# Oppure includi direttamente
include_directories(external/italoccult/astdyn/tools)

add_executable(my_program src/main.cpp)
target_link_libraries(my_program m)  # Link math library
```

### Passo 3: Nel codice
```cpp
#define ASTDYN_NO_MAIN
#include "astdyn_propagator.cpp"

// ... usa come normale
```

---

## API Essenziale

### Strutture Dati

```cpp
// Elementi orbitali kepleriani
struct OrbitalElements {
    std::string name;
    double epoch;   // JD
    double a;       // AU
    double e;       // eccentricità
    double i;       // inclinazione (gradi)
    double Omega;   // longitudine nodo (gradi)
    double omega;   // argomento perielio (gradi)
    double M;       // anomalia media (gradi)
};

// Elementi equinoziali (formato AstDyS)
struct EquinoctialElements {
    std::string name;
    double epoch_mjd;  // MJD
    double a;          // AU
    double h, k;       // e*sin(ϖ), e*cos(ϖ)
    double p, q;       // tan(i/2)*sin(Ω), tan(i/2)*cos(Ω)
    double lambda;     // longitudine media
};

// Coordinate equatoriali
struct EquatorialCoords {
    double ra;   // gradi
    double dec;  // gradi
    double dist; // AU
};
```

### Metodi Principali

```cpp
class AstDynPropagator {
    // Costruttore
    AstDynPropagator(double tolerance = 1e-12);
    
    // Configurazione
    void setTolerance(double tol);
    void usePlanets(bool enable);
    void useAST17(bool enable);
    void useRelativity(bool enable);
    
    // Conversioni
    OrbitalElements equinoctialToKeplerian(const EquinoctialElements& eq);
    State elementsToState(const OrbitalElements& elem);
    
    // Propagazione
    State propagate(const State& y0, double jd0, double jd1);
    EquatorialCoords propagateElements(const OrbitalElements& elem, 
                                      double target_jd);
    
    // Coordinate
    EquatorialCoords getEquatorialCoords(const State& state, double jd);
};
```

---

## Esempio Completo: Calcolo Occultazione

```cpp
#define ASTDYN_NO_MAIN
#include "astdyn_propagator.cpp"

using namespace astdyn;

int main() {
    // Elementi asteroide 17030
    OrbitalElements elem;
    elem.name = "17030";
    elem.epoch = 2458193.5;  // 2018-03-16
    elem.a = 3.173489964321051;
    elem.e = 0.04796607451625862;
    elem.i = 2.904309538190326;
    elem.Omega = 104.1845838362649;
    elem.omega = 102.1497438064497;
    elem.M = 99.03517819281583;
    
    // Setup propagatore
    AstDynPropagator prop(1e-12);
    prop.usePlanets(true);
    prop.useRelativity(true);
    
    // Coordinate stella
    double star_ra = 73.41610815;   // gradi
    double star_dec = 20.33166161;
    
    // Propaga per finestra temporale
    double jd_start = 2460277.5;  // 28 Nov 2025 00:00
    double min_dist = 1e10;
    double best_jd = jd_start;
    
    for (int min = 0; min <= 60; min += 5) {
        double jd = jd_start + min / 1440.0;
        
        EquatorialCoords ast = prop.propagateElements(elem, jd);
        
        // Distanza angolare
        double dra = (ast.ra - star_ra) * cos(star_dec * M_PI/180);
        double ddec = ast.dec - star_dec;
        double dist_deg = sqrt(dra*dra + ddec*ddec);
        double dist_arcsec = dist_deg * 3600;
        
        if (dist_deg < min_dist) {
            min_dist = dist_deg;
            best_jd = jd;
        }
        
        printf("%02d:%02d UTC - Distance: %.2f arcsec\n", 
               min/60, min%60, dist_arcsec);
    }
    
    printf("\nMinimum distance: %.2f arcsec at JD %.6f\n", 
           min_dist * 3600, best_jd);
    
    return 0;
}
```

Compila:
```bash
g++ -std=c++17 -O2 occultation.cpp -o occultation -lm
./occultation
```

---

## Dipendenze

**Richieste**:
- C++17 compiler (g++ 7+, clang 5+)
- Math library (`-lm`)

**Opzionali**:
- OpenMP per parallelizzazione (`-fopenmp`)
- Eigen3 per algebra lineare avanzata

---

## Troubleshooting

### Errore: "redefinition of main"
```cpp
#define ASTDYN_NO_MAIN
#include "astdyn_propagator.cpp"
```

### Errore: "undefined reference to sqrt"
```bash
g++ ... -lm  # Aggiungi -lm alla fine
```

### Warning: "unused variable"
```bash
g++ -Wno-unused-variable ...
```

---

## Performance

- **Compilazione**: ~2-3 secondi (O2)
- **Propagazione 8 anni**: ~0.5 secondi (1e-12 tol)
- **Memoria**: ~10 MB
- **Thread-safe**: Sì (ogni istanza indipendente)

---

## Licenza e Crediti

AstDynPropagator fa parte di **ITALOccultLibrary**
- Repository: https://github.com/manvalan/ITALOccultLibrary
- Autore: Michele Bigi
- Licenza: [specificare]

---

## 🚀 Novità: astdyn_trajectory_export (Tool ad Alta Precisione)

Per casi di ricerca scientifica che richiedono integratori impliciti o simplettici (Gauss, SABA4, IAS15), utilizzare il nuovo tool `astdyn_trajectory_export.cpp` situato in `astdyn/tools/`.

### Esempio: Esportare traiettoria ad alta fedeltà
```bash
# Esempio per Ceres usando IAS15 (Radau di ordine 15)
./bin/astdyn_trajectory_export --asteroid 1 --t0 60000 --tf 61000 --integrator IAS15 --tolerance 1e-12 --output ceres_ias15.csv
```

---

**© 2026 AstDyS Team - Dipartimento di Matematica, Università di Pisa**
*Ultima modifica: 17 marzo 2026*
