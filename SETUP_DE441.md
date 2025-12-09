# Setup JPL DE441 con CSPICE

## Prerequisiti

Hai già il file `de441.bsp`. Serve solo installare CSPICE.

## 1. Installazione CSPICE

### Download
```bash
cd ~/Downloads
wget https://naif.jpl.nasa.gov/pub/naif/toolkit/C/MacIntel_OSX_AppleC_64bit/packages/cspice.tar.Z
```

### Estrazione e Compilazione
```bash
tar xzf cspice.tar.Z
cd cspice
./makeall.csh
```

### Installazione (opzionale, in cartella comune)
```bash
# Crea directory comune per SPICE
sudo mkdir -p /usr/local/spice
sudo cp -r include /usr/local/spice/
sudo cp -r lib /usr/local/spice/
```

## 2. Posizionamento DE441.bsp

### Opzione A: Cartella comune (raccomandato)
```bash
# Crea directory dati SPICE
mkdir -p ~/data/spice
# Sposta de441.bsp (se non è già lì)
mv /path/to/de441.bsp ~/data/spice/

# Imposta variabile ambiente (aggiungi a ~/.zshrc)
echo 'export SPICE_DATA=~/data/spice' >> ~/.zshrc
source ~/.zshrc
```

### Opzione B: Nel progetto
```bash
mkdir -p astdyn/data
cp /path/to/de441.bsp astdyn/data/
```

## 3. Compilazione con DE441

### CMakeLists.txt
Aggiungi al tuo `CMakeLists.txt`:

```cmake
# Find CSPICE
find_library(CSPICE_LIB cspice 
    PATHS /usr/local/spice/lib ~/cspice/lib
    REQUIRED)

find_path(CSPICE_INCLUDE SpiceUsr.h
    PATHS /usr/local/spice/include ~/cspice/include
    REQUIRED)

# Link con astdyn
target_include_directories(astdyn PUBLIC ${CSPICE_INCLUDE})
target_link_libraries(astdyn ${CSPICE_LIB})
```

### Compilazione Manuale
```bash
g++ -std=c++17 -O2 \
    -I./astdyn/include \
    -I/usr/local/spice/include \
    -I/opt/homebrew/include/eigen3 \
    test_de441.cpp \
    astdyn/src/ephemeris/DE441Provider.cpp \
    -L/usr/local/spice/lib \
    -lcspice \
    -o test_de441
```

## 4. Uso

### Codice Esempio
```cpp
#include "astdyn/ephemeris/DE441Provider.hpp"
#include "astdyn/ephemeris/VSOP87Provider.hpp"

using namespace astdyn::ephemeris;

int main() {
    // Opzione 1: VSOP87 (built-in, sempre disponibile)
    VSOP87Provider vsop87;
    auto pos_vsop = vsop87.getPosition(CelestialBody::EARTH, 2451545.0);
    
    // Opzione 2: DE441 (massima precisione)
    std::string bsp_path = std::getenv("SPICE_DATA") 
        ? std::string(std::getenv("SPICE_DATA")) + "/de441.bsp"
        : "astdyn/data/de441.bsp";
    
    DE441Provider de441(bsp_path);
    auto pos_de441 = de441.getPosition(CelestialBody::EARTH, 2451545.0);
    
    // Confronto
    double diff_km = (pos_de441 - pos_vsop).norm() * 149597870.691;
    std::cout << "Differenza VSOP87 vs DE441: " << diff_km << " km\n";
    
    return 0;
}
```

## 5. Test Rapido

Crea `test_de441_quick.cpp`:

```cpp
#include <iostream>
#include "astdyn/ephemeris/DE441Provider.hpp"

int main() {
    try {
        // Prova a caricare DE441
        astdyn::ephemeris::DE441Provider de441(
            std::getenv("SPICE_DATA") 
                ? std::string(std::getenv("SPICE_DATA")) + "/de441.bsp"
                : "de441.bsp"
        );
        
        // Test: posizione Terra a J2000.0
        auto pos = de441.getPosition(
            astdyn::ephemeris::CelestialBody::EARTH, 
            2451545.0
        );
        
        std::cout << "✓ DE441 funziona!\n";
        std::cout << "Terra a J2000.0: " 
                  << pos.transpose() << " AU\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "✗ Errore: " << e.what() << "\n";
        return 1;
    }
}
```

Compila e testa:
```bash
g++ -std=c++17 -O2 \
    -I./astdyn/include \
    -I/usr/local/spice/include \
    -I/opt/homebrew/include/eigen3 \
    test_de441_quick.cpp \
    astdyn/src/ephemeris/DE441Provider.cpp \
    -L/usr/local/spice/lib \
    -lcspice \
    -o test_de441_quick

./test_de441_quick
```

## 6. Troubleshooting

### Errore: "Cannot find cspice"
```bash
# Verifica installazione
ls -la /usr/local/spice/lib/libcspice.a
# O
ls -la ~/cspice/lib/libcspice.a
```

### Errore: "Cannot load de441.bsp"
```bash
# Verifica file
ls -lh ~/data/spice/de441.bsp
# Dovrebbe essere ~3.3 GB

# Verifica variabile ambiente
echo $SPICE_DATA
```

### Errore di linking
```bash
# Aggiungi percorso library
export DYLD_LIBRARY_PATH=/usr/local/spice/lib:$DYLD_LIBRARY_PATH
```

## 7. Performance

**VSOP87:**
- Velocità: ⭐⭐⭐⭐⭐ (istantaneo)
- Accuratezza: ~20 arcsec
- Dimensione: 0 (built-in)

**DE441:**
- Velocità: ⭐⭐⭐⭐ (interpolazione da file)
- Accuratezza: ~0.001 arcsec (cm-level)
- Dimensione: 3.3 GB

**Raccomandazione:**
- Usa VSOP87 per default
- Usa DE441 solo se serve precisione ultra-alta

## 8. Prossimi Passi

Una volta installato CSPICE:

1. ✅ Compila test_de441_quick
2. ✅ Verifica funzionamento
3. ✅ Integra in AstDyn
4. ✅ Crea test di validazione vs VSOP87

---

**Note:**
- DE441 copre 1550-2650
- Per date fuori range, usa VSOP87
- CSPICE è thread-safe dopo inizializzazione
