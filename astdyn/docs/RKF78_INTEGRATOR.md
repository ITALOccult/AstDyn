
# AstDynPropagator - Documentazione Completa

## Propagatore Orbitale ad Alta Precisione con RKF78

**Autore**: AstDyS Team  
**Data**: 29 Novembre 2025  
**Versione**: 1.1

---

## 1. Panoramica

**AstDynPropagator** è un propagatore orbitale ad alta precisione basato sull'integratore **RKF78** (Runge-Kutta-Fehlberg di ordine 7 con stima dell'errore di ordine 8).

### Caratteristiche Principali

| Proprietà | Valore |
|-----------|--------|
| **Ordine di propagazione** | 7 |
| **Ordine stima errore** | 8 |
| **Numero di stadi** | 13 |
| **Tipo** | Embedded, adattivo |
| **Valutazioni per passo** | 13 |
| **Frame di riferimento** | ICRF (equatoriale J2000) |
| **Supporto formati** | Kepleriano, Equinoziale (AstDyS) |
| **Riferimento** | Fehlberg (1968) NASA TR R-287 |

### Validazione

| Sorgente Elementi | Intervallo | Errore vs JPL |
|-------------------|------------|---------------|
| JPL Horizons | 6.9 anni | ~11" |
| AstDyS OEF2.0 | ±30 giorni | ~11" |
| Round-trip | ±30 giorni | < 1 metro |

---

## 2. Formulazione Matematica

### 2.1 Schema Generale

Per un sistema $\dot{y} = f(t, y)$, il metodo calcola:

$$y_{n+1} = y_n + h \sum_{i=1}^{13} b_i k_i$$

dove gli stadi sono:

$$k_i = f\left(t_n + c_i h, \, y_n + h \sum_{j=1}^{i-1} a_{ij} k_j\right)$$

### 2.2 Stima dell'Errore

L'errore locale è stimato come differenza tra soluzioni di ordine 7 e 8:

$$\epsilon = h \sum_{i=1}^{13} (b_i^{(7)} - b_i^{(8)}) k_i$$

### 2.3 Controllo Adattivo del Passo

Il nuovo passo è calcolato come:

$$h_{new} = h \cdot \min\left(5, \max\left(0.1, 0.9 \cdot \left(\frac{\text{tol}}{|\epsilon|}\right)^{1/8}\right)\right)$$

---

## 3. Coefficienti di Fehlberg

### 3.1 Nodi $c_i$

```cpp
static constexpr double c[13] = {
    0.0,
    2.0/27.0,
    1.0/9.0,
    1.0/6.0,
    5.0/12.0,
    1.0/2.0,
    5.0/6.0,
    1.0/6.0,
    2.0/3.0,
    1.0/3.0,
    1.0,
    0.0,
    1.0
};
```

### 3.2 Pesi Ordine 7 ($b^{(7)}$)

```cpp
static constexpr double b7[13] = {
    41.0/840.0,
    0.0,
    0.0,
    0.0,
    0.0,
    34.0/105.0,
    9.0/35.0,
    9.0/35.0,
    9.0/280.0,
    9.0/280.0,
    41.0/840.0,
    0.0,
    0.0
};
```

### 3.3 Pesi Ordine 8 ($b^{(8)}$)

```cpp
static constexpr double b8[13] = {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    34.0/105.0,
    9.0/35.0,
    9.0/35.0,
    9.0/280.0,
    9.0/280.0,
    0.0,
    41.0/840.0,
    41.0/840.0
};
```

---

## 4. API dell'Integratore

### 4.1 Classe `RKF78Integrator`

```cpp
class RKF78Integrator {
public:
    // Costruttore con tolleranza (default 1e-12)
    explicit RKF78Integrator(double tol = 1e-12);
    
    // Imposta tolleranza
    void setTolerance(double tol);
    
    // Imposta limiti passo
    void setStepLimits(double h_min, double h_max);
    
    // Integrazione principale
    State integrate(
        const State& y0,           // Stato iniziale [r, v]
        double t0,                 // Tempo iniziale (JD)
        double t1,                 // Tempo finale (JD)
        AccelFunction accel,       // Funzione accelerazione
        IntegrationStats* stats    // Statistiche (opzionale)
    );
    
    // Singolo passo adattivo
    State step_adaptive(
        const State& y,
        double t,
        double& h,                 // Passo (in/out)
        AccelFunction accel,
        bool& accepted             // Passo accettato?
    );
};
```

### 4.2 Struttura `State`

```cpp
struct State {
    std::array<double, 3> r;  // Posizione [AU]
    std::array<double, 3> v;  // Velocità [AU/day]
    
    // Operatori aritmetici
    State operator+(const State& other) const;
    State operator*(double scalar) const;
    
    // Norma posizione
    double norm_r() const;
};
```

### 4.3 Struttura `IntegrationStats`

```cpp
struct IntegrationStats {
    int steps_accepted;     // Passi accettati
    int steps_rejected;     // Passi rifiutati
    int func_evaluations;   // Valutazioni f(t,y)
    double h_min;           // Passo minimo usato
    double h_max;           // Passo massimo usato
};
```

### 4.4 Tipo `AccelFunction`

```cpp
using AccelFunction = std::function<
    std::array<double, 3>(double t, const std::array<double, 3>& r, 
                          const std::array<double, 3>& v)
>;
```

---

## 5. Modello Dinamico

### 5.1 Equazione del Moto

L'accelerazione totale è:

$$\ddot{\mathbf{r}} = \ddot{\mathbf{r}}_\odot + \ddot{\mathbf{r}}_{\text{pianeti}} + \ddot{\mathbf{r}}_{\text{AST17}} + \ddot{\mathbf{r}}_{\text{rel}}$$

### 5.2 Termine Kepleriano (Sole)

$$\ddot{\mathbf{r}}_\odot = -\frac{GM_\odot}{r^3} \mathbf{r}$$

con $GM_\odot = 0.01720209895^2$ AU³/day² (costante di Gauss).

### 5.3 Perturbazioni Planetarie

$$\ddot{\mathbf{r}}_{\text{pianeti}} = \sum_{i=1}^{8} GM_i \left( \frac{\mathbf{r}_i - \mathbf{r}}{|\mathbf{r}_i - \mathbf{r}|^3} - \frac{\mathbf{r}_i}{r_i^3} \right)$$

**Pianeti inclusi** (effemeridi Simon et al. 1994):

| Pianeta | GM [AU³/day²] | Precisione |
|---------|---------------|------------|
| Mercurio | 4.9125e-11 | ~1" |
| Venere | 7.2435e-10 | ~1" |
| Terra-Luna | 8.9970e-10 | ~1" |
| Marte | 9.5495e-11 | ~1" |
| Giove | 2.8253e-07 | ~1" |
| Saturno | 8.4597e-08 | ~1" |
| Urano | 1.2920e-08 | ~5" |
| Nettuno | 1.5244e-08 | ~10" |

### 5.4 Perturbazioni AST17

$$\ddot{\mathbf{r}}_{\text{AST17}} = \sum_{j=1}^{16} GM_j \left( \frac{\mathbf{r}_j - \mathbf{r}}{|\mathbf{r}_j - \mathbf{r}|^3} - \frac{\mathbf{r}_j}{r_j^3} \right)$$

**Asteroidi AST17**:

| # | Nome | GM [AU³/day²] | Massa relativa |
|---|------|---------------|----------------|
| 1 | Ceres | 1.392e-13 | 1.000 |
| 2 | Pallas | 3.036e-14 | 0.218 |
| 4 | Vesta | 3.803e-14 | 0.273 |
| 10 | Hygiea | 1.231e-14 | 0.088 |
| 704 | Interamnia | 5.467e-15 | 0.039 |
| 511 | Davida | 5.467e-15 | 0.039 |
| 52 | Europa | 3.893e-15 | 0.028 |
| 15 | Eunomia | 4.565e-15 | 0.033 |
| 16 | Psyche | 3.389e-15 | 0.024 |
| 3 | Juno | 3.893e-15 | 0.028 |
| 87 | Sylvia | 2.155e-15 | 0.015 |
| 88 | Thisbe | 2.491e-15 | 0.018 |
| 31 | Euphrosyne | 2.323e-15 | 0.017 |
| 324 | Bamberga | 1.515e-15 | 0.011 |
| 451 | Patientia | 1.851e-15 | 0.013 |
| 65 | Cybele | 1.683e-15 | 0.012 |

### 5.5 Correzione Relativistica (Schwarzschild)

$$\ddot{\mathbf{r}}_{\text{rel}} = \frac{GM_\odot}{c^2 r^3} \left[ \left(4\frac{GM_\odot}{r} - v^2\right)\mathbf{r} + 4(\mathbf{r}\cdot\mathbf{v})\mathbf{v} \right]$$

con $c = 173.1446$ AU/day.

---

## 6. Esempio d'Uso

### 6.1 Propagazione Base

```cpp
#include "rkf78_integrator.hpp"

int main() {
    // Stato iniziale (asteroide Sierks)
    State y0;
    y0.r = {1.082716234178542, 3.086569667750047, -0.091582857678774};
    y0.v = {-0.008923488396027, 0.002807706967566, 0.000404145908519};
    
    // Epoche
    double jd0 = 2461000.5;    // Epoca elementi
    double jd1 = 2461008.0913; // Target
    
    // Crea integratore
    RKF78Integrator rk78(1e-12);  // Tolleranza 1e-12
    
    // Funzione accelerazione (include pianeti, AST17, relatività)
    auto accel = [](double t, const auto& r, const auto& v) {
        return compute_acceleration(t, r, v);
    };
    
    // Integra
    IntegrationStats stats;
    State y1 = rk78.integrate(y0, jd0, jd1, accel, &stats);
    
    // Risultati
    std::cout << "Posizione finale: " << y1.r[0] << ", " 
              << y1.r[1] << ", " << y1.r[2] << " AU\n";
    std::cout << "Passi: " << stats.steps_accepted << "\n";
    std::cout << "Valutazioni: " << stats.func_evaluations << "\n";
    
    return 0;
}
```

### 6.2 Output Tipico

```
================================================================
  TEST INTEGRATORE RKF78 (Fehlberg 1968)
================================================================

ASTEROIDE: (17030) Sierks
Epoca:  JD 2461000.5000
Target: JD 2461008.0913
Δt = 7.5913 giorni

Passi accettati: 12
Passi rifiutati: 0
Valutazioni f:   157
Passo min:       7.591e-02 giorni
Passo max:       7.591e-01 giorni

POSIZIONE FINALE:
  RA  = 04 53 11.597
  Dec = +20 19 26.23

CONFRONTO CON JPL HORIZONS:
  Errore: ΔRA*cos(δ)=4.88", ΔDec=0.43" → Tot=4.89"

ROUND-TRIP TEST:
  Errore posizione: 4.4e-16 AU = 0.07 mm ✅
  Errore velocità:  7.2e-18 AU/day = 0.01 nm/s ✅
================================================================
```

---

## 7. Confronto con Altri Integratori

| Integratore | Ordine | Stadi | Adattivo | Precisione | Efficienza |
|-------------|--------|-------|----------|------------|------------|
| RK4 | 4 | 4 | No | ★★★ | ★★★★ |
| RKF45 | 4(5) | 6 | Sì | ★★★★ | ★★★★ |
| **RKF78** | **7(8)** | **13** | **Sì** | **★★★★★** | **★★★★** |
| Bulirsch-Stoer | Variabile | Variabile | Sì | ★★★★★ | ★★★ |

### Vantaggi RKF78

1. **Alta precisione**: Ordine 7 con stima errore ordine 8
2. **Efficienza**: Solo 12 passi per propagazioni di settimane
3. **Robustezza**: Controllo adattivo del passo
4. **Reversibilità**: Errore round-trip a livello di macchina

---

## 8. Limitazioni e Avvertenze

### 8.1 Effemeridi Planetarie

Le posizioni planetarie usano le formule analitiche di **Simon et al. (1994)**:
- Precisione: 1-20" a seconda del pianeta
- Per massima precisione, usare effemeridi JPL DE441

### 8.2 Elementi Orbitali

L'accuratezza finale dipende dalla qualità degli elementi orbitali:
- Elementi MPC: precisione tipica ~5-10"
- Per sub-arcsecond, usare elementi JPL o fit differenziale

### 8.3 Perturbazioni Aggiuntive

Non incluse in questa versione:
- Yarkovsky effect
- Radiazione solare
- Maree terrestri
- Oblateness planetaria (J2)

---

## 9. Riferimenti

1. **Fehlberg, E.** (1968). "Classical fifth-, sixth-, seventh-, and eighth-order Runge-Kutta formulas with stepsize control." NASA TR R-287.

2. **Simon, J.L., et al.** (1994). "Numerical expressions for precession formulae and mean elements for the Moon and the planets." Astronomy & Astrophysics, 282, 663-683.

3. **Standish, E.M.** (1998). "JPL Planetary and Lunar Ephemerides, DE405/LE405." JPL IOM 312.F-98-048.

4. **Farnocchia, D., et al.** (2015). "Asteroid masses from the analysis of their close encounters with other asteroids." Icarus, 245, 94-111.

---

## 10. Supporto Formati AstDyS

### 10.1 Elementi Equinoziali

Il propagatore supporta nativamente il formato **OEF2.0** di AstDyS:

```cpp
// Struttura elementi equinoziali
struct EquinoctialElements {
    double a;       // Semiasse maggiore [AU]
    double h;       // e*sin(ϖ)
    double k;       // e*cos(ϖ)
    double p;       // tan(i/2)*sin(Ω)
    double q;       // tan(i/2)*cos(Ω)
    double lambda;  // Longitudine media [deg]
    double epoch;   // Epoca [MJD]
    bool is_mjd;    // true se MJD
};

// Uso
EquinoctialElements eq;
eq.a = 2.6808535916678031;
eq.h = 0.032872036471001;
eq.k = 0.036254405825130;
eq.p = 0.103391596538937;
eq.q = -0.042907901689093;
eq.lambda = 235.8395861037268;
eq.epoch = 61000.0;
eq.is_mjd = true;

AstDynPropagator prop;
State y0 = prop.equinoctialToState(eq);
```

### 10.2 Conversione Elementi

```cpp
// Equinoziali → Kepleriani
OrbitalElements kep = prop.equinoctialToKeplerian(eq);

// Formule di conversione:
// e = sqrt(h² + k²)
// i = 2·atan(sqrt(p² + q²))
// Ω = atan2(p, q)
// ω = atan2(h, k) - Ω
// M = λ - atan2(h, k)
```

### 10.3 Download Dati AstDyS

```bash
# Elementi orbitali
curl -s "https://newton.spacedys.com/~astdys2/epoch/numbered/11/11234.eq1"

# Osservazioni con residui
curl -s "https://newton.spacedys.com/~astdys2/mpcobs/numbered/11/11234.rwo"
```

Per dettagli completi, vedere **ASTDYS_FORMAT.md**.

---

## 11. Guida all'Integrazione

Questa sezione spiega come integrare `AstDynPropagator` nel proprio progetto.

### 11.1 Requisiti

- **Compilatore**: C++17 o superiore
- **Librerie**: Solo standard library (no dipendenze esterne)
- **Sistema**: Linux, macOS, Windows

### 11.2 Installazione

Il propagatore è un **header-only** file singolo. Per usarlo:

```bash
# Copia il file nel tuo progetto
cp astdyn/tools/astdyn_propagator.cpp my_project/

# Oppure includi come header (rinomina in .hpp)
cp astdyn/tools/astdyn_propagator.cpp my_project/astdyn_propagator.hpp
```

### 11.3 Uso Base: Propagazione da Elementi Kepleriani

```cpp
#include <iostream>
#include "astdyn_propagator.cpp"  // Include tutto il codice

int main() {
    using namespace astdyn;
    
    // 1. Definisci elementi orbitali
    OrbitalElements elem;
    elem.name = "MyAsteroid";
    elem.a = 2.5;        // Semiasse maggiore [AU]
    elem.e = 0.15;       // Eccentricità
    elem.i = 10.5;       // Inclinazione [deg]
    elem.Omega = 120.0;  // Longitudine nodo [deg]
    elem.omega = 45.0;   // Argomento perielio [deg]
    elem.M = 180.0;      // Anomalia media [deg]
    elem.epoch = 2460000.5;  // Epoca [JD]
    
    // 2. Crea propagatore
    AstDynPropagator prop(1e-12);  // Tolleranza
    
    // 3. Converti in stato cartesiano ICRF
    State y0 = prop.elementsToState(elem);
    
    // 4. Propaga a data target
    double jd_target = 2460030.5;  // +30 giorni
    PropagationStats stats;
    State y1 = prop.propagate(y0, elem.epoch, jd_target, &stats);
    
    // 5. Calcola coordinate equatoriali
    EquatorialCoords coords = prop.getEquatorialCoords(y1, jd_target);
    
    // 6. Output
    std::cout << "RA  = " << coords.formatRA() << "\n";
    std::cout << "Dec = " << coords.formatDec() << "\n";
    std::cout << "Δ   = " << coords.dist << " AU\n";
    std::cout << "Passi: " << stats.steps_accepted << "\n";
    
    return 0;
}
```

**Compilazione:**
```bash
g++ -std=c++17 -O3 -o my_propagator my_propagator.cpp
```

### 11.4 Uso con Elementi AstDyS (Equinoziali)

```cpp
#include "astdyn_propagator.cpp"

int main() {
    using namespace astdyn;
    
    // Elementi equinoziali da file AstDyS .eq1
    EquinoctialElements eq;
    eq.a = 2.6808535916678031;      // a [AU]
    eq.h = 0.032872036471001;       // e*sin(ϖ)
    eq.k = 0.036254405825130;       // e*cos(ϖ)
    eq.p = 0.103391596538937;       // tan(i/2)*sin(Ω)
    eq.q = -0.042907901689093;      // tan(i/2)*cos(Ω)
    eq.lambda = 235.8395861037268;  // Longitudine media [deg]
    eq.epoch = 61000.0;             // MJD
    eq.is_mjd = true;
    eq.name = "(11234) 1999 JS82";
    
    // Crea propagatore e converti
    AstDynPropagator prop(1e-12);
    
    // Metodo 1: Conversione diretta a stato
    State y0 = prop.equinoctialToState(eq);
    
    // Metodo 2: Prima converti a Kepleriano, poi a stato
    OrbitalElements kep = prop.equinoctialToKeplerian(eq);
    State y0_alt = prop.elementsToState(kep);
    
    // Propaga
    double jd0 = eq.getJD();  // Converte MJD → JD
    double jd1 = jd0 + 30;    // +30 giorni
    
    State y1 = prop.propagate(y0, jd0, jd1);
    EquatorialCoords coords = prop.getEquatorialCoords(y1, jd1);
    
    std::cout << "RA  = " << coords.formatRA() << "\n";
    std::cout << "Dec = " << coords.formatDec() << "\n";
    
    return 0;
}
```

### 11.5 Configurazione Avanzata

```cpp
AstDynPropagator prop(1e-12);  // Tolleranza base

// Configura perturbazioni
prop.usePlanets(true);     // Perturbazioni 8 pianeti (default: true)
prop.useAST17(true);       // Perturbazioni 16 asteroidi (default: true)
prop.useRelativity(true);  // Correzione Schwarzschild (default: true)

// Configura limiti passo
prop.setStepLimits(0.001, 10.0);  // h_min, h_max [giorni]

// Configura tolleranza
prop.setTolerance(1e-14);  // Per massima precisione
```

### 11.6 Propagazione Multi-Epoca (Effemeridi)

```cpp
#include <vector>

// Genera effemeridi per un mese
std::vector<EquatorialCoords> generateEphemeris(
    const OrbitalElements& elem,
    double jd_start,
    double jd_end,
    double step_days)
{
    AstDynPropagator prop(1e-12);
    State y = prop.elementsToState(elem);
    double t = elem.epoch;
    
    std::vector<EquatorialCoords> ephem;
    
    for (double jd = jd_start; jd <= jd_end; jd += step_days) {
        // Propaga allo step corrente
        y = prop.propagate(y, t, jd);
        t = jd;
        
        // Calcola coordinate
        EquatorialCoords coords = prop.getEquatorialCoords(y, jd);
        ephem.push_back(coords);
    }
    
    return ephem;
}

int main() {
    OrbitalElements elem = {...};
    
    auto ephem = generateEphemeris(elem, 2461000.5, 2461030.5, 1.0);
    
    for (size_t i = 0; i < ephem.size(); i++) {
        std::cout << "Day " << i << ": "
                  << ephem[i].formatRA() << " "
                  << ephem[i].formatDec() << "\n";
    }
    
    return 0;
}
```

### 11.7 Lettura File AstDyS

```cpp
#include <fstream>
#include <sstream>

// Parser semplice per file .eq1
EquinoctialElements parseEq1File(const std::string& filename) {
    EquinoctialElements eq;
    std::ifstream file(filename);
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.find("KEP") != std::string::npos) {
            // Leggi i 6 elementi dalla riga successiva
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> eq.a >> eq.h >> eq.k >> eq.p >> eq.q >> eq.lambda;
        }
        if (line.find("MJD") != std::string::npos) {
            size_t pos = line.find("MJD");
            eq.epoch = std::stod(line.substr(pos + 3));
            eq.is_mjd = true;
        }
    }
    
    return eq;
}

int main() {
    auto eq = parseEq1File("11234.eq1");
    
    AstDynPropagator prop;
    State y0 = prop.equinoctialToState(eq);
    // ...
}
```

### 11.8 Integrazione con Python (ctypes)

Per usare il propagatore da Python:

```cpp
// propagator_wrapper.cpp
extern "C" {
    #include "astdyn_propagator.cpp"
    
    void propagate_asteroid(
        double a, double e, double i, double Omega, double omega, double M,
        double epoch, double target_jd,
        double* ra_out, double* dec_out, double* dist_out)
    {
        using namespace astdyn;
        
        OrbitalElements elem;
        elem.a = a; elem.e = e; elem.i = i;
        elem.Omega = Omega; elem.omega = omega;
        elem.M = M; elem.epoch = epoch;
        
        AstDynPropagator prop(1e-12);
        State y0 = prop.elementsToState(elem);
        State y1 = prop.propagate(y0, epoch, target_jd);
        EquatorialCoords coords = prop.getEquatorialCoords(y1, target_jd);
        
        *ra_out = coords.ra;
        *dec_out = coords.dec;
        *dist_out = coords.dist;
    }
}
```

```python
# propagator.py
import ctypes

lib = ctypes.CDLL('./libpropagator.so')
lib.propagate_asteroid.argtypes = [
    ctypes.c_double] * 8 + [ctypes.POINTER(ctypes.c_double)] * 3

def propagate(a, e, i, Omega, omega, M, epoch, target):
    ra = ctypes.c_double()
    dec = ctypes.c_double()
    dist = ctypes.c_double()
    
    lib.propagate_asteroid(
        a, e, i, Omega, omega, M, epoch, target,
        ctypes.byref(ra), ctypes.byref(dec), ctypes.byref(dist))
    
    return ra.value, dec.value, dist.value

# Uso
ra, dec, dist = propagate(3.175, 0.045, 2.9, 104.2, 100.5, 229.8,
                          2461000.5, 2461030.5)
print(f"RA={ra:.6f} rad, Dec={dec:.6f} rad, Dist={dist:.4f} AU")
```

### 11.9 Best Practices

1. **Tolleranza**: Usa `1e-12` per precisione sub-arcsecond, `1e-10` per velocità
2. **Perturbazioni**: Disabilita AST17 per propagazioni molto brevi (< 1 giorno)
3. **Frame**: Lo stato è sempre in ICRF - non mescolare con eclittico
4. **Epoca**: Usa sempre JD (non MJD) nelle funzioni propagate
5. **Round-trip**: Testa sempre la reversibilità per validare i risultati

### 11.10 Troubleshooting

| Problema | Causa | Soluzione |
|----------|-------|-----------|
| Errore > 1' | Elementi sbagliati | Verifica frame (ECLM vs ICRF) |
| Round-trip > 1m | Tolleranza bassa | Aumenta a 1e-12 o 1e-14 |
| Crash/NaN | Eccentricità > 1 | Verifica elementi orbitali |
| Lentezza | Tolleranza troppo stretta | Usa 1e-10 per preview |
| Offset costante | Posizione Terra errata | Verifica effemeridi Terra |

---

## 12. File Sorgente

| File | Descrizione |
|------|-------------|
| `tools/astdyn_propagator.cpp` | Classe AstDynPropagator completa |
| `tools/test_rkf78.cpp` | Implementazione standalone RKF78 |
| `tools/test_jpl_full.cpp` | Test validazione vs JPL (7 anni) |
| `tools/test_astdys_full.cpp` | Test con elementi AstDyS |
| `data/11234.eq1` | Elementi esempio (11234) |
| `data/11234.rwo` | Osservazioni esempio |

---

**© 2025 AstDyS Team - Dipartimento di Matematica, Università di Pisa**
