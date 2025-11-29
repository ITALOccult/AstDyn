# Formato File AstDyS - Guida Completa

## Panoramica

Il database **AstDyS** (Asteroids Dynamic Site) dell'Università di Pisa fornisce elementi orbitali
ad alta precisione per asteroidi e oggetti NEO. Questo documento descrive i formati supportati
e le conversioni implementate nel propagatore AstDyn.

## Sorgenti Dati

### URL per Download

```
# Elementi orbitali (formato equinoziale)
https://newton.spacedys.com/~astdys2/epoch/numbered/{XX}/{NNNNN}.eq1

# Osservazioni con residui (formato RWO)
https://newton.spacedys.com/~astdys2/mpcobs/numbered/{XX}/{NNNNN}.rwo
```

Dove:
- `{XX}` = prime 2 cifre del numero (es. `11` per 11234)
- `{NNNNN}` = numero completo dell'asteroide

### Esempio

```bash
# Asteroide (11234)
curl -s "https://newton.spacedys.com/~astdys2/epoch/numbered/11/11234.eq1" > 11234.eq1
curl -s "https://newton.spacedys.com/~astdys2/mpcobs/numbered/11/11234.rwo" > 11234.rwo
```

---

## Formato OEF2.0 (Elementi Equinoziali)

### Struttura File .eq1

```
format  = 'OEF2.0'       ! file format
rectype = 'ML'           ! record type (1L/ML)
refsys  = ECLM J2000     ! default reference system
END_OF_HEADER
11234
! Equinoctial elements: a, e*sin(LP), e*cos(LP), tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.
 EQU   2.6808535916678031E+00   0.032872036471001   0.036254405825130   0.103391596538937   -0.042907901689093   235.8395861037268
 MJD     61000.000000000 TDT
 MAG  12.874  0.150
```

### Elementi Equinoziali

| Simbolo | Nome | Formula | Unità |
|---------|------|---------|-------|
| `a` | Semiasse maggiore | - | AU |
| `h` | Componente eccentricità | e·sin(ϖ) | - |
| `k` | Componente eccentricità | e·cos(ϖ) | - |
| `p` | Componente inclinazione | tan(i/2)·sin(Ω) | - |
| `q` | Componente inclinazione | tan(i/2)·cos(Ω) | - |
| `λ` | Longitudine media | Ω + ω + M | gradi |

Dove:
- `ϖ = Ω + ω` = longitudine del perielio
- `Ω` = longitudine del nodo ascendente
- `ω` = argomento del perielio
- `M` = anomalia media

### Frame di Riferimento

- **ECLM J2000**: Eclittica media J2000.0
- **Epoca**: MJD (Modified Julian Date) in TDT (Terrestrial Dynamical Time)

---

## Conversione Elementi

### Equinoziali → Kepleriani

```cpp
// Input: a, h, k, p, q, λ (AstDyS)
// Output: a, e, i, Ω, ω, M (Kepler)

e = sqrt(h² + k²)
i = 2·atan(sqrt(p² + q²))
Ω = atan2(p, q)
ϖ = atan2(h, k)           // longitudine del perielio
ω = ϖ - Ω                 // argomento del perielio
M = λ - ϖ                 // anomalia media
```

### Esempio (11234)

**Input AstDyS:**
```
a = 2.6808535916678031 AU
h = 0.032872036471001
k = 0.036254405825130
p = 0.103391596538937
q = -0.042907901689093
λ = 235.8395861037268°
```

**Output Kepleriani:**
```
a = 2.6808536 AU
e = 0.0489383
i = 12.7744°
Ω = 112.5386°
ω = 289.6601°
M = 193.6408°
```

---

## Conversione a Stato Cartesiano

### Kepleriani → Eclittico

1. Risolvi equazione di Keplero: `E - e·sin(E) = M`
2. Calcola anomalia vera: `ν = 2·atan2(sqrt(1+e)·sin(E/2), sqrt(1-e)·cos(E/2))`
3. Calcola distanza: `r = a·(1 - e·cos(E))`
4. Posizione nel piano orbitale:
   ```
   x_orb = r·cos(ν)
   y_orb = r·sin(ν)
   ```
5. Velocità nel piano orbitale:
   ```
   vx_orb = -n·a·sin(E) / (1 - e·cos(E))
   vy_orb = n·a·sqrt(1-e²)·cos(E) / (1 - e·cos(E))
   ```
6. Rotazione al frame eclittico (matrice P, Q, W)

### Eclittico → ICRF (Equatoriale)

```cpp
// Obliquità J2000
ε = 23.4392911°

// Rotazione
x_icrf = x_ecl
y_icrf = cos(ε)·y_ecl - sin(ε)·z_ecl
z_icrf = sin(ε)·y_ecl + cos(ε)·z_ecl
```

---

## Formato RWO (Osservazioni)

### Struttura File

```
! Object 11234 - (11234) 1999 JS82
! Obs.: 6844 from 1989-10-17.16014 to 2025-11-26.24920
! Format: [designation] [date] [RA] [Dec] [mag] [observatory] [rms_ra] [rms_dec] [residuals]

     K99J82S  1989 10 17.16014 1 02 58 05.43 +00 44 57.2          16.5 Ru~05AZ 691  0.30  0.30  +0.16 -0.33
     K99J82S  1989 10 17.19653 1 02 58 06.02 +00 45 00.9          16.7 Ru~05AZ 691  0.30  0.30  +0.18 -0.35
```

### Campi Principali

| Campo | Descrizione |
|-------|-------------|
| Designation | Designazione provvisoria (K99J82S = 1999 JS82) |
| Date | Data UTC (YYYY MM DD.ddddd) |
| RA | Ascensione retta (HH MM SS.ss) |
| Dec | Declinazione (±DD MM SS.s) |
| Mag | Magnitudine |
| Observatory | Codice osservatorio MPC |
| RMS RA/Dec | RMS stimato [arcsec] |
| Residui | O-C in RA e Dec [arcsec] |

---

## Uso nel Propagatore

### Caricamento Elementi AstDyS

```cpp
#include "astdyn_propagator.cpp"

using namespace astdyn;

// Definisci elementi equinoziali
EquinoctialElements eq;
eq.a = 2.6808535916678031;
eq.h = 0.032872036471001;
eq.k = 0.036254405825130;
eq.p = 0.103391596538937;
eq.q = -0.042907901689093;
eq.lambda = 235.8395861037268;
eq.epoch = 61000.0;
eq.is_mjd = true;
eq.name = "(11234)";

// Crea propagatore
AstDynPropagator prop(1e-12);

// Converti a stato ICRF
State y0 = prop.equinoctialToState(eq);

// Propaga
double jd_target = 2461030.5;  // 2025-Dec-21
State y1 = prop.propagate(y0, eq.getJD(), jd_target);

// Ottieni coordinate
EquatorialCoords coords = prop.getEquatorialCoords(y1, jd_target);
std::cout << "RA:  " << coords.formatRA() << std::endl;
std::cout << "Dec: " << coords.formatDec() << std::endl;
```

### Conversione Diretta

```cpp
// Converti elementi equinoziali a kepleriani
OrbitalElements kep = prop.equinoctialToKeplerian(eq);

std::cout << "e = " << kep.e << std::endl;
std::cout << "i = " << kep.i << "°" << std::endl;
std::cout << "Ω = " << kep.Omega << "°" << std::endl;
```

---

## Validazione

### Confronto con JPL Horizons

| Epoca | Intervallo | Errore |
|-------|------------|--------|
| 2025-Oct-22 | -30 giorni | 11.4" |
| 2025-Dec-21 | +30 giorni | 11.4" |

### Test Round-Trip

```
Propagazione: epoca → +30d → -30d → epoca
Errore: 15 metri
```

---

## Note Tecniche

### Precisione

- Gli elementi AstDyS sono derivati da osservazioni ottiche e hanno precisione tipica di ~0.1" nel residuo
- La propagazione introduce errori aggiuntivi dovuti a:
  - Effemeridi planetarie approssimate (~10")
  - Mancanza di perturbatori minori
  - Effetti non-gravitazionali per NEO

### Frame di Riferimento

- AstDyS usa **ECLM J2000** (eclittica media)
- JPL Horizons usa **ICRF** (equatoriale)
- Il propagatore converte automaticamente da eclittica a ICRF

### Scala Temporale

- AstDyS: TDT (Terrestrial Dynamical Time)
- JPL: TDB (Barycentric Dynamical Time)
- Differenza: ~1.7 ms (trascurabile per propagazione)

---

## Riferimenti

1. **AstDyS**: https://newton.spacedys.com/astdys/
2. **NEODyS**: https://newton.spacedys.com/neodys/
3. **OEF Format**: AstDyS documentation
4. **JPL Horizons**: https://ssd.jpl.nasa.gov/horizons/

---

*Documento aggiornato: 29 Novembre 2025*
*Versione: 1.0*
