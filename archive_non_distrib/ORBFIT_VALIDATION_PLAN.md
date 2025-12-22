# Piano di Validazione: AstDyn vs OrbFit

## Obiettivo
Verificare che il propagatore e l'integratore di AstDyn producano risultati equivalenti a quelli di OrbFit, lo standard de facto per la determinazione orbitale asteroidale.

## Metodologia

### 1. Test di Propagazione Orbitale
**Input:** Elementi orbitali da file `.eq1` (formato AstDyS/OrbFit)
**Output:** Stati cartesiani a vari istanti temporali

**Confronto:**
- Posizione (X, Y, Z) in AU
- Velocità (VX, VY, VZ) in AU/day
- Errore relativo e assoluto

**Configurazione:**
- Integratore AstDyn: RKF78 (tolleranza 1e-12)
- Integratore OrbFit: Radau15 o RKG (configurazione standard)
- Perturbazioni: Tutti i pianeti + relatività
- Frame: ECLM J2000 (mean ecliptic J2000)

### 2. Test di Residui Astrometrici
**Input:** 
- Elementi orbitali
- Osservazioni astrometriche da file `.rwo` (formato OrbFit)

**Output:**
- Residui O-C (Observed - Computed) in RA e Dec
- RMS dei residui
- Chi-quadro

**Confronto:**
- Residui individuali per ogni osservazione
- Statistiche globali (RMS, chi²)

### 3. Configurazione Test

**Asteroide di Test:** 17030 Sierks
- File elementi: `astdyn/tools/203_astdys.eq1` (alias per 17030)
- File osservazioni: `astdyn/tools/203.rwo`
- Periodo: 1990-2025 (3765 osservazioni)

**Epoche di Propagazione:**
- Epoca iniziale: MJD 61000.0 (25 Nov 2025)
- Intervallo: 30 giorni
- Step: 1 giorno

## Implementazione

### Fase 1: Wrapper Fortran per OrbFit
Creare un programma Fortran che:
1. Legge elementi da `.eq1`
2. Propaga con il propagatore OrbFit
3. Scrive stati cartesiani in formato CSV

```fortran
PROGRAM orbfit_propagator_test
  ! Legge 203_astdys.eq1
  ! Propaga da MJD 61000 a 61030 (step 1 giorno)
  ! Output: orbfit_states.csv
END PROGRAM
```

### Fase 2: Test C++ AstDyn
Programma C++ che:
1. Legge stessi elementi da `.eq1`
2. Propaga con AstDyn RKF78
3. Confronta con output OrbFit

```cpp
// test_astdyn_vs_orbfit.cpp
// 1. Load orbital elements
// 2. Propagate with AstDyn
// 3. Read OrbFit reference data
// 4. Compute differences
// 5. Report statistics
```

### Fase 3: Test Residui
Confronto calcolo residui:
1. OrbFit: usa routine native `pred_obs` + `compute_residuals`
2. AstDyn: usa `ResidualCalculator` 
3. Confronto residuo per residuo

## Criteri di Successo

| Metrica | Soglia Accettazione | Note |
|:--------|:-------------------:|:-----|
| **Errore Posizione RMS** | < 1 km | Su 30 giorni |
| **Errore Posizione Max** | < 10 km | Singolo punto |
| **Errore Velocità RMS** | < 1 mm/s | Equivalente |
| **Differenza Residui** | < 0.01 arcsec | Per singola osservazione |
| **Differenza RMS Totale** | < 0.001 arcsec | Statistica globale |

## File da Creare

1. `orbfit_propagator_wrapper.f90` - Wrapper Fortran per propagazione
2. `test_astdyn_vs_orbfit.cpp` - Test C++ di confronto
3. `compile_orbfit_test.sh` - Script compilazione Fortran
4. `run_comparison.sh` - Pipeline completa
5. `ORBFIT_COMPARISON_REPORT.md` - Report risultati

## Dipendenze

**OrbFit:**
- Compilatore Fortran (gfortran)
- Librerie OrbFit (già presenti in `astdyn/tools/`)
- File dati: `.eq1`, `.rwo`

**AstDyn:**
- Già compilato e testato
- Parser `.eq1` funzionante
- Propagatore RKF78 validato vs JPL

## Note Implementative

### Sistemi di Coordinate
**Attenzione:** OrbFit usa ECLM J2000 (eclittica media J2000) come frame di riferimento standard per gli elementi orbitali. AstDyn deve:
1. Leggere elementi in ECLM J2000
2. Propagare in ECLM J2000
3. Confrontare stati nello stesso frame

### Conversioni Temporali
Entrambi usano:
- TDB (Barycentric Dynamical Time) per la propagazione
- UTC per le osservazioni
- Conversioni UTC↔TDB devono essere identiche

### Perturbazioni
Verificare che entrambi includano:
- Perturbazioni planetarie (8 pianeti)
- Relatività generale (correzione Schwarzschild)
- Eventualmente: asteroidi massivi (opzionale)

## Timeline

| Fase | Durata | Deliverable |
|:-----|:------:|:------------|
| 1. Wrapper Fortran | 2h | `orbfit_propagator_wrapper.f90` compilato |
| 2. Test C++ | 3h | `test_astdyn_vs_orbfit.cpp` funzionante |
| 3. Esecuzione e analisi | 1h | CSV con confronti |
| 4. Report | 1h | `ORBFIT_COMPARISON_REPORT.md` |
| **Totale** | **~7h** | Validazione completa |

## Prossimi Passi

1. ✅ Identificare file di test (203_astdys.eq1, 203.rwo)
2. ⏳ Creare wrapper Fortran per propagazione OrbFit
3. ⏳ Implementare test C++ di confronto
4. ⏳ Eseguire confronto e generare report
5. ⏳ Certificare compatibilità AstDyn-OrbFit
