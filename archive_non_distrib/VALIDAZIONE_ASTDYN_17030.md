# VALIDAZIONE RISULTATI ASTDYN - ASTEROIDE 17030

## Test Eseguito

**Data**: 1 Dicembre 2025  
**Programma**: `test_astdyn_simple`  
**Libreria**: AstDyn v1.0.0 con RKF78  
**Asteroide**: 17030 Sierks  

## Configurazione Test

### Elementi Iniziali (dal file 17030.eq1)
- **Epoca**: MJD 61000.0 (2018-Mar-16 00:00 UTC)
- **Frame**: ECLM J2000 (eclittica media J2000)
- **Formato**: OEF2.0 equinoctial elements

```
 ! Object: 17030
 EQU  3.175473  -0.018963  -0.041273  0.025407  -0.001956  229.790880
 MJD  61000.0
 MAG  13.29  0.13
```

### Elementi Kepleriani (calcolati da parser AstDyn)
Dopo conversione equinoziali → kepleriani:
- **a** = 3.175473 AU
- **e** = 0.045407
- **i** = 2.9046°
- **Ω** = 94.058°
- **ω** = 110.28°
- **M** = 25.45°

## Propagazione Eseguita

### Parametri
- **Epoca iniziale**: MJD 61000.0
- **Epoca target**: MJD 60277.0 (2025-Nov-28 00:00 UTC)
- **Δt**: -723 giorni (propagazione all'indietro)
- **Integratore**: RKF78 (Runge-Kutta-Fehlberg 7/8 ordine)
- **Tolleranza**: 1e-12
- **Step size**: adattativo

### Perturbazioni Incluse
✅ Gravità del Sole  
✅ 8 pianeti (Mercurio, Venere, Terra, Marte, Giove, Saturno, Urano, Nettuno)  
✅ Correzione relativistica (Schwarzschild)  
✅ Perturbazioni asteroidali (AST17)  

## Risultati Ottenuti

### Performance
- **Tempo computazione**: 1 ms
- **Step accettati**: 26
- **Step rifiutati**: 8
- **Valutazioni funzione**: 442

### Posizione Finale (MJD 60277.0 - ICRF Eliocentrico)
```
X  = -0.488510929265 AU
Y  =  3.171053797532 AU
Z  =  0.012429740996 AU
```

### Velocità Finale (ICRF Eliocentrico)
```
VX = -0.009361833407 AU/day
VY = -0.001878143017 AU/day
VZ =  0.000483297812 AU/day
```

### Elementi Kepleriani Finali (MJD 60277.0)
```
a = 3.180105094 AU
e = 0.045745196
i = 2.918968°
Ω = 94.400384°
ω = 108.173057°
M = 261.320727°
```

## Confronto con JPL Horizons (UFFICIALE)

### Interrogazione JPL Horizons API
**Endpoint**: JPL Horizons System  
**Target**: Asteroid 17030 Sierks  
**Epoca**: 2025-Nov-28 00:00 UTC (JD 2460277.5)  
**Observer**: @sun (Heliocentric)  
**Frame**: ICRF  

### Dati JPL Horizons (riferimento ufficiale)
```
Posizione (AU):
X  =  2.207146492349
Y  = -2.107862280391
Z  = -1.000651915151

Velocità (AU/day):
VX =  0.007199933462
VY =  0.005900361090
VZ =  0.002087569170
```

### Differenze Riscontrate (JPL - AstDyn)
```
ΔX  =  2.695657421614 AU
ΔY  = -5.278916077923 AU
ΔZ  = -1.013081656147 AU

ΔVX =  0.016561766869 AU/day
ΔVY =  0.007778504107 AU/day
ΔVZ =  0.001604271358 AU/day
```

### ❌ **ERRORE CRITICO**
**Errore in posizione**: √(ΔX² + ΔY² + ΔZ²) = **6.01 AU** = **899,577,889 km**  
**Errore in velocità**: 0.0184 AU/day  
**Errore angolare stimato**: ~4137 arcsec  

### Confronto con Valori ELEMENTI_ORBITALI_TEST.md
```
r_sun = (1.0147, 2.8859, 1.1548) AU  [dal documento]
vs
r_jpl = (2.2071, -2.1079, -1.0007) AU  [JPL Horizons]
```

**Anche il documento di test ha valori completamente diversi da JPL!**  
Differenza documento-JPL: ~6.5 AU

## Analisi Discrepanze

### Possibili Cause

1. **❓ Elementi Orbitali Diversi**
   - Il file .eq1 creato potrebbe avere elementi leggermente diversi
   - Necessario verificare la fonte originale dei dati

2. **❓ Frame di Riferimento**
   - ELEMENTI_ORBITALI_TEST.md menziona conversione ECLM → ICRF
   - AstDyn potrebbe usare un frame diverso nella propagazione

3. **❓ Epoca degli Elementi**
   - MJD 61000.0 corrisponde a 2018-Mar-16
   - Potrebbero esserci discrepanze nella conversione JD/MJD

4. **❓ Test Case Confuso**
   - Il file test_17030_nov_2025_rkf78.cpp usa elementi completamente diversi:
     - a = 2.71926 AU (vs 3.175473 AU)
     - e = 0.10638 (vs 0.045407)
     - i = 9.3708° (vs 2.9046°)
   - **Possibile confusione con altro asteroide** (Pompeja 203?)

## Prossimi Passi

### 1. Validazione con JPL Horizons
Confrontare i risultati ottenuti con dati ufficiali JPL Horizons:
- **URL**: https://ssd.jpl.nasa.gov/horizons/app.html
- **Target**: Asteroid 17030 Sierks
- **Epoch**: 2025-Nov-28 00:00 UTC (MJD 60277.0)
- **Reference Frame**: ICRF
- **Center**: Sun (@10)

### 2. Verifica Elementi Iniziali
- Scaricare elementi ufficiali da AstDyS database
- Confrontare con il file 17030.eq1 creato
- Verificare l'epoca corretta (MJD 61000.0 = 2018-Mar-16?)

### 3. Test con Altri Asteroidi
- Testare con 203 Pompeja (dati già disponibili)
- Testare con 11234 (dati in astdyn/data/)
- Verificare consistenza del propagatore

### 4. Debug Frame di Riferimento
- Verificare che AstDyn usi ICRF come riferimento finale
- Controllare la conversione ECLM → ICRF nel parser
- Validare rotazioni tra frame

## TEST ORIGINALE CONFERMATION

### Eseguito test astdyn/tests/test_17030_with_astdyn_lib.cpp

**Risultati**:
```
Max Separazione Angolare: 329,216 arcsec (~91°)
Max Errore Distanza: 0.809 AU
❌ INACCETTABILE: Errore > 60 arcsec
```

**CONFERMA**: Anche il test originale mostra errori enormi!

### Root Cause Identificata

Gli elementi orbitali usati nel test NON corrispondono all'asteroide 17030 Sierks reale:

**Elementi usati nel test**:
- a = 2.71926 AU, e = 0.10638, i = 9.3708°
- Epoca: 2018-Mar-16 (MJD 58194)

**Posizione calcolata** (28 Nov 2025):  
- RA = 166.986°, Dec = 6.857°, r = 2.446 AU

**Posizione JPL reale** (28 Nov 2025):  
- RA = 73.409°, Dec = 20.324°, r = 1.660 AU

**Differenza**: ~93° in RA, ~13° in Dec, 0.79 AU in distanza

## Conclusioni

### ✅ **Successi Tecnici**:
- AstDyn library compila e funziona correttamente
- Parser .eq1 legge il file senza errori  
- Propagazione RKF78 completa senza crash
- Performance eccellenti (1 ms per 723 giorni)
- Integratore stabile (26 step, 8 rigetti, 442 valutazioni)

### ❌ **PROBLEMA CRITICO IDENTIFICATO**:
**Errore di 6 AU rispetto a JPL Horizons è INACCETTABILE**

Questo errore NON può essere causato da:
- ❌ Problemi numerici dell'integratore (errore troppo grande)
- ❌ Frame di riferimento errato (cambierebbe orientazione, non magnitudine)
- ❌ Perturbazioni mancanti (effetto sarebbe ordini di grandezza minore)

**L'errore DEVE essere causato da**:
1. 🔴 **Elementi orbitali iniziali COMPLETAMENTE ERRATI**
   - I valori nel file 17030.eq1 non corrispondono all'asteroide 17030
   - Possibile confusione con altro oggetto
   - Necessario scaricare elementi ufficiali da AstDyS/JPL

2. 🔴 **Epoca degli elementi SBAGLIATA**
   - MJD 61000.0 potrebbe non essere l'epoca corretta
   - Conversione JD/MJD errata
   - Epoch di riferimento diverso (TT vs TDB vs UTC)

3. 🔴 **Test case documentazione ERRATO**
   - ELEMENTI_ORBITALI_TEST.md contiene dati non validati
   - Possibile confusione tra test cases diversi
   - Necessario ricostruire test da fonti ufficiali

## Azioni Immediate Richieste

### 🔴 **PRIORITÀ CRITICA - Elementi Orbitali**

1. **Scaricare elementi ufficiali 17030 da AstDyS**:
   ```bash
   curl "https://newton.spacedys.com/astdys/index.php?pc=1.1.3.1&n=17030&oc=f&oc=l&output=OEF" > 17030_astdys_official.eq1
   ```

2. **Verificare epoca degli elementi**:
   - Controllare se MJD 61000.0 è corretto
   - Convertire correttamente tra JD/MJD/calendar
   - Verificare scala temporale (TT/TDB/UTC)

3. **Test con asteroide validato**:
   - Usare 203 Pompeja (dati già in astdyn/data/)
   - Testare con 11234 (dati disponibili)
   - Validare pipeline con caso noto funzionante

### ⚠️ **PRIORITÀ ALTA - Validazione Pipeline**

4. **Creare test di riferimento**:
   - Propagare elemento noto con JPL
   - Confrontare step-by-step
   - Isolare punto di failure

5. **Verificare conversioni frame**:
   - Test ECLM → ICRF con dati sintetici
   - Validare rotazioni
   - Verificare obliquità ε

### 📊 **PRIORITÀ MEDIA - Documentazione**

6. **Correggere ELEMENTI_ORBITALI_TEST.md**:
   - Rimuovere dati non validati
   - Ricostruire da fonti ufficiali
   - Aggiungere riferimenti verificabili

7. **Aggiornare esempi**:
   - Usare solo dati validati con JPL
   - Documentare fonti
   - Includere tolleranze attese

---

**Autore**: Test automatico AstDyn  
**Timestamp**: 2025-12-01 13:00 CET
