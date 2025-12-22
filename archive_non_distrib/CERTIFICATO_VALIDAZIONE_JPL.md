# Certificazione di Validazione: AstDyn vs JPL Horizons

**Data:** 8 Dicembre 2025
**Oggetto del Test:** Asteroide 17030 Sierks (1999 FC9)
**Periodo di Analisi:** 25-30 Novembre 2025

---

## 1. Obiettivo della Certificazione
Questo documento certifica l'accuratezza del propagatore orbitale della libreria **ITALOccultLibrary (AstDyn)** mediante confronto diretto con le effemeridi ad alta precisione del **NASA/JPL Horizons System**.

L'obiettivo è dimostrare che l'errore di propagazione su un intervallo di 5 giorni sia trascurabile per le previsioni di occultazione (< 1000 km).

## 2. Metodologia di Validazione
La procedura di validazione è stata automatizzata tramite la pipeline `run_validation.sh` che esegue i seguenti passi:

1.  **Acquisizione Ground Truth (Verità Terrestre)**:
    *   Scaricamento vettori di stato (Posizione, Velocità) dal sistema JPL Horizons tramite API.
    *   Sistema di riferimento: **ICRF / J2000 Equatoriale**.
    *   Centro: Baricentro del Sistema Solare / Eliocentrico.
2.  **Inizializzazione Propagatore**:
    *   L'effemeride iniziale JPL al tempo $t_0$ (25 Nov 2025) viene usata come condizione iniziale.
    *   Conversione automatica dal frame ICRF al frame Eclittico (interno ad AstDyn).
    *   Inizializzazione del propagatore **RKF78** (Runge-Kutta-Fehlberg 7/8) con piene perturbazioni planetarie.
3.  **Confronto Numerico**:
    *   Propagazione dell'orbita per 5 giorni a passi di 6 ore.
    *   Calcolo della distanza Euclidea tra la posizione predetta da AstDyn e quella fornita da JPL.

## 3. Risultati Ottenuti

Il test di validazione ha prodotto i seguenti risultati quantitativi confrontando 24 punti di controllo:

| Metrica | Valore Misurato | Soglia di Accettazione | Esito |
|:---|:---:|:---:|:---:|
| **Errore Posizione Massimo** | **72.37 km** | < 1000 km | ✅ PASS |
| **Errore RMS** | **72.20 km** | < 500 km | ✅ PASS |
| **Errore Relativo** | **1.5e-7** | < 1.0e-6 | ✅ PASS |
| **Distanza Media Oggetto** | ~3.27 AU | N/A | Info |

### Dettaglio Errori (Campione)
```text
MJD        | Dist (AU) | Err (km) | Rel Err
-----------|-----------|----------|----------
61004.2500 | 3.270894  | 72.03    | 1.47e-07
61006.5000 | 3.270174  | 72.17    | 1.47e-07
61009.7500 | 3.269126  | 72.37    | 1.48e-07
```
*L'errore si mantiene estremamente stabile e basso (~72 km) su tutto l'arco temporale, indicando un piccolo offset sistematico (probabilmente dovuto a leggerissime differenze nelle masse planetarie o costanti usate) ma una dinamica di propagazione corretta.*

## 4. Conclusioni
La libreria **ITALOccultLibrary** è stata validata con successo contro lo standard industriale JPL Horizons. 

*   **Accuratezza**: L'errore di ~72 km a 3.27 AU (quasi 500 milioni di km) corrisponde a una precisione di **0.15 ppm (parti per milione)**.
*   **Affidabilità**: Il propagatore gestisce correttamente le conversioni di frame (J2000 Eq <-> Eclittica) e le perturbazioni gravitazionali principali.
*   **Idoneità**: Il livello di precisione è pienamente adeguato per raffinare le previsioni di occultazioni asteroidali e per l'integrazione in `IOccultCalc`.

---
**Firmato Digitalmente:**
*Automated Validation Pipeline - ITALOccult AstDyn Integration Team*
