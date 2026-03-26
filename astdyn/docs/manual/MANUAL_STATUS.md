# AstDyn Scientific Manual - Status Report

**Last Update**: 25 March 2026
**Current State**: 220 pages, 19/26 chapters complete (73%)
**Current PDF**: `docs/build/manual/main.pdf` (910 KB)

## 🎯 Progress Overview

- ✅ **Part I: Theoretical Foundations** (Chapters 1-7) - COMPLETE
- ✅ **Part II: Numerical Methods** (Chapters 8-11) - COMPLETE (Updated with Uncertainty Propagation)
- ✅ **Part III: Orbit Determination** (Chapters 12-15) - COMPLETE
- ⏳ **Part IV: Software Implementation** (Chapters 16-20) - Next
- ✅ **Part V: Applications** (Chapters 21-25) - COMPLETE
- ⏳ **Chapter 26: Future Developments** - UPDATED (Revised Roadmap)

## 📚 Detailed Chapter Status

Ho creato la struttura completa per il **manuale scientifico della libreria AstDyn** in due lingue:

### Versione Inglese (`en/`)
- ✅ **Frontespizio** completo con figura TikZ
- ✅ **Prefazione** (motivazioni, struttura, pubblico, notazioni)
- ✅ **Capitolo 1: Introduzione** (23 pagine di contenuto)
  - Cos'è la meccanica celeste
  - Problema dei due corpi e N-corpi
  - Panoramica della libreria AstDyn
  - Architettura software con diagrammi
  - Applicazioni pratiche
  - Esempio di codice funzionante
- ✅ **Capitolo 2: Sistemi di Tempo** (contenuto completo)
  - Julian Day e MJD
  - UTC, UT1, TAI, TT, TDB
  - Conversioni e leap seconds
  - Diagrammi e grafici
  - Esempi di implementazione
- ✅ **Capitolo 3: Sistemi di Coordinate** (contenuto completo)
  - Coordinate equatoriali ed eclittiche
  - Trasformazioni tra sistemi
  - Frame di riferimento J2000.0
  - Precessione
  - Diagrammi TikZ illustrativi
- ✅ **Capitolo 4: Sistemi di Riferimento** (~15 pagine)
  - ICRS (International Celestial Reference System)
  - J2000.0 equatoriale (più comune per asteroidi)
  - Sistemi eclittici
  - Trasformazioni con matrici di rotazione
  - Precessione (50.3 arcsec/anno)
  - FK5, invariable plane, body-centric frames
  - 2 diagrammi TikZ (assi J2000.0, precessione)
- ✅ **Capitolo 5: Elementi Orbitali** (~18 pagine)
  - Elementi kepleriani classici (a, e, i, Ω, ω, M)
  - Terza legge di Keplero
  - Singolarità (Ω indefinito a i=0, ω indefinito a e=0)
  - Vettori di stato cartesiani
  - Elementi equinoziali (evitano singolarità, usati da JPL)
  - Elementi di Delaunay (variabili canoniche)
  - Algoritmi di conversione (Kep↔Cart, accuratezza ~10⁻¹⁵ AU)
  - Diagramma TikZ geometria orbitale 3D
- ✅ **Capitolo 6: Problema Due Corpi** (~23 pagine)
  - Leggi di conservazione (energia, momento angolare, LRL)
  - Leggi di Keplero (derivazione completa)
  - Sezioni coniche (ellisse, parabola, iperbole)
  - Equazione di Keplero (M = E - e sin E)
  - Risoluzione Newton-Raphson (3-5 iterazioni)
  - Relazioni tra anomalie (vera, eccentrica, media)
  - Equazione vis-viva
  - Orbite paraboliche e iperboliche
  - Coefficienti di Lagrange per propagazione
  - 3 diagrammi TikZ (coniche, anomalie, geometria)
- ✅ **Capitolo 7: Perturbazioni** (~15 pagine) - COMPLETE
- ✅ **Capitolo 8: Integrazione Numerica** (~12 pagine) - COMPLETE
- ✅ **Capitolo 9: Propagazione Orbitale** (~12 pagine) - COMPLETE
- ✅ **Capitolo 10: Matrice Transizione di Stato** (~14 pagine) - COMPLETE
- ✅ **Capitolo 11: Effemeridi** (~12 pagine) - COMPLETE
- ✅ **Capitolo 12: Osservazioni** (~16 pagine) - COMPLETE
- ✅ **Capitolo 13: Determinazione Orbita Iniziale** (~12 pagine) - COMPLETE
- ✅ **Capitolo 14: Correzione Differenziale** (~12 pagine) - COMPLETE
- ✅ **Capitolo 15: Analisi Residui** (~12 pagine) - COMPLETE
- ⏳ **Capitoli 16-26**: Struttura pronta, contenuto da scrivere

### Versione Italiana (`it/`)
- ✅ **Frontespizio** completo (tradotto)
- ✅ **Capitoli 4-15**: Traduzioni complete
  - Cap 4: Sistemi di Riferimento
  - Cap 5: Elementi Orbitali
  - Cap 6: Problema Due Corpi
  - Cap 7: Perturbazioni
  - Cap 8: Integrazione Numerica
  - Cap 9: Propagazione Orbitale
  - Cap 10: Matrice Transizione di Stato
  - Cap 11: Effemeridi
  - Cap 12: Osservazioni
  - Cap 13: Determinazione Orbita Iniziale
  - Cap 14: Correzione Differenziale
  - Cap 15: Analisi Residui
- ⚠️ **Capitoli 1-3**: Stub minimali (necessitano traduzione completa)
- 🔄 **Struttura parallela** alla versione inglese

## 🎨 Caratteristiche Tipografiche

**Come configurato (aggiornato dopo test compilazione):**
- ✅ Font: **Palatino** (mathpazo package - sostituto per URW Classico non disponibile)
- ✅ Dimensione: **12pt**
- ✅ Interlinea: **1.3** (setstretch{1.3})
- ✅ Formato: A4, fronte-retro
- ✅ Margini: 3.5cm sinistro, 2.5cm destro
- ✅ Unicode: Tutti i simboli convertiti a LaTeX (°→^\circ, ≈→~~, ₁→_1)

## 📐 Caratteristiche Grafiche

- **Diagrammi TikZ**: Orbite, sistemi di coordinate, architettura software
- **Grafici PGFPlots**: Leap seconds, funzioni matematiche
- **Syntax highlighting**: Codice C++ colorato
- **Hyperlink**: Collegamenti interni ed esterni colorati
- **Ambienti matematici**: Teoremi, definizioni, esempi

## 🛠️ Sistema di Build

### Makefile Incluso
```bash
cd astdyn/docs/manual/en
make          # Compila versione inglese
make quick    # Test rapido
make clean    # Pulisce file ausiliari
```

### Script di Test
```bash
./compile_test.sh en    # Testa versione inglese
./compile_test.sh it    # Testa versione italiana
```

## 📖 Contenuto Completo (già scritto)

### Capitolo 1: Introduzione (~23 pagine)
1. **Cos'è la Meccanica Celeste**
   - Definizione e ambito
   - Storia (Kepler, Newton)
   - Problema dei due corpi
   - Problema degli N corpi

2. **Panoramica Libreria AstDyn**
   - Filosofia di design
   - Funzionalità principali
   - Architettura (con diagramma a strati)
   - Dipendenze (Eigen, Boost)

3. **Applicazioni**
   - Determinazione orbite asteroidi
   - Analisi traiettorie spaziali
   - Evoluzione orbitale a lungo termine
   - Strumento educativo

4. **Validazione**
   - Confronti con OrbFit
   - JPL Horizons
   - Soluzioni analitiche

5. **Getting Started**
   - Installazione
   - Esempio minimale funzionante
   - Organizzazione capitoli

### Capitolo 2: Sistemi di Tempo (~15 pagine)
- Perché servono sistemi multipli
- Julian Day e MJD
- UT, UTC, UT1
- TAI, TT, TDB
- Leap seconds (con grafico storico)
- Trasformazioni tra scale temporali
- $\Delta T$ e approssimazioni
- Esempi di codice implementativo
- Precisione richiesta

### Capitolo 3: Coordinate (~12 pagine)
- Frame inerziali vs rotanti
- Sistema equatoriale (RA, Dec)
- Sistema eclittico ($\lambda$, $\beta$)
- Diagrammi illustrativi
- Trasformazioni matriciali
- Obliquità
- J2000.0 (epoch vs equinox)
- Precessione
- Esempi pratici

## 📦 Struttura delle Parti

### Parte I: Fondamenti Teorici (Capitoli 1-7)
- ✅ Introduzione
- ✅ Sistemi di tempo
- ✅ Sistemi di coordinate
- ✅ Sistemi di riferimento
- ✅ Elementi orbitali
- ✅ Problema dei due corpi
- 🔄 Perturbazioni (PROSSIMO)

### Parte II: Metodi Numerici (Capitoli 8-11)
- 🔄 Integrazione numerica
- 🔄 Propagazione orbite
- 🔄 Matrice di transizione
- 🔄 Effemeridi

### Parte III: Determinazione Orbitale (Capitoli 12-15)
- 🔄 Osservazioni
- 🔄 Orbita iniziale (Gauss)
- 🔄 Correzione differenziale
- 🔄 Analisi residui

### Parte IV: Implementazione (Capitoli 16-20)
- 🔄 Architettura
- 🔄 Moduli core
- 🔄 **Parser configurabili** (da documentare!)
- 🔄 API reference
- 🔄 Esempi

### Parte V: Validazione (Capitoli 21-23)
- 🔄 **Caso Pompeja** (da documentare!)
- 🔄 Casi studio
- 🔄 Performance

## 🎯 Prossimi Passi

### Priorità Alta
1. **Completare Cap. 7 (Perturbazioni)**: N-body, J2, SRP, relatività (~20 pagine) - PROSSIMO
2. **Completare Cap. 18 (Parser)**: Documentare il sistema configurabile appena creato
3. **Completare Cap. 22 (Pompeja)**: Caso studio completo con grafici dei residui

### Priorità Media
4. **Cap. 8-11**: Metodi numerici (integrazione, propagazione)
5. **Cap. 12-15**: Orbit determination (Gauss, diff. correction)
6. **Cap. 19-20**: API ed esempi

### Priorità Bassa
7. **Traduzione italiana**: Parallela man mano che si completa l'inglese
8. **Figure aggiuntive**: Più diagrammi e grafici
9. **Appendici**: Tabelle di costanti, formule di riferimento

## 📝 Note per l'Espansione

### Template per Nuovi Capitoli
Ogni capitolo dovrebbe includere:
1. **Introduzione** (motivazione, contesto)
2. **Teoria matematica** (derivazioni, equazioni)
3. **Implementazione** (algoritmi, codice)
4. **Esempi** (casi d'uso pratici)
5. **Validazione** (test, confronti)

### Stile di Scrittura
- Livello: Laurea magistrale/PhD
- Rigoroso ma accessibile
- Equazioni numerate
- Riferimenti bibliografici
- Codice commentato

## 🔧 Requisiti di Compilazione

### LaTeX Packages Necessari
- `classico` (font URW Classico)
- `tikz`, `pgfplots` (grafici)
- `listings` (codice)
- `hyperref` (collegamenti)
- `amsmath`, `amssymb` (matematica)

### Test di Compilazione
```bash
cd astdyn/docs/manual/en
pdflatex -interaction=nonstopmode main.tex
# Dovrebbe produrre main.pdf
```

## 📊 Statistiche Attuali (26 novembre 2025)

- **Pagine scritte**: 98 (versione inglese compilata)
- **Dimensione PDF**: 369 KB
- **Capitoli completi**: 6/26 (23%)
  - EN: Capitoli 1-6 completi
  - IT: Capitoli 4-6 tradotti (1-3 in stub)
- **Figure**: 8 diagrammi TikZ
- **Esempi di codice**: ~12
- **Equazioni**: ~80+
- **Compilazione**: ✅ Successo (EN), IT struttura pronta

## 🎓 Riferimenti Bibliografici

Da aggiungere nel Cap. 24:
- Bate, Mueller, White: *Fundamentals of Astrodynamics*
- Curtis: *Orbital Mechanics for Engineering Students*
- Danby: *Fundamentals of Celestial Mechanics*
- OrbFit documentation
- IAU SOFA documentation

---

**Stato**: Struttura completa creata, primi 3 capitoli scritti, pronto per espansione incrementale.

**Tempo stimato per completamento**: 2-3 settimane a tempo pieno (data la complessità tecnica e matematica).
