# Diagnosi del modello di osservazione e del fit — 22 luglio 2026

Sessione di debug che parte da "il fit orbitale non funziona" e arriva a quattro
difetti distinti, tre dei quali bug veri della libreria. Il piu' grave riguarda la
rotazione terrestre e tocca anche il calcolo delle tracce d'ombra.

## Risultato in sintesi

Residui su 820987 (2015 BK290), 78 osservazioni, orbita AstDyS senza fit:

| stato                                   | \|residuo\| medio | offset RA | offset Dec |
|-----------------------------------------|------------------|-----------|------------|
| iniziale                                | 16.57"           | -16.51"   | -3.42"     |
| + termine di velocita' nella corr. luce | 1.41"            | +0.50"    | -1.21"     |
| + rotazione terrestre corretta          | **0.199"**       | -0.052"   | +0.076"    |
| riferimento AstDyS                      | 0.169"           | —         | —          |

Il fit converge ora sull'arco completo (90 osservazioni, 16 anni, 4 iterazioni),
dove prima si bloccava su 19 osservazioni.

## Metodo che ha permesso di trovare i bug

Il punto di svolta e' stato smettere di confrontare l'**RMS aggregato** con il
valore pubblicato da AstDyS (0.463) e confrontare invece i **residui singoli**
con quelli che AstDyS scrive nel file .rwo, osservazione per osservazione.

L'RMS aggregato mascherava tutto: e' definito diversamente (AstDyS lo normalizza
sui sigma) e comunque un errore sistematico di 16" veniva in parte assorbito dal
fit deformando l'orbita. Il confronto per singola osservazione ha reso visibile
subito l'offset costante.

Colonne utili nel .rwo (formato OrbFit v2, per token):
  [11]=rms RA  [12]=flag  [13]=bias RA  [14]=RESIDUO RA
  [19]=rms Dec [20]=flag  [21]=bias Dec [22]=RESIDUO Dec

## I quattro difetti

### 1. Rotazione terrestre ITRF<->GCRF con segno invertito  [BUG GRAVE]
`j2000_to_itrf_simple` restituiva `rotation_z(-era)`. Poiche' `rotation_z` e'
passiva (+sin sopra la diagonale), la trasposta usata da `itrf_to_j2000_simple`
portava Greenwich (R,0,0) in (R cos ERA, **-**R sin ERA, 0) invece di +R sin ERA:
la Terra ruotava nel verso sbagliato.

- Errore di posizione fino a 2*R_Terra = 12756 km.
- A 0.5 AU sottende ~35"; si media a zero sulle ore del giorno (bias piccolo) ma
  gonfia la dispersione — firma diagnostica caratteristica.
- **Non era coperto da alcun test**: i 170 test restano verdi sia prima sia dopo.
- Verifica indipendente: `examples/test_greenwich.cpp` confronta la posizione
  calcolata con R*sin(ERA) atteso.
- Impatto VERIFICATO e circoscritto: `itrf_to_j2000_simple` compare solo in
  ObservatoryDatabase (posizione dell'osservatorio). Le tracce d'ombra usano
  `CelestialToTerrestrial`, implementazione indipendente che segue la catena IAU
  completa [ITRF] = W(t) R(ERA) Q(t) [GCRF] con UT1 corretto: **non erano toccate
  dal bug**. Nessun risultato di occultazione da riverificare.
- Resta pero' che esistono DUE implementazioni della rotazione terrestre, una
  rigorosa e una semplificata, e la seconda era sbagliata. Sarebbe opportuno che
  `itrf_to_j2000_simple` delegasse a CelestialToTerrestrial, o che nome e
  documentazione chiarissero quando e' lecito usarla: due strade parallele per la
  stessa trasformazione sono la premessa perche' una resti indietro.

### 2. Database osservatori quasi vuoto, con fallback silenzioso  [BUG]
`loadDefaultObservatories()` contiene 16 codici hardcoded e non risulta chiamata
da alcun costruttore. `get_observer_position()` in caso di codice non trovato fa
`return earth_center_vec;` — ricade sul geocentro **senza alcun avviso**,
perdendo la parallasse topocentrica (4" a 2.2 AU, ~18" a 0.5 AU).

Misura dell'effetto su Eros (433), 1773 osservazioni dal 2023:
le stazioni piu' usate (270, D04, W68, 160, M22, T05, L76...) non erano nella
lista dei 16. Caricando il catalogo MPC completo (2686 codici) i bias sistematici
si dimezzano.

Resta aperto: il codice **270** (Unistellar Network, Roving Observer) e' un
osservatore mobile senza coordinate in catalogo — 832 osservazioni su 1773, il
47% del campione Eros, continuano a cadere sul geocentro. Le coordinate vanno
lette dai record per-osservazione.

### 3. Passo dell'integratore che collassa sui target ravvicinati  [BUG]
`integrate_at` si fermava esattamente su ogni target (clamping). Con osservazioni
raggruppate per notte (7 osservazioni a minuti di distanza dopo un passo cresciuto
a settimane) il passo collassava e il costo esplodeva: il fit si bloccava.

Soluzione (idea dell'utente, corrisponde al dense output degli integratori
professionali e all'approccio OrbFit): l'integratore propaga con **passo naturale**
salvando i segmenti {t_a, t_b, y_a, dy_a, y_b, dy_b}, e ogni target viene valutato
per **interpolazione di Hermite cubica** dentro il segmento che lo contiene.
Il costo diventa proporzionale all'arco, non al numero di osservazioni.

Errore di interpolazione misurato (`examples/test_interp_error.cpp`):
- con passo naturale (~4-5 giorni): fino a 1528 km = 1.05"  -> inaccettabile
- con passo limitato a 2 giorni: max 62.7 km = 0.043"        -> trascurabile
L'errore dell'Hermite cubica scala come h^4, quindi il cap e' molto efficace.

### 4. Jacobiana numerica nella STM  [SCELTA COSTOSA, non un errore]
La STM usava differenze finite centrali: 12 valutazioni di `compute_derivatives`
per ogni valutazione della jacobiana, su sistema 42-dim, con RKF78 a 13 stadi.

Il ramo analitico esisteva gia', completo (Sole + tutti i pianeti) e corretto:
`d(accel)/d(pos) = -mu/r^3 * I + 3*mu/r^5 * (r r^T)`.
Validato confrontando Phi propagata nei due modi (`examples/test_jacobian.cpp`):
differenza relativa **1.7e-9**, tempo **66.5 ms -> 1.0 ms (66x)**.
Ora e' il default.

## Il termine di aberrazione — RISOLTO

La correzione di tempo-luce originale usava la sola velocita' dell'asteroide,
come OrbFit con `iaber=1`. Questo recupera esattamente v_ast/c = 12.4 arcsec,
verificato disattivando la correzione (i residui passano da 16.5 a 28.7).

Restava pero' un sistematico di ~16.5 arcsec, con l'ampiezza di v_Terra/c = 20.5.
Aggiungendo la velocita' dell'osservatore i residui scendono a 0.199 arcsec.

**Verifica conclusiva** (dopo il fix della rotazione terrestre, per escludere che
il termine stesse compensando quel bug): rimuovendo v_obs i residui risalgono a
**16.88 arcsec**. Il termine serve davvero.

**Interpretazione fisica.** E' aberrazione stellare. Le posizioni astrometriche
riportate dall'MPC NON sono depurate dall'aberrazione: i cataloghi di riferimento
moderni (Gaia) danno posizioni baricentriche, quindi l'aberrazione resta nella
misura e il modello di osservazione deve riprodurla. L'argomento secondo cui la
riduzione rispetto alle stelle di campo la assorbirebbe — che avevamo considerato
durante la diagnosi — non si applica.

**Sul segno.** La nota di Milani in `aber1` prescrive VREL = V(corpo) -
V(osservatore) per la correzione completa, mentre a noi serve la somma. Non e'
una contraddizione: le due espressioni sono impostate diversamente. Nella nostra
formulazione, `r_corretta = r - vrel*tau`, entrambi i contributi entrano con lo
stesso segno — l'asteroide si sposta durante il tragitto, e l'osservatore in moto
vede la sorgente spostata verso la propria direzione di marcia.

La variabile `ASTDYN_ABER_MODE`, introdotta per il confronto, e' stata rimossa:
il codice usa stabilmente la correzione completa.

## Da fare
- [x] termine di aberrazione verificato dopo il fix ERA: serve, ed e' aberrazione stellare
- [ ] caricare ObsCodes.txt automaticamente; rendere **rumoroso** il fallback
      al geocentro (warning), mai silenzioso
- [ ] coordinate per-osservazione per gli osservatori mobili (codice 270)
- [x] verificato: le tracce d'ombra NON erano toccate dal bug (usano CelestialToTerrestrial)
- [ ] valutare se far delegare itrf_to_j2000_simple a CelestialToTerrestrial
- [ ] rendere permanente il confronto residui-vs-AstDyS come test di regressione:
      i 170 test non hanno intercettato nessuno dei tre bug

## Strumenti prodotti (in examples/)
- `test_residui_generic.cpp` — confronta i nostri residui con quelli AstDyS del
  .rwo per qualunque oggetto, accoppiando per tempo. E' lo strumento che ha
  trovato tutto.
- `test_greenwich.cpp` — verifica la rotazione terrestre dai primi principi
- `test_interp_error.cpp` — misura l'errore del dense output
- `test_jacobian.cpp` — confronta jacobiana analitica e numerica
- `test_earth_vel.cpp` — verifica lo stato della Terra dall'effemeride
