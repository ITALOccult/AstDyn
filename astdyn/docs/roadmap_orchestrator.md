# Roadmap — Orchestratore di campagne per ITALOccult

## Visione

Un orchestratore Python sta *davanti* a `ioccultcalc` e prepara tutto: espande
la lista oggetti, procura gli elementi orbitali (da fonti intercambiabili), li
depone in una cartella convenzionale, invoca il motore, raccoglie i risultati.

`ioccultcalc` resta piccolo nel ruolo: legge elementi da una cartella e scrive
risultati. Non sa nulla di `allnum.db`, di NEODyS, di query SQL, di "tutti gli
asteroidi". La complessità di selezione e provenienza vive in Python, dove è
ispezionabile e facile da tenere in testa.


## Stato (aggiornato 2026-07-22)

- [x] **Fase 0 — Elementi da file nel motore** (C++). Gancio `objects.elements_dir`
      in ioccultcalc: legge `<id>.eq1` locale via la catena tipizzata
      (read_eq1 -> equinoctial_to_keplerian), fallback a Horizons. Corretti 3 bug
      reali lungo il percorso (doppia conversione a*AU; lambda gradi->rad; memoria
      DE441). read_eq1 reso tollerante agli .eq1 senza covarianza. Validato su
      BK290 (ellisse esatta 0.117721" x 0.0644153").
- [x] **Fase 1 — Orchestratore minimo** (Python, tools/orchestrator.py). Config
      UNICA (astdyn_base obbligatoria); fonti db|user|jpl-horizons|astdys-neodys;
      KEP->EQU validata contro l'oracolo BK290 (.eq1 dal DB bit-identico ad AstDyS);
      guardia frame (solo ECLIPTIC_J2000); output XML Occult4 + mappe in results/.
- [x] **Fase 2 — `*` con sistema di filtri SQL** (Python). objects.asteroids='*'
      seleziona dal DB con un blocco `filters`: range [min,max] su qualsiasi
      colonna in whitelist (a,e,i,H,G,...) piu' `diameter_km` derivato (H via
      albedo). Filtri combinabili; guardie ('*' senza filtri rifiutato, avviso+
      conferma sopra soglia). Nota: l'albedo assunto (0.14) sovrastima il diametro
      dei corpi scuri (TNO) -> per la fascia principale combinare con a:[2.0,3.5].
- [x] **Fase 3+4 — Ricerca a due stadi** (Python + C++). L'intento reale che
      unifica le fasi: screening veloce (force field leggero) su tutti i corpi,
      poi raffinamento completo (force field pesante) SOLO sui corpi che hanno
      dato occultazione. Profili `physics_profiles` liberi (light/full o altri),
      applicati via `first_pass`/`second_stage`; default con WARNING rosso se non
      esplicitati; retrocompatibile (senza i blocchi -> passata singola). Output
      tracciabile: screening in results/pass1/, raffinamento in results/.
      Validato: 12 giganti -> screening -> 2 positivi -> raffinamento solo su 2.
      Supporto C++: nuovo OccultationJSONIO (--json-output) per il formato interno
      ricco; i positivi si estraggono da object.number nel JSON (robusto).
      Nota: 'fonte online astdys/neodys' (scope originale Fase 3) resta da fare
      (il container non raggiunge newton.spacedys.com; da implementare/testare sul
      sistema con rete). Nota scientifica: lo screening leggero rischia falsi
      negativi -> renderlo generoso; la differenza light/full pesa su archi lunghi
      o corpi perturbati, mentre su archi brevi il vantaggio e' soprattutto velocita'.

- [x] **Fase 3 (online) - Sorgente AstDyS + ellissi per-asteroide in campagna**.
      source: astdys scarica gli .eq1 individuali da newton.spacedys.com
      (raggruppamento numero//1000), CON covarianza, con cache locale e download
      di cortesia (User-Agent, throttle). La covarianza e' un PARAMETRO DI CALCOLO
      (physics.covariance: true): il gancio elements_dir cattura orb.covariance, la
      trasforma (jacobiana) e la salva per-asteroide (stored_covariances[id]); il
      loop applica a ciascun evento la covarianza del suo corpo -> ellissi reali in
      campagna, in una passata. Osservazioni .rwo scaricate on-demand per i positivi
      (fit, Fase 5). Validato: astdys su BK290 -> ellisse 0.117721 x 0.0644153 @
      71.39 (identica ad AstDyS), vs source:db che da' 0. Modello: covarianza/fit
      sono parametri di calcolo, on nel raffinamento, off nello screening.

- [ ] Fase 5 — Fit orbita on-demand (prima riparare il fit in libreria)
- [ ] Fase 6 — Asteroidi multipli (satelliti/binari)
- [ ] Fase 7 — Installer e distribuzione (catalogo incluso)
- [ ] **Manuale utente finale** (MD organico dell'intero sistema) — ultimo punto.

## Nota trasversale: tipizzare gli scalari

Ogni bug scovato nelle sessioni Fase 0/1 nasceva da uno scalare non tipizzato a
un confine (a in AU vs km, lambda in gradi vs rad). La libreria e' tipizzata
sugli stati e sui frame (KeplerianStateTyped<Frame>), ma scende a `double` /
`Eigen::Vector3d` nudi nelle conversioni, ed e' li' che i bug di unita' si
annidano silenziosamente. Direzione futura: tipizzare anche distanze e angoli
(Distance, Angle) nelle conversioni, cosi' un mismatch di unita' non compili
invece di produrre numeri sbagliati a runtime.

## Struttura di lavoro (una campagna = una cartella)

```
base/
  config.yaml            # la config della campagna (usata da tool e motore)
  elements/              # elementi orbitali preparati dal tool (un file per corpo)
  obs/                   # osservazioni per il fit (solo se richiesto)
  results/
    xml/                 # predizioni Occult4/OWC
    map/                 # mappe
    log/                 # log per corpo / di campagna
```

Il tool prepara `elements/` (e `obs/` se serve); il motore legge di lì e scrive
in `results/`. La cartella è un artefatto riproducibile e ispezionabile.

## Il vincolo tecnico che orienta tutto

`ioccultcalc` oggi prende gli elementi SOLO da JPL Horizons (online, un corpo
per volta): `horizons.query_elements(id, epoch)`. Il `.eq1` che gia' legge
(`CovarianceIO::read_eq1`) porta elementi + covarianza, ma il flusso lo usa
solo per la covarianza, a valle. Quindi:

- Il **lettore di elementi da file esiste gia'** (`read_eq1`).
- Manca il **canale** che lo usi al posto di Horizons per la propagazione.

Aprire quel canale e' un intervento C++ piccolo e localizzato (non tocca la
mole del file): "se esiste un file elementi locale per il corpo, usalo; altrimenti
Horizons". E' la prima pietra: senza, il tool Python non ha come passare gli
elementi al motore.

---

## Fasi

### Fase 0 — Fondamenta: elementi da file nel motore  [C++, piccolo]
Obiettivo: `ioccultcalc` sa caricare elementi da una cartella locale.
- Nuova opzione (config/CLI): `elements_dir` (o `objects.elements_dir`).
- Nel loop di caricamento: se esiste `elements_dir/<id>.eq1` (o formato scelto),
  usa `read_eq1` per elementi (e covarianza) invece di `query_elements`.
- Fallback a Horizons se il file non c'e' (comportamento attuale invariato).
- Validazione: BK290 con elemento locale da' la stessa ellisse di oggi.
Dipendenze: nessuna. E' il prerequisito di tutto il resto.

### Fase 1 — Orchestratore minimo (lista esplicita)  [Python]
Obiettivo: il tool legge la config, prepara `elements/` per una lista di corpi
data (es. "1,4,704"), invoca il motore, raccoglie i risultati.
- Legge la config (la stessa che il motore usa).
- Per ogni corpo: estrae gli elementi dal nostro DB e scrive `elements/<id>.eq1`.
- Invoca `ioccultcalc` con `elements_dir` e i parametri di config.
- Struttura le cartelle `results/`.
Dipendenze: Fase 0.

### Fase 2 — Espansione "*" con filtro diametro  [Python]
Obiettivo: `objects: "*"` diventa "tutti i numerati nella fascia diametro".
- Query su `allnum.db`: fascia diametro -> fascia H (formula, albedo config,
  con margine) -> lista dei numeri. Selezione ispezionabile (scrivibile come
  lista prima di lanciare).
- Guardia: `*` senza filtro diametro (o senza tetto) e' rifiutato, per non
  enumerare 1.5M corpi per errore.
Dipendenze: Fase 1.

### Fase 3 — Fonti intercambiabili degli elementi  [Python]
Obiettivo: gli elementi in `elements/` possono venire dal nostro DB, o essere
scaricati da AstDyS/NEODyS, o forniti dall'utente.
- Astrazione "source": db locale | astdys online | neodys online | file utente.
- Config sceglie la fonte; il tool prepara `elements/` di conseguenza.
- Il motore non cambia: legge sempre da `elements/`.
Dipendenze: Fase 1 (idealmente dopo Fase 2).

### Fase 4 — Multipassata / second-stage  [Python + forse C++]
Obiettivo: rielaborare solo gli eventi positivi con parametri diversi.
- Config: blocco `second_stage` con `only_positive: true` e i propri parametri.
- Prima passata: predizione su tutti i corpi selezionati.
- Seconda passata: solo i corpi che hanno prodotto un'occultazione, ri-processati
  con la config del secondo stage (es. finestra piu' fine, filtri diversi).
- Il fit orbita (vedi Fase 5) si innesta qui, solo se richiesto.
Dipendenze: Fase 1.

### Fase 5 — Fit orbita on-demand  [C++ da sistemare + Python]
Obiettivo: fit orbitale esplicito, attivato solo quando richiesto, tipicamente
nel second-stage sui pochi eventi positivi.
- Il fit nella libreria e' attualmente buggato: prima va sistemato (lavoro C++
  a se', con validazione su un caso noto).
- Attivazione esplicita in config (`fit: {enabled: true, ...}`), default off.
- Il tool prepara `obs/` (osservazioni da AstDyS/MPC) per i corpi da fittare.
- Gira solo sui positivi (poche unita'), non sull'intera selezione.
Dipendenze: Fase 4 per l'innesto; il debug del fit e' indipendente e puo'
iniziare quando si vuole.

### Fase 6 — Asteroidi multipli (satelliti/binari)  [C++ + Python]
Obiettivo: predizione di occultazioni per corpi con satelliti.
- Il motore ha gia' un percorso "system-ids"/BSP per corpi multipli: va
  verificato ed esteso.
- Il tool prepara gli elementi del sistema (primario + satelliti).
Dipendenze: Fase 1; ortogonale alle altre, si puo' collocare quando serve.

### Fase 7 — Installer e distribuzione  [Python]
Obiettivo: un tool che porta da zero a sistema funzionante, catalogo incluso.
- Passi: costruzione/So download del catalogo -> DB iniziale (via
  update_asteroids.py) -> build dell'eseguibile -> verifica.
- Distribuisce anche il catalogo (troppo grande per GitHub) da una sorgente
  definita.
Dipendenze: le fasi che definiscono la struttura finale (idealmente ultima,
ma l'ossatura si puo' preparare prima).

---

## Ordine consigliato

0 -> 1 -> 2  danno subito valore (campagne su selezione o su "*" filtrato, con
elementi dal nostro DB, senza Horizons come collo di bottiglia).

Poi 3 (fonti) e 4 (multipassata) a seconda della priorita'. 5 (fit) quando il
bug del fit e' risolto. 6 (multipli) e 7 (installer) quando l'ossatura e' stabile.

## Principio guida

Ogni fase e' un pezzo che sta in testa, con un confine netto. Il C++ si tocca
il meno possibile e solo dove serve un canale (Fase 0, forse 4/5/6). Tutta la
logica di selezione, provenienza e orchestrazione vive in Python, ispezionabile.
