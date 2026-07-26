# Piano di intervento — Fase 5, nodi aperti

Ipotesi guida: **diversi sintomi apparentemente distinti sono lo stesso
problema visto da angolazioni diverse.** Verificato in apertura: ne'
ioccultcalc.cpp ne' test_fit_bk290.cpp caricano il catalogo osservatori, quindi
il fit colloca gli osservatori nel geocentro perdendo la parallasse topocentrica.

I sintomi da spiegare:
- residui 0.199" partendo dall'orbita AstDyS, ma RMS 1.59" **dopo** il fit —
  un fit non puo' peggiorare cio' che minimizza, quindi le due misure non
  misurano la stessa cosa;
- chi quadro ridotto ~10 (residui tre volte i sigma dichiarati);
- ellisse del secondo stadio identica a sei cifre a quella del primo.

---

## Intervento 1 — Caricare il catalogo osservatori  [PRIORITA' ALTA, costo basso]

**Diagnosi.** `test_residui_generic` carica ObsCodes.txt (2686 codici) e ottiene
0.199"; `test_fit_bk290` non lo carica e ottiene 1.59". E' l'unica differenza
sostanziale fra i due percorsi. Senza catalogo, `get_observer_position` ricade
sul geocentro **in silenzio** (`return earth_center_vec`), perdendo la parallasse
topocentrica: ~4" a 2.2 AU, molto di piu' sui singoli casi.

**Da fare.**
1. `ioccultcalc` deve caricare il catalogo all'avvio, da una chiave di
   configurazione (`observatories.file`) con ricerca in posizioni note come
   ricaduta (~/.ioccultcalc/observatories/ObsCodes.txt).
2. Il fallback al geocentro deve **emettere un avviso**, non restare silenzioso.
   E' la lezione ricorrente di questa giornata.
3. Stessa cosa in test_fit_bk290, per poter confrontare i due percorsi.

**Verifica.** L'RMS del fit su BK290 deve scendere da 1.59" verso ~0.2",
allineandosi al confronto sui residui singoli. Se accade, i tre sintomi elencati
sopra si spiegano insieme e il chi quadro torna vicino a 1.

**Se NON accade**, l'ipotesi e' sbagliata e si passa all'intervento 2 — ma
sapendo qualcosa in piu'.

---

## Intervento 2 — Perturbatori asteroidali  [se serve dopo l'1]

**Diagnosi.** `include_asteroids` e' true di default ma `asteroid_ephemeris_file`
e' vuoto: i perturbatori non vengono applicati benche' sb441-n16.bsp (i 16 corpi
massivi di JPL) sia gia' sul disco. AstDyS li include sempre. Su un arco di 16
anni l'effetto e' dell'ordine dell'arcsec.

**Da fare.** Impostare il file nella configurazione di campagna e nel profilo
`full`, poi misurare l'effetto sui residui **prima** di trarre conclusioni: nella
prova del 2026-07-22 l'attivazione non aveva cambiato nulla, segno che il file
non veniva davvero caricato. Va verificato che l'effemeride sia aperta,
non solo che il parametro sia impostato.

**Verifica.** Confronto dei residui con e senza perturbatori sullo stesso caso.
Una differenza nulla a quattro cifre significa che non vengono applicati.

---

## Intervento 3 — Debiasing dei cataloghi stellari  [valore atteso 0.2-0.3"]

**Diagnosi.** Il .rwo contiene le correzioni di bias per osservazione (colonne
bias RA/Dec, tipicamente -0.199 / +0.278 arcsec) che AstDyS applica e noi no.

**Da fare.** Leggerle nel RWOReader e sottrarle ai residui.

**Verifica.** I bias medi residui su BK290 devono avvicinarsi a zero.

---

## Intervento 4 — Covarianza dal fit  [dipende da 1-3]

Solo dopo che il chi quadro e' vicino a 1 ha senso collegare `(A^T W A)^-1`
alle ellissi. Con chi quadro ~10 la covarianza formale sarebbe ottimistica di
un fattore 3, e in una predizione di occultazione sottostimare l'ellisse
significa mancare l'evento.

Se dopo gli interventi 1-3 il chi quadro resta lontano da 1, applicare il
fattore di gonfiamento sqrt(chi2/dof) — pratica standard — documentandolo.

**Verifica.** Ellisse dal fit confrontata con quella AstDyS: accordo entro
qualche percento conferma l'intera catena.

---

## Intervento 5 — Manuale del configuratore  [indipendente]

Documentare le chiavi introdotte: `objects.observations_dir`,
`diffcorr.enabled`, `diffcorr.obs_years`, `diffcorr.tolerance`,
`second_stage.fit`, e la nuova `observatories.file`.

Segnalare che `diffcorr.tolerance` e' misurato come **inefficace** (quattro
ordini di grandezza senza effetto sul costo): meglio dirlo che lasciar credere
il contrario.

---

## Ordine e razionale

1 per primo: costo basso, e se l'ipotesi regge spiega tutti i sintomi insieme.
2 e 3 solo se dopo l'1 resta uno scarto da spiegare.
4 quando i residui sono sotto controllo — non prima, perche' la covarianza
erediterebbe ogni errore residuo.
5 in qualsiasi momento, e' indipendente.

## Criterio trasversale

Ogni intervento si chiude con un confronto contro un riferimento **esterno**
(AstDyS, JPL), non con un controllo interno di coerenza. I quattro bug del
2026-07-22 erano tutti invisibili ai 177 test e sono emersi solo cosi'.
