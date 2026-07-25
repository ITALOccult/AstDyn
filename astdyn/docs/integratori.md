# Gli integratori di AstDyn

Capitolo di manuale. Descrive i metodi implementati, cosa li distingue, e in quali
casi conviene usarli. I valori misurati provengono dal banco di prova
`examples/test_validate_integrators.cpp`, che confronta ogni metodo con l'oracolo
JPL Horizons su archi crescenti.

> Misure del 2026-07-22 su 820987 (2015 BK290), banco di prova
> `examples/test_integrators_bench.cpp`. Manca RADAU, ancora da caratterizzare.

## Come leggere le misure

Ogni integratore viene caratterizzato da tre grandezze:

- **accuratezza** — scarto dalla posizione JPL dopo la propagazione, in AU.
  Riferimento utile: 1e-6 AU = 150 km, che a 2 AU sottendono ~0.1 arcsec.
- **costo** — tempo di calcolo per arco, a parita' di modello di forze.
- **stabilita' sull'arco lungo** — come cresce l'errore allungando la
  propagazione. Un metodo puo' essere accurato su un anno e divergere su dieci.

Le prove usano lo stesso oggetto, la stessa epoca e lo stesso modello di forze,
cosi' che le differenze siano attribuibili solo al metodo.

## Metodi implementati

### RK4 — Runge-Kutta classico, passo fisso
Il metodo di riferimento didattico: quarto ordine, quattro valutazioni per passo,
nessun controllo dell'errore.

Il passo fisso e' il suo limite fondamentale in astrodinamica: il passo va scelto
per il tratto piu' impegnativo dell'orbita (il perielio) e resta inutilmente
piccolo altrove. Su orbite eccentriche il costo diventa proibitivo, e non esiste
alcun segnale quando la precisione degrada.

*Uso consigliato:* verifiche rapide, confronti didattici, archi brevi con orbite
quasi circolari. **Sconsigliato come default** — in AstDyn lo e' ancora, ed e' una
trappola documentata: propagazioni lunghe risultano lentissime senza spiegazione.

### RKF78 — Runge-Kutta-Fehlberg 7(8), passo adattivo
Coppia di formule di ordine 7 e 8 valutate con gli stessi tredici stadi: la loro
differenza stima l'errore locale e guida il passo. E' il cavallo di battaglia
della libreria.

L'adattivita' e' decisiva: il passo si allarga nei tratti lisci e si stringe dove
serve, senza intervento dell'utente. L'ordine elevato lo rende molto efficiente
per la precisione richiesta dall'astrometria.

*Misurato (BK290, tolleranza 1e-12):* errore ~2e-6 AU su tutti gli archi provati
fino a 10 anni; 28 ms per 10 anni. Miglior rapporto precisione/costo fra tutti.

*Uso consigliato:* **default per propagazione, fit orbitale e STM.** E' il metodo
rispetto al quale vengono validati gli altri.

### AAS — simplettico adattivo (integratore proprietario AstDyn)
Metodo simplettico con controllo del passo basato sull'Hessiana, sviluppato per
la libreria e oggetto di un articolo dedicato.

I metodi simplettici conservano la struttura geometrica del problema hamiltoniano:
l'errore sull'energia resta limitato invece di accumularsi secolarmente, proprieta'
preziosa su integrazioni molto lunghe. Il passo adattivo, di norma incompatibile
con la simpletticita', e' qui recuperato tramite il controllo hessiano.

*Misurato (BK290, tolleranza 1e-12):* accuratezza pari a RKF78 (~2e-6 AU) su tutti
gli archi, ma **circa 40 volte piu' lento** (1190 ms contro 28 ms su 10 anni).

*Uso consigliato:* studi di stabilita' su tempi lunghi, dove conta il comportamento
qualitativo dell'energia piu' della singola posizione. Non conviene per il fit, dove
gli archi sono di decenni e la precisione richiesta e' gia' raggiunta da RKF78.

### SABA4 — NON SUPPORTATO
L'implementazione presente e' difettosa: errore di **1.60 AU gia' su un arco di
0.1 giorni**, dove gli altri metodi danno ~1e-15, e 79 secondi per arco. Dal
2026-07-22 la costruzione fallisce con un messaggio esplicito invece di
restituire numeri sbagliati.

Diagnosi — tre difetti concorrenti:

1. **Il drift e' rettilineo** (`q += h*p`) invece che kepleriano esatto. I metodi
   SABA propriamente detti (Laskar & Robutel) separano l'hamiltoniana in
   kepleriano + perturbazione e risolvono *esattamente* la parte kepleriana a
   ogni sottopasso: e' da li' che viene la loro efficienza. Con il drift
   rettilineo lo schema e' un semplice simplettico alla Yoshida, che su un'orbita
   curva richiede passi minuscoli — e viene costruito con passo iniziale di mezza
   giornata, in cui BK290 percorre ~780.000 km.
2. **La stima d'errore divide per il passo**:
   `err = |y_saba4 - y_saba2| / h`. Accorciando il passo l'errore apparente
   *cresce*, innescando un meccanismo auto-alimentato. La normalizzazione
   corretta e' rispetto alla grandezza dello stato.
3. **Il passo resta bloccato al minimo**: dopo l'accettazione forzata a `h_min_`
   (1e-6 giorni), `adapt_step_size` con errore enorme satura a scala 0.2 e
   riporta il passo a `h_min_`. Per coprire 0.1 giorni servirebbero 100.000
   passi; il ciclo esaurisce il tetto dei 5.000.000 e restituisce uno stato
   quasi immobile — da cui l'errore di 1.6 AU, che e' semplicemente la distanza
   non percorsa.

*Per riabilitarlo* servirebbe riscrivere il drift come propagazione kepleriana
esatta (esiste gia' `kepler_propagator.hpp`), applicare nel kick la sola
perturbazione e non la forza totale, e rifare la stima d'errore. Il codice resta
nel repository. In alternativa, chiamarlo con il suo nome — simplettico di
Yoshida — e correggere solo il controllo del passo, accettandone la lentezza.

### Gauss-Jackson — multistep
Predittore-correttore multistep, storicamente il metodo di elezione per orbite
quasi circolari, ancora molto usato nella meccanica orbitale operativa. Riutilizza
le valutazioni dei passi precedenti, quindi e' economico per passo, ma richiede
un avviamento con un metodo a passo singolo e mal sopporta i cambi di passo.

*Misurato:* accurato (3.7e-8 AU sul decennio) ma lento: 248.8 ms, oltre dieci
volte RKF78.

*Uso consigliato:* nessun vantaggio rispetto a RKF78 nel nostro impiego.

### Radau-IIA — implicito (Everhart RA15)
Metodo di collocazione implicito, A-stabile, adatto ai problemi rigidi e agli
incontri ravvicinati. E' il metodo che OrbFit usa come opzione `imet=3`.

Il carattere implicito richiede la soluzione iterativa di un sistema a ogni passo:
molto costoso per passo, ma con passi molto piu' lunghi e ottima stabilita'.

*Misurato:* troppo lento per uso pratico nella nostra implementazione (prova
interrotta). Da profilare: puo' trattarsi di un problema di implementazione
piu' che del metodo.

*Uso consigliato:* da riconsiderare dopo l'ottimizzazione. Naturale candidato per
gli incontri planetari ravvicinati.

### GRKN64 — Runge-Kutta-Nystrom generalizzato 6(4)
I metodi di Nystrom sfruttano il fatto che l'equazione del moto e' del
second'ordine con la derivata seconda dipendente dalla sola posizione, evitando di
propagare la velocita' come variabile indipendente: meno valutazioni a parita' di
ordine.

*Misurato:* il piu' efficiente sugli archi lunghi — 19.1 ms e 1.9e-7 AU sul
decennio, meglio di RKF78 su entrambi i fronti. Da verificare il trattamento
della correzione relativistica (dipendente dalla velocita') prima di promuoverlo.

*Uso consigliato:* candidato a diventare il default per la propagazione, dopo la
verifica di cui sopra.

## Tabella riassuntiva

Errore rispetto a RKF78 @ 1e-13, propagazione all'indietro dall'epoca.
Tempi su MacBook Air, modello di forze con perturbazioni planetarie.

| metodo | 10 giorni | 100 g   | 1000 g  | 3650 g  | tempo 3650 g | andata-ritorno | energia   |
|--------|-----------|---------|---------|---------|--------------|----------------|-----------|
| RKF78  | 5.3e-14   | 2.9e-10 | 2.1e-08 | 7.9e-07 | **21.9 ms**  | 9.0e-11        | 5.1e-10   |
| GRKN64 | 1.1e-11   | 1.1e-09 | 4.2e-08 | 1.9e-07 | **19.1 ms**  | 2.7e-13        | 2.5e-14   |
| GAUSS  | 1.1e-13   | 2.8e-11 | 4.0e-10 | 3.7e-08 | 248.8 ms     | 1.4e-13        | 8.6e-15   |
| RK4    | 1.7e-12   | 4.3e-11 | 2.8e-10 | 3.7e-08 | 308.8 ms     | 3.5e-15 (*)    | 1.5e-15   |
| AAS    | 2.8e-11   | 1.6e-09 | 1.8e-07 | 4.7e-07 | 1125.4 ms    | 2.2e-09        | 1.1e-08   |
| SABA4  | —         | —       | —       | —       | 79.5 s       | —              | —         |
| RADAU  | da misurare | | | | | | |

Errori in AU. Riferimento utile: 1e-6 AU = 150 km ~ 0.1 arcsec a 2 AU.

(*) Il valore di RK4 nell'andata-ritorno e' ingannevole: a passo fisso i passi
all'indietro ripercorrono esattamente quelli in avanti e l'errore si cancella.
Non indica accuratezza.

### Cosa dicono queste misure

**GRKN64 e' il piu' efficiente sugli archi lunghi**: sul decennio e' insieme il
piu' rapido (19.1 ms) e il piu' accurato (1.9e-7 AU). Coerente con la teoria —
i metodi di Nystrom sfruttano la struttura del second'ordine. Da verificare
prima di promuoverlo a default: il modello include la correzione relativistica,
che dipende dalla velocita', mentre i Nystrom puri assumono accelerazione
funzione della sola posizione. Il "generalized" del nome dovrebbe coprire il
caso, ma va confermato leggendo l'implementazione.

**I cinque metodi funzionanti concordano fra loro entro ~1e-7 AU**, mentre nella
prova contro JPL Horizons RKF78 mostrava uno scarto di ~2e-6 AU. Se tutti i
metodi concordano fra loro dieci volte meglio di quanto concordino con JPL, lo
scarto residuo non e' dell'integrazione ma del **modello di forze** (perturbatori
asteroidali assenti, o altro). E' l'indicazione piu' utile emersa dal banco.

**GAUSS e RK4 sono accurati ma lenti**: buona precisione anche sul decennio, a
un costo dieci volte superiore a RKF78.

**AAS paga la simpletticita'**: precisione paragonabile agli altri ma ~50 volte
piu' lento di RKF78. La deriva d'energia (1.1e-8) e' superiore a quella degli
altri metodi su questo arco: il vantaggio dei simplettici si manifesta su
integrazioni molto piu' lunghe di un decennio.

## Come scegliere

- **Propagazione, fit, STM, uso generale** -> RKF78
- **Stabilita' su tempi lunghissimi, studi dinamici** -> AAS
- **Verifiche didattiche, archi brevi** -> RK4
- Gli altri: non usare finche' non validati.

## Nota metodologica sulla validazione

La validazione usa JPL Horizons come oracolo, interrogato a **MJD interi esatti**.
Un disallineamento di mezzo giorno fra l'epoca degli elementi e quella
dell'oracolo introduce un errore spurio di ~5e-3 AU, sufficiente a far apparire
rotti integratori perfettamente funzionanti. Errore commesso e corretto durante
le prove del 2026-07-22.

Conversione: JD = MJD + 2400000.5

## Esecuzione dei test di validazione esterna

I test `ValidazioneEsterna.*` confrontano i risultati con riferimenti
indipendenti. Richiedono dati locali e **vengono saltati** se mancano:

    export ASTDYN_EPHEMERIS_PATH=~/.ioccultcalc/ephemerides/de440s.bsp
    export ASTDYN_OBSCODES=~/.ioccultcalc/observatories/ObsCodes.txt
    export ASTDYN_TEST_DATA=<repo>/astdyn/tests/data
    ctest --test-dir build -R ValidazioneEsterna --output-on-failure

Un test saltato non e' un test passato: senza queste variabili la copertura
sulla catena osservativa e' assente.

## Esecuzione dei test di validazione esterna

I test `ValidazioneEsterna.*` confrontano i risultati con riferimenti
indipendenti. Richiedono dati locali e **vengono saltati** se mancano:

    export ASTDYN_EPHEMERIS_PATH=~/.ioccultcalc/ephemerides/de440s.bsp
    export ASTDYN_OBSCODES=~/.ioccultcalc/observatories/ObsCodes.txt
    export ASTDYN_TEST_DATA=<repo>/astdyn/tests/data
    ctest --test-dir build -R ValidazioneEsterna --output-on-failure

Un test saltato non e' un test passato: senza queste variabili la copertura
sulla catena osservativa e' assente.
