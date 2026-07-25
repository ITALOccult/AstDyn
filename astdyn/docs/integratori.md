# Gli integratori di AstDyn

Capitolo di manuale. Descrive i metodi implementati, cosa li distingue, e in quali
casi conviene usarli. I valori misurati provengono dal banco di prova
`examples/test_validate_integrators.cpp`, che confronta ogni metodo con l'oracolo
JPL Horizons su archi crescenti.

> STATO: le sezioni teoriche sono complete; le tabelle di misura vanno riempite
> eseguendo il banco di prova. I valori gia' presenti sono quelli raccolti il
> 2026-07-22 su 820987 (2015 BK290).

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

### SABA4 — simplettico di Yoshida
Metodo simplettico a composizione, pensato per hamiltoniane separabili in una parte
kepleriana dominante e una perturbativa piccola — esattamente la struttura del
problema asteroidale.

*Misurato:* **non funzionante.** Errore di ~5 AU (divergenza completa) e ~80
secondi per arco. Da considerare rotto finche' non verificato.

*Uso consigliato:* nessuno, allo stato attuale.

### Gauss-Jackson — multistep
Predittore-correttore multistep, storicamente il metodo di elezione per orbite
quasi circolari, ancora molto usato nella meccanica orbitale operativa. Riutilizza
le valutazioni dei passi precedenti, quindi e' economico per passo, ma richiede
un avviamento con un metodo a passo singolo e mal sopporta i cambi di passo.

*Misurato:* da caratterizzare.

*Uso consigliato:* da definire dopo le misure.

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

*Misurato:* da caratterizzare.

*Uso consigliato:* da definire dopo le misure.

## Tabella riassuntiva

| metodo    | tipo                  | passo     | accuratezza | costo    | stato       |
|-----------|-----------------------|-----------|-------------|----------|-------------|
| RK4       | Runge-Kutta 4         | fisso     | da misurare | —        | funzionante |
| RKF78     | RK-Fehlberg 7(8)      | adattivo  | ~2e-6 AU    | 28 ms/10a| **default** |
| AAS       | simplettico adattivo  | adattivo  | ~2e-6 AU    | 1190 ms/10a | funzionante |
| SABA4     | simplettico Yoshida   | fisso     | ~5 AU       | 80 s/arco| **rotto**   |
| GAUSS     | multistep             | quasi-fisso | da misurare | —      | non provato |
| RADAU     | implicito IIA         | adattivo  | da misurare | molto alto | da profilare |
| GRKN64    | Nystrom 6(4)          | adattivo  | da misurare | —        | non provato |

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
