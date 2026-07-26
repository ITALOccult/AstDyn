# ioccultcalc-setup — installer e aggiornatore dei dati esterni

## Il problema

Chi clona il repository oggi non ha `ObsCodes.txt`: il motore riduce le
osservazioni dal geocentro e il fit peggiora di un fattore sei (RMS 1.59"
invece di 0.26" su BK290). L'avviso c'e', ma puo' sfuggire.

Piu' in generale il sistema dipende da sette categorie di dati esterni, alcune
pesanti, alcune con procedure di reperimento non ovvie, nessuna installata
automaticamente. Il lavoro di oggi sulla Fase 6 si e' fermato proprio per questo:
il motore multi-corpo funziona (validato su Haumea con Hi'iaka e Namaka) ma
l'effemeride dei satelliti non e' piu' sul disco e non c'e' modo documentato di
riprocurarsela.

## Convenzioni esistenti da rispettare

Struttura in `~/.ioccultcalc/`:

    cache/          (vuota)
    config/         (vuota)
    data/
    database/       allnum.db e altri
    elements/       .eq1 AstDyS
    ephemerides/    de440s.bsp, sb441-n16.bsp, de441.bsp, ...
    gaia/           gaia_cache/
    iers/           (vuota)
    logs/           (vuota)
    observatories/  ObsCodes.txt
    obscodes_extended.json      <- in radice, formato diverso

Da chiarire: `obscodes_extended.json` (radice) e `observatories/ObsCodes.txt`
sono lo stesso dato in due formati e due posti. Il motore legge il secondo.
Il tool dovrebbe unificare o almeno documentare la differenza.

## Modello da seguire: update_asteroids.py

Ha gia' le qualita' giuste e va riusato, non sostituito:
- lavora su una COPIA temporanea, verifica, poi swap atomico con backup .bak;
- nessun DELETE: cio' che le fonti non coprono resta intatto;
- accetta un file locale al posto del download (`--neodys-cat`);
- `--dry-run` che riporta senza modificare.

Il tool nuovo estende questa logica alle altre categorie e per gli asteroidi
invoca lo script esistente.

## Le sette categorie

| dato | dove | dimensione | sorgente | priorita' |
|------|------|-----------|----------|-----------|
| codici osservatorio | observatories/ObsCodes.txt | ~200 KB | MPC (pagina HTML) | **alta** |
| effemeridi planetarie | ephemerides/de440s.bsp | 31 MB | JPL NAIF | **alta** |
| perturbatori asteroidali | ephemerides/sb441-n16.bsp | 645 MB | JPL NAIF | media |
| catalogo asteroidi | database/allnum.db | — | AstDyS + NEODyS | **alta** |
| catalogo stellare Gaia | gaia/…db | grande | da definire | **alta** |
| effemeridi satelliti | ephemerides/ | varie | JPL Horizons, caso per caso | bassa |
| elementi/osservazioni | elements/, obs/ | piccole | AstDyS, per oggetto | (orchestrator) |

### Note per categoria

**Codici osservatorio.** La sorgente e' una pagina HTML
(`minorplanetcenter.net/iau/lists/ObsCodes.html`) con il contenuto a colonne
fisse dentro `<pre>`: va ripulita dai tag. Verifica di integrita': attese ~2700
righe, e devono comparire codici noti (500, 691, F51, G96).

**Effemeridi planetarie.** `de440s.bsp` copre 1849-2150 ed e' sufficiente.
ATTENZIONE: `de441.bsp` (3 GB) blocca la macchina — non scaricarlo per default.

**Effemeridi satelliti.** Per i sistemi ben studiati (Plutone, Haumea) esistono
kernel SPK JPL; per la gran parte dei binari noti **non esistono affatto** e le
orbite vanno da modelli pubblicati caso per caso. Il tool puo' scaricare i
kernel disponibili, non risolvere il problema generale.

**Catalogo Gaia.** Referenziato come
`~/.catalog/crossreference/gaia_dr3_occult_pro.db` — fuori da `~/.ioccultcalc/`.
Da chiarire come viene prodotto e se sia distribuibile.

## Interfaccia

    ioccultcalc-setup [--all] [--obscodes] [--ephemerides] [--catalog]
                      [--gaia] [--check] [--force] [--dry-run]

- `--check` (default se nessuna opzione): verifica cosa c'e' e cosa manca,
  senza scaricare nulla. E' il comando che un nuovo utente lancia per primo.
- `--all`: installa il minimo necessario a un sistema funzionante.
- `--force`: riscarica anche se presente e recente.

## Requisiti

- **Cache con timestamp**: non riscaricare se recente, salvo `--force`.
- **Verifica dopo il download**: righe attese, formato, dimensione minima. Un
  file HTML di errore salvato come .bsp deve essere riconosciuto subito.
- **Cortesia**: User-Agent identificabile e throttle, come gia' fatto per AstDyS.
- **Nessun degrado silenzioso**: se un dato manca, dirlo chiaramente. E' la
  lezione ricorrente — il fallback al geocentro e' costato mezza giornata di
  diagnosi.
- **Idempotenza**: rilanciarlo non deve rompere nulla.

## Verifica

`--check` su un sistema vuoto deve elencare tutto cio' che manca; `--all` deve
portare da zero a un fit funzionante su BK290 con RMS ~0.26".
