# Catalogo stellare: migrazione al formato multifile

Nota del 27 luglio 2026. Vale la pena farlo, ma non subito: serve spazio disco
che al momento non c'e'.

---

## Il problema

Il catalogo che distribuiamo — SQLite, estratto di Gaia DR3 fino a magnitudine
14 — contiene **una sola banda fotometrica**, la G. Lo si vede dallo schema:

```sql
CREATE TABLE stars (
    sid INTEGER PRIMARY KEY, ra REAL, dec REAL,
    pmra REAL, pmdec REAL, plx REAL, mag REAL, ruwe REAL,
    name TEXT, bayer TEXT, flam INTEGER, const TEXT,
    hd INTEGER, hip INTEGER, sao INTEGER
);
```

Mancano BP e RP. Conseguenze sull'output consegnato agli osservatori:

- i campi B e R di `<Star>` nell'XML escono a `0.00`, che si legge come
  "stella di magnitudine zero";
- il **calo di magnitudine in banda R non e' calcolabile**: con `mag_r = 0` il
  calcolo dava 21.11 su 820987, cioe' un calo di ventun magnitudini — un numero
  privo di senso in un file destinato a chi va a osservare.

Il calo in **V** si calcola correttamente (10.69 su quel caso, coincidente con
la predizione di riferimento) perche' la G di Gaia e' un buon surrogato della V.

## Cosa esiste gia'

`external/ioc_gaialib` implementa un secondo formato, **multifile_v2**, molto
piu' ricco:

| | SQLite (in uso) | multifile_v2 |
|---|---|---|
| stelle | 78 milioni | 231-303 milioni |
| limite | G ≤ 16 | G ≤ 18 |
| bande | G | G, BP, RP |
| indice | R-tree | HEALPix |
| dimensione | 14.7 GB | 14-19 GB |
| cono 0.5° | — | 0.001 ms |

Il codice per leggerlo c'e' ed e' compilato
(`src/catalog/gaia/concurrent_multifile_catalog_v2.cpp`), legge le tre bande
(righe 402-403) e si configura cosi':

```json
{
  "catalog_type": "multifile_v2",
  "multifile_directory": "~/.catalog/gaia_mag18_v2_multifile",
  "max_cached_chunks": 100
}
```

**Manca solo il catalogo**: la directory e' stata cancellata, e
`~/.catalog/gaia_mag18_v2.mag18v2` e' oggi un collegamento rotto.

## Come ricostruirlo

Procedura documentata in `external/ioc_gaialib/CATALOG_DOWNLOADS.md`:

```bash
# 1. GRAPPA3E dal FTP dell'IMCCE — 146 GB
cd ~/catalogs && mkdir -p GRAPPA3E && cd GRAPPA3E
wget -r -np -nH --cut-dirs=3 ftp://ftp.imcce.fr/pub/catalogs/GRAPPA3E/

# 2. conversione — circa 75 minuti
./build_mag18_catalog_v2 ~/catalogs/GRAPPA3E ~/catalogs/gaia_mag18_v2.cat
```

**Requisiti**: ~160 GB liberi durante la lavorazione (146 di sorgente piu' 14 di
risultato). Al 27 luglio 2026 il disco ne ha 34.

## Cosa guadagneremmo

- **BP e RP**, quindi i campi B e R dell'XML valorizzati e il calo in R
  calcolabile — informazione che serve all'osservatore per scegliere il filtro;
- **magnitudine 18 invece di 16**: piu' stelle deboli, quindi piu' eventi
  candidati per chi ha strumentazione adeguata;
- **quattro volte le stelle** (303 contro 78 milioni);
- **ricerca per cono molto piu' rapida** grazie all'indice HEALPix, il che conta
  sulle campagne massive dove la scansione del corridoio domina il tempo.

## Cosa non guadagneremmo

Nulla sul **calcolo dell'evento**: istante, geometria dell'ombra, incertezza e
calo in V dipendono dalla posizione, dal moto proprio e dalla banda G, che gia'
abbiamo. Il multifile aggiunge completezza e informazione fotometrica, non
accuratezza.

## Nota sulla distribuzione

Il pacchetto attuale (`gaia_mag14.db.xz`) e' 1.2 GB compressi da 3.2. Il
multifile a magnitudine 18 sarebbe 14 GB, forse 5 compressi: un ordine di
grandezza in piu' da ospitare e da far scaricare.

Se si fara', conviene produrre **estratti per magnitudine** anche da quel
formato, come gia' fatto per lo SQLite: un `mag14` con tre bande resterebbe
intorno al gigabyte e coprirebbe la gran parte degli usi.

## Nel frattempo

I campi non disponibili vanno **dichiarati tali**, non riempiti di zeri:
- B e R nell'XML: il valore convenzionale per "assente", non `0.00`;
- calo in R: non scritto quando `mag_r` manca, invece del `21.11` attuale;
- nel JSON, che e' formato nostro, si puo' usare `null` senza ambiguita'.

Un campo vuoto dice la verita'; un campo con dentro un numero sbagliato no.
