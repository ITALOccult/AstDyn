# Validazione delle tracce d'ombra

Aggiornato al 29 luglio 2026. Strumenti: `tools/confronta_tracce.py` per il
confronto fra predizioni, `tools/verifica_corde.py` per il confronto con le
osservazioni.

---

## Sintesi

| verifica | prima | dopo |
|----------|-------|------|
| (316) Goberta, contro due corde osservate | 70-83 km | **3-13 km** |
| (1216) Askania, contro JPL, elementi AstDyS | 332 km | **235 km** |
| (1216) Askania, contro JPL, con fit locale | 169 km | **82 km** |
| posizione dell'asteroide, contro Horizons | 0.052 arcsec | **0.010 arcsec** |

Il "prima" e' lo stato del 27 luglio; il "dopo" include la correzione
dell'obliquita' del frame ECLIPJ2000 e la riduzione astrometrica completa della
stella.

---

## Il difetto trovato: l'obliquita' del frame

### Come e' emerso

Il confronto fra predizioni non lo mostrava. Su (1216) Askania i coefficienti
besseliani coincidevano con quelli di riferimento alla quarta cifra
significativa, la geometria era identica, la stella la stessa.

Il difetto e' emerso confrontando con due corde **realmente osservate** di
(316) Goberta, 13 luglio 2026: Faustino Garcia (La Vara, Spagna) e Marek Harman
(Vartovka, Slovacchia), entrambe positive con tempi GPS e NTP, durate 3.41 e
3.48 secondi.

La traccia predetta passava **70-83 km** fuori posto. Su un'ombra larga 56 km,
un osservatore posizionato secondo quella predizione manca l'evento.

### La causa

La conversione da coordinate eclittiche a equatoriali usava l'obliquita' IAU
2006 (84381.406 arcsec), mentre gli elementi orbitali di AstDyS, NEODyS e JPL
sono riferiti a un piano eclittico costruito con la convenzione IAU 1976
(84381.448).

I 42 mas di differenza sono una rotazione spuria fra i due sistemi: 55 km sulla
traccia di un oggetto a 2.7 AU, 22 km a 1.07 AU.

Non e' una questione di accuratezza. Per ruotare un vettore fra due sistemi
serve l'angolo con cui quei sistemi sono stati definiti; usarne un altro
introduce un errore indipendentemente da quale valore sia fisicamente migliore.

**Occult4 mantiene la stessa distinzione**: nel suo codice due costanti separate,
`Ecliptic2000_deg_1976` e `Ecliptic2000_deg_2006`, e la prima nel modulo delle
occultazioni (`PlanetOccultationElements.cs:2187`).

### Cosa era stato escluso prima di trovarlo

Ogni anello della catena, verificato contro riferimenti esterni:

| anello | metodo | esito |
|--------|--------|-------|
| elementi orbitali | JPL Horizons, stessa epoca | identici |
| kepleriano -> cartesiano | calcolo indipendente | 12 cifre |
| effemeride planetaria | vettori Terra e Marte contro Horizons | 12 cifre |
| propagazione | confronto con moto rettilineo | 0.5 km su 1363 s |
| tempo-luce | valore e convergenza contro Horizons | 9 ms |
| posizione stellare | contro il campo `<Star>` del riferimento | 0.0004 arcsec |
| geometria besseliana | formula contro il sorgente di Occult4 | identica |
| perturbatori asteroidali | stima analitica | < 1 km |
| deflessione, parallasse | calcolo per la geometria reale | 0.002 arcsec |
| tolleranza dell'integratore | 1e-12 contro 1e-14 | invariato |

Una nota metodologica: il confronto con Horizons in modalita' OBSERVER si e'
rivelato **fuorviante**. Quelle coordinate includono la deflessione
relativistica, mentre il vettore geocentrico della modalita' VECTORS no: le due
uscite differiscono di 0.2 arcsec, e confrontare la nostra posizione
astrometrica con le prime ha portato fuori strada per un'ora. Il confronto va
fatto sui vettori.

---

## (316) Goberta contro due corde osservate

Evento del 13 luglio 2026, stella TYC 6321-01104-1 di magnitudine 12.25, ombra
larga 56 km, durata massima 3.99 s.

| osservatore | distanza dall'asse predetta | dalla corda | scarto sull'istante |
|-------------|------------------------------|-------------|---------------------|
| F. Garcia (ES) | 12.2 km | 14.5 km | −3.0 s |
| M. Harman (SK) | 0.7 km | 13.7 km | −3.3 s |

La predizione di riferimento (JPL#90+INTG) sugli stessi dati da' 5.5 e 10.7 km,
con scarti sull'istante di +0.1 e −0.2 s.

Lo scarto residuo sull'istante corrisponde ai 19 km di residuo sulla posizione,
visti lungo la traccia invece che trasversalmente.

---

## (1216) Askania contro la predizione JPL

Evento del 15 agosto 2026, occultazione di J175336.90−214720.9.

| predizione | scarto da JPL#67 |
|------------|------------------|
| elementi AstDyS, come pubblicati | 235 km |
| elementi AstDyS, raffinati sulle osservazioni | **82 km** |

**Il fit locale riduce lo scarto di un fattore tre.** Il fit converge su 4015
osservazioni degli ultimi quindici anni con RMS 0.233 arcsec, tre iterazioni, un
solo scarto.

Buona parte della discrepanza fra i due centri non sta quindi nei dati ma negli
elementi pubblicati: la loro epoca, il modello di forze, o il momento
dell'ultimo aggiornamento.

Il confronto non afferma che una soluzione sia migliore dell'altra. JPL e Pisa
usano insiemi di osservazioni diversi, debiasing diverso, e per i NEO JPL
include il radar. Quel che si misura e' quanto due determinazioni indipendenti e
accurate finiscano distanti sul terreno.

### Le incertezze formali sono ottimistiche

I corridoi a 1 sigma valgono ±20 km (nostro), ±22 km (nostro con fit) e ±32 km
(JPL). Gli scarti sono di centinaia di chilometri: le tre determinazioni si
escludono a vicenda.

Le covarianze formali assumono che i pesi delle osservazioni siano corretti e il
modello di forze completo. Un osservatore che si fidasse del corridoio a 1 sigma
si posizionerebbe con sicurezza nel posto sbagliato.

---

## Il diametro conta quanto la posizione

Il diametro decide la larghezza dell'ombra e la durata massima. Fino al 28
luglio veniva stimato dalla magnitudine assoluta con albedo fissa 0.10.

Su (316) Goberta questo dava **31.9 km** contro i **56.1 ± 0.6** misurati: la
durata prevista era 2.27 s, e una corda osservata di 3.4 s sarebbe stata piu'
lunga dell'intera ombra — cioe' geometricamente impossibile.

`tools/importa_diametri.py` porta ora nel catalogo locale 135 475 diametri
misurati da JPL SBDB (IRAS, WISE, Akari, radar, occultazioni) con la loro
incertezza. Copertura: 8.8% degli asteroidi numerati, ma comprende quasi tutti
quelli che producono occultazioni osservabili.

| oggetto | nostro | riferimento |
|---------|--------|-------------|
| (316) Goberta | 56.07 ± 0.64 | 54.6 |
| (1216) Askania | 10.53 ± 0.09 | 10.66 |
| (4769) Castalia | 1.400 | 1.400 |
| (65803) Didymos | 0.780 ± 0.03 | 0.78 |

---

## Nota sulla precisione degli elementi trasmessi

Il campo `<Orbit>` del formato occelmnt riporta l'anomalia media con quattro
decimali. La specifica lo dichiara *"low-precision, used to plot the motion of
the object on star charts"*.

Su (65803) Didymos, a 1.27 AU, un'unita' sull'ultima cifra vale 429 km lungo
l'orbita: riprodurre una predizione ACROSS con quegli elementi da' 730 km di
scarto, contro i 46 km ottenuti su Askania a 2.23 AU con lo stesso metodo.

Quel campo non e' adatto a ricostruire predizioni su oggetti vicini.

---

## Metodo

Ogni verifica confronta con qualcosa di esterno: un'osservazione, una predizione
pubblicata, un calcolo indipendente. In questo progetto il confronto del codice
con se stesso non ha mai trovato nulla.

L'errore dell'obliquita' e' stato invisibile per mesi a 177 test automatici e a
tutti i confronti fra predizioni. Lo hanno rivelato due corde registrate da
osservatori che non sapevano di star facendo debug di un software.

## Dati

Corde: SODIS (IOTA/ES), inserimenti 14414 e 14415.
Predizioni di riferimento: OccultWatcher Cloud, ACROSS (Observatoire de la
Côte d'Azur).
Elementi orbitali: AstDyS e NEODyS (<https://newton.spacedys.com>), JPL Horizons.
Diametri: JPL Small-Body Database.
