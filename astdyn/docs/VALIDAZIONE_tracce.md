# Validazione delle tracce d'ombra

Aggiornato all'11 ottobre 2026. Strumenti: `tools/confronta_tracce.py` per il
confronto fra predizioni, `tools/verifica_corde.py` per il confronto con le
osservazioni.

## 11 ottobre 2026: frame, tempo-luce, integratore, riduzione del fit

Campagna su (3830) Trelleborg, ottobre 2026 - dicembre 2027. Ogni difetto e'
stato misurato contro un riferimento esterno prima di essere corretto.

| difetto | riferimento | effetto | dopo |
|---------|-------------|---------|------|
| frame bias FK5->ICRS dentro ECLIPJ2000 | vettori Horizons ICRF/eclittica | 47 km sulla posizione, fino a 52 km sulla traccia | < 1 m |
| RKF78 a tolleranza 1e-11 (default) | RADAU a 1e-13 | 55 km in 565 giorni | 0.07 km (default 1e-13) |
| tempo-luce eliocentrico (moto del Sole ignorato) | RA/Dec astrometrica Horizons | 1-2 mas | 0.05 mas |
| deflessione piena applicata all'asteroide nel fit | RA/Dec astrometrica Horizons | fino a 12 mas | 0.1 mas |
| sito MPC con doppio schiacciamento | residui AstDyS per osservatorio | ~15 mas, segno dell'emisfero | < 5 mas |
| precessione-nutazione assente sul sito | catena IAU completa | 11-16 km di sito, 2-5 mas | — |
| debiasing FCCT14 del .rwo ignorato | residui AstDyS | 0.10" di dispersione | 0.014" RA, 0.005" Dec |
| osservazioni da satellite (TESS, WISE) scartate | residui AstDyS | 2521 osservazioni perse | lette, concordi entro 8 mas |
| pesi indipendenti per lotti notturni | soluzione JPL#66 | orbita a 180-265 km da JPL | 68-85 km (Veres 2017) |
| perturbatori asteroidali dallo SPK non ruotati in eclittica | vettori Horizons 1949-2024 di (3930) | 145 km nel 2015, 2718 km nel 1949 | 2.9 km, 213 km |
| TAI-UTC fisso a 37 s per ogni data | storia IERS/USNO, Delta T | 1-2 s su Gaia 2014-16, 30-40 s sulle lastre | per data |
| astrometria Gaia FPR non usata | soluzione JPL#64 di (3930) | fit a 37 km da JPL | <= 9 km con le 472 osservazioni FPR |
| aberrazione annua assente dal punto sub-stellare e dalla proiezione | elementi Occult4 di (316) Goberta | 20" sulla stella apparente, 1.3 s di rotazione | 0.005" |
| elementi besseliani esportati su assi J2000 invece che della data | elementi Occult4, corde di Goberta | 0.14 gradi, 14 km lungo la traccia, -1.0 s | -0.2/+0.1 s come JPL |
| riga satellite .rwo letta per token (nota in colonna 15) | residui TESS 2023 di Goberta | 836 osservazioni a 50" | lette per colonna |
| date prima di mezzogiorno UT (32 ottobre) | calendario | 6 date impossibili su 205 | test di andata e ritorno |
| filtro diurno valutato fuori dalla Terra | campionamento della traccia | 54 eventi a 2-5 gradi dal Sole | 95 eventi osservabili su 205 |
| occelmnt: campi vuoti letti con double.Parse | lettore di Occult4 (sorgente) | file rifiutato da Occult4 | test sulle regole del lettore |

Lo scarto diagnostico fra la covarianza propagata (STT) e la traccia, che era
di 32-52 km su tutti gli eventi, e' sceso a 1-10 km: i primi due difetti ne
spiegavano quasi tutto, il resto veniva dal modello di forze ridotto della STT,
ora allineato a quello della traccia.

### (316) Goberta contro le corde, catena corretta (12 ottobre 2026)

Corde di M. Harman (Vartovka, SK, 48 43 04.4 N, 19 09 13.7 E, 570 m) e
F. Garcia (La Vara, Valdes, ES, 43 28 21.6 N, 6 23 26.7 W, 219 m), da SODIS.

| predizione | Harman | Garcia | trasversale |
|------------|--------|--------|-------------|
| Occult4, JPL#90+INTG (riferimento) | -0.2 s | +0.1 s | 10.7 / 5.5 km |
| nostra catena, orbita JPL | -0.2 s | +0.1 s | 10.6 / 5.4 km |
| nostra catena, nostro fit (rwo + Gaia FPR) | -2.2 s | -1.9 s | 7.8 / 8.3 km |

Con la stessa orbita la catena riproduce Occult4 e le corde. Il nostro fit di
Goberta resta 16 mas avanti lungo la traccia (2 s): sui dati da terra 2020-2026
le due orbite differiscono di ~10 mas in RA, dentro quanto quei dati
vincolano (~7 mas), e la scelta la fanno pesi e debiasing (JPL: Eggl 2020 e
Veres 2017; noi: FCCT14 di AstDyS). Le corde danno ragione a JPL. L'astrometria
TESS di Goberta ha un sistematico di -0.20" in RA e +0.17" in Dec rispetto
all'orbita JPL, ma col peso per notte non sposta il fit.

### Accordo con JPL, a catena corretta

Su (3930) Vasilev, con perturbatori, scala dei tempi e astrometria Gaia FPR, il
nostro fit e JPL#64 coincidono entro 3.6 mas (<= 9 km sulla traccia) dal 2014 al
2027, con gli stessi residui: 2.5 mas lungo la scansione Gaia, 0.305" da terra.
La discrepanza descritta sotto per (3830) era dovuta a questi difetti, non ai
dati.

### Le soluzioni orbitali non concordano entro le loro incertezze (superato)

Anche corretta la catena, la nostra orbita e quella di JPL distano 25-55 mas in
RA, 68-85 km sulla traccia, mentre JPL dichiara 3-5 mas a 3 sigma. Le due
orbite adattano le osservazioni da terra 2018-2026 con lo stesso RMS (0.445" e
0.450"): i dati non bastano a scegliere. Come su (1216) Askania, le incertezze
formali sono ottimistiche: per (3830), diametro 17.8 km, l'incertezza realistica
e' di alcune larghezze della traccia, e un aggiornamento astrometrico prima
dell'evento resta necessario.

Restano da ripetere, con il codice corretto, il confronto con le corde di
(316) Goberta e quello con JPL su (1216) Askania.

---

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

Una nota metodologica, corretta l'11 ottobre 2026: la RA/Dec "astrometric" di
Horizons in modalita' OBSERVER **non** include la deflessione relativistica —
la documentazione di Horizons la descrive corretta per il solo tempo-luce — ed
e' un riferimento valido. Con gli elementi osculanti JPL di (3830) la nostra
catena la riproduce entro 0.05 mas su tre date del 2027. I 0.2 arcsec osservati
a luglio venivano da difetti allora presenti (frame bias, sotto).

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
