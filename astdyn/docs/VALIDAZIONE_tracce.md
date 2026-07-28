# Validazione delle tracce d'ombra contro predizioni pubblicate

Misure del 28 luglio 2026. Strumento: `tools/confronta_tracce.py`, che legge
file occelmnt e ne ricostruisce il percorso dagli elementi besseliani.

---

## Sintesi

| caso | oggetto | distanza | confronto | scarto |
|------|---------|----------|-----------|--------|
| A | (1216) Askania | 1.07 AU | stessi elementi, contro Occult4 | **46 km** |
| B | (1216) Askania | 1.07 AU | elementi AstDyS, contro JPL | 332 km |
| C | (1216) Askania | 1.07 AU | AstDyS + nostro fit, contro JPL | 169 km |
| D | (65803) Didymos | 1.27 AU | stessi elementi, contro ACROSS/JPL | **730 km** |

Il caso A stabilisce che il calcolo geometrico e' corretto. Gli altri misurano
quanto pesano, rispettivamente, la differenza fra soluzioni orbitali (B), il
raffinamento locale dell'orbita (C) e la **precisione con cui gli elementi
vengono trasmessi** (D).

---

## A. Askania con gli stessi elementi: 46 km

Riproducendo l'evento del 15 agosto 2026 con gli elementi del campo `<Orbit>`
della predizione di riferimento — anziche' con quelli AstDyS — la nostra traccia
coincide con quella pubblicata entro **46 km al massimo avvicinamento**.

Questo isola il calcolo dalla sorgente orbitale: geometria besseliana,
propagazione, riduzione astrometrica e proiezione sull'ellissoide sono corrette.
I 46 km residui sono compatibili con la precisione degli elementi ricopiati
(quattro decimali sull'anomalia media, che a 2.23 AU valgono gia' ~400 km lungo
l'orbita).

I termini besseliani di ordine superiore coincidono alla quarta cifra
significativa: `d2x = 0.0039880` contro `0.0039883` del riferimento.

## B e C. Askania: soluzioni orbitali diverse

Quattro predizioni dello stesso evento, misurate al massimo avvicinamento:

| predizione | scarto da JPL#67 |
|------------|------------------|
| AstDyS, senza fit | 332 km |
| AstDyS + fit sulle osservazioni | **169 km** |
| Occult4 (stessi elementi) | 46 km |

**Il fit locale dimezza lo scarto.** Partendo da elementi AstDyS distanti 332 km
dalla soluzione JPL, il raffinamento sulle osservazioni astrometriche porta a
169 km — pur usando l'insieme osservativo e la pesatura di AstDyS.

Buona parte della discrepanza fra i due centri non sta quindi nei dati ma negli
elementi pubblicati: nella loro epoca, nel modello di forze, o nel momento
dell'ultimo aggiornamento.

**Le incertezze formali sono ottimistiche.** I corridoi a 1 sigma valgono
+/-20 km (nostro), +/-22 km (nostro fittato) e +/-32 km (JPL), mentre gli scarti
sono di centinaia di chilometri. Tre determinazioni indipendenti si escludono a
vicenda di oltre cinque deviazioni standard.

## D. Didymos: la precisione degli elementi diventa dominante

(65803) Didymos, evento del 27 aprile 2026, predizione ACROSS
(Observatoire de la Côte d'Azur, finanziato ESA) da soluzione JPL#240+INTG.

Riproducendo con **gli stessi elementi** del campo `<Orbit>`:

| | ACROSS | nostro |
|---|---|---|
| x, y (raggi terrestri) | 0.3177, −0.9344 | 0.2794, −0.8274 |
| modulo | 0.9869 | 0.8733 |
| UT massimo avvicinamento | 5.2094 | 5.0056 |
| d2x | −0.0038989 | −0.0038998 |
| d3x | 0.0000008 | 0.0000009 |

I termini di ordine superiore coincidono, ma la **posizione** differisce di 0.11
raggi terrestri — 730 km — e l'istante di 12 minuti.

La causa e' la precisione degli elementi trasmessi. Il campo `<Orbit>` riporta
l'anomalia media con quattro decimali, e su Didymos 0.0001 gradi valgono 429 km
lungo l'orbita. A 1.27 AU dalla Terra, con moto apparente tre volte piu' rapido
di un asteroide di fascia principale, l'arrotondamento si amplifica: cinque
unita' sull'ultima cifra bastano a produrre lo scarto osservato.

**La specifica del formato occelmnt lo dice**: quel campo e' *"low-precision,
used to plot the motion of the object on star charts"*. Non e' pensato per
ricostruire una predizione, e su un oggetto vicino non lo consente.

### Corollario: due orbite eccellenti, un evento incerto

NEODyS e AstDyS pubblicano per Didymos soluzioni che coincidono a **otto cifre
significative** sul semiasse. Riportate alla stessa epoca, la loro anomalia media
differisce dalla soluzione JPL usata da ACROSS di 0.0027 gradi — che a 1.27 AU
sottendono **12.6 arcsec** sul cielo.

Su un evento radente come questo, dove il centro dell'ombra passa a 0.987 raggi
terrestri dal centro della Terra, quella differenza sposta il candidato: con
l'orbita NEODyS la stella UCAC4 338-104450 non risulta piu' occultata, e diventa
tale un'altra stella distante due primi d'arco.

Non e' un difetto del calcolo. E' la ragione per cui le campagne su NEA
mobilitano decine di osservatori distribuiti su fasce larghe, e per cui ogni
corda positiva raccolta migliora sensibilmente le predizioni successive.

---

## Nota metodologica

Lo scarto fra tracce si misura in tre modi, che danno risultati diversi:

- **al massimo avvicinamento** — la geometria e' definita senza ambiguita', ed e'
  il valore confrontabile con l'ellisse d'incertezza;
- **perpendicolarmente**, come minimo della distanza da tutta la traccia di
  riferimento — e' di quanto un osservatore deve spostarsi;
- fra punti **contemporanei** — sovrastima molto sulle code, dove le tracce
  divergono in direzione, e non va usato.

Su tracce parallele i primi due coincidono, ed e' il caso di tutti i confronti
qui riportati.

## Dati e riproducibilita'

Predizioni di riferimento: OccultWatcher Cloud per Askania; ACROSS
(<https://www.oca.eu/fr/prediction>) per Didymos, file
`Didymos_April-May.xml`, 375 eventi.

Elementi orbitali: AstDyS e NEODyS (<https://newton.spacedys.com>).

Configurazioni di campagna in `tools/campagna_*.yaml`.
