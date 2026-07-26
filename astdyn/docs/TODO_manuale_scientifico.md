# Manuale scientifico — aggiornamenti necessari

Stato al 25 marzo 2026: 220 pagine, 19/26 capitoli completi (73%).
Le correzioni del 22 e 26 luglio 2026 rendono **inesatti** alcuni capitoli gia'
scritti e ne arricchiscono altri con misure nuove.

Ordinato per gravita': prima cio' che oggi il manuale afferma di sbagliato,
poi cio' che manca.

---

## A. Contenuti da CORREGGERE (il manuale dice cose non piu' vere)

### Cap. 4 — Reference Frames
La rotazione terrestre ITRF↔GCRF aveva il **segno invertito**: `rotation_z(-ERA)`
invece di `rotation_z(+ERA)` in `j2000_to_itrf_simple`. Portava Greenwich in
(R cos ERA, **-**R sin ERA, 0) invece di +R sin ERA — errore fino a 2 raggi
terrestri (12756 km).

Da verificare che il capitolo non riporti la formula sbagliata. Da aggiungere:
la distinzione fra le DUE implementazioni presenti in libreria — quella
semplificata (`itrf_to_j2000_simple`, solo ERA) e quella rigorosa
(`CelestialToTerrestrial`, catena IAU completa W(t)·R(ERA)·Q(t) con UT1) — e
quando usare l'una o l'altra.

### Cap. 8 — Numerical Integration
- **SABA4 non e' piu' supportato**: implementazione difettosa (drift rettilineo
  invece che kepleriano, stima d'errore proporzionale ad h). Se il capitolo lo
  presenta come metodo disponibile, va corretto.
- **RADAU**: i coefficienti erano corretti ma la stima d'errore era improvvisata
  (`b_hat_` con l'ultimo peso azzerato, dichiarato "crude" nel codice), e il
  metodo non compiva nemmeno un passo. Ora usa estrapolazione di Richardson.
  Vale la pena spiegare perche' i Radau IIA **non ammettono** una coppia
  incorporata alla maniera dei Runge-Kutta espliciti.
- Sostituire le eventuali tabelle di confronto con quelle misurate (vedi sotto).

### Cap. 14 — Differential Correction
- La **reiezione outlier** aveva la soglia legata a `iter/max_iterations`: con
  fit che convergono in 2-3 iterazioni su 20 restava a ~9.3 sigma, praticamente
  inattiva. Su 433 Eros scartava 0 osservazioni su 6155.
- La **semantica di convergenza**: un line search che non trova passi migliori
  segnala un minimo raggiunto, non un fallimento.
- Il **line search** provava fino a 14 dimezzamenti di alpha, ciascuno con una
  propagazione completa: ora ha un tetto.

### Cap. 15 — Residuals
Il modello di osservazione richiede **entrambe** le aberrazioni:
- planetaria (v_asteroide/c ≈ 12.4 arcsec per un corpo di fascia principale);
- stellare (v_Terra/c = 20.5 arcsec).

Verificato sperimentalmente: senza il termine dell'osservatore i residui su
820987 passano da 0.199 a 16.88 arcsec. Le posizioni MPC **non** sono depurate
dall'aberrazione stellare — i cataloghi di riferimento (Gaia) danno posizioni
baricentriche.

Sul segno: la nota di Milani in OrbFit (`aber1`) prescrive V(corpo) −
V(osservatore), ma in un'espressione impostata diversamente; nella nostra
formulazione `r − v·tau` i due contributi entrano con lo stesso segno.

### Cap. 12 — Observations
Manca del tutto la **parallasse topocentrica** come fattore critico. Senza il
catalogo MPC degli osservatori le osservazioni vengono ridotte dal geocentro:
~4 arcsec a 2.2 AU, molto di piu' sui NEO. Misurato sul fit di 820987: RMS 1.59
arcsec senza catalogo contro 0.26 con.

Da trattare anche il caso degli **osservatori mobili** (codice 270, Unistellar
Network): non hanno coordinate in catalogo e le portano nei record della singola
osservazione. Su 433 Eros sono 832 osservazioni su 16324.

---

## B. Misure NUOVE da inserire

### Cap. 8 / 23 — confronto degli integratori
Propagazione di 820987 dall'epoca MJD 61200, riferimento RKF78 @ 1e-13,
modello con perturbazioni planetarie:

| metodo | 10 g | 100 g | 1000 g | 3650 g | tempo 3650 g |
|--------|------|-------|--------|--------|--------------|
| RKF78  | 5.3e-14 | 2.9e-10 | 2.1e-08 | 7.9e-07 | 21.9 ms |
| GRKN64 | 1.1e-11 | 1.1e-09 | 4.2e-08 | 1.9e-07 | 19.1 ms |
| RADAU  | 7.7e-13 | 2.5e-10 | 7.5e-09 | **2.1e-08** | 177.8 ms |
| GAUSS  | 1.1e-13 | 2.8e-11 | 4.0e-10 | 3.7e-08 | 248.8 ms |
| RK4    | 1.7e-12 | 4.3e-11 | 2.8e-10 | 3.7e-08 | 308.8 ms |
| AAS    | 2.8e-11 | 1.6e-09 | 1.8e-07 | 4.7e-07 | 1125.4 ms |

Errori in AU. Osservazione di rilievo: **i cinque metodi funzionanti concordano
fra loro entro ~1e-7 AU**, mentre nel confronto con JPL Horizons RKF78 mostra
uno scarto di ~2e-6. Se i metodi concordano fra loro dieci volte meglio di
quanto concordino con l'oracolo, lo scarto residuo non e' dell'integrazione ma
del **modello di forze**.

### Cap. 10 — State Transition Matrix
Jacobiana **analitica** contro numerica, confronto della matrice Phi propagata:
differenza relativa 1.7e-9, tempo 66.5 ms contro 1.0 ms (**66x**). L'analitica e'
ora il default. Forma chiusa: `d(accel)/d(pos) = -mu/r^3 · I + 3mu/r^5 · (r r^T)`,
sommata su Sole e perturbatori planetari.

### Cap. 9 / 23 — dense output
L'interpolazione di Hermite cubica per valutare la soluzione a istanti arbitrari
senza fermare l'integratore. Errore misurato: 1528 km (1.05 arcsec) con passo
naturale, 62.7 km (0.043 arcsec) limitando il passo a 2 giorni — scala come h^4.

Motivazione: le osservazioni sono raggruppate per notte (fino a 7 misure a
minuti di distanza), e fermare l'integratore su ciascuna faceva collassare il
passo. Il costo diventa proporzionale all'arco, non al numero di osservazioni.

### Cap. 18 / 21 — risultati di validazione
Residui su 820987 (2015 BK290), 78 osservazioni, confronto con quelli calcolati
da AstDyS nello stesso file .rwo:

| stato | residuo medio | offset RA | offset Dec |
|-------|--------------|-----------|------------|
| iniziale | 16.57" | -16.51" | -3.42" |
| + aberrazione stellare | 1.41" | +0.50" | -1.21" |
| + rotazione terrestre corretta | **0.199"** | -0.052" | +0.076" |
| riferimento AstDyS | 0.169" | — | — |

Fit orbitale completo: RMS 0.26 arcsec su 90 osservazioni (16 anni di arco),
convergenza in 2 iterazioni.

Covarianza dal fit contro quella AstDyS, ellisse di predizione per l'evento del
27 luglio 2026:

| | semiasse magg. | semiasse min. | PA | cross-track |
|---|---|---|---|---|
| AstDyS | 0.117721" | 0.064415" | 71.3895 | 131.80 km |
| dal fit | 0.107151" | 0.063296" | 71.3871 | 120.40 km |

Accordo del 9% sul semiasse maggiore, 2% sul minore, orientamento coincidente.

---

## C. Capitoli MANCANTI o da estendere

### Parte IV — Software Implementation (cap. 16-20), l'unica incompleta
Da scrivere. Materiale disponibile oggi:
- architettura a due stadi delle campagne (screening leggero, raffinamento);
- sistema di frame tipizzato C++20 (gli errori di frame sono impossibili a
  compile-time — e' la scelta di progetto piu' caratterizzante);
- covarianza per-asteroide dagli .eq1 AstDyS;
- il fit come parametro di calcolo, non come passo obbligato.

### Cap. 26 — Future Developments
Aggiornare con lo stato reale: Fase 5 completata (fit operativo in campagna con
covarianza propria), Fase 6 in attesa delle effemeridi di satelliti, Fase 7
(installer) avviata.

---

## D. Nota metodologica che meriterebbe un posto nel manuale

Tutti i difetti trovati nel luglio 2026 erano **silenziosi** e **invisibili ai
177 test unitari**: rotazione terrestre invertita, database osservatori vuoto con
ricaduta muta sul geocentro, SABA4 che restituiva spazzatura, RADAU che non
compiva un passo, reiezione outlier inattiva.

Sono emersi solo confrontando i risultati con riferimenti **esterni** — AstDyS
per i residui, JPL Horizons per le posizioni, i primi principi per la rotazione
terrestre — e osservazione per osservazione, non sull'RMS aggregato, che
maschera i sistematici.

E' una lezione metodologica che vale oltre questo software: la coerenza interna
non e' correttezza.
