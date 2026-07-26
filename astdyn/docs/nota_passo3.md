# Passo 3 (covarianza dal fit) — stato e nodo da sciogliere

## Cosa e' gia' verificato

`compute_covariance_internal` calcola `(A^T W A)^-1` e la propaga fino a
`OrbitDeterminationResult::covariance`. Misurato su BK290 (90 osservazioni):

- **unita': AU e AU/giorno**, le stesse della covarianza da .eq1.
  Sigma di posizione 120, 100, 77 km. L'interpretazione in SI darebbe
  nanometri, quindi non c'e' ambiguita' e **non serve alcuna conversione**.
- **frame ed epoca**: ECLIPJ2000 cartesiano all'epoca degli elementi, come
  `stored_covariances`.
- **struttura sensata**: correlazione 0.86 fra le componenti di posizione nel
  piano orbitale, -0.963 fra posizione e velocita' lungo x, 0.952 lungo y.
  Sono le correlazioni forti attese in un'orbita ben determinata.
- **ordine di grandezza compatibile** con AstDyS: incertezza cross-track
  dell'evento 131.80 km, nostre sigma 77-120 km.

Il collegamento e' quindi di poche righe: nel punto di ioccultcalc.cpp dove il
fit sostituisce l'orbita, catturare anche `esito.covariance` e assegnarla a
`stored_covariances[id]`.

## Il nodo da sciogliere PRIMA di collegarla

`(A^T W A)^-1` e' la covarianza **formale**: presuppone che i pesi siano
corretti, cioe' che i sigma delle osservazioni descrivano la dispersione reale.

Ma i numeri dicono altro:
- RMS del nostro fit: **1.59 arcsec**
- sigma tipici nel .rwo: **0.4-0.5 arcsec**
- rapporto ~3, quindi chi quadro ridotto intorno a **10**

I residui sono tre volte piu' grandi di quanto i pesi prevedano. In questi casi
la pratica standard e' gonfiare la covarianza di sqrt(chi2/dof): con un fattore
3 le nostre sigma passerebbero da ~120 a ~360 km, e l'ellisse risulterebbe circa
**tre volte piu' larga** di quella AstDyS.

Collegare la covarianza formale senza affrontare questo significherebbe
dichiarare un'incertezza ottimistica — l'opposto di cio' che serve in una
predizione di occultazione, dove sottostimare l'ellisse porta a mancare l'evento.

## La domanda vera

Perche' il nostro RMS e' tre volte i sigma delle osservazioni? Due ipotesi:

1. **Modello dinamico incompleto su 16 anni.** I perturbatori asteroidali non
   vengono applicati (asteroid_ephemeris_file vuoto per default, benche'
   sb441-n16.bsp sia disponibile). AstDyS li include. Su archi lunghi l'effetto
   e' dell'ordine giusto.
2. **Pesi da rivedere.** Il .rwo porta anche le correzioni di debiasing dei
   cataloghi stellari (colonne bias RA/Dec), che non applichiamo. Valgono
   0.2-0.3 arcsec: non spiegano da soli un fattore 3, ma contribuiscono.

Nota: AstDyS pubblica RMSast = 0.463, ma e' un RMS **normalizzato sui sigma**,
quindi non confrontabile direttamente con il nostro. Il confronto onesto e'
quello sui residui singoli, che dopo le correzioni del 2026-07-22 danno 0.199
arcsec contro 0.169 di AstDyS — molto migliore del rapporto suggerito dall'RMS
del fit. La discrepanza fra i due confronti merita essa stessa chiarimento.

## Come procedere

1. capire il fattore 3 (perturbatori asteroidali per primi: sono gia' scaricati);
2. decidere se applicare il fattore di gonfiamento sqrt(chi2/dof);
3. collegare la covarianza (poche righe);
4. **verificare l'ellisse risultante contro AstDyS** — e' il criterio di
   accettazione: un accordo entro qualche percento conferma la catena, uno
   scarto grande indica che qualcosa nella normalizzazione non torna.

Strumento pronto: `examples/test_covarianza.cpp`.
