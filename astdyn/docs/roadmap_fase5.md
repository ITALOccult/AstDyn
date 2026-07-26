# Fase 5 — Fit orbitale on-demand: roadmap

Stato al 2026-07-22: il nodo tecnico e' risolto. Il fit converge sull'arco pieno
di BK290 (90 osservazioni, 16 anni, 4 iterazioni) e il modello di osservazione e'
validato contro AstDyS (residui 0.199 vs 0.169 arcsec). Resta il lavoro di
integrazione: rendere il fit un parametro di calcolo utilizzabile in campagna.

## Cosa e' gia' fatto

- modello di osservazione corretto e validato (tempo-luce + aberrazione stellare,
  parallasse topocentrica, rotazione terrestre);
- `integrate_at` con dense output: le osservazioni raggruppate per notte non
  fermano piu' l'integratore;
- jacobiana analitica nella STM (66x piu' veloce della numerica);
- scarico dei `.rwo` per i positivi, gia' implementato nella Fase 3;
- test di regressione contro riferimenti esterni.

---

## Passo 1 — Parametri del fit nella configurazione  [piccolo]

Oggi il fit ha valori fissi nel codice. Servono come parametri:

| parametro          | significato                                  | default |
|--------------------|----------------------------------------------|---------|
| `fit.enabled`      | attiva il fit                                | false   |
| `fit.obs_years`    | usa solo le osservazioni degli ultimi N anni | 0 (tutte) |
| `fit.tolerance`    | tolleranza dell'integratore durante il fit   | 1e-10   |
| `fit.max_iter`     | iterazioni massime del correttore            | 20      |
| `fit.outlier_sigma`| soglia di reiezione                          | 3.0     |

`obs_years` e' un compromesso da documentare: arco corto significa fit veloce ma
orbita meno vincolata, perche' l'accuratezza dipende dalla baseline temporale.
Misurato su BK290: 19 osservazioni su 3 anni danno RMS 1.46 arcsec, 90
osservazioni su 16 anni danno 2.03 (non pesato).

*Verifica:* una configurazione con `fit.enabled: true` produce un fit riproducibile
sul caso BK290.

---

## Passo 2 — Innesto nell'orchestrator  [medio]

Il fit gira **solo sui positivi** del secondo stadio, che sono poche unita'.

Sequenza per ciascun positivo:
1. scarico del `.rwo` (gia' implementato, Fase 3);
2. filtro delle osservazioni secondo `fit.obs_years`;
3. `fit_orbit()` a partire dagli elementi AstDyS;
4. se converge -> si usa l'orbita fittata per la predizione;
   se non converge -> si mantiene l'orbita AstDyS e si segnala nel report.

Il punto 4 e' importante: un fit che non converge non deve far fallire la
campagna ne' — peggio — sostituire silenziosamente un'orbita buona con una
peggiore. Il criterio di accettazione va esplicitato: `converged == true` e RMS
non superiore a quello di partenza.

*Verifica:* campagna su BK290 con `fit.enabled: true` che completa e riporta
l'esito del fit per ciascun positivo.

---

## Passo 3 — Covarianza dal fit  [medio, alto valore]

Oggi la covarianza viene scaricata da AstDyS (Fase 3). Il fit puo' produrla
direttamente: `compute_covariance_internal` esiste gia' nel DifferentialCorrector.

Questo chiude il cerchio: orbita **e** incertezza calcolate da noi, con
l'ellisse di predizione derivata dalla nostra soluzione invece che importata.

*Verifica:* l'ellisse ottenuta dal fit va confrontata con quella AstDyS sullo
stesso oggetto. Un accordo entro qualche percento conferma tutta la catena;
uno scarto grande indica un problema nei pesi o nella normalizzazione.

Attenzione: e' proprio il tipo di confronto che oggi ha rivelato i bug. Va fatto
prima di dichiarare il passo concluso, non dopo.

---

## Passo 4 — Prestazioni in campagna  [da misurare]

Il fit su 90 osservazioni richiede alcuni secondi. Su una campagna con decine di
positivi il costo va misurato, non stimato.

Leve disponibili, in ordine di efficacia attesa:
- `fit.obs_years` per accorciare l'arco;
- tolleranza del fit meno stretta di quella di propagazione;
- parallelizzazione sui positivi (sono indipendenti fra loro).

*Verifica:* tempo totale di una campagna reale con e senza fit.

---

## Passo 5 — Documentazione  [piccolo]

Capitolo di manuale sul fit: quando serve, cosa aspettarsi, come leggere RMS e
numero di osservazioni rigettate, il compromesso di `obs_years`.

Da chiarire esplicitamente che il nostro `rms_total` NON e' confrontabile con
l'`RMSast` pubblicato da AstDyS: il primo e' la radice della media dei residui in
arcsec, il secondo e' normalizzato sui sigma delle osservazioni. Confrontarli
direttamente e' fuorviante — errore commesso durante la diagnosi del 2026-07-22.

---

## Ordine consigliato

1 -> 2 danno subito un fit utilizzabile in campagna.
3 e' il passo di maggior valore scientifico (indipendenza da AstDyS per
l'incertezza) e conviene affrontarlo quando 1 e 2 sono stabili.
4 e 5 a chiusura.

## Rischi noti

- **Fit che peggiora l'orbita.** Con pochi dati o arco corto il fit puo' produrre
  una soluzione peggiore di quella AstDyS. Il criterio di accettazione del passo 2
  e' la protezione; va verificato su casi reali, non solo su BK290.
- **Osservatori mobili.** Il codice 270 (Unistellar) e altri non hanno coordinate
  in catalogo e cadono sul geocentro: sugli oggetti vicini l'effetto e' rilevante
  (misurato su Eros: residui 4.05 arcsec contro 0.42 di AstDyS). Per gli asteroidi
  di fascia principale l'impatto e' minore, ma va tenuto presente.
- **Oggetti con pochi dati.** Il comportamento del fit su archi di poche settimane
  o poche decine di osservazioni non e' stato provato.
