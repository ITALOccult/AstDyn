# TODO consolidato — Fase 5 (fit) e lavori aperti (2026-07-22)

## FASE 5 — FIT ORBITALE (il fardello: da sistemare prima o poi)

### Stato: il fit NON e' rotto nella logica, ma ha un collo di bottiglia strutturale
Diagnosi completa di questa sessione (3 cause distinte, tutte identificate):

1. [RISOLTO - era errore di test] Il "fit rotto" iniziale era in gran parte un bug
   del test di diagnosi: from_traditional(epoch, a, e, i, node, omega, M) vuole gli
   angoli in GRADI, non radianti. Passando radianti, l'orbita era sbagliata -> RMS
   3508 arcsec. Corretto -> il modello di osservazione predice giusto:
   computed_ra=145.36 == JPL == osservazione. Il fit e' SANO.

2. [RISOLTO - committato 0190e6e] integrate_at (propagazione batch, target multipli)
   aveva un bug: il passo di CENTRAGGIO di un target, se rifiutato da adaptive_step,
   non veniva ridotto -> loop di 1000 rifiuti -> eccezione. Fix v3: il centraggio
   rifiutato viene ridotto e riusato piu' piccolo, senza degradare il passo di lavoro.
   Validato: caso C (30 target) passa; 170/170 test verdi.

3. [DA RIVEDERE COMPLETAMENTE] La STM (StateTransitionMatrix) e' il collo di bottiglia
   che ancora BLOCCA il fit su archi lunghi (>~1-2 anni):
   - propaga un sistema 42-dim (6 stato + 36 STM) con RKF78 @ tolleranza 1e-13
     (strettissima, ingiustificata per le equazioni variazionali).
   - la jacobiana e' NUMERICA per default (use_numerical_jacobian_=true, "OrbFit style"):
     12 compute_derivatives per ogni valutazione della jacobiana.
   - RKF78 fa 13 valutazioni/step -> ~169 compute_derivatives per step sul sistema 42-dim.
   - su archi lunghi il costo esplode: il fit "si blocca" (macina, non va in loop).
   - PROVATO senza successo: abbassare tolleranza STM a 1e-9 (non basta); forzare
     jacobiana analitica (use_numerical_jacobian_=false) (non basta da solo).
   - CONCLUSIONE: la STM va RIPENSATA, non e' un parametro. Opzioni da valutare:
     a) jacobiana analitica completa (gravita' + perturbazioni), non solo gravita';
     b) propagare STM e stato in UN unico integratore condiviso (un solo passo),
        invece di due propagazioni separate;
     c) studiare come fa OrbFit (clonato in /home/claude/orbfit come oracolo):
        propag.f90, pred_obs.f90 - OrbFit fitta archi di decenni velocemente,
        quindi NON ricalcola jacobiana numerica su tutto l'arco a ogni step cosi';
     d) integratore dedicato per il sistema variazionale (piu' adatto del RKF78 esplicito).

### Quando la STM sara' sistemata, completare la Fase 5:
- collegare fit_orbit all'orchestrator come parametro di calcolo (fit: true nel
  second_stage sui positivi);
- filtro osservazioni "ultimi N anni" (obs_years, opzionale) - gia' predisposto nel
  test test_fit_bk290. Tradeoff documentato: arco corto = fit veloce ma orbita meno
  precisa (l'accuratezza dipende dalla baseline temporale delle osservazioni).
- tolleranza fit configurabile (fit_tolerance);
- scarico osservazioni .rwo per i positivi (gia' implementato in orchestrator Fase 3);
- validare: fit da orbita AstDyS -> RMS ~0.46 arcsec (RMSast dell'header .rwo),
  converged=YES. Oracolo JPL disponibile (vedi test_validate_integrators).

## TODO — REVISIONE INTEGRATORI (validazione vs JPL)
Test pronto: examples/test_validate_integrators.cpp (BK290 dall'epoca vs JPL, 4 archi).
Risultati sessione (oracolo JPL a MJD interi esatti):
- RKF78 @1e-12: [OK] tutti gli archi (err ~2e-6 AU), VELOCISSIMO (28ms/10anni). SCELTO.
- AAS @1e-12: [OK] tutti gli archi, ma ~40x piu' lento (1190ms/10anni). Riserva robusta.
- RKF78 @1e-10: [OK]/[~], degrada un po' su archi lunghi ma velocissimo.
- SABA4: ROTTO. err ~5 AU (diverge) + lentissimo (~80s/arco). PRIORITA' ALTA.
- RADAU: troppo lento (interrotto). Da profilare o marcare "solo stiff".
- GAUSS, GRKN64, RK4: NON ancora testati. Completare.
Migliorare il test: aggiungere propagazione in AVANTI e round-trip (reversibilita');
escludere SABA4/RADAU dai run rapidi o dare timeout.
LEZIONE: usare oracolo JPL a MJD INTERI esatti (mezzo giorno di disallineamento
introduce ~5e-3 AU di errore spurio - ci ha depistati a lungo).

## TODO — PAPER JOSS (Journal of Open Source Software) a fine lavoro
- JOSS accetta paper.md in Markdown con metadati YAML (formato preferito).
- Gratuito, open-access, peer review sul software. DOI CrossRef all'accettazione.
- Descrive il SOFTWARE (AstDyn + ioccultcalc), non i risultati scientifici
  (quelli -> paper RASTI sull'integratore AAS). Due paper complementari.
- Requisiti: OSI license, repo pubblico (ITALOccult/AstDyn OK), documentazione,
  test automatici (170 test + validazione vs JPL/AstDyS OK). 750-1750 parole.
- Include AI usage disclosure (richiesto da JOSS).
- Sezioni: Summary, Statement of need, State of the field, Software design,
  Research impact, References. Evidenziare: two-stage search, covarianza
  per-asteroide -> ellissi, sistema frame tipizzato C++20, validazione vs AstDyS/JPL.
