Per rispondere alla richiesta [M3] del Referee ("No numerical experiments") e validare l'integratore AAS su Celestial Mechanics and Dynamical Astronomy (CeMDA), l'agente deve generare dati strutturati per quattro grafici fondamentali. Ecco le istruzioni dettagliate da dare all'agente (Claude Code in Antigravity) per preparare la suite di test e i file di output.

**NOTA DI AGGIORNAMENTO (17/03/2026)**: Abbiamo risolto il problema del "tailspin" numerico negli integratori Gauss e SABA4 negli scenari con perturbazioni planetarie. Il benchmark ora può essere eseguito in modo robusto su archi temporali lunghi senza che l'integratore tenti di "correggere" la non-conservazione fisica dell'energia tramite riduzioni eccessive del passo.

Istruzioni per l'Agente: Protocollo di Validazione Numerica [M3]
Obiettivo: Generare file CSV strutturati per la creazione dei grafici di validazione scientifica dell'integratore AAS.

1. Test di Convergenza (Convergence Plot)
Target: Un'orbita di prova (es. Ceres per 1 periodo orbitale).
Metodo: Eseguire la propagazione variando il parametro precision_ ($\varepsilon$) in un range logaritmico: $[10^{-2}, 10^{-3}, 10^{-4}, 10^{-5}, 10^{-6}, 10^{-7}, 10^{-8}]$.
Dato da salvare: Errore di posizione relativo finale rispetto a un riferimento ad alta precisione (es. JPL Horizons o un'integrazione IAS15 a tolleranza $10^{-16}$).
Output CSV: convergence_test.csv (Colonne: epsilon, relative_pos_error, nfe [Number of Force Evaluations]).

2. Confronto Efficienza (Work-Precision Diagram)
Target: Due scenari: (A) Ceres ($e=0.07$) e (B) 1566 Icarus ($e=0.82$).
Metodo: Confrontare AAS con RK4 (baseline) e IAS15 (se disponibile).
Dato da salvare: Errore di posizione finale in funzione del numero totale di valutazioni della funzione forza (stats_.num_function_evals).
Output CSV: work_precision.csv (Colonne: integrator_name, nfe, final_error, eccentricity).

3. Conservazione dell'Energia a Lungo Termine (Long-term Energy)
Target: Ceres o un asteroide della fascia principale per $10.000$ orbite.
Metodo: Utilizzare una precisione fissa (es. $\varepsilon = 10^{-5}$).
Dato da salvare: Ad ogni step (o ogni $N$ step), registrare il tempo, l'errore dell'Hamiltoniana ($\delta H$) e l'errore dell'Hamiltoniana Ombra ($\delta \tilde{H}$).
Output CSV: energy_conservation.csv (Colonne: time_years, delta_H, delta_H_shadow). Nota: Ci aspettiamo che $\delta H$ oscilli e $\delta \tilde{H}$ rimanga piatto.

4. Validazione della STM (Symplecticity & Accuracy)
Target: Qualsiasi asteroide NEA per 1 anno.
Metodo: Calcolare la STM $\Phi(t)$ con AAS.
Dati da salvare:
L'errore del determinante: $|\det(\Phi) - 1|$.
L'errore rispetto a una STM calcolata con differenze finite (FD) sullo stato iniziale.
Output CSV: stm_validation.csv (Colonne: time, det_error, diff_with_FD).

Formato dei Risultati (Standard richiesto)
L'agente deve salvare i file in una cartella /benchmarks con la seguente struttura:
Header chiaro: Ogni CSV deve avere un header con i parametri fisici usati ($\mu, J_2$).
Unità di misura: Tempo in Giorni Giuliani o Anni, distanze in Metri o AU.
Statistiche extra: Includere sempre num_steps e num_function_evals per ogni riga di test.

Script di Visualizzazione (Python/Matplotlib)
Chiedi all'agente di generare anche uno script plot_results.py che:
Usi la scala Log-Log per i grafici di convergenza e work-precision.
Usi la scala Semi-Log (Y logaritmica) per il grafico dell'energia.
Verifichi che la pendenza del grafico di convergenza sia esattamente -4.