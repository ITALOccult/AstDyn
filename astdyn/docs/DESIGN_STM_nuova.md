# Proposta di redesign della STM (ispirata a OrbFit)

## Cosa fa OrbFit (dall'esame dell'oracolo /home/claude/orbfit)

1. **Jacobiana ANALITICA, non numerica.**
   In force9d.f90 la matrice delle derivate seconde del potenziale d2udx(3,3) e'
   calcolata con formule esplicite chiuse (righe 347-358: termini in x2,y2,z2,cj2r9...),
   incluse gravita' + J2 + relativita'. Zero differenze finite. La variazionale usa
   direttamente d2U/dx2 analitica.

2. **Stato + variazionali in UN vettore, UN integratore, UN passo.**
   Tutto e' impacchettato in y2(nvar) con nvar=6*(nbod-1+na+nvz); reord/stack
   gestiscono il layout. Lo stato dell'asteroide e le sue equazioni variazionali sono
   propagati INSIEME, dallo stesso integratore, nello stesso passo (una sola CALL propag
   restituisce xastr E dxdpar). Non due propagazioni separate.

3. **Integratore scelto da CONFIG, non fisso.**
   imet in propag_state.f90: =1 multistep, =2 runge-kutta, =3 everhart (RA15).
   Il metodo di integrazione e' un parametro, non hardcoded.

## I nostri 3 problemi attuali (rispetto a OrbFit)

| Aspetto            | Nostro (rotto)                          | OrbFit (riferimento)          |
|--------------------|-----------------------------------------|-------------------------------|
| Jacobiana          | NUMERICA (12 compute_derivatives/valut.)| ANALITICA (formule chiuse)    |
| Integratore STM    | RKF78 dedicato HARDCODED @ 1e-13        | stesso integratore, da config |
| Struttura          | 2 propagazioni separate (stato / STM)   | 1 sistema integrato insieme   |

Effetto combinato: su archi lunghi il costo esplode (sistema 42-dim x jacobiana FD),
il fit "si blocca". E' un problema STRUTTURALE, non di parametro.

## Proposta di design nuovo

### Principio 1 — Jacobiana analitica come default
- Completare/usare il ramo analitico gia' presente (use_numerical_jacobian_=false):
  A = d(f)/d(x,v) con
    - blocco d(accel)/d(pos) = -GM/r^3 * (I - 3 rr^T/r^2)  (two-body, forma chiusa)
    - + contributi analitici dei perturbatori planetari principali
    - + J2 se attivo (come d2udx di OrbFit)
- Verificare il ramo analitico esistente contro la jacobiana numerica su un caso
  noto (devono coincidere entro ~1e-8). Se il ramo analitico e' incompleto (solo
  two-body), completarlo con i perturbatori dominanti.
- La jacobiana numerica resta disponibile come fallback/validazione (flag), ma NON
  e' piu' il default sul percorso caldo.

### Principio 2 — Integratore della STM da CONFIG (tua indicazione)
- RIMUOVERE l'RKF78 hardcoded @1e-13 dal costruttore di StateTransitionMatrix.
- La STM deve usare l'integratore scelto in AstDynConfig (integrator_type, tolerance),
  lo STESSO del propagatore, coerente con il resto del sistema.
- Cosi' se l'utente sceglie AAS o RKF78 @1e-10, la STM segue. Niente due tolleranze
  scollegate. La tolleranza per le variazionali puo' avere un DEFAULT piu' morbido
  (le derivate parziali non richiedono la precisione dello stato) ma resta derivata
  dalla config, non un magic number 1e-13.

### Principio 3 — Stato + variazionali propagati INSIEME (gia' parzialmente fatto)
- Gia' oggi propaghiamo il sistema aumentato 42-dim con f_augmented in un solo
  integrate_at: questo aspetto e' gia' allineato a OrbFit. Bene.
- Il problema NON e' qui: e' la jacobiana numerica dentro f_augmented + la tolleranza
  1e-13 hardcoded. Risolti i punti 1 e 2, questo terzo e' gia' a posto.

## Passi concreti proposti
1. Verificare il ramo analitico esistente di compute_jacobian (correttezza vs numerica).
2. Se incompleto, completarlo (two-body chiuso + perturbatori principali; J2 se attivo).
3. Rendere l'integratore della STM guidato da AstDynConfig (togliere RKF78 @1e-13 fisso);
   tolleranza variazionali derivata dalla config, con default ragionevole.
4. Ri-test: fit BK290 arco pieno (~15 anni di osservazioni) -> deve completare e
   convergere a RMS ~0.46 arcsec (RMSast). Confronto con AstDyS/JPL.
5. 170 test verdi.
6. Poi innesto Python (Fase 5): fit come parametro nel second-stage sui positivi.

## Nota
OrbFit resta l'oracolo: force9d.f90 (jacobiana analitica d2udx), propag_state.f90
(scelta integratore da config imet). Consultare per le formule esatte dei perturbatori.
