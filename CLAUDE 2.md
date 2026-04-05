## Stack tecnologico
- **Linguaggio**: C++20  ← non C++17

## Principi CTFYH (Code That Fits in Your Head)

### Regole strutturali
- Max 7 oggetti per scope cognitivo (Miller's Law) — funzioni, variabili locali,
  parametri: se superi 7, decomponi
- Funzioni max ~24 righe (una schermata senza scroll)
- Classi con una sola responsabilità misurabile — se il nome contiene "And", splitta
- Max 3 livelli di indentazione — oltre: estrai funzione o usa early return
- Parametri funzione: max 3-4; oltre, usa un struct o builder

### Regole di leggibilità
- Il nome della funzione deve descrivere COSA fa, non COME
- Nessun commento che spiega il "come" — il codice lo mostra già;
  i commenti spiegano solo il "perché" (decisioni non ovvie)
- No abbreviazioni nei nomi: `sep` → `separation_arcsec`, `dyn` → `dynamics`
- Variabili booleane con prefisso: `is_`, `has_`, `should_`

### Regole per AstDyn specificamente
- Ogni funzione che tocca unità fisiche deve avere il suffisso dell'unità
  nel nome della variabile: `dist_au`, `vel_km_s`, `gm_km3_s2`
- Il funtore Dynamics non supera 60 righe — se cresce, estrai metodi privati
  `solar_acceleration()`, `nbody_acceleration()`, `relativistic_correction()`
- Nessun magic number: `86400.0` → `constants::SECONDS_PER_DAY`
- Ogni nuovo force model viene aggiunto sia in ForceField che in
  MultiBodyPropagator::Dynamics (checklist nel PR)



### Checklist pre-commit
- [ ] La funzione ha un nome che descrive il risultato, non il processo?
- [ ] Il funtore Dynamics usa ForceField o replica logica già esistente?
- [ ] Le unità fisiche sono esplicite nel nome delle variabili locali?
- [ ] Ci sono magic numbers da sostituire con costanti?
- [ ] Il diff tocca MultiBodyPropagator? → verifica anche RelativeMultiBodyPropagator
```

Con questo setup, ogni volta che apri Claude Code nella directory del progetto, ha automaticamente:
- il contesto architetturale
- i bug noti aperti
- i principi CTFYH come vincoli di stile attivi

Puoi poi rinforzarli nella sessione con un prompt iniziale tipo:
```
Rispetta sempre i principi CTFYH definiti in CLAUDE.md.
Prima di proporre qualsiasi modifica, verifica che il codice risultante
rispetti i limiti di complessità cognitiva definiti lì.sicuro che sia in c++17 o. in c++2322:58Ottima domanda. Guardo il CMakeLists.txt.Eseguiti 8 comandiEseguiti 8 comandiRisposta certa: C++20, non C++17 né C++23.
La prova è in frame_rotation.hpp:
cpp#include <concepts>                          // header C++20
concept IsRotatable = ...;                   // syntax C++20
requires IsRotatable<Frame, TargetFrame>     // requires-clause C++20