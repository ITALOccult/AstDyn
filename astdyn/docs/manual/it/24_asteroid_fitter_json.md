# AsteroidFitter – Configurazione JSON

## Introduzione
Il nuovo **AsteroidFitter** supporta la configurazione tramite file JSON, consentendo di specificare in modo flessibile:
- Percorsi locali o URL per i file OrbFit (`.eq1` e `.rwo`).
- Parametri opzionali di configurazione (`.oop`).
- Elementi orbitali e epoche di osservazione **in‑memoria** quando non sono disponibili file.

Questa documentazione descrive lo **schema JSON**, le funzioni di caricamento e l’uso dell’API `fitFromConfig`.

---

## Schema JSON
Il file di configurazione è un oggetto JSON con i seguenti campi (tutti opzionali, tranne quelli richiesti per il caso d’uso):

```json
{
  "eq1_file": "path/to/file.eq1",          // Percorso locale (stringa, opzionale)
  "rwo_file": "path/to/file.rwo",          // Percorso locale (stringa, opzionale)
  "oop_file": "path/to/config.oop",        // Configurazione AstDyn (opzionale)
  "eq1_url": "https://example.com/file.eq1", // URL da scaricare (stringa, opzionale)
  "rwo_url": "https://example.com/file.rwo", // URL da scaricare (stringa, opzionale)
  "orbit": {                                 // Elementi orbitali (obbligatori se non si usano file)
    "a": 2.5,                               // Semiasse maggiore (AU)
    "e": 0.1,                               // Eccentricità
    "i": 0.05,                              // Inclinazione (rad)
    "Omega": 1.0,                           // Nodo ascendente (rad)
    "omega": 0.5,                           // Argomento del pericentro (rad)
    "M": 0.0                                 // Anomalia media (rad)
  },
  "mjd_observations": [61000.0, 61001.5, 61003.2], // Epoche di osservazione (MJD, UTC)
  "outputEquatorial": true                 // true = Equatorial J2000, false = Ecliptic J2000
}
```

- **eq1_file / rwo_file**: se specificati, il fitting avviene usando i file locali.
- **eq1_url / rwo_url**: se i percorsi locali sono vuoti, il file viene scaricato in una directory temporanea tramite `curl`.
- **orbit + mjd_observations**: usati quando non sono presenti file; il metodo `fitFromConfig` chiama direttamente `computeFromMemory`.
- **outputEquatorial**: definisce il frame di output delle posizioni.

---

## Funzione di caricamento
Il file `AsteroidFitConfig.hpp` fornisce la funzione:

```cpp
astdyn::ephemeris::AsteroidFitConfig
loadAsteroidFitConfig(const std::string& json_path);
```
Questa funzione legge il file JSON, popola la struttura `AsteroidFitConfig` e gestisce i valori di default (stringhe vuote, `outputEquatorial = true`).

---

## Utilizzo dell’API `fitFromConfig`
```cpp
#include "astdyn/ephemeris/AsteroidFitConfig.hpp"
#include "astdyn/ephemeris/AsteroidFitter.hpp"

int main(){
    // 1. Carica la configurazione JSON
    auto cfg = astdyn::ephemeris::loadAsteroidFitConfig("config.json");

    // 2. Esegue il fitting (file, URL o in‑memory)
    auto result = astdyn::ephemeris::AsteroidFitter::fitFromConfig(cfg);

    // 3. Verifica il risultato
    if(result.success){
        std::cout << "Fit completato con successo!" << std::endl;
        // accesso a result.fitted_positions, result.fitted_orbit, ecc.
    } else {
        std::cerr << "Fit fallito: " << result.message << std::endl;
    }
    return 0;
}
```
Il metodo gestisce automaticamente:
1. **Download** dei file se sono forniti URL.
2. **Fallback** a calcolo in‑memory quando i file non sono disponibili.
3. **Propagazione** delle posizioni con `PositionCalculator`.
4. **Restituzione** di un `AsteroidFitResult` completo di posizioni, RMS e messaggi.

---

## Note sul download dei file
- Il download è effettuato con `curl -L -s -o <temp_path> <url>`.
- Il file temporaneo viene creato in `std::filesystem::temp_directory_path()` con un nome univoco basato su hash dell’URL.
- In caso di errore (`curl` restituisce codice diverso da 0) viene lanciata un’eccezione `std::runtime_error` e il risultato `AsteroidFitResult` con `success = false` e `message` contenente il motivo.

---

## Gestione degli errori
- **File non trovati** (locale o download fallito): il risultato contiene `success = false` e un messaggio descrittivo.
- **JSON malformato**: la funzione `loadAsteroidFitConfig` lancia `std::runtime_error` con il motivo.
- **Mancanza di dati**: se né file né elementi orbitali sono forniti, il risultato è fallito con messaggio “No input data provided”.

---

## Test consigliati
1. **File locale** – fornire `eq1_file` e `rwo_file` validi e verificare che `fitFromConfig` chiami `fit`.
2. **URL** – impostare `eq1_url` e/o `rwo_url` a risorse pubbliche e controllare che il download avvenga e il fitting sia eseguito.
3. **In‑memory** – omettere tutti i percorsi e fornire `orbit` + `mjd_observations`; verificare che le posizioni vengano calcolate con `computeFromMemory`.
4. **Error handling** – fornire URL non raggiungibili o JSON incompleti e verificare che il messaggio di errore sia chiaro.

---

## Integrazione nella documentazione del progetto
Il file dovrebbe essere inserito nella sezione **Manuale Italiano** sotto `astdyn/docs/manual/it/`. Si consiglia di aggiungere il riferimento in `main_it.tex`:

```tex
\include{24_asteroid_fitter_json}
```

---

## Conclusioni
Questa documentazione fornisce tutti gli elementi necessari per utilizzare la nuova API basata su JSON, garantendo flessibilità nella gestione dei dati di input e semplificando l’integrazione nei workflow di AstDyn.
