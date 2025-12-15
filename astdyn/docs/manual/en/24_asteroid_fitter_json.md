# AsteroidFitter – JSON Configuration

## Introduction
The new **AsteroidFitter** supports configuration via a JSON file, allowing flexible specification of:
- Local paths or URLs for OrbFit files (`.eq1` and `.rwo`).
- Optional engine configuration (`.oop`).
- Orbital elements and observation epochs **in‑memory** when no files are provided.

This documentation describes the JSON schema, loading function, and usage of the `fitFromConfig` API.

---

## JSON Schema
The configuration file is a JSON object with the following fields (all optional, except those required for the chosen use case):

```json
{
  "eq1_file": "path/to/file.eq1",          // Local path (string, optional)
  "rwo_file": "path/to/file.rwo",          // Local path (string, optional)
  "oop_file": "path/to/config.oop",        // AstDyn configuration (optional)
  "eq1_url": "https://example.com/file.eq1", // URL to download (string, optional)
  "rwo_url": "https://example.com/file.rwo", // URL to download (string, optional)
  "orbit": {                                 // Orbital elements (required if no files are used)
    "a": 2.5,                               // Semi‑major axis (AU)
    "e": 0.1,                               // Eccentricity
    "i": 0.05,                              // Inclination (rad)
    "Omega": 1.0,                           // Longitude of ascending node (rad)
    "omega": 0.5,                           // Argument of periapsis (rad)
    "M": 0.0                                 // Mean anomaly (rad)
  },
  "mjd_observations": [61000.0, 61001.5, 61003.2], // Observation epochs (MJD, UTC)
  "outputEquatorial": true                 // true = Equatorial J2000, false = Ecliptic J2000
}
```

- **eq1_file / rwo_file**: If provided, fitting is performed using the local files.
- **eq1_url / rwo_url**: If the local paths are empty, the file is downloaded to a temporary location using `curl`.
- **orbit + mjd_observations**: Used when no files are supplied; `fitFromConfig` calls `computeFromMemory` directly.
- **outputEquatorial**: Determines the output reference frame for the propagated positions.

---

## Loading Function
`AsteroidFitConfig.hpp` provides the helper:

```cpp
astdyn::ephemeris::AsteroidFitConfig
loadAsteroidFitConfig(const std::string& json_path);
```
It parses the JSON file, fills the `AsteroidFitConfig` structure, and applies sensible defaults (empty strings, `outputEquatorial = true`).

---

## Using the `fitFromConfig` API
```cpp
#include "astdyn/ephemeris/AsteroidFitConfig.hpp"
#include "astdyn/ephemeris/AsteroidFitter.hpp"

int main(){
    // 1. Load the JSON configuration
    auto cfg = astdyn::ephemeris::loadAsteroidFitConfig("config.json");

    // 2. Run the fitting (file, URL or in‑memory)
    auto result = astdyn::ephemeris::AsteroidFitter::fitFromConfig(cfg);

    // 3. Check the result
    if(result.success){
        std::cout << "Fit completed successfully!" << std::endl;
        // Access result.fitted_positions, result.fitted_orbit, etc.
    } else {
        std::cerr << "Fit failed: " << result.message << std::endl;
    }
    return 0;
}
```
The method automatically:
1. **Downloads** files when URLs are supplied.
2. **Falls back** to in‑memory computation when no files are present.
3. **Propagates** positions using `PositionCalculator`.
4. **Returns** an `AsteroidFitResult` containing positions, RMS values, and status messages.

---

## Download Details
- Performed with `curl -L -s -o <temp_path> <url>`.
- The temporary file is created in `std::filesystem::temp_directory_path()` with a unique name derived from a hash of the URL.
- If `curl` returns a non‑zero exit code, a `std::runtime_error` is thrown and the result contains `success = false` with an explanatory message.

---

## Error Handling
- **Missing files / download failure**: `AsteroidFitResult.success = false` and `message` explains the problem.
- **Malformed JSON**: `loadAsteroidFitConfig` throws `std::runtime_error` with details.
- **No input data**: If neither files nor orbital elements are provided, the result is failed with the message “No input data provided”.

---

## Recommended Tests
1. **Local files** – Provide valid `eq1_file` and `rwo_file` and verify that `fitFromConfig` invokes the file‑based `fit`.
2. **URL download** – Set `eq1_url` and/or `rwo_url` to reachable resources and confirm that the download occurs and fitting proceeds.
3. **In‑memory** – Omit all file fields, supply `orbit` and `mjd_observations`, and check that positions are computed via `computeFromMemory`.
4. **Error cases** – Use unreachable URLs or malformed JSON and verify that clear error messages are returned.

---

## Integration into the Project Manual
Place this file under `astdyn/docs/manual/en/` and reference it in the English main LaTeX file (`main_en.tex`):

```tex
\include{24_asteroid_fitter_json}
```

---

## Conclusion
The documentation above provides all the information needed to use the new JSON‑based configuration API for `AsteroidFitter`, offering flexibility in data input and simplifying integration into user workflows.
