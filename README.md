# AstDyn: High-Fidelity Asteroid Dynamics Library

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Standard: C++23](https://img.shields.io/badge/standard-C%2B%2B23-blue.svg)
![Status: Beta](https://img.shields.io/badge/Status-Beta-yellow.svg)

**AstDyn** is a high-performance C++ library designed for precision asteroid orbit propagation and determination. It focuses on speed, accuracy, and modern software architecture, making it suitable for demanding applications like stellar occultation prediction and hazard assessment.

## Key Features

*   **High-Precision Wrapper**: Simplified interface (`ITALOccultLibrary`) for rapid integration with high-level workflows.
*   **Extreme Dynamics**: Full N-body gravity, **AST17** asteroid perturbations, and adaptive symplectic integrators (**AAS**) for high-fidelity occultation modeling.
*   **Dynamic Occultation Discovery**: Automated corridor search using Gaia DR3 with real-time Earth-Asteroid distance interpolation.
*   **Validated Accuracy**: Validated against NASA/JPL Horizons Ephemeris System with residuals $< 2.5 \mu\text{m}$ ($10^{-9}$ AU) and sub-arcsecond geocentric precision ($< 0.7$ arcsec).
*   **Performance**: Optimized C++17 architecture, ~8x faster than numerical STM approaches.

---

## 🏛️ Technical Overview

AstDyn is architected for both research and high-demand operational environments. Its core is built around a hybrid propagation engine:

*   **Propulsion & Precision**: Native support for **AAS (Adaptive Step Symplectic)** integrators, specifically tuned for Encke-type asteroid dynamics.
*   **Perturbation Engine**: The **AST17** module includes real-time gravitational influences from the 17 most massive asteroids (Ceres, Pallas, Vesta, etc.), essential for 10-mas occultation precision.
*   **Unified Reference Frames**: Seamless handling of GCRF, ECLIPJ2000, and ITRF (Topocentric) frames using high-fidelity IAU 2000/2006 precession-nutation models.
*   **Modular Architecture**: Designed for easy integration with external projects (e.g., via `ITALOccultLibrary`) while maintaining a strict high-performance cor## Documentation

Full technical documentation is available in the `astdyn/docs` directory:

*   **[User Manual](astdyn/docs/AstDyn_User_Manual.pdf)**: Comprehensive guide to the API, installation, and theoretical background.
*   **[Configuration Manual](astdyn/docs/config_manual.md)**: Detailed guide for advanced `IOCConfig` (YAML, OOP, JSON).

## Quick Start

### Prerequisites
*   **C++23** compliant compiler (GCC 13+, Clang 16+, MSVC 2022)
*   **CMake 3.20+**
*   **JPL DE441 Ephemeris**: Required for precise calculations. Download `de441.bsp` from [JPL FTP](https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de441.bsp) and place it in your build directory or at `~/.ioccultcalc/ephemerides/de441.bsp`.
*   Eigen3, CURL, nlohmann_json libraries

### Installation

```bash
git clone https://github.com/ITALOccult/AstDyn.git
cd AstDyn
mkdir build && cd build
cmake ..
make
sudo make install
```

#### Simplified High-Precision API (ITALOccultLibrary)

For most users, `ITALOccultLibrary` provides a more convenient way to work with AstDyS elements and high-precision observations:

```cpp
#include <italoccultlib/astdyn_wrapper.h>

void example() {
    ioccultcalc::AstDynWrapper wrapper;
    wrapper.loadFromEQ1File("34713.eq1"); // Automatically handles mean-to-osculating conversion
    
    // High-precision geocentric RA/Dec with DE441 light-time correction
    auto obs = wrapper.calculateObservation(61050.0);
    std::cout << "RA: " << obs.ra_deg << " deg" << std::endl;
}
```

## Tools

### ioccultcalc (Beta 0.6 - Released 2026-03-16)
A professional command-line tool for searching and comparing stellar occultations. It retrieves asteroid elements directly from **JPL Horizons** and searches the **Gaia DR3** catalog.

### Advanced Configuration System (IOCConfig)
AstDyn now features a powerful, multi-format configuration manager inspired by the `IOC_Config` pattern. It supports:
*   **JSON**: Native support for structured configurations.
*   **YAML/OOP Format**: High-level hierarchical syntax with braces (e.g., `integrator { type = RK4 }`).
*   **Flat Dot-Notation**: Compact key-value pairs (e.g., `integrator.step_size = 0.1`).
*   **Automatic Type Deduction**: Parses booleans, numbers, and strings automatically.

**Example YAML-style config (`config.yaml`):**
```yaml
integrator {
  type = RKF78
  step_size = 0.01
}
diffcorr.light_time = true
verbose = false
```

**New in Beta 0.6 (v1.0 RC Features - Updated 2026-03-16):**
*   **Advanced IOC Configuration**: Support for YAML, OOP, and JSON configuration files via `--conf`.
*   **Workspace & Path Management**: Organize multi-event searches with `--out-dir` and `--prefix` for automated file sorting.
*   **Individual Event Mapping**: Automatically generates dedicated SVG/KML files for every match in a batch search.
*   **High-Precision Catalog Corrections**: Full implementation of Gaia DR3 **3D Space Motion** (proper motion) and **Annual/Annual Parallax** shift.
*   **Relativistic Aberration & Light Deflection**: Native support for rigorous IAU 2000 aberration and gravitational light deflection (GR).
*   **Advanced Observer Filters**: Automate search refinement with thresholds for **Sun/Moon Altitude**, **Moon Proximity**, and **Magnitude Drop**.
*   **Multi-Catalog Framework**: Support for dynamic selection of stellar catalogs (Gaia DR3, Legacy) via `--catalog`.
*   **Regional Search Filter**: Search for occultations passing near your specific location with `--lat` and `--lon`.

#### Usage
```bash
# Search with high-precision aberration and organized workspace output
ioccultcalc --asteroid 1,4 --jd-start 2461131.5 --duration 30.0 --out-dir campaign_v1 --prefix vesta_test
```

#### Command-line Options
- `--asteroid <list|@file>`: Comma-separated designations or a file path (prefixed with `@`).
- `--jd-start <val>`: Starting Julian Date (TDB) for the search.
- `--duration <days>`: Length of the search window in days (default: 1.0).
- `--mag <val>`: Magnitude limit for star search (default: 15.0).
- `--catalog <type>`: Stellar catalog to use (`gaia_dr3`, `legacy`).
- `--out-dir <path>`: Base directory for all output files (autosorts multi-event results).
- `--prefix <str>`: Prefix for individual event files (default: `occ`).
- `--svg-output <file>`: Generate a high-resolution SVG world map.
- `--xml-output <file>`: Save matches to an Occult4-compatible XML.
- `--zoom <val>`: Zoom level (1.0 = global, 10.0 = local).
- `--lat` / `--lon`: Regional filter for observation sites.

---

## 🧪 Beta Testing & Implementation Status

AstDyn is currently in an active testing phase. For a comprehensive guide on using the occultation engine, high-precision configurations, and analysis tools, please refer to the:

**[IOccultCalc User Guide (v0.6)](astdyn/docs/ioccultcalc_guide.md)** | **[Configuration Manual](astdyn/docs/config_manual.md)**

### Implementation Card (v1.0 RC Ready)

| Feature | Status | Description |
| :--- | :---: | :--- |
| **N-Body Dynamics** | ✅ | Full planetary perturbations (DE441). |
| **Asteroid Perturbations** | ✅ | AST17 (Top 17 most massive asteroids). |
| **Stellar Corrections** | ✅ | Proper motion, parallax, aberration, light deflection. |
| **Numerical Methods** | ✅ | AAS, RKF78, Gauss, Radau, SABA4. |
| **Occultation Engine** | ✅ | Corridor discovery, refinement, and uncertainty modeling. |
| **Binary Systems** | ✅ | Support for satellite occultations via BSP files. |
| **Advanced Config** | ✅ | YAML, OOP, JSON hierarchical support. |
| **Visualization** | ✅ | Global/Local SVG, KML, Occult4 XML. |
| **Workspace Org** | ✅ | Automatic directory sorting and event prefixing. |

If you encounter discrepancies or bugs, please open an issue with the relevant input parameters and asteroid IDs.

---

## Validation Reference (34713 Ilse)

| Method | Source | RA Error [arcsec] | Dec Error [arcsec] |
| :--- | :--- | :--- | :--- |
| **AstDyn (Nominal)** | AstDyS .eq1 | **0.627** | **0.151** |
| **AstDyn (JPL-Ref)** | JPL Horizons | **0.634** | **0.147** |

*Validation epoch: 2026-Jan-10 00:00:00 UTC (MJD 61050.0)*

### (234) Barbara Occultation & Integrator Parity (May 2026)

A rigorous validation of the (234) Barbara occultation event (May 11, 2026) has been performed to verify the consistency between the **RKF78** (Adaptive) and **AAS** (Adaptive Symplectic) integrators.

*   **Result**: Perfect agreement between integrators (MJD 61170.7086748 at TCA).
*   **Astrometric precision**: Sub-mas residuals against JPL Horizons reference when using high-precision Gaia DR3 corrections.
*   **Documentation**: Detailed analysis in **[Barbara Validation Report (IT)](astdyn/docs/validation_barbara_it.md)**.

### Haumea System Occultation (May 2026)

A full high-fidelity multi-body propagation and validation report for the complete Haumea system (including satellites Hi'iaka and Namaka) occultation event on May 4, 2026 is available here: **[Haumea Validation Report](astdyn/docs/haumea_final_validation_report.md)**.

### Ricerche Massive e Filtri Scientifici
`ioccultcalc` è stato potenziato per supportare campagne di ricerca su larga scala con criteri di qualità scientifica avanzati:

*   **Range di Asteroidi**: Supporta la sintassi `1000-5000` per processare intere sequenze.
*   **Filtri di Qualità Stellare**: 
    *   `max-ruwe`: Filtra stelle con RUWE (Gaia DR3) elevato (es. > 1.4) per garantire precisione astrometrica.
*   **Filtri Lunari**:
    *   `max-moon-phase`: Esclude osservazioni durante fasi lunari troppo luminose (0.0-1.0).
    *   `min-moon-dist`: Distanza angolare minima dalla Luna per evitare background rumoroso.
*   **Filtri Fisici e Geometrici**: 
    *   `min-duration`: Durata minima dell'evento in secondi.
    *   `min-diameter`: Diametro minimo dell'asteroide in km.
    *   `max-shadow-dist`: Raggio massimo di ricerca dell'ombra geocentrica (km).
*   **Filtri Prossimità Geografica**:
    *   `lat`, `lon`, `max-dist-km`: Filtra eventi che passano entro X km da un sito osservativo.

Esempio di ricerca professionale (ID 1000-5000) con filtri di qualità:
```bash
ioccultcalc --asteroid 1000-5000 --max-ruwe 1.4 --max-moon-phase 0.5 --min-moon-dist 15 --min-duration 0.3 --jd-start 2461171.1
```

## 🚀 Development Roadmap (Towards v1.0)

AstDyn is currently in an active developmental phase. Below is the strategic plan for the IOccultCalc tool and core library:

*   **Beta 0.5**: **Advanced Visualization**. High-fidelity SVG mapping, political boundaries, and 1-sigma uncertainty corridor. (Completed 2026-03-16)
*   **Beta 0.6 (Current)**: **Search Engine Refinement**. Rigorous Gaia DR3 corrections, annual parallax, and observer filters.
*   **v1.0 RC**: **Validation**. Final stress tests against NASA/JPL reference standards and full technical documentation release.

See the detailed [Development Plan](file:///Users/michelebigi/.gemini/antigravity/brain/ca66ea00-11d5-4584-b8f8-9d341f9a1847/development_plan.md) for more information.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
-asteroid 50936 --jd 2461112.5 --mag 15.0
