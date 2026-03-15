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

## Documentation

Full technical documentation is available in the `astdyn/docs` directory:

*   **[User Manual](astdyn/docs/AstDyn_User_Manual.pdf)**: Comprehensive guide to the API, installation, and theoretical background.

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

### ioccultcalc (v2.0 - Released 2026-03-15 16:10 CET)
A professional command-line tool for searching and comparing stellar occultations. It retrieves asteroid elements directly from **JPL Horizons** and searches the **Gaia DR3** catalog online.

**New in v2.0:**
*   **Multi-Asteroid Batch Processing**: Search multiple bodies in a single run.
*   **Flexible Input**: Supports comma-separated IDs or list files (e.g., `--asteroid @list.txt`).
*   **Chebyshev Acceleration**: High-performance daily segmented polynomials for ultra-fast discovery over long durations.
*   **Analytical Derivatives**: Exact shadow velocity and trend calculation via polynomial derivatives.

#### Usage
```bash
# Search for Ceres and Vesta over a 10-day window
ioccultcalc --asteroid 1,4 --jd-start 2461131.5 --duration 10.0 --mag 14.0
```

#### Command-line Options
- `--asteroid <list|@file>`: Comma-separated designations or a file path (prefixed with `@`).
- `--jd-start <val>`: Starting Julian Date (TDB) for the search.
- `--duration <days>`: Length of the search window in days (default: 1.0).
- `--mag <val>`: Magnitude limit for star search (default: 15.0).
- `--xml-output <file>`: Save matches to an Occult4-compatible XML.
- `--kml <file>`: Generate a KML path for Google Earth (first match).

---

## 🧪 Beta Testing
AstDyn is currently in an active testing phase. For a comprehensive guide on using the occultation engine, high-precision configurations, and analysis tools, please refer to the:

**[IOccultCalc User Guide](astdyn/docs/ioccultcalc_guide.md)**

If you encounter discrepancies or bugs, please open an issue with the relevant input parameters and asteroid IDs.

---

## Validation Reference (34713 Ilse)

| Method | Source | RA Error [arcsec] | Dec Error [arcsec] |
| :--- | :--- | :--- | :--- |
| **AstDyn (Nominal)** | AstDyS .eq1 | **0.627** | **0.151** |
| **AstDyn (JPL-Ref)** | JPL Horizons | **0.634** | **0.147** |

*Validation epoch: 2026-Jan-10 00:00:00 UTC (MJD 61050.0)*



## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
