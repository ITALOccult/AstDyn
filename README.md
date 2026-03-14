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

Full scientific and technical documentation is available in the `astdyn/docs` directory:

*   **[Scientific Paper (RASTI Submission)](astdyn/docs/AstDyn_RASTI_Paper.pdf)**: Detailed methodology, algorithms, and validation results.
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

### ioccultcalc
A professional command-line tool for searching and comparing stellar occultations. It retrieves asteroid elements directly from **JPL Horizons** and searches the **Gaia DR3** catalog online.

#### Usage
Run the tool directly from the command line:
```bash
ioccultcalc --asteroid 704 --jd 2461131.61 --mag 15.0
```

#### Command-line Options
- `--asteroid <num>`: Asteroid number or designation.
- `--jd <val>`: Julian Date (TDB) for the search.
- `--mag <val>`: Magnitude limit for star search (default: 15.0).
- `--xml-output <file>`: Save matches to an Occult4-compatible XML.
- `--kml <file>`: Generate a KML path for Google Earth (first match).
- `--xml-check <file>`: Compare results against a reference XML.

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
