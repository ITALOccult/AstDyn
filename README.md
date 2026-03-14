# AstDyn: High-Fidelity Asteroid Dynamics Library

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Standard: C++17](https://img.shields.io/badge/standard-C%2B%2B17-blue.svg)
![Status: Stable](https://img.shields.io/badge/Status-Stable-green.svg)

**AstDyn** is a high-performance C++ library designed for precision asteroid orbit propagation and determination. It focuses on speed, accuracy, and modern software architecture, making it suitable for demanding applications like stellar occultation prediction and hazard assessment.

## Key Features

*   **Precision Integrators**: Implements adaptive Runge-Kutta-Fehlberg 7(8), RK4, and Symplectic Gauss-Legendre methods.
*   **Perturbation Modeling**: Full N-body solar system gravity (JPL DE440/441), Relativistic effects (1PN), and Solar Radiation Pressure.
*   **Orbit Determination**: Robust differential correction (least squares) with analytic State Transition Matrix (STM) computation.
*   **High-Precision Wrapper**: Simplified interface (`ITALOccultLibrary`) for rapid integration with high-level workflows.
*   **Validated Accuracy**: Validated against NASA/JPL Horizons Ephemeris System with residuals $< 2.5 \mu\text{m}$ ($10^{-9}$ AU) and sub-arcsecond geocentric precision ($< 0.7$ arcsec).
*   **Performance**: Optimized C++17 architecture, ~8x faster than numerical STM approaches.

## Documentation

*   **[User Manual](astdyn/docs/AstDyn_User_Manual.pdf)**: Comprehensive guide to the API, installation, and theoretical background.

## Quick Start

### Prerequisites
*   C++17 compliant compiler (GCC 9+, Clang 10+, MSVC 2019+)
*   CMake 3.15+
*   Eigen3 library

### Installation

```bash
git clone https://github.com/your-username/astdyn.git
cd astdyn
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

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

1. Requisiti Fondamentali
Aggiungi o aggiorna la sezione dei requisiti per essere molto chiaro sul file JPL:

markdown
### Prerequisites
- **C++20 Compiler** (GCC 11+, Clang 13+, or MSVC 2022)
- **CMake 3.20+**
- **JPL DE441 Ephemeris**: The tools require the `de441.bsp` file. 
  - Download it from [JPL FTP](https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de441.bsp).
  - Place it in the executable directory or at `~/.ioccultcalc/ephemerides/de441.bsp`.
2. Sezione "High-Precision Occultations"
Metti in risalto le nuove capacità dinamiche:

markdown
## 💎 High-Precision Dynamics
AstDyn 3.0 introduces an extreme-precision mode for occultation discovery, featuring:
- **AAS Integrator**: Adaptive step-size Symplectic integrator for long-term stability.
- **Full N-Body**: Perturbations from all planets, the Moon, and the **AST17** (17 most massive asteroids).
- **Relativistic Corrections**: Post-Newtonian (PPN) gravity model.
- **Dynamic Distance**: Real-time interpolation of Earth-Asteroid distance during the event refinement.
3. Guida e Beta Testing
Invita i tester a leggere la nuova documentazione:

markdown
## 🧪 Beta Testing
We are currently in active beta. For a detailed walkthrough on how to use the search engine, please refer to our **[IOccultCalc User Guide](docs/ioccultcalc_guide.md)**.
### Reporting Issues
If you find a discrepancy with Occult4 or JPL Horizons, please open an issue providing:
1. The asteroid name/number.
2. The exact Julian Date used.
3. Your `ioccultcalc_precision.json` configuration.
4. Aggiornamento ioccultcalc
Assicurati che l'esempio nel README sia quello più aggiornato:

markdown
### ioccultcalc Quick Start
To find occultations for a Jupiter Trojan (e.g., 50936 Nireus) with full perturbations:
```bash
# Compile
cmake -B build
cmake --build build --target ioccultcalc
# Run search for March 13, 2026
./build/astdyn/tools/bin/ioccultcalc --asteroid 50936 --jd 2461112.5 --mag 15.0
