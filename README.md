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

Full scientific and technical documentation is available in the `astdyn/docs` directory:

*   **[Scientific Paper (RASTI Submission)](astdyn/docs/AstDyn_RASTI_Paper.pdf)**: Detailed methodology, algorithms, and validation results.
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

## Validation Reference (34713 Ilse)

| Method | Source | RA Error [arcsec] | Dec Error [arcsec] |
| :--- | :--- | :--- | :--- |
| **AstDyn (Nominal)** | AstDyS .eq1 | **0.627** | **0.151** |
| **AstDyn (JPL-Ref)** | JPL Horizons | **0.634** | **0.147** |

*Validation epoch: 2026-Jan-10 00:00:00 UTC (MJD 61050.0)*

## Citation

If you use AstDyn in your research, please cite the accompanying paper (in review):

> Bigi, M., et al. (2025). "AstDyn: A High-Fidelity C++ Library for Asteroid Dynamics and Occultation Prediction." Submitted to RASTI.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
