# AstDyn: High-Fidelity Asteroid Dynamics Library

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Standard: C++17](https://img.shields.io/badge/standard-C%2B%2B17-blue.svg)
![Status: Stable](https://img.shields.io/badge/Status-Stable-green.svg)

**AstDyn** is a high-performance C++ library designed for precision asteroid orbit propagation and determination. It focuses on speed, accuracy, and modern software architecture, making it suitable for demanding applications like stellar occultation prediction and hazard assessment.

## Key Features

*   **Precision Integrators**: Implements adaptive Runge-Kutta-Fehlberg 7(8), RK4, and Symplectic Gauss-Legendre methods.
*   **Perturbation Modeling**: Full N-body solar system gravity (JPL DE440/441), Relativistic effects (1PN), and Solar Radiation Pressure.
*   **Orbit Determination**: Robust differential correction (least squares) with analytic State Transition Matrix (STM) computation.
*   **Validated Accuracy**: Validated against NASA/JPL Horizons Ephemeris System with residuals $< 2.5 \mu\text{m}$ ($10^{-9}$ AU) over 10-year propagation arcs.
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

### Example Usage

```cpp
#include <astdyn/propagation/Propagator.hpp>
#include <iostream>

int main() {
    using namespace astdyn;

    // 1. Setup Force Model (Sun + Planets)
    ForceModel forces;
    forces.enable_planets({"Earth", "Jupiter", "Mars"});

    // 2. Initial State (Position/Velocity in AU, AU/day)
    Vector6d y0; 
    y0 << 2.5, 0.0, 0.0, 0.0, 0.01, 0.0; 
    
    // 3. Propagate for 30 days
    Propagator prop(forces);
    auto result = prop.propagate(y0, 0.0, 30.0);

    std::cout << "Final State: " << result.state.transpose() << std::endl;
    return 0;
}
```

## Citation

If you use AstDyn in your research, please cite the accompanying paper (in review):

> Bigi, M., et al. (2025). "AstDyn: A High-Fidelity C++ Library for Asteroid Dynamics and Occultation Prediction." Submitted to RASTI.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
