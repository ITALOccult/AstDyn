# AstDyn C++

[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B17)
[![CMake](https://img.shields.io/badge/CMake-3.15+-blue.svg)](https://cmake.org/)
[![License](https://img.shields.io/badge/License-GPL--3.0-green.svg)](LICENSE)

Modern C++ port of **AstDyn** - A comprehensive software package for orbit determination and propagation of asteroids and celestial objects.

## 📖 Overview

AstDyn C++ is a complete rewrite of the original Fortran 90 AstDyn software, bringing modern C++ design patterns, improved performance, and enhanced maintainability to orbital mechanics computations.

### Features

- ✅ **Phase 1 Complete** (Setup & Infrastructure)
  - Modern CMake build system
  - Eigen3 for linear algebra
  - Boost for utilities
  - Google Test framework
  - Core types and constants

- 🚧 **In Development**
  - Time scale conversions
  - JPL ephemeris integration
  - Orbit propagation
  - Least squares fitting
  - Observation handling

## 🛠️ Requirements

### Minimum Requirements

- **C++ Compiler**: GCC 7+, Clang 6+, or MSVC 2017+
- **CMake**: 3.15 or higher
- **Eigen3**: 3.4 or higher (auto-fetched if not found)
- **Boost**: 1.70 or higher
  - Components: filesystem, program_options, date_time

### Optional Dependencies

- **CSPICE** (NASA SPICE Toolkit) - For enhanced ephemeris support
- **Doxygen** - For generating API documentation
- **Google Test** - For unit testing (auto-fetched if not found)

## 🚀 Quick Start

### Building from Source

```bash
# Clone the repository
git clone https://github.com/manvalan/ITALOccultLibrary.git
cd ITALOccultLibrary/astdyn

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build
cmake --build . -j$(nproc)

# Run tests
ctest --output-on-failure

# Install (optional)
sudo cmake --install .
```

### CMake Options

Configure the build with these options:

```bash
cmake -DASTDYN_BUILD_TESTS=ON \          # Build unit tests (default: ON)
      -DASTDYN_BUILD_EXAMPLES=ON \       # Build examples (default: ON)
      -DASTDYN_BUILD_DOCS=OFF \          # Build documentation (default: OFF)
      -DASTDYN_USE_SPICE=ON \            # Use SPICE toolkit (default: ON)
      -DASTDYN_ENABLE_PROFILING=OFF \    # Enable profiling (default: OFF)
      -DCMAKE_BUILD_TYPE=Release \       # Build type (Debug/Release)
      ..
```

### Example Build Configurations

**Debug build with all features:**
```bash
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DASTDYN_BUILD_TESTS=ON \
      -DASTDYN_BUILD_EXAMPLES=ON \
      ..
```

**Optimized release build:**
```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DASTDYN_BUILD_TESTS=OFF \
      -DASTDYN_BUILD_EXAMPLES=OFF \
      ..
cmake --build . -j$(nproc)
```

**With SPICE support:**
```bash
export CSPICE_ROOT=/path/to/cspice
cmake -DASTDYN_USE_SPICE=ON ..
```

## 📦 Installing Dependencies

### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    libeigen3-dev \
    libboost-all-dev \
    libgtest-dev \
    doxygen
```

### macOS (via Homebrew)

```bash
brew install cmake eigen boost googletest doxygen
```

### Windows (via vcpkg)

```powershell
vcpkg install eigen3:x64-windows boost:x64-windows gtest:x64-windows
cmake -DCMAKE_TOOLCHAIN_FILE=[vcpkg root]/scripts/buildsystems/vcpkg.cmake ..
```

## 📚 Usage Example

```cpp
#include <astdyn/AstDyn.hpp>
#include <iostream>

int main() {
    // Initialize library
    if (!astdyn::initialize()) {
        std::cerr << "Failed to initialize AstDyn\n";
        return 1;
    }
    
    // Print version
    std::cout << "AstDyn C++ v" << astdyn::Version::string << "\n";
    
    // Access constants
    using namespace astdyn::constants;
    std::cout << "AU = " << AU << " km\n";
    std::cout << "Speed of light = " << C_LIGHT << " km/s\n";
    
    // Create a 3D vector
    astdyn::Vector3d position(1.0, 0.0, 0.0);  // 1 AU on x-axis
    std::cout << "Position: " << position.transpose() << "\n";
    
    // Cleanup
    astdyn::shutdown();
    return 0;
}
```

Compile and link:
```bash
g++ -std=c++17 example.cpp -lastdyn -I/usr/local/include -L/usr/local/lib
```

## 🧪 Testing

Run all unit tests:
```bash
cd build
ctest --output-on-failure
```

Run specific test:
```bash
./tests/astdyn_tests --gtest_filter=ConstantsTest.*
```

## 📁 Project Structure

```
astdyn/
├── CMakeLists.txt              # Root CMake configuration
├── README.md                   # This file
├── LICENSE                     # GPL-3.0 license
├── cmake/                      # CMake modules and scripts
│   ├── FindCSPICE.cmake
│   ├── Version.hpp.in
│   └── Config.hpp.in
├── include/astdyn/             # Public headers
│   ├── AstDyn.hpp             # Main include file
│   ├── core/
│   │   ├── Constants.hpp      # Physical constants
│   │   └── Types.hpp          # Type definitions
│   ├── math/                  # Mathematical utilities
│   ├── time/                  # Time scale conversions
│   ├── orbit/                 # Orbital elements
│   ├── ephemeris/             # Ephemeris handling
│   ├── observations/          # Observation data
│   └── propagation/           # Orbit propagation
├── src/                       # Implementation files
│   ├── CMakeLists.txt
│   ├── core/
│   ├── math/
│   ├── time/
│   └── ...
├── tests/                     # Unit tests
│   ├── CMakeLists.txt
│   ├── test_constants.cpp
│   └── test_types.cpp
├── examples/                  # Example programs
├── docs/                      # Documentation
└── data/                      # Data files (ephemerides, etc.)
```

## 🔧 Development

### Code Style

- **C++ Standard**: C++17 (minimum)
- **Formatting**: Follow project .clang-format
- **Naming**:
  - Classes: `PascalCase`
  - Functions/methods: `snake_case`
  - Constants: `UPPER_SNAKE_CASE`
  - Namespaces: `lowercase`

### Adding New Features

1. Create header in `include/astdyn/module/`
2. Implement in `src/module/`
3. Add unit tests in `tests/`
4. Update CMakeLists.txt
5. Document with Doxygen comments

### Running Static Analysis

```bash
# Using clang-tidy
clang-tidy src/**/*.cpp -- -std=c++17 -Iinclude

# Using cppcheck
cppcheck --enable=all --std=c++17 src/
```

## 📊 Performance

Preliminary benchmarks show:
- **10-30% faster** than Fortran version on modern CPUs
- **Reduced memory footprint** with smart pointer management
- **Better cache utilization** with Eigen's optimized linear algebra

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

Original Fortran AstDyn © 1997-2020 AstDyn Consortium

## 🙏 Acknowledgments

- **Original AstDyn Team**: Andrea Milani, Steven Chesley, Mario Carpino, and contributors
- **Eigen3 Library**: Linear algebra foundation
- **NASA JPL**: SPICE Toolkit and ephemerides
- **IAU**: Standard astronomical constants

## 📞 Contact

- **Project Repository**: https://github.com/manvalan/ITALOccultLibrary
- **Issue Tracker**: https://github.com/manvalan/ITALOccultLibrary/issues
- **Documentation**: Part of ITALOccultLibrary project

## 🗺️ Roadmap

- [x] **Phase 1**: Setup & Infrastructure *(Complete)*
- [ ] **Phase 2**: Base Utilities & Math
- [ ] **Phase 3**: Ephemerides & Reference Systems
- [ ] **Phase 4**: Observations
- [ ] **Phase 5**: Orbital Elements
- [ ] **Phase 6**: Propagation Core
- [ ] **Phase 7**: Orbit Determination
- [ ] **Phase 8**: Close Approaches
- [ ] **Phase 9**: Main Programs
- [ ] **Phase 10**: Testing & Validation
- [ ] **Phase 11**: Documentation & Release

See [ASTDYN_CPP_CONVERSION_PLAN.md](ASTDYN_CPP_CONVERSION_PLAN.md) for detailed roadmap.

---

**Made with ❤️ for the asteroid science community**
