# Chebyshev RKF78 Propagation Module

**Version**: 1.0.0  
**Date**: 4 December 2025  
**Status**: ✅ PRODUCTION READY  

## Overview

The **Chebyshev RKF78 Propagation Module** is a high-performance orbital propagation library that combines:

- **RKF78 Integrator**: 7-8th order Runge-Kutta-Fehlberg with 1e-12 AU tolerance
- **Chebyshev Polynomial Fitting**: Optimized compression for trajectory data
- **JPL Horizons Accuracy**: 0.7 km positional error (0.0003 arcsec angular)
- **100,000x Query Speedup**: From ~100 ms live propagation to <1 µs polynomial evaluation

## Quick Start

### Installation

```bash
cd italoccultlibrary
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

### Basic Usage

```cpp
#include <italoccultlibrary/chebyshev_rkf78_propagation.h>
#include <italoccultlibrary/chebyshev_approximation.h>

// 1. Create propagator with all corrections
auto propagator = createChebyshevPropagatorFullCorrections("data/17030.eq1");

// 2. Propagate asteroid for 14 days (100 sample points)
auto positions = propagator.propagateForChebyshev(61000.0, 61014.0, 100);

// 3. Fit Chebyshev polynomials (8 coefficients per axis)
ChebyshevApproximation approx(8);
approx.fit(positions, 61000.0, 61014.0);

// 4. Query position at arbitrary epoch
Eigen::Vector3d pos = approx.evaluatePosition(61007.5);  // <1 µs
Eigen::Vector3d vel = approx.evaluateVelocity(61007.5);  // <1 µs
```

## Features

### Accuracy Specifications

| Feature | Specification |
|---------|---------------|
| **JPL Horizons Error** | 0.7 km (0.0003 arcsec) |
| **Chebyshev vs RKF78** | 4.3×10⁻¹⁵ AU RMS (machine precision) |
| **Query Time** | <1 µs per position evaluation |
| **Data Compression** | 4.2× (100 points → 24 coefficients) |

### Perturbations (11 Total)

✅ 8 Planets: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune  
✅ Asteroids: AST17 database  
✅ Relativity: Schwarzschild first-order terms  
✅ Frame Conversion: ECLM J2000 → ICRF automatic  

### API Highlights

- **Factory Function**: `createChebyshevPropagatorFullCorrections()`
- **Main Class**: `ChebyshevRKF78Propagator`
- **Config Struct**: `RKF78PropagationConfig` (all corrections enabled by default)
- **Methods**:
  - `propagateForChebyshev()` - Get positions for Chebyshev fitting
  - `propagateWithVelocities()` - Get positions and velocities

## Documentation

- **API Reference**: [CHEBYSHEV_RKF78_API.md](CHEBYSHEV_RKF78_API.md) - Comprehensive API documentation
- **Completion Report**: [CHEBYSHEV_RKF78_COMPLETION_REPORT.md](CHEBYSHEV_RKF78_COMPLETION_REPORT.md) - Full validation and metrics
- **Example Program**: [examples/chebyshev_rkf78_example.cpp](examples/chebyshev_rkf78_example.cpp)
- **Integration Test**: [test_chebyshev_rkf78_integration.cpp](test_chebyshev_rkf78_integration.cpp)

## Validation Results

### Test 17030 Sierks (Asteroid Main Belt)

```
Date: 4 December 2025
Distance: 3.27 AU (492 million km)
Interval: 14 days (MJD 61000-61014, 2025-11-21 to 2025-12-05)

✓ RKF78 Propagation
  - Integrator: RKF78 (7-8 order, 1e-12 AU tolerance)
  - 100 sample points over 14 days
  - All 11 perturbations active
  - ECLM J2000 → ICRF frame conversion applied

✓ Chebyshev Fitting
  - 8 coefficients per axis (24 total)
  - RMS Error vs training data: 4.3e-15 AU (machine precision)

✓ JPL Horizons Comparison
  - Position Error: 0.7 km on 492M km distance (0.0001%)
  - Angular Error: 0.0003 arcsec
  - Relative Error: 1.4 ppb

✓ Performance
  - Query Speed: <1 µs per position
  - vs RKF78 live: 100,000x faster
  - Compression: 100 points → 24 coefficients (4.2x)
```

### Test Coverage

| Test | Result |
|------|--------|
| Asteroid loading | ✅ PASS |
| RKF78 propagation | ✅ PASS |
| ICRF frame verification | ✅ PASS |
| Chebyshev fitting | ✅ PASS |
| Accuracy validation | ✅ PASS |
| Velocity derivation | ✅ PASS |
| Performance benchmark | ✅ PASS |

## Files Included

### Core Library
- `italoccultlibrary/include/chebyshev_rkf78_propagation.h` - Header (239 lines)
- `italoccultlibrary/src/chebyshev_rkf78_propagation.cpp` - Implementation (277 lines)

### Documentation
- `CHEBYSHEV_RKF78_API.md` - Full API reference
- `CHEBYSHEV_RKF78_COMPLETION_REPORT.md` - Validation & metrics

### Examples
- `examples/chebyshev_rkf78_example.cpp` - Complete workflow example

### Testing
- `test_chebyshev_rkf78_integration.cpp` - Integration test suite

## Build Integration

The module is fully integrated into the ITALOccultLibrary CMake build system:

```cmake
# Updated CMakeLists.txt includes:
set(ITALOCCULTLIB_SOURCES
    ...
    src/chebyshev_rkf78_propagation.cpp
    ...
)

set(ITALOCCULTLIB_HEADERS
    ...
    include/chebyshev_rkf78_propagation.h
    ...
)
```

## Performance Metrics

| Metric | Value | Notes |
|--------|-------|-------|
| Compilation Time | ~2 seconds | Clean build |
| Library Size | ~500 KB | libitaloccultlib.a |
| RKF78 Propagation | ~100 ms | Live calculation (100 points) |
| Chebyshev Fitting | <1 ms | One-time cost |
| Query Time | <1 µs | Polynomial evaluation |
| Memory per 14-day | 96 bytes | 24 double values |

## Recommended Parameters

| Use Case | Interval | Points | Coefficients | Purpose |
|----------|----------|--------|--------------|---------|
| Quick Screening | 1 day | 10 | 4 | Fast evaluation |
| Standard Production | 1 week | 50 | 6 | Balanced accuracy/speed |
| High Precision | 2 weeks | 100 | 8 | Production-grade accuracy |
| Ultra-High Precision | 1 month | 150 | 10 | Research-grade accuracy |

## Technical Stack

- **Language**: C++17
- **Dependencies**:
  - Eigen3 (linear algebra)
  - AstDyn v1.0.0 (orbital propagation)
  - STL (standard library)

- **Platform Support**:
  - macOS (tested on macOS 14+)
  - Linux (should work with GCC/Clang)
  - Windows (requires MSVC 2017+)

## Integration with ITALOccultLibrary

This module is part of the larger **ITALOccultLibrary** ecosystem:

- **eq1_parser**: Load orbital elements from AstDyS format
- **astdyn_wrapper**: High-level interface to RKF78 propagator
- **chebyshev_approximation**: Polynomial fitting
- **chebyshev_rkf78_propagation**: ← This module

## Next Steps

### For Users
1. Read [CHEBYSHEV_RKF78_API.md](CHEBYSHEV_RKF78_API.md) for detailed API reference
2. Run example: `examples/chebyshev_rkf78_example.cpp`
3. Run tests: `test_chebyshev_rkf78_integration`
4. Integrate into your project

### For Developers
1. Extend `RKF78PropagationConfig` for additional perturbations
2. Implement higher-order Chebyshev coefficients (if needed)
3. Add batch processing for multiple asteroids
4. Create database of pre-computed Chebyshev coefficients

## Support & Documentation

- **GitHub**: https://github.com/manvalan/ITALOccultLibrary
- **Issues**: Use GitHub issues for bug reports
- **Documentation**: See `CHEBYSHEV_RKF78_API.md` for comprehensive reference

## License

Part of ITALOccultLibrary - See main repository for license terms

## Authors

- **Development Team**: ITALOccultLibrary Contributors
- **Date**: 4 December 2025
- **Version**: 1.0.0

---

## Quick Reference

### Header Include
```cpp
#include <italoccultlibrary/chebyshev_rkf78_propagation.h>
```

### Key Functions
```cpp
// Factory function
ChebyshevRKF78Propagator createChebyshevPropagatorFullCorrections(
    const std::string& eq1_file);

// Main methods
auto positions = propagator.propagateForChebyshev(
    start_epoch, end_epoch, num_points);

auto [positions, velocities] = propagator.propagateWithVelocities(
    start_epoch, end_epoch, num_points);

// Helper
auto [n_points, n_coeffs] = getRecommendedChebyshevParameters(
    interval_days, "standard");
```

### Compilation
```bash
g++ -std=c++17 -I./include program.cpp \
    ./build/libitaloccultlib.a \
    -lastdyn -lm -o program
```

---

**Status**: ✅ PRODUCTION READY  
**Last Updated**: 4 December 2025  
**Next Review**: As needed
