# Ephemeris Support & Native SPK Reader

## Overview

AstDyn C++ provides a high-precision, thread-safe ephemeris engine designed for parallel computing. It has transitioned away from the external CSPICE toolkit to a **native C++ SPK (stateless) implementation**, which allows massive multi-core scaling for occultation discovery and batch orbit propagation.

### Available Sources

1. **JPL DE441/440** (Standard) - Numerical integration kernels from JPL (mm-level precision for planets).
2. **VSOP87** (Analytical) - Built-in analytical theory for planetary positions (fallback when kernels are missing).
3. **AST17 / AST30** - Perturbation models for the most massive asteroids.
4. **Custom SPK** - Support for any SPICE `.bsp` kernel via the native reader.

## 🚀 Native SPK Reader (Thread-Safe)

Unlike the traditional CSPICE toolkit, which uses global buffers and is not thread-safe, the AstDyn **Native SPK Reader** (`astdyn::io::SPKReader`) is completely stateless.

### Key Benefits:
- **Zero Lock Contention**: Multiple integrators can query the same ephemeris file simultaneously.
- **Low Memory Overhead**: Direct DAF (Double Precision Array File) parsing with optimized caching.
- **Stateless Design**: Perfect for OpenMP, MPI, or SIMD-accelerated applications.

### Usage Example

```cpp
#include <astdyn/io/SPKReader.hpp>
#include <astdyn/ephemeris/DE441Provider.hpp>

// 1. Point to your JPL DE441 file
auto ephem = std::make_shared<DE441Provider>("/path/to/de441.bsp");

// 2. Query states directly (Thread-safe)
// Midnight UTC 2026-03-22 -> TDB conversion
double et = ...; 
Eigen::VectorXd state = ephem->getState(CelestialBody::EARTH, et);
```

## ☄️ Asteroid Perturbations

AstDyn includes a flexible model for asteroid-induced perturbations.

### Default Sets
- **AST17**: 16 most massive asteroids (Ceres, Pallas, Vesta, etc.) + Pluto.
- **BC405**: Extended set for ultra-high precision research.

### Dynamic Loading
You can load asteroid masses and states directly from a dedicated SPK file:

```cpp
#include <astdyn/ephemeris/AsteroidPerturbations.hpp>

AsteroidPerturbations asteroids;
// Loads masses and uses the native reader for positions
asteroids.loadSPK("/path/to/sb441-n16.bsp"); 

// Compute acceleration in m/s²
Eigen::Vector3d accel = asteroids.computePerturbation(asteroid_pos, mjd_tdb);
```

## 🛠️ Configuration

Ephemeris settings are controlled via the `ephemeris` block in the config manual:

```yaml
ephemeris {
  type = DE441
  file = "/path/to/de441.bsp"
  asteroid_file = "/path/to/sb441-n16.bsp"
}
```

## 📊 Performance & Accuracy

| Metric | VSOP87 (Analytical) | JPL DE441 (Native Reader) |
| :--- | :--- | :--- |
| **Precision** | 1-20 arcsec | < 1 mm (Planetary) |
| **Thread-Safety** | Yes | Yes (Stateless) |
| **IO Style** | Memory-resident | Random access (cached) |
| **Startup** | Instant | < 10ms (Index loading) |

## References
1. **Park et al. (2021)**: "The JPL Planetary and Lunar Ephemerides DE440 and DE441".
2. **Bretagnon (1988)**: "Planetary theories in rectangular and spherical variables (VSOP87)".
3. **NASA NAIF**: DAF and SPK file format specifications.

---
*AstDyn Documentation - 2026*
