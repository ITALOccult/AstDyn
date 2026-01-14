# ITALOccultLibrary Integration Guide

This guide describes how to integrate the high-precision asteroid dynamics and frame conversion modules into larger projects, such as **IOccultCalc**.

## Core Components

The library consists of three main modules:

1.  **`EQ1Parser`**: Parses AstDyS `.eq1` files (OEF2.0 format).
2.  **`OrbitalConversions`**: Handles transformations between Equinoctial, Keplerian, and Cartesian states, including the critical **ECLM J2000 ↔ ICRF** frame conversion.
3.  **`AstDynWrapper`**: A high-level interface to the `AstDyn` core library, managing complex propagation settings and providing JPL-grade accuracy.
4.  **`ChebyshevApproximation`**: Provides fast interpolation for massive occultation screening using Chebyshev polynomials.
5.  **`OrbitContext`**: Metadata structure defining reference frame, origin, and element type (e.g., Mean vs Osculating).
6.  **`OrbitFitAPI`**: Simplified workflow for performing JPL-grade orbit fits from OrbFit files.

## Integration Steps

### 1. CMake Setup

Add `ITALOccultLibrary` as a dependency in your project's `CMakeLists.txt`:

```cmake
find_package(ITALOccultLibrary REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(AstDyn REQUIRED)

target_link_libraries(your_target 
    PUBLIC 
        ITALOccultLibrary::italoccultlib
        Eigen3::Eigen
        AstDyn::astdyn
)
```

### 2. Basic Propagation Flow

To propagate an asteroid from an `.eq1` file to a target epoch in the ICRF frame:

```cpp
#include <italoccultlib/astdyn_wrapper.h>
#include <iostream>

using namespace ioccultcalc;

void propagateAsteroid(const std::string& eq1_path, double target_mjd) {
    // 1. Initialize wrapper with high accuracy settings
    AstDynWrapper wrapper(PropagationSettings::highAccuracy());
    
    // 2. Load elements
    if (!wrapper.loadFromEQ1File(eq1_path)) {
        std::cerr << "Failed to load: " << eq1_path << std::endl;
        return;
    }
    
    // 3. Propagate (automatically returns ICRF state)
    try {
        CartesianStateICRF state = wrapper.propagateToEpoch(target_mjd);
        
        std::cout << "Position (AU): " << state.position.transpose() << std::endl;
        std::cout << "Velocity (AU/day): " << state.velocity.transpose() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}
```

### 3. Database Integration (Allnum.cat / SQLite)

If your application loads elements from a database (e.g., `allnum.cat` or an SQLite DB), you often have **Mean Equinoctial Elements** (AstDyS standard) but no `.eq1` file. Use `setMeanEquinoctialElements` to handle this programmatically.

The library will automatically convert these Mean elements to Osculating ICRF before propagation.

```cpp
#include <italoccultlib/astdyn_wrapper.h>

// Example: Function called after fetching a row from SQLite
void processAsteroidFromDB(const DBRow& row, double epoch_mjd, double target_mjd) {
    // 1. Initialize Wrapper
    ioccultcalc::AstDynWrapper wrapper(ioccultcalc::PropagationSettings::highAccuracy());
    
    // 2. Pass Raw Mean Elements from DB (AstDyS format)
    //    Note: Check if angles are Deg or Rad in your DB. Wrapper expects RADIANS.
    //    Allnum.cat usually has Mean Anomaly/Longitude in Degrees.
    wrapper.setMeanEquinoctialElements(
        row.a,          // Semi-major axis [AU]
        row.h,          // e * sin(w + Omega)
        row.k,          // e * cos(w + Omega)
        row.p,          // tan(i/2) * sin(Omega)
        row.q,          // tan(i/2) * cos(Omega)
        row.lambda_rad, // Mean Longitude [rad] (convert from deg if needed!)
        epoch_mjd,      // Epoch of elements [MJD TDB]
        row.name        // Asteroid Name/ID
    );
    
    // 3. Calculate Observation (High Precision)
    auto result = wrapper.calculateObservation(target_mjd);
    
    std::cout << "Asteroid: " << row.name << std::endl;
    std::cout << "  RA: " << result.ra_deg << " deg" << std::endl;
    std::cout << "  Dec: " << result.dec_deg << " deg" << std::endl;
    std::cout << "  Mag: " << result.visual_magnitude << std::endl;
}
```

### 4. Advanced High-Precision & Fitting

For high-precision work, it is critical to distinguish between **Mean** and **Osculating** elements. `AstDyn` uses `OrbitContext` to track this metadata.

```cpp
#include <astdyn/api/OrbitFitAPI.hpp>
#include <astdyn/propagation/HighPrecisionPropagator.hpp>

void highPrecisionWorkflow() {
    // 1. Convert AstDyS Mean elements to Osculating Equatorial
    auto equ_mean = astdyn::api::OrbitFitAPI::parse_eq1("34713.eq1");
    auto kep_osc = astdyn::api::OrbitFitAPI::convert_mean_equinoctial_to_osculating(equ_mean);
    
    // 2. Propagate using N-Body integrator with DE441
    astdyn::propagation::HighPrecisionPropagator::Config config;
    config.de441_path = "path/to/de441.bsp";
    astdyn::propagation::HighPrecisionPropagator prop(config);
    
    // 3. Track context
    astdyn::OrbitContext ctx;
    ctx.frame = astdyn::ReferenceFrame::EQUATORIAL_J2000;
    ctx.model = astdyn::OrbitModel::OSCULATING;
    
    auto obs = prop.calculateGeocentricObservation(kep_osc, target_jd, astdyn::propagation::HighPrecisionPropagator::InputFrame::EQUATORIAL);
    std::cout << "RA: " << obs.ra_deg << " arcsec" << std::endl;
}
```

## Precision Validation (Asteroid 34713)

The following results were obtained for Asteroid 34713 at epoch MJD 61050.0, comparing `ITALOccultLibrary` against JPL Horizons truth data.

| Case | Starting Elements | Model | RA Diff [arcsec] | Dec Diff [arcsec] |
| :--- | :--- | :--- | :--- | :--- |
| **Nominal**| AstDyS .eq1 | Mean -> Osc | **0.627** | **0.151** |
| **JPL** | JPL Horizons | Osculating | **0.634** | **0.147** |

> [!IMPORTANT]
> A difference $< 1$ arcsecond over a multi-year propagation arc confirms the library's high fidelity in handling planetary perturbations, light-time correction, and frame transformations.

## Best Practices

*   **Always use `propagateToEpoch`**: This ensures the frame is correctly converted to ICRF, which is necessary for comparison with JPL Horizons and star catalogues.
*   **Coordinate Systems**: Internal calculations in `AstDyn` are often in Ecliptic J2000. `ITALOccultLibrary` serves as the bridge to standard equatorial ICRF.
*   **Tolerances**: For final occultation predictions, keep tolerance at `1e-12`. For wide-area screening, `1e-9` is acceptable and faster.

## Troubleshooting

*   **Error vs JPL**: If errors are > 100 km, check if you accidentally skipped the `eclipticToICRF` step.
*   **Performance**: Ensure you are building in `Release` mode (`-DCMAKE_BUILD_TYPE=Release`).
