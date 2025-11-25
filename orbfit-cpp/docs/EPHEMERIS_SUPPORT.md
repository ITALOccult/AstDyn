# Enhanced Ephemeris Support

## Overview

OrbFit C++ now supports multiple ephemeris sources and asteroid perturbations:

1. **VSOP87** - Built-in analytical planetary ephemerides (1-20 arcsec accuracy)
2. **JPL DE405/DE441** - High-precision numerical ephemerides via CSPICE
3. **AST17 Asteroids** - Perturbations from 16 most massive asteroids

## Installation

### CSPICE Toolkit (for JPL DE support)

1. Download CSPICE from NASA NAIF:
   ```bash
   # macOS
   wget https://naif.jpl.nasa.gov/pub/naif/toolkit/C/MacIntel_OSX_AppleC_64bit/packages/cspice.tar.Z
   
   # Linux
   wget https://naif.jpl.nasa.gov/pub/naif/toolkit/C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
   ```

2. Extract and build:
   ```bash
   uncompress cspice.tar.Z
   tar xf cspice.tar
   cd cspice
   # Add to your CMakeLists.txt:
   # find_library(CSPICE_LIBRARY cspice PATHS /path/to/cspice/lib)
   # include_directories(/path/to/cspice/include)
   ```

3. Download JPL ephemeris files:
   ```bash
   # DE441 (1550-2650, ~100 MB)
   wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441.bsp
   
   # DE405 (1600-2200, ~60 MB)
   wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de405.bsp
   ```

## Usage

### 1. Using VSOP87 (Default)

```cpp
#include <orbfit/ephemeris/PlanetaryEphemeris.hpp>

// Automatic - uses VSOP87
auto pos = PlanetaryEphemeris::getPosition(CelestialBody::EARTH, jd_tdb);
```

### 2. Using JPL DE Ephemerides

```cpp
#include <orbfit/ephemeris/JPLDEProvider.hpp>
#include <orbfit/ephemeris/EphemerisFactory.hpp>

// Create JPL DE441 provider
auto provider = std::make_unique<JPLDEProvider>(
    "path/to/de441.bsp", 
    EphemerisSource::JPL_DE441
);

// Get Earth position
auto earth_pos = provider->getPosition(CelestialBody::EARTH, jd_tdb);

// Get Mars state (position + velocity)
auto mars_state = provider->getState(CelestialBody::MARS, jd_tdb);
```

### 3. Asteroid Perturbations

```cpp
#include <orbfit/ephemeris/AsteroidPerturbations.hpp>

// Load default AST17 asteroids (16 most massive)
AsteroidPerturbations asteroids;

// Compute perturbation at spacecraft position
Eigen::Vector3d position = ...; // spacecraft position [AU]
double mjd_tdb = 60000.0;
Eigen::Vector3d accel = asteroids.computePerturbation(position, mjd_tdb);

// Enable/disable specific asteroids
asteroids.setAsteroidEnabled(1, true);   // Ceres
asteroids.setAsteroidEnabled(4, false);  // Vesta disabled

// Get total mass
double total_mass_msun = asteroids.getTotalMass();
```

### 4. Custom Asteroid Data

Create a CSV file `my_asteroids.csv`:
```csv
number,name,GM,a,e,i,omega,Omega,M0,epoch_mjd,n
1,Ceres,62.6284,2.7675,0.0760,10.593,73.597,80.3932,113.410,51544.5,0.2141
16,Psyche,1.8,2.9216,0.1339,3.096,227.305,150.2873,179.942,51544.5,0.1990
```

Load:
```cpp
AsteroidPerturbations custom_asteroids("my_asteroids.csv");
```

### 5. Integration with Propagator

```cpp
#include <orbfit/propagation/Propagator.hpp>

// Create propagator with asteroid perturbations
PropagatorSettings settings;
settings.include_asteroids = true;
settings.num_asteroids = 16;  // Use all AST17 asteroids

Propagator propagator(integrator, ephemeris, settings);

// Propagation will now include asteroid perturbations
auto final_state = propagator.propagate_cartesian(initial_state, target_mjd);
```

## Performance Considerations

### Ephemeris Sources

| Source | Accuracy | Speed | Memory |
|--------|----------|-------|--------|
| VSOP87 | 1-20" | Fast (analytical) | ~1 KB |
| JPL DE405 | ~1 km | Medium (interpolation) | ~60 MB |
| JPL DE441 | ~cm | Medium (interpolation) | ~100 MB |

### Asteroid Perturbations

- **16 asteroids**: ~16 extra force evaluations per step
- **Impact on step size**: Minimal (~5% reduction)
- **When to include**:
  - Inner solar system (< 3 AU): Use Ceres, Pallas, Vesta
  - Outer solar system: Usually negligible
  - High precision orbits: Include all 16

### Optimization Tips

```cpp
// Only include nearby asteroids
AsteroidPerturbations asteroids;
for (int i = 5; i <= 16; ++i) {
    asteroids.setAsteroidEnabled(i, false);  // Disable small ones
}

// Use VSOP87 for propagation, JPL DE for final comparison
auto vsop_provider = EphemerisFactory::createDefault();
auto jpl_provider = EphemerisFactory::create(
    EphemerisSource::JPL_DE441, "de441.bsp"
);
```

## Accuracy Comparison

### Planetary Positions (2000-2050)

| Body | VSOP87 vs DE441 |
|------|-----------------|
| Mercury | ~10 arcsec |
| Venus | ~5 arcsec |
| Earth | ~1 arcsec |
| Mars | ~5 arcsec |
| Jupiter | ~1 arcsec |
| Saturn | ~2 arcsec |
| Uranus | ~10 arcsec |
| Neptune | ~20 arcsec |

### Asteroid Perturbations

| Object | Max perturbation (1 AU) |
|--------|-------------------------|
| Ceres | ~0.3 m/s² |
| Pallas | ~0.1 m/s² |
| Vesta | ~0.1 m/s² |
| All 16 asteroids | ~1 m/s² |

**Rule of thumb**: Include asteroids for propagations > 1 year in inner solar system.

## References

1. **VSOP87**: Bretagnon, P., & Francou, G. (1988). "Planetary theories in rectangular and spherical variables"
2. **JPL DE441**: Park et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440 and DE441"
3. **AST17**: Baer, J., et al. (2011). "Astrometric masses of 21 asteroids"
4. **CSPICE**: https://naif.jpl.nasa.gov/naif/toolkit.html

## Troubleshooting

### CSPICE not found
```
CMake Error: Could not find CSPICE library
```
Solution: Install CSPICE and set `CSPICE_ROOT` in CMakeLists.txt

### Ephemeris file error
```
SPICE(SPKINSUFFDATA) -- Insufficient ephemeris data
```
Solution: Download correct DE file for your time range

### Asteroid perturbation too large
If asteroid perturbations seem unrealistic, check:
- Units: GM should be in km³/s²
- Positions: Should be in AU
- Time: Should be MJD TDB

## Future Enhancements

- [ ] Moon perturbations (ELP/MPP02)
- [ ] More asteroids (top 100)
- [ ] INPOP ephemerides support
- [ ] Automatic ephemeris selection based on time range
- [ ] Cached asteroid positions for efficiency
