# 🌌 IOccultCalc: The Ultimate Guide to Occultation Discovery

`ioccultcalc` is a high-precision command-line engine built on top of **AstDyn**. It integrates JPL Horizons ephemeris, Gaia DR3 stellar catalogs, and advanced numerical integration (AAS/RKF78) to predict and analyze stellar occultations.

---

## 🚀 Quick Start

To run a basic search with the internal high-precision defaults:

```bash
# Search for Ceres and Vesta over 5 days starting from JD 2461112.5
ioccultcalc --asteroid 1,4 --jd-start 2461112.5 --duration 5.0 --mag 15.0
```

### Command Line Arguments

| Argument | Description |
| :--- | :--- |
| `--asteroid <list|@f>` | JPL Target ID(s) or file (e.g., `1,4` or `@target_list.txt`). |
| `--jd-start <value>` | Start Julian Date (TDB) for the search window. |
| `--duration <days>` | Duration of the search in days (Default: 1.0). |
| `--mag <limit>` | Magnitude limit for the Gaia DR3 query (Default: 15.0). |
| `--conf <file>` | Path to a custom JSON configuration file. |
| `--xml-output <f>` | Save found events to an Occult-compatible XML file. |
| `--kml <file>` | Export the ground track of the first event to Google Earth. |

---

## 🛠️ Precision Configuration (AAS)

For maximum accuracy, especially with perturbed bodies like Jupiter Trojans (e.g., Nireus), use the **AAS (Adaptive Step size)** integrator.

### `ioccultcalc_precision.json`

```json
{
    "ephemeris_file": "/Users/michelebigi/.ioccultcalc/ephemerides/de441.bsp",
    "ephemeris_type": "DE441",
    "integrator_type": "AAS",
    "aas_precision": 1e-04,
    "initial_step_size": 0.01,
    "propagator_settings": {
        "include_planets": true,
        "include_moon": true,
        "include_relativity": true,
        "include_asteroids": true
    },
    "light_time_correction": true,
    "aberration_correction": true
}
```

To use it:
```bash
ioccultcalc --conf ioccultcalc_precision.json --asteroid 50936 --jd 2461112.5
```

---

## 🔍 Discovery Logic

When you run `ioccultcalc`, the engine performs the following steps:

1.  **Preparation**: Fetches orbital states from **JPL Horizons** for all listed asteroids.
2.  **Chebyshev Pre-calculation**: Fits high-degree polynomials for each body across the entire search window (optimized for duration).
3.  **Corridor Scan**: Iterates daily through the window, searching for stars (DR3) within a **10 arcminute** corridor for each asteroid.
4.  **Refinement**: For every candidate impact, it performs analytical refinement of the **Closest Approach (CA)** using polynomial derivatives.
5.  **Filtering**: Events are reported if the **Impact Parameter** is within the discovery threshold (Default: 50,000 km).

---

## 📊 Interpreting Results

The output list shows:
- **TCA (MJD)**: Time of Closest Approach.
- **Impact**: Distance from the Earth's center to the shadow axis (km). If `< 6371 km`, the shadow hits the Earth!
- **Vel**: Velocity of the shadow on the fundamental plane (km/s).

### Ground Track Visualization
Use the `--kml` flag to generate a file you can open in **Google Earth** or **Marble** to see exactly where the shadow will pass:

```bash
ioccultcalc --asteroid 50936 --jd 2461112.5 --kml track.kml
```

---

> [!TIP]
> If you are looking for a specific star and it's not appearing, try increasing the `--mag` limit or check if the asteroid orbit in JPL is up-to-date with the latest astrometry.
