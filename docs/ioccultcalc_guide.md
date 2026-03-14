# 🌌 IOccultCalc: The Ultimate Guide to Occultation Discovery

`ioccultcalc` is a high-precision command-line engine built on top of **AstDyn**. It integrates JPL Horizons ephemeris, Gaia DR3 stellar catalogs, and advanced numerical integration (AAS/RKF78) to predict and analyze stellar occultations.

---

## 🚀 Quick Start

To run a basic search with the internal high-precision defaults:

```bash
ioccultcalc --asteroid 50936 --jd 2461112.5 --mag 15.0
```

### Command Line Arguments

| Argument | Description |
| :--- | :--- |
| `--asteroid <id>` | JPL Target ID (e.g., `1` for Ceres, `50936` for Nireus). |
| `--jd <value>` | Julian Date for the center of the search window. |
| `--mag <limit>` | Magnitude limit for the Gaia DR3 query (Default: 15.0). |
| `--conf <file>` | Path to a custom JSON configuration file. |
| `--xml-output <f>` | Save found events to an Occult-compatible XML file. |
| `--kml <file>` | Export the ground track of the first event to Google Earth. |
| `--xml-check <f>` | Compare results against a reference Occult4 XML file. |

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

1.  **Element Retrieval**: Fetches the state of the asteroid from **JPL Horizons** for the exact epoch.
2.  **Corridor Scan**: Computes the asteroid's apparent path over a 24-hour window using a high-degree Chebyshev fit.
3.  **Gaia Query**: Downloads all stars (DR3) within a **10 arcminute** corridor of the path.
4.  **Refinement**: For each candidate, it performs a numerical refinement of the **Closest Approach (CA)**.
5.  **Filtering**: Events are reported if the **Impact Parameter** is less than **100,000 km** (approx. 15 Earth radii).

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
