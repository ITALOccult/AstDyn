# 🌌 IOccultCalc: The Ultimate Guide to Occultation Discovery (Beta 0.6)

`ioccultcalc` is a high-precision command-line engine built on top of **AstDyn**. It integrates JPL Horizons ephemeris, Gaia DR3 stellar catalogs, and advanced numerical integration (AAS/RKF78) to predict and analyze stellar occultations.

---

## 🚀 Quick Start (v0.6 Features)

To run a basic search with organized workspace output:

```bash
# Search for Ceres and Vesta over 30 days using all CPU cores
ioccultcalc --asteroid 1,4 --jd-start 2461131.5 --duration 30.0 --out-dir campaign_2026 --prefix v_test
```

### v0.7 Features (Latest)
- **OpenMP Acceleration**: Automatically scales search tasks across all CPU cores.
- **CSPICE-Free Engine**: No external dependencies required for high-precision ephemeris.
- **Bessel-FP Refinement**: Validated occultation logic with milli-arcsecond accuracy.

### Essential Command Line Arguments

| Argument | Description |
| :--- | :--- |
| `--asteroid <list|@f>` | JPL Target ID(s) or file (e.g., `1,4` or `@targets.txt`). |
| `--jd-start <value>` | Start Julian Date (TDB) for the search window. |
| `--duration <days>` | Duration of the search in days (Default: 1.0). |
| `--mag <limit>` | Magnitude limit for stars (Default: 15.0). |
| `--conf <file>` | Path to a custom YAML, OOP, or JSON configuration file. |
| `--out-dir <path>` | Base directory for all output files. |
| `--prefix <str>` | Prefix for individual match files (default: `occ`). |
| `--xml-output <f>` | Save summary results to an Occult4-compatible XML. |
| `--svg-output <f>` | Generate a high-resolution Global SVG map. |

---

## 📂 Workspace Management & Batch Processing

Starting from **Beta 0.6**, `ioccultcalc` treats output sets as cohesive "workspaces".

### Directory Sorting
Using `--out-dir <path>` ensures that:
1.  The specified folder is created if it doesn't exist.
2.  The global XML (`events.xml`) and KML files are placed inside.
3.  **Individual SVG Maps** are automatically generated for *every* occultation found, centered on the event path.

### Filename Prefixing
Individual maps follow the naming convention:
`{prefix}_{body}_{star_id}_{mjd}.svg`

Example: `v_test_Vesta_12345678_61135.svg`

---

## 🛠️ Advanced Configuration (IOCConfig)

`ioccultcalc` now uses the unified `IOCConfig` system for all engine and tool parameters. You can pass a comprehensive setup file instead of long CLI strings.

**Example `occult_setup.yaml`:**
```yaml
integrator {
  type = RKF78
  step_size = 0.05
}
# Results will go here
out-dir = vesta_summer_campaign
prefix = vesta
jd-start = 2461131.5
duration = 60.0
asteroid = 4
```

To run:
```bash
ioccultcalc --conf occult_setup.yaml
```

---

## 🔍 Discovery & Precision Logic

When you run `ioccultcalc`, the engine performs:

1.  **Preparation**: Fetches orbital states from **JPL Horizons** or local **SPK** files.
2.  **Stellar Correctors**: Applies Gaia DR3 proper motion, **annual parallax**, and relativistic **aberration**.
3.  **Corridor Scan**: Iterates through the window searching for stars within a discovery corridor.
4.  **Refinement**: Performs analytical refinement of the **Closest Approach (CA)** using polynomial derivatives and high-precision Bessel proejction.
5.  **Multi-core Scaling**: Search tasks are distributed via **OpenMP** for maximum throughput.
6.  **Uncertainty Modeling**: If a `--covariance` file is provided, it calculates 1-sigma uncertainty cross-tracks.

---

## 📊 Outputs & Visualization

### Occult4 XML
The summary XML output is fully compatible with **Occult4**, including the necessary asteroid elements and star properties for further refinement.

### Global & Local SVG
Use `--svg-output` for a world-wide overview. For regional views, use `--zoom <v>` and `--map-lat`/`--map-lon` to center the camera.

---

> [!IMPORTANT]
> Always verify that your `de441.bsp` path is correctly set in your configuration file or points to the default `~/.ioccultcalc/ephemerides/` location.
