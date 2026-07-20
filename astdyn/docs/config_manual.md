# AstDyn Configuration Manual

This manual describes the configuration parameters for the AstDyn library and
the `ioccultcalc` tool.

## 1. Configuration Syntax

AstDyn reads configuration through the `IOCConfig` system, which accepts two
formats, selected by file extension. Both are parsed into the same internal
JSON representation, so any parameter can be written in either.

### 1.1 YAML (`.yaml` / `.yml`)

Best for structured, hierarchical configurations. Parsed by the bundled fkYAML
and converted internally to JSON.

```yaml
integrator:
  type: RKF78
  step_size: 0.05
  tolerance: 1.0e-11
```

### 1.2 JSON (`.json`)

Standard for automated workflows.

```json
{
  "integrator": { "type": "RKF78", "step_size": 0.05, "tolerance": 1e-11 }
}
```

> **Note.** `IOCConfig::load()` accepts only `.json` and `.yaml`/`.yml`. The
> earlier OrbFit-style braced format (`integrator { type = RKF78 }`) is no
> longer supported — such files are rejected rather than silently ignored.

Values are addressed by dot-path (`integrator.tolerance`, `physics.asteroids.enabled`),
which is how they are named throughout this manual.

---

## 2. Core Library Parameters

These affect `AstDynEngine`: integration, ephemerides, force model, and orbit
fitting.

### 2.1 `integrator`

Controls numerical integration of asteroid orbits.

- `type` (string): `RK4` (default), `RKF78`, `GAUSS`, `RADAU`, `AAS`, `SABA4`, `GRKN64`.
- `step_size` (double): initial step in days (default `0.1`).
- `tolerance` (double): relative error tolerance for adaptive integrators (default `1e-11`).
- `aas_precision` (double): step-control metric for the `AAS` integrator (default `1e-4`).

> **On tolerance.** The default is `1e-11`, not `1e-12`. For RKF78 in AU/day a
> relative tolerance of `1e-12` is below the arithmetic noise of the method's
> thirteen stages: it is never met, the step collapses to its minimum, and
> integration becomes extremely slow. `1e-11` is accurate and reachable; set a
> tighter value explicitly only if a specific case needs it.

### 2.2 `ephemeris`

Solar-system ephemeris sources.

- `type` (string): `Analytical` (low precision) or `DE441` (high-precision JPL, **default**).
- `file` (string): path to the JPL `.bsp` file (native stateless reader). `~` is expanded.
- `asteroid_file` (string): path to a `.bsp` with asteroid ephemerides.

### 2.3 `physics`

The high-precision force model. Planets can be toggled globally or per body;
asteroid perturbers can be a default set, a 30-body set, or an explicit list.

- `planets` (bool): all major planets on/off (default `true`).
- `mercury` … `neptune` (bool): per-planet toggles (each default `true`).
- `moon` (bool): include the Moon as a distinct perturber (default `true`).
- `relativity` (bool): general-relativistic PPN correction (default `true`).
- `ppn_beta`, `ppn_gamma` (double): PPN parameters (default `1.0`).
- `sun_j2` (bool): Sun J2 oblateness (default `true`).
- `earth_j2` (bool): Earth J2 oblateness (default `true`).
- `asteroids`:
  - `enabled` (bool): asteroid perturbations on/off (default `true`).
  - `use_default_17` (bool): AstDyn default set of 17 massive asteroids + Pluto (default `true`).
  - `use_30_set` (bool): top-30 massive asteroids, BC405 set (default `false`).
  - `list` (int array): explicit perturber list, e.g. `[1, 4, 10]`. When present and non-empty it **becomes** the perturber set — it enables asteroids and disables the default sets. "The list you give is the list you use."
- `yarkovsky`:
  - `enabled` (bool): Yarkovsky acceleration (default `false`).
  - `a2` (double): Yarkovsky A2 coefficient.

> An arbitrary perturber *range* (e.g. `"1-300"`) is not yet supported: it needs
> a GM mass table for asteroids beyond the tabulated sets.

### 2.4 `diffcorr`

Differential-correction (least-squares) orbit fitting.

- `max_iter` (int): maximum iterations (default `10`).
- `convergence` (double): threshold in AU on the state vector (default `1e-6`).
- `outlier_threshold` (double): sigma-clipping threshold for rejected observations.
- `light_time` (bool): light-time delay correction.
- `aberration` (bool): annual stellar aberration.
- `light_deflection` (bool): gravitational light deflection (GR).

### 2.5 `occultation` (engine-level)

Discovery and refinement logic for occultation candidates.

- `min_sun_alt` (double): maximum Sun altitude for visibility (default `-12.0`).
- `min_obj_alt` (double): minimum asteroid altitude at the centre line (default `10.0`).
- `min_moon_dist` (double): minimum angular distance from the Moon, deg (default `5.0`).
- `min_mag_drop` (double): minimum stellar magnitude drop for detectability (default `0.05`).
- `max_mag_star` (double): faint magnitude limit for the search (default `16.0`).
- `max_shadow_dist` (double): maximum shadow search distance, km (default `10000`).
- `min_duration` (double): minimum event duration, s (default `0.0`).
- `filter_daylight` (bool): skip events in local daylight (default `true`).
- `use_proper_motion` (bool): apply Gaia DR3 proper motion to star positions (default `true`).
- `use_parallax` (bool): apply annual parallax to star positions (default `true`).

---

## 3. Dynamic Configuration for `ioccultcalc`

`ioccultcalc` reads a single file that fixes an entire prediction campaign. It
adds the blocks below **on top of** the core parameters of §2, and provides a
more readable syntax for object selection, the time window, and filters.

**Precedence.** For each setting: a value in one of the blocks below wins; then
the matching CLI flag (`--asteroid`, `--jd-start`, `--duration`, `--mag`, …);
then the legacy flat key (`asteroid`, `jd-start`, …). Existing command lines and
older config files keep working unchanged.

### 3.1 Object selection — `objects`

- `objects.asteroids` (string):
  - `"100-34244"` — inclusive numeric range
  - `"1,5,203,820987"` — explicit list
  - `"@file.txt"` — one designation per line from a file

(`"*"` for all numbered asteroids is not yet implemented: it needs a defined
enumeration source.)

### 3.2 Time window — `time`

`start` and `end` accept three interchangeable forms:

| Form | Example |
|------|---------|
| Calendar `YYYY-MM-DD` | `"2026-07-27"` |
| Modified Julian Date | `"mjd:61248.5"` |
| Julian Date | `"jd:2461248.5"` |

The window end is either `time.end` (same three forms; duration is `end − start`)
or `time.duration_days`. A bare number is read as JD if ≥ 2 400 000, else MJD.

```yaml
time:
  start: "2026-07-27"
  end:   "2026-08-15"      # or: duration_days: 19
```

### 3.3 Filters — `filters`

The diameter filter runs **before** path/cone generation, so a body outside the
bounds is skipped without paying the search cost; the others act on the events
found.

- `filters.diameter_min_km`, `filters.diameter_max_km` (double): asteroid diameter bounds. `diameter_max_km` also excludes bodies whose long shadow paths would generate many query cones.
- `filters.star_mag_max` (double): faint magnitude limit for the star.
- `filters.event_duration_s_min` (double): minimum event duration, s.
- `filters.drop_mag_min` (double): minimum light drop, mag.
- `filters.max_gaia_ruwe` (double): reject stars above this RUWE.

> A body whose diameter is unknown is treated as the default (100 km), so
> `diameter_max_km` also excludes it — appropriate for a "small bodies only"
> campaign.

### 3.4 System / satellite bodies

- `bsp` (string): SPK file with secondary/satellite ephemerides.
- `system-ids` (string): comma-separated NAIF IDs (e.g. `101,201`).

### 3.5 Uncertainty & analysis

- `covariance` (string): path to a `.cor`/`.csv` 6×6 covariance matrix.
- `clones` (int): Monte-Carlo clones for probability analysis (experimental).

### 3.6 Observer & regional filtering

- `lat`, `lon` (double): observer geocentric latitude/longitude, deg.
- `alt` (double): observer altitude, m.

### 3.7 Visualization & mapping

- `svg-output` (string): world-map SVG filename.
- `kml` (string): Google Earth path filename.
- `zoom` (double): 1.0 = global, >10 = local.
- `map-lat` / `map-lon` (double): map centre.

### 3.8 Output — `output`, catalogue, ephemeris

- `xml-output` (string): Occult4/OWC-compatible XML results.
- `out-dir` (string): base directory for output files.
- `output.write_empty` (bool): write `<Occultations></Occultations>` even with
  zero events (default `false`). In a batch this distinguishes "ran, found
  nothing" from "run failed".
- `catalog_config` (string): inline JSON selecting the star catalogue, e.g.
  `'{"catalog_type":"sqlite_dr3","sqlite_file_path":"~/.catalog/crossreference/gaia_dr3_occult_pro.db"}'`.
  Default is the local SQLite DR3 catalogue; `~` is expanded.
- `verbose` (bool): detailed logging.

### 3.9 Batch robustness

The per-asteroid loop is wrapped in try/catch: a body that fails to load (no
Horizons elements) or to integrate is logged and skipped, and the batch
continues. Query-cone generation is capped (200): a path that would exceed it is
re-generated with a wider search radius, covering the same path with fewer,
larger queries — no events lost.

---

## 4. Example: full campaign config

```yaml
# Comprehensive occultation campaign — YAML
verbose: true

objects:
  asteroids: "100-34244"

time:
  start: "2026-07-27"
  duration_days: 10

filters:
  diameter_max_km: 5.0
  star_mag_max: 14.5

integrator:
  type: RKF78
  step_size: 0.05
  tolerance: 1.0e-11

ephemeris:
  type: DE441
  file: "~/.ioccultcalc/ephemerides/de441.bsp"
  asteroid_file: "~/.ioccultcalc/ephemerides/sb441-n16.bsp"

physics:
  sun_j2: true
  earth_j2: true
  relativity: false
  asteroids:
    enabled: true
    use_default_17: true

occultation:
  min_sun_alt: -18.0
  use_proper_motion: true
  use_parallax: true

output:
  write_empty: true

xml-output: "results.xml"
svg-output: "campaign_map.svg"
lat: 41.9
lon: 12.5
zoom: 5.0
```

The same configuration in JSON is accepted identically; only the syntax differs.
