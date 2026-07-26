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

- `type` (string): `RK4` (default), `RKF78`, `GAUSS`, `RADAU`, `AAS`, `GRKN64`.
  `SABA4` is **no longer supported**: its step control is defective and produces
  meaningless results (1.60 AU error over a 0.1-day arc, against ~1e-15 for every
  other method). Selecting it raises an exception rather than returning wrong
  numbers silently. See `docs/integratori.md` for the diagnosis and for measured
  accuracy/cost figures of the working methods.
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

- `enabled` (bool): run the orbit fit (default `false`). When on, `ioccultcalc`
  looks for `<observations_dir>/<id>.rwo` for each body and refines its orbit
  before predicting. **The fitted orbit replaces the starting one only if the
  fit converges**; a failed fit leaves the original orbit untouched and says so.
- `obs_years` (double): use only observations from the last N years relative to
  the element epoch (default `0` = all). Shorter arcs fit faster but constrain
  the orbit less: accuracy depends on the temporal baseline. Measured on 820987,
  all 90 observations give RMS 0.26″, while a 3-year cut leaves 19 observations.
  The fit is refused below six observations — six orbital elements cannot be
  determined with fewer constraints than unknowns.
- `tolerance` (double): integrator tolerance during the fit (default `0` = keep
  the current one). **Measured as ineffective**: four orders of magnitude produce
  no measurable change in cost, because the expense lies in the state-transition
  matrix rather than in state propagation. Kept for completeness.
- `compute_covariance` (bool): compute the formal covariance `(AᵀWA)⁻¹` from the
  fit (default `true`). When available it **replaces** the `.eq1` covariance for
  the prediction ellipse, so orbit and uncertainty both come from our own
  solution. On 820987 this gives 0.107″ × 0.063″ @ PA 71.39 (cross-track
  120.4 km) against AstDyS 0.118″ × 0.064″ @ PA 71.39 (131.8 km).
- `max_iter` (int): maximum iterations (default `10`).
- `convergence` (double): threshold in AU on the state vector (default `1e-6`).
- `outlier_sigma` (double): sigma-clipping threshold (default `3.0`).
- `outlier_max_sigma`, `outlier_min_sigma` (double): start and end thresholds for
  progressive outlier rejection (defaults `10.0` and `3.0`).
- `light_time` (bool): light-time delay correction.
- `aberration` (bool): annual stellar aberration.
- `light_deflection` (bool): gravitational light deflection (GR).

> **Correction.** Earlier revisions of this manual listed `outlier_threshold`.
> That key is not read by the code; the correct name is `outlier_sigma`.

> **On convergence.** A fit that stops because the line search finds no better
> step has reached a minimum, and is reported as *converged*. This is the normal
> outcome when the starting orbit is already good.

### 2.5 `observatories`

The MPC observatory code catalogue, needed to place each observing site on the
Earth's surface.

- `file` (string): path to MPC `ObsCodes.txt`. If unset, these are tried in
  order: `~/.ioccultcalc/observatories/ObsCodes.txt`,
  `~/.ioccultcalc/ObsCodes.txt`, `./ObsCodes.txt`.

> **Why this matters.** Without the catalogue every observation is reduced from
> the **geocentre**, losing topocentric parallax — about 4″ at 2.2 AU and far
> more for nearby objects. Measured effect on the 820987 fit: RMS **1.59″
> without** the catalogue against **0.26″ with** it. When the catalogue is
> missing the tool now says so explicitly instead of degrading quietly.
>
> Get the file from
> `https://www.minorplanetcenter.net/iau/lists/ObsCodes.html` (an HTML page:
> strip the tags, the fixed-column content is inside `<pre>`).
>
> Mobile observers — code `270` (Unistellar Network) and similar — carry no
> catalogue coordinates and still fall back to the geocentre. On 433 Eros they
> account for 47% of recent observations.

### 2.6 `occultation` (engine-level)

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

- `objects.elements_dir` (string): directory of `<id>.eq1` element files
  (AstDyS format). When set, elements are read from disk instead of Horizons.
- `objects.observations_dir` (string): directory of `<id>.rwo` observation files,
  used by the orbit fit (`diffcorr.enabled`). The campaign orchestrator downloads
  them for second-stage positives only.

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

observatories:
  file: "~/.ioccultcalc/observatories/ObsCodes.txt"

# Orbit fit (optional): refines each orbit on its own observations and derives
# the prediction ellipse from the resulting covariance.
diffcorr:
  enabled: true
  obs_years: 0          # 0 = use all observations
  compute_covariance: true

objects:
  observations_dir: "~/campaigns/obs"

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
