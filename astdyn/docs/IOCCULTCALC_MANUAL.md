# ioccultcalc — User Manual

A high-precision engine for predicting stellar occultations by asteroids,
built on the AstDyn celestial-mechanics library.

This manual supersedes the earlier separate documents `ioccultcalc_guide.md`,
`config_manual.md`, `EPHEMERIS_SUPPORT.md`, `ASTDYS_FORMAT.md`,
`OOP_CONFIG_SHORT.md` and `integratori.md`.

---

## Table of contents

1. [What ioccultcalc does](#1-what-ioccultcalc-does)
2. [Installation and external data](#2-installation-and-external-data)
3. [Quick start](#3-quick-start)
4. [Configuration files](#4-configuration-files)
5. [Selecting objects and the time window](#5-selecting-objects-and-the-time-window)
6. [The physical model](#6-the-physical-model)
7. [Numerical integrators](#7-numerical-integrators)
8. [Observations and orbit fitting](#8-observations-and-orbit-fitting)
8bis. [Binary asteroids and their satellites](#8bis-binary-asteroids-and-their-satellites)
9. [How events are discovered](#9-how-events-are-discovered)
10. [Output and visualisation](#10-output-and-visualisation)
11. [Two-stage campaigns](#11-two-stage-campaigns)
12. [Data formats](#12-data-formats)
13. [Appendix A — complete key reference](#appendix-a--complete-key-reference)
14. [Appendix B — troubleshooting](#appendix-b--troubleshooting)

---

## 1. What ioccultcalc does

Given a set of asteroids and a time window, `ioccultcalc` finds the stars that
each asteroid will occult, computes when and where the shadow crosses the Earth,
and estimates how uncertain that prediction is.

The chain of computation is:

1. orbital elements are read from a local catalogue, from AstDyS files, or from
   JPL Horizons;
2. each orbit is propagated through the search window with a full force model;
3. the resulting sky path is searched against a Gaia DR3 catalogue for stars
   that fall within a discovery corridor;
4. candidate events are refined to sub-millisecond closest-approach timing;
5. the shadow track is projected onto the Earth's surface;
6. the positional uncertainty is propagated into an error ellipse on the sky
   plane.

Optionally the orbit can be **refined on its own astrometric observations**
before prediction, in which case both the orbit and its uncertainty are computed
locally rather than imported (see §8).

---

## 2. Installation and external data

`ioccultcalc` needs data that is not distributed with the source code:
observatory codes, planetary ephemerides, an asteroid catalogue and a stellar
catalogue. **The program runs without them, but produces worse results — and
not always visibly.**

The clearest example: without the MPC observatory catalogue every observation is
reduced from the geocentre, losing topocentric parallax. On a test fit of
asteroid 820987 this changed the RMS of the residuals from 0.26″ to 1.59″ — a
factor of six, with no error message beyond a warning that is easy to miss.

### 2.1 Checking what you have

    python3 tools/ioccultcalc_setup.py

This downloads nothing. It reports, for each dataset, whether it is present, how
large and how old it is, and — when missing — what the consequences are and how
to install it.

### 2.2 Installing the essentials

    python3 tools/ioccultcalc_setup.py --essential

This fetches the observatory codes, the planetary ephemeris and the asteroid
catalogue. Individual datasets can be requested with `--obscodes`,
`--ephemerides`, `--perturbers` or `--catalog`; `--force` re-downloads something
already present.

Downloads go to a temporary file, are verified, and only then installed — a
previous copy is kept as `.bak`. An interrupted download never damages data in
use.

### 2.3 The datasets

| dataset | location | size | source |
|---------|----------|------|--------|
| observatory codes | `~/.ioccultcalc/observatories/ObsCodes.txt` | 150 KB | Minor Planet Center |
| planetary ephemeris | `~/.ioccultcalc/ephemerides/de440s.bsp` | 31 MB | JPL NAIF |
| asteroid perturbers | `~/.ioccultcalc/ephemerides/sb441-n16.bsp` | 645 MB | JPL |
| asteroid catalogue | `~/.ioccultcalc/database/allnum.db` | 295 MB | AstDyS + NEODyS |
| stellar catalogue | `~/.catalog/crossreference/gaia_dr3_occult_pro.db` | 14.7 GB | prepared package |

The stellar catalogue is not available from a public source: it is a processed
extract of Gaia DR3 carrying a spatial index, proper motions, RUWE and
cross-identifications with the classical catalogues. Packages limited by
magnitude are in preparation.

`data/manifest.json` in the repository lists all sources with their expected
sizes and verification criteria; the setup tool reads it, so a source can move
without requiring a new release.

### 2.4 A note on de441

`de441.bsp` covers a far wider time span than `de440s.bsp`, but it is 3 GB and
will stall machines with limited memory. It is deliberately **not** among the
essential datasets. `de440s.bsp` covers 1849–2150, which is ample for
occultation work.

---

## 3. Quick start

A single asteroid over a short window:

```bash
ioccultcalc --asteroid 4 --jd-start 2461131.5 --duration 30.0 \
            --out-dir campaign_2026 --prefix vesta
```

The same through a configuration file, which is the preferred way for anything
beyond a quick test:

```yaml
# vesta.yaml
objects:
  asteroids: "4"

time:
  start: "2026-07-01"
  duration_days: 30

filters:
  star_mag_max: 14.5

out-dir: "campaign_2026"
prefix: "vesta"
```

```bash
ioccultcalc --conf vesta.yaml
```

### Command-line arguments

| argument | meaning |
|----------|---------|
| `--conf <file>` | configuration file (YAML or JSON) |
| `--asteroid <list>` | target IDs, e.g. `1,4` or `@targets.txt` |
| `--jd-start <jd>` | start of the search window, Julian Date (TDB) |
| `--duration <days>` | length of the window (default 1.0) |
| `--mag <limit>` | faint magnitude limit for stars |
| `--out-dir <path>` | base directory for output |
| `--prefix <str>` | prefix for per-event files (default `occ`) |
| `--xml-output <f>` | Occult4-compatible XML summary |
| `--svg-output <f>` | world map |
| `--bsp <file>` | SPK file with satellite ephemerides |
| `--system-ids <ids>` | NAIF IDs of system bodies, e.g. `101,201` |

Where a setting exists in both places, the configuration file takes precedence
over the command line, which in turn takes precedence over the legacy flat keys.

---

## 4. Configuration files

Configuration is read through the `IOCConfig` system, which accepts **YAML**
(`.yaml`, `.yml`) and **JSON** (`.json`). Both are parsed into the same internal
representation, so any setting can be written in either.

```yaml
integrator:
  type: RKF78
  step_size: 0.05
  tolerance: 1.0e-11
```

```json
{ "integrator": { "type": "RKF78", "step_size": 0.05, "tolerance": 1e-11 } }
```

> **The OrbFit-style braced format is no longer supported.** Files written as
> `integrator { type = RKF78 }` are rejected rather than silently ignored.
> Older documentation describing `.oop` files is obsolete.

Settings are addressed by dot-path — `integrator.tolerance`,
`physics.asteroids.enabled` — which is how they are named throughout this
manual.

---

## 5. Selecting objects and the time window

### 5.1 Objects

```yaml
objects:
  asteroids: "100-34244"    # inclusive numeric range
  # asteroids: "1,5,203,820987"   # explicit list
  # asteroids: "@targets.txt"     # one designation per line
```

Two further keys control where orbital elements come from:

- `objects.elements_dir` — directory of `<id>.eq1` files in AstDyS format. When
  set, elements are read from disk instead of being fetched from Horizons.
- `objects.observations_dir` — directory of `<id>.rwo` observation files, used
  by the orbit fit (§8).

`"*"` for all numbered asteroids is not implemented: it would need a defined
enumeration source.

### 5.2 Time window

`start` and `end` accept three interchangeable forms:

| form | example |
|------|---------|
| calendar date | `"2026-07-27"` |
| Modified Julian Date | `"mjd:61248.5"` |
| Julian Date | `"jd:2461248.5"` |

```yaml
time:
  start: "2026-07-27"
  end:   "2026-08-15"      # or: duration_days: 19
```

A bare number is read as JD if ≥ 2 400 000, otherwise as MJD.

### 5.3 Filters

The diameter filter is applied **before** path generation, so a body outside the
bounds costs nothing; the others act on the events found.

```yaml
filters:
  diameter_min_km: 0.0
  diameter_max_km: 5.0
  star_mag_max: 14.5
  event_duration_s_min: 0.5
  drop_mag_min: 0.2
  max_gaia_ruwe: 1.4
```

> A body whose diameter is unknown is treated as the default 100 km, so
> `diameter_max_km` excludes it too — which is appropriate for a "small bodies
> only" campaign, but worth knowing.

RUWE deserves a note: it measures how well a star's astrometry fits a
single-star model. Values above about 1.4 suggest an unresolved binary, whose
photocentre may not coincide with the occulting body's target — a real source of
negative observations.

---

## 6. The physical model

### 6.1 Ephemerides

```yaml
ephemeris:
  type: DE441                                    # or Analytical
  file: "~/.ioccultcalc/ephemerides/de440s.bsp"
  asteroid_file: "~/.ioccultcalc/ephemerides/sb441-n16.bsp"
```

Planetary positions come from a native, stateless SPK reader rather than from
CSPICE. Being stateless it carries no global buffers, so several integrators can
query the same kernel concurrently without locking — which is what makes the
parallel corridor scan possible.

`Analytical` selects a built-in VSOP87 theory as a fallback when no kernel is
available. It is accurate to a few arcseconds, against sub-millimetre for DE440,
so it is suitable for screening but not for prediction.

### 6.2 Force model

```yaml
physics:
  planets: true             # all major planets
  moon: true                # the Moon as a distinct perturber
  relativity: true          # PPN correction
  sun_j2: true
  earth_j2: true
  asteroids:
    enabled: true
    use_default_17: true    # 16 massive asteroids + Pluto
    # use_30_set: true      # BC405 extended set
    # list: [1, 4, 10]      # explicit set
```

Individual planets can be switched off with `physics.mercury` … `physics.neptune`.

When `physics.asteroids.list` is present and non-empty it **becomes** the
perturber set: it enables asteroid perturbations and disables the default sets.
The list you give is the list you use.

> Asteroid perturbations require `ephemeris.asteroid_file` to point at a kernel
> containing those bodies. With `enabled: true` but no file, the perturbations
> are simply not applied.

---

## 7. Numerical integrators

```yaml
integrator:
  type: RKF78
  step_size: 0.05          # initial step, days
  tolerance: 1.0e-11
```

### 7.1 Choosing one

| purpose | method |
|---------|--------|
| propagation, fitting, general use | **RKF78** |
| maximum accuracy over long arcs, close encounters | **RADAU** |
| efficiency over long arcs | **GRKN64** |
| long-term stability studies | **AAS** |
| teaching, short arcs | **RK4** |

### 7.2 Measured behaviour

Asteroid 820987 propagated backwards from MJD 61200, compared against RKF78 at
tolerance 1e-13, with planetary perturbations. Errors in AU; 1e-6 AU is 150 km,
about 0.1″ at 2 AU.

| method | 10 d | 100 d | 1000 d | 3650 d | time (3650 d) |
|--------|------|-------|--------|--------|---------------|
| RKF78  | 5.3e-14 | 2.9e-10 | 2.1e-08 | 7.9e-07 | 21.9 ms |
| GRKN64 | 1.1e-11 | 1.1e-09 | 4.2e-08 | 1.9e-07 | 19.1 ms |
| RADAU  | 7.7e-13 | 2.5e-10 | 7.5e-09 | **2.1e-08** | 177.8 ms |
| GAUSS  | 1.1e-13 | 2.8e-11 | 4.0e-10 | 3.7e-08 | 248.8 ms |
| RK4    | 1.7e-12 | 4.3e-11 | 2.8e-10 | 3.7e-08 | 308.8 ms |
| AAS    | 2.8e-11 | 1.6e-09 | 1.8e-07 | 4.7e-07 | 1125.4 ms |

Two observations worth drawing from this table.

**RADAU is the most accurate over a decade** — ten times better than RKF78 — at
eight times the cost. It is the implicit collocation method Everhart introduced
as RA15 and which OrbFit offers as `imet=3`; when accuracy matters more than
speed, it is the right choice.

**The five working methods agree with each other to about 1e-7 AU**, while
RKF78 differs from JPL Horizons by about 2e-6. When independent methods agree
among themselves ten times better than they agree with an external oracle, the
residual discrepancy is not in the integration but in the **force model**.

### 7.3 The methods

**RK4** — fixed-step fourth-order Runge-Kutta. The step must be chosen for the
most demanding part of the orbit and stays needlessly small elsewhere; on
eccentric orbits the cost becomes prohibitive, and nothing signals when accuracy
degrades. It is still the compiled-in default, which is a trap for long
propagations.

**RKF78** — embedded 7(8) Runge-Kutta-Fehlberg pair sharing thirteen stages;
their difference estimates the local error and drives the step. Best
accuracy-to-cost ratio, and the reference against which the others are validated.

**RADAU** — three-stage implicit Radau IIA, order 5, A-stable. Each step
requires solving a nonlinear system, so the per-step cost is high, but the steps
are long and the accuracy excellent.

**GRKN64** — generalised Runge-Kutta-Nyström. Nyström methods exploit the fact
that the equation of motion is second-order, saving evaluations at a given
order. Fastest over long arcs in these tests.

**GAUSS** — Gauss-Jackson multistep. Accurate but roughly ten times slower than
RKF78 here, with no compensating advantage in our use.

**AAS** — adaptive symplectic method with Hessian-based step control, developed
for this library. Symplectic integrators preserve the geometric structure of the
Hamiltonian problem, so energy error stays bounded rather than accumulating —
valuable over very long integrations. About fifty times slower than RKF78 for
comparable accuracy over a decade.

**SABA4** — **not supported.** The implementation is defective: 1.60 AU of error
over an arc of 0.1 days, where every other method gives about 1e-15. Selecting
it raises an exception rather than returning meaningless numbers.

### 7.4 On tolerance

The default is `1e-11`, not `1e-12`. For RKF78 working in AU/day a relative
tolerance of 1e-12 sits below the arithmetic noise of the method's thirteen
stages: it is never met, the step collapses to its minimum and integration
becomes extremely slow. Set something tighter only when a specific case needs it.

---

## 8. Observations and orbit fitting

By default `ioccultcalc` predicts from the orbital elements it is given. It can
instead **refine each orbit on its own astrometric observations** first, and
derive the prediction uncertainty from that fit.

```yaml
objects:
  observations_dir: "~/campaigns/obs"

diffcorr:
  enabled: true
  obs_years: 0              # 0 = all observations
  compute_covariance: true
```

For each body the engine looks for `<observations_dir>/<id>.rwo`. If found, the
orbit is refined by weighted least squares; **the fitted orbit replaces the
starting one only if the fit converges.** A failed fit leaves the original orbit
untouched and says so, so a bad fit can never silently degrade a good orbit.

### 8.1 What the fit does

Each iteration computes residuals between observed and predicted positions,
builds the design matrix from the state-transition matrix, solves the normal
equations, and applies the correction through a line search that only accepts
steps which reduce the RMS.

The observation model includes light-time delay, both planetary and stellar
aberration, topocentric parallax from the observatory position, and relativistic
light deflection.

A fit that stops because the line search finds no better step has reached a
minimum and is reported as **converged** — this is the normal outcome when the
starting orbit is already good.

### 8.2 Choosing the arc

`obs_years` restricts the fit to observations from the last N years relative to
the element epoch. Shorter arcs fit faster but constrain the orbit less:
accuracy depends on the temporal baseline, not on the number of points alone.

Measured on 820987, whose observations span 16 years:

| arc | observations | RMS |
|-----|--------------|-----|
| all | 90 | 0.26″ |
| 3 years | 19 | 0.32″ |
| 1 year | 3 | *refused* |

The fit is refused below six observations: six orbital elements cannot be
determined with fewer constraints than unknowns.

### 8.3 Outlier rejection

Observations whose residuals exceed a threshold are rejected. The threshold
tightens progressively from `outlier_max_sigma` to `outlier_min_sigma` over
`outlier_ramp_steps` iterations.

On asteroid 433 Eros, with 16 324 observations going back to 1893, this removes
around 700 aberrant measurements regardless of the arc used — and the RMS
improves as the arc lengthens, which is the expected behaviour of an orbit
determination:

| arc | used | rejected | RMS |
|-----|------|----------|-----|
| 1 year | 2899 | 679 | 0.855″ |
| 3 years | 4047 | 695 | 0.776″ |
| 10 years | 6155 | 708 | 0.634″ |

### 8.4 Covariance and the error ellipse

With `compute_covariance: true` the formal covariance `(AᵀWA)⁻¹` from the fit
**replaces** the covariance imported from the AstDyS elements. Orbit and
uncertainty then both come from the same solution.

The two stay coupled: a converged fit contributes both, a failed fit neither.
Using a fitted orbit with an imported covariance would describe an uncertainty
that does not belong to the orbit actually used.

For the 27 July 2026 event of 820987:

| | semi-major | semi-minor | PA | cross-track |
|---|-----------|------------|-----|-------------|
| AstDyS covariance | 0.117721″ | 0.064415″ | 71.3895 | 131.80 km |
| from our fit | 0.107151″ | 0.063296″ | 71.3871 | 120.40 km |

The orientations coincide and the axes agree to 9% and 2% — the kind of
agreement expected between two independent determinations of the same quantity.

> A formal covariance assumes the observation weights are correct. If the RMS of
> the fit greatly exceeds the quoted sigmas, the resulting ellipse will be
> optimistic. Compare the two before trusting it on an unfamiliar object.

### 8.5 Observatory positions

Topocentric parallax matters: about 4″ at 2.2 AU, far more for near-Earth
objects. It requires the MPC observatory catalogue (§2).

Some codes carry no catalogue coordinates — `270` (Unistellar Network) and other
mobile observers give their position per observation. Those still fall back to
the geocentre. On 433 Eros they account for 832 of 16 324 observations.

---

## 8bis. Binary asteroids and their satellites

`ioccultcalc` can predict occultations by the satellites of binary asteroids,
producing a separate shadow track for each body of the system.

Two ways of describing a satellite are supported, and they differ in what data
they need.

### From an SPK kernel

If a kernel with the satellite's ephemeris exists, it is the better source: it
carries positions that already include the measured perturbations.

```yaml
bsp: "~/.ioccultcalc/ephemerides/haumea_system.bsp"
system-ids: "136108,20136108,30136108"
```

NAIF numbers satellites of asteroids as `20000000 + number` for the first and
`30000000 + number` for the second, so Hi'iaka is `20136108` and Namaka
`30136108`. Kernels can be generated through the JPL Horizons SPK service.

Such kernels exist only for a handful of well-studied systems.

### From the mutual orbit

For most known binaries no kernel exists, but the literature publishes the
orbit of the satellite **around the primary**: semi-major axis, eccentricity,
inclination and period. The satellite's position is then derived by propagating
that orbit and adding the result to the primary's heliocentric position.

```yaml
objects:
  asteroids: "22"

system_bodies:
  "22":
    name: "Kalliope"
    diameter_km: 150.0
    H: 6.45
  "22-S1":
    name: "Linus"
    primary: "22"              # explicit: which body it orbits
    diameter_km: 28.0
    H: 10.4
    orbit:
      a_km: 1099.0             # semi-major axis of the mutual orbit
      e: 0.005
      i_deg: 103.7
      node_deg: 100.6
      peri_deg: 0.0
      M_deg: 0.0
      period_days: 3.5951      # required: the system's GM comes from this
      epoch_mjd: 61200.0
      plane: equatorial        # equatorial | ecliptic
```

The primary must be among the campaign's asteroids: its elements are reused.

### Two settings that decide whether the result is meaningful

**`plane`** — published orbital angles are usually referred to the **J2000
equator**, while the library works in the ecliptic. Getting this wrong rotates
the mutual orbit by 23.4°, which places the satellite's track elsewhere in a way
that still looks plausible. The default is `equatorial`, matching the common
convention, but check your source.

**`diameter_km`** — the diameter sets the shadow width and the event duration.
Without it a placeholder of 10 km is used, with a warning: for a body like
Hi'iaka (~320 km) that would make the track thirty-two times too narrow.

### Where the system's GM comes from

The gravitational parameter is derived from the **period**, through Kepler's
third law μ = 4π²a³/P², rather than from published masses. Periods are measured
well from lightcurves; binary masses are often uncertain by 20–30%.

This gives a useful consistency check: the implied mass, divided by the volume
of a sphere of the primary's diameter, should give a plausible asteroid density.
For the three systems used in testing:

| system | implied mass | implied density |
|--------|--------------|-----------------|
| (22) Kalliope | 8.14e18 kg | 4605 kg/m³ (M-type, metallic — expected high) |
| (45) Eugenia | 5.89e18 kg | 1365 kg/m³ (C-type) |
| (87) Sylvia | 1.48e19 kg | 1423 kg/m³ (C-type) |

A density outside roughly 500–6000 kg/m³ means the semi-major axis, the period
or the diameter do not agree with each other.

### Limits worth knowing

**The mutual orbit is treated as Keplerian.** Real satellite orbits are
perturbed by the primary's non-sphericity — the J2 of an irregular body is
large — by tides and by the Sun. Over long spans the node and pericentre precess
appreciably. The prediction is reliable near the element epoch and degrades away
from it. Where an SPK kernel exists, prefer it.

**The satellite's position is usually less certain than the system's
heliocentric orbit.** The prediction ellipse does not currently account for this,
so the uncertainty on a satellite track is understated.

## 9. How events are discovered

1. **Preparation** — orbital states are obtained from the local catalogue, from
   AstDyS elements or from JPL Horizons.
2. **Stellar corrections** — Gaia DR3 proper motion, annual parallax and
   relativistic aberration are applied to catalogue positions.
3. **Corridor scan** — the sky path is swept for stars within a discovery
   corridor, in parallel across cores.
4. **Refinement** — closest approach is refined analytically to sub-millisecond
   timing.
5. **Shadow projection** — the track is projected onto the Earth's surface
   through the full IAU transformation chain, using UT1.
6. **Uncertainty** — the covariance is mapped onto the sky plane to give the
   1-sigma ellipse and the cross-track uncertainty in kilometres.

Occultation-specific thresholds:

```yaml
occultation:
  min_sun_alt: -12.0        # Sun below this altitude
  min_obj_alt: 10.0         # asteroid above this altitude
  min_moon_dist: 5.0        # degrees from the Moon
  min_mag_drop: 0.05        # detectable light drop
  max_mag_star: 16.0
  max_shadow_dist: 10000    # km
  min_duration: 0.0         # seconds
  filter_daylight: true
  use_proper_motion: true
  use_parallax: true
```

### Batch robustness

The per-asteroid loop is wrapped so that a body which fails to load or to
integrate is logged and skipped, and the run continues. Query-cone generation is
capped at 200: a path that would exceed it is regenerated with a wider radius,
covering the same ground with fewer, larger queries — no events are lost.

---

## 10. Output and visualisation

```yaml
out-dir: "campaign_2026"
prefix: "occ"
xml-output: "xml/occultations.xml"
json-output: "json/occultations.json"
svg-output: "map/worldmap.svg"
kml: "map/paths.kml"

output:
  write_empty: true
```

**XML** is compatible with Occult4, carrying the asteroid elements and star
properties needed for further refinement.

**JSON** is the native format, richer than the XML, and is what the orchestrator
reads to identify positive detections.

**SVG** maps are produced globally, and individually for every event found,
named `{prefix}_{body}_{star_id}_{mjd}.svg`. For regional views use `zoom` with
`map-lat` and `map-lon`.

`output.write_empty` decides whether an empty result file is written when
nothing is found. In a batch this distinguishes "ran, found nothing" from "run
failed" — worth keeping on.

### Observer-specific output

```yaml
lat: 41.9      # geocentric latitude, degrees
lon: 12.5      # longitude, degrees
alt: 100.0     # altitude, metres
```

---

## 11. Two-stage campaigns

`tools/orchestrator.py` runs a campaign in two passes: a cheap screening over
many bodies, then a full-physics refinement over the few that showed an event.

```yaml
astdyn_base: ~/campaigns/july2026
ioccultcalc: ~/.../build/tools/bin/ioccultcalc
ephemeris_file: ~/.ioccultcalc/ephemerides/de440s.bsp

source: astdys            # db | astdys | user | jpl-horizons
objects:
  asteroids: "100-34244"

time:
  start: "2026-07-01"
  duration_days: 31

first_pass:
  profile: light

second_stage:
  profile: full
  fit: true               # refine orbits on their observations
```

```bash
python3 tools/orchestrator.py campaign.yaml --dry-run   # check first
python3 tools/orchestrator.py campaign.yaml --yes
```

The screening pass uses a light force model; bodies with a detected event are
extracted from the JSON output and re-run with full physics. With `fit: true`
the orchestrator downloads `.rwo` observation files for the positives only —
they are large, and are needed only where fitting happens.

`source` selects where elements come from: `db` from the local catalogue,
`astdys` downloaded per object with covariance, `user` from files you place in
`elements/`, `jpl-horizons` from the network at run time.

---

## 12. Data formats

### 12.1 AstDyS elements (`.eq1`)

```
format  = 'OEF2.0'
rectype = 'ML'
refsys  = ECLM J2000
END_OF_HEADER
11234
! Equinoctial elements: a, e*sin(LP), e*cos(LP), tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.
 EQU   2.6808535916678031E+00   0.032872036471001   0.036254405825130   0.103391596538937   -0.042907901689093   235.8395861037268
 MJD     61000.000000000 TDT
 MAG  12.874  0.150
```

Equinoctial elements avoid the singularities of the Keplerian set (Ω undefined
at zero inclination, ω undefined at zero eccentricity):

| symbol | quantity | definition | unit |
|--------|----------|------------|------|
| a | semi-major axis | — | AU |
| h | eccentricity component | e·sin(ϖ) | — |
| k | eccentricity component | e·cos(ϖ) | — |
| p | inclination component | tan(i/2)·sin(Ω) | — |
| q | inclination component | tan(i/2)·cos(Ω) | — |
| λ | mean longitude | Ω + ω + M | degrees |

with ϖ = Ω + ω the longitude of perihelion. Conversion back:

    e = √(h² + k²)
    i = 2·atan(√(p² + q²))
    Ω = atan2(p, q)
    ϖ = atan2(h, k)
    ω = ϖ − Ω
    M = λ − ϖ

The reference system is the mean ecliptic of J2000.0 and the epoch is MJD in TDT.

Files carrying `COV` records also contain the 6×6 covariance matrix, which
ioccultcalc uses for the prediction ellipse unless a fit supersedes it.

### 12.2 AstDyS observations (`.rwo`)

Astrometric observations with AstDyS's own residuals and weights. Each optical
line carries, for right ascension and declination separately: the assumed sigma,
a rejection flag, a catalogue debiasing correction, and the residual from
AstDyS's own fit.

Those residuals are useful as an external reference: comparing them with ours,
observation by observation, is far more revealing than comparing aggregate RMS
values — which are not even defined the same way. AstDyS publishes an RMS
normalised by the sigmas; ours is the plain arcsecond RMS of the residuals.

### 12.3 Download URLs

    https://newton.spacedys.com/~astdys2/epoch/numbered/{group}/{number}.eq1
    https://newton.spacedys.com/~astdys2/mpcobs/numbered/{group}/{number}.rwo

where `{group}` is the asteroid number divided by 1000 (so 820987 → 820,
11234 → 11, 704 → 0). The orchestrator handles this automatically, with caching
and a courtesy delay between requests.

---

## Appendix A — complete key reference

### integrator

| key | type | default | meaning |
|-----|------|---------|---------|
| `type` | string | `RK4` | `RK4`, `RKF78`, `GAUSS`, `RADAU`, `AAS`, `GRKN64` |
| `step_size` | double | 0.1 | initial step, days |
| `tolerance` | double | 1e-11 | relative error tolerance |
| `aas_precision` | double | 1e-4 | step-control metric for AAS |

### ephemeris

| key | type | default | meaning |
|-----|------|---------|---------|
| `type` | string | `DE441` | `DE441` or `Analytical` |
| `file` | string | — | path to the SPK kernel |
| `asteroid_file` | string | — | SPK with asteroid perturbers |

### observatories

| key | type | meaning |
|-----|------|---------|
| `file` | string | MPC `ObsCodes.txt`; if unset, standard locations are tried |

### physics

| key | type | default |
|-----|------|---------|
| `planets` | bool | true |
| `mercury` … `neptune` | bool | true |
| `moon` | bool | true |
| `relativity` | bool | true |
| `ppn_beta`, `ppn_gamma` | double | 1.0 |
| `sun_j2`, `earth_j2` | bool | true |
| `asteroids.enabled` | bool | true |
| `asteroids.use_default_17` | bool | true |
| `asteroids.use_30_set` | bool | false |
| `asteroids.list` | int array | — |
| `yarkovsky.enabled` | bool | false |
| `yarkovsky.a2` | double | — |

### diffcorr

| key | type | default | meaning |
|-----|------|---------|---------|
| `enabled` | bool | false | run the orbit fit |
| `obs_years` | double | 0 | limit to the last N years (0 = all) |
| `tolerance` | double | 0 | integrator tolerance during the fit (0 = unchanged) |
| `compute_covariance` | bool | true | derive the ellipse from the fit |
| `max_iter` | int | 10 | maximum iterations |
| `convergence` | double | 1e-6 | threshold in AU on the state |
| `outlier_sigma` | double | 3.0 | rejection threshold |
| `outlier_max_sigma` | double | 10.0 | initial threshold of the ramp |
| `outlier_min_sigma` | double | 3.0 | final threshold of the ramp |
| `outlier_ramp_steps` | int | 3 | iterations over which it tightens |
| `light_time` | bool | true | light-time delay |
| `aberration` | bool | — | annual stellar aberration |
| `light_deflection` | bool | — | gravitational light deflection |

> `diffcorr.tolerance` is measured to have **no effect** on cost: four orders of
> magnitude produce no measurable difference, because the expense lies in the
> state-transition matrix rather than in state propagation. It is kept for
> completeness.

### system_bodies

Physical properties and, optionally, mutual orbit of the bodies of a multiple
system. Keys are the body identifiers.

| key | type | meaning |
|-----|------|---------|
| `<id>.name` | string | display name |
| `<id>.diameter_km` | double | diameter; sets shadow width and duration |
| `<id>.H` | double | absolute magnitude; sets the light drop |
| `<id>.primary` | string | for a satellite: which body it orbits |
| `<id>.orbit.a_km` | double | semi-major axis of the mutual orbit |
| `<id>.orbit.e` | double | eccentricity |
| `<id>.orbit.i_deg` | double | inclination |
| `<id>.orbit.node_deg` | double | longitude of the ascending node |
| `<id>.orbit.peri_deg` | double | argument of pericentre |
| `<id>.orbit.M_deg` | double | mean anomaly at epoch |
| `<id>.orbit.period_days` | double | orbital period — **required**, the GM derives from it |
| `<id>.orbit.epoch_mjd` | double | epoch of the elements |
| `<id>.orbit.plane` | string | `equatorial` (default) or `ecliptic` |

### objects, time, filters, occultation, output

See §5, §9 and §10; the keys are listed there with their defaults.

---

## Appendix B — troubleshooting

**Residuals of the fit are several arcseconds, not sub-arcsecond.**
The most likely cause is a missing observatory catalogue: check with
`ioccultcalc_setup.py`. Without it every observation is reduced from the
geocentre. Measured effect: 1.59″ instead of 0.26″.

**A propagation is extremely slow.**
Check `integrator.type`. The compiled-in default is `RK4`, which is fixed-step:
over long arcs it needs an enormous number of steps. Use `RKF78`.

**Selecting SABA4 raises an exception.**
It is deliberate. That implementation is defective and would return meaningless
results.

**No events found, but the object should show one.**
Check `filters.diameter_max_km`: a body whose diameter is unknown is treated as
100 km and will be excluded by a small limit. Check also `star_mag_max` and
`occultation.min_obj_alt`.

**The prediction ellipse does not change when the fit is enabled.**
If `diffcorr.compute_covariance` is false, the ellipse keeps coming from the
imported covariance and the fit affects only the geometry — which typically
moves it by a negligible amount.

**A satellite's track is far from the primary's, or missing.**
Check `system_bodies.<id>.primary`: without it the body cannot be placed and is
skipped with a message. Check `orbit.plane` as well — published angles are
usually equatorial, and reading them as ecliptic rotates the orbit by 23.4°.

**A `.oop` configuration file is rejected.**
That format was removed. Convert to YAML or JSON; the settings map one to one.

**The engine reports fewer observatory codes than expected.**
It is falling back to a built-in list of sixteen codes. Install the full MPC
catalogue with `ioccultcalc_setup.py --obscodes`.

---

*ioccultcalc — ITALOccult / AstDyn. This manual documents the state of the
software as of July 2026.*
