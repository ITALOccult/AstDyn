# AstDyn — asteroid dynamics and stellar occultation prediction

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Standard: C++23](https://img.shields.io/badge/standard-C%2B%2B23-blue.svg)
![Status: Beta](https://img.shields.io/badge/Status-Beta-yellow.svg)

This repository contains two things:

- **AstDyn**, a C++23 library for high-precision asteroid orbit propagation and
  determination;
- **ioccultcalc**, a command-line engine built on it that predicts stellar
  occultations by asteroids — when they happen, where the shadow falls, and how
  uncertain the prediction is.

---

## Using ioccultcalc

### 1. Install the external data

The engine needs data that is not in the repository: observatory codes,
planetary ephemerides, an asteroid catalogue. It runs without them but produces
worse results, and not always visibly — without the MPC observatory catalogue
every observation is reduced from the geocentre, which on a test fit changed the
residual RMS from 0.26″ to 1.59″.

```bash
python3 astdyn/tools/ioccultcalc_setup.py            # what is present, what is missing
python3 astdyn/tools/ioccultcalc_setup.py --essential # install the minimum
```

### 2. Run a search

```bash
ioccultcalc --asteroid 1,4 --jd-start 2461131.5 --duration 30.0 \
            --out-dir campaign --prefix test
```

Or through a configuration file, which is preferable for anything beyond a quick
test:

```yaml
objects:
  asteroids: "4"
time:
  start: "2026-07-01"
  duration_days: 30
filters:
  star_mag_max: 14.5
out-dir: "campaign"
```

```bash
ioccultcalc --conf campaign.yaml
```

### 3. Read the manual

**[ioccultcalc User Manual](astdyn/docs/IOCCULTCALC_MANUAL.md)** — installation,
configuration, physical model, integrators, orbit fitting, campaigns, data
formats, complete key reference and troubleshooting.

---

## What it does

**Occultation discovery.** Orbits are propagated through the search window,
their sky path is swept against the Gaia DR3 catalogue for stars inside a
discovery corridor, candidate events are refined to sub-millisecond closest
approach, and the shadow track is projected onto the Earth's surface through the
full IAU transformation chain.

**Orbit fitting.** Orbits can be refined on their own astrometric observations
before prediction, so that both the orbit and its uncertainty are computed
locally rather than imported. On asteroid 820987 the fit reaches an RMS of
0.26″ over 90 observations spanning 16 years, and the resulting prediction
ellipse agrees with the AstDyS one to 9% on the major axis with coincident
orientation.

**Binary asteroids.** Satellites can be described either by an SPK kernel, where
one exists, or by the parameters of their **mutual orbit** — which is what the
literature publishes for most known binaries. Each body of the system gets its
own shadow track.

**Two-stage campaigns.** `tools/orchestrator.py` screens many bodies with a
light force model, then re-runs the few that showed an event with full physics
and, optionally, orbit fitting.

---

## Using the library

AstDyn is designed so that frame and unit errors are impossible at compile time.
A position in GCRF and one in ECLIPJ2000 are different types; so are a time in
TT and one in UTC. This is the choice that most shapes the codebase.

```cpp
#include <astdyn/AstDynEngine.hpp>

AstDynEngine engine;
auto cfg = engine.config();
cfg.integrator_type = IntegratorType::RKF78;
engine.set_config(cfg);

auto orbit = physics::KeplerianStateTyped<core::ECLIPJ2000>::from_traditional(
    time::EpochTDB::from_mjd(61200.0),
    2.688, 0.103, 11.85, 253.16, 98.14, 333.02);   // angles in degrees
engine.set_initial_orbit(orbit);

auto propagated = engine.propagate_to(time::EpochTDB::from_mjd(60300.0));
```

Further documentation is in `astdyn/docs/`: the
[API reference](astdyn/docs/API_REFERENCE.md), the
[integrator comparison](astdyn/docs/integratori.md), and the
[Scientific Reference Manual](astdyn/docs/manual/) (LaTeX, in progress).

---

## Numerical integrators

Measured on asteroid 820987 propagated backwards from MJD 61200, against RKF78
at tolerance 1e-13, with planetary perturbations. Errors in AU; 1e-6 AU is
150 km, about 0.1″ at 2 AU.

| method | 10 d | 100 d | 1000 d | 3650 d | time (3650 d) |
|--------|------|-------|--------|--------|---------------|
| RKF78  | 5.3e-14 | 2.9e-10 | 2.1e-08 | 7.9e-07 | 21.9 ms |
| GRKN64 | 1.1e-11 | 1.1e-09 | 4.2e-08 | 1.9e-07 | 19.1 ms |
| RADAU  | 7.7e-13 | 2.5e-10 | 7.5e-09 | **2.1e-08** | 177.8 ms |
| GAUSS  | 1.1e-13 | 2.8e-11 | 4.0e-10 | 3.7e-08 | 248.8 ms |
| RK4    | 1.7e-12 | 4.3e-11 | 2.8e-10 | 3.7e-08 | 308.8 ms |
| AAS    | 2.8e-11 | 1.6e-09 | 1.8e-07 | 4.7e-07 | 1125.4 ms |

**RKF78** is the default choice; **RADAU** when accuracy over long arcs matters
more than speed; **AAS**, an adaptive symplectic method developed for this
library, for long-term stability studies. **SABA4 is not supported** — its
implementation is defective and raises an exception rather than returning
meaningless numbers.

A note worth drawing from the table: the five working methods agree with each
other to about 1e-7 AU, while RKF78 differs from JPL Horizons by about 2e-6.
When independent methods agree among themselves ten times better than they agree
with an external oracle, the residual discrepancy is in the force model, not in
the integration.

---

## Validation

The library is checked against external references rather than only against
itself. `ctest` runs 177 tests, of which the `ValidazioneEsterna` group compares
results with AstDyS, with JPL Horizons, and with first principles; those tests
skip with an explanatory message when the required data is absent.

**Residuals against AstDyS** — asteroid 820987, 78 observations, comparing our
residuals with those AstDyS records in the same `.rwo` file:

| | mean \|residual\| | RA bias | Dec bias |
|---|---|---|---|
| ours | 0.199″ | −0.052″ | +0.076″ |
| AstDyS | 0.169″ | — | — |

**Positions against JPL Horizons** — better than 1e-5 AU over arcs up to ten
years.

**Earth rotation** — verified from first principles: a point at Greenwich must
land where the Earth Rotation Angle says it does.

Older validation reports are in `astdyn/docs/`. Note that
`haumea_final_validation_report.md` predates corrections made in July 2026 to
the SPK ephemeris path and should be redone.

---

## Building

**Requirements**: a C++23 compiler (GCC 13+, Clang 16+, MSVC 2022), CMake 3.20+,
Eigen3, CURL, nlohmann_json.

```bash
git clone https://github.com/ITALOccult/AstDyn.git
cd AstDyn
cmake -S astdyn -B build
cmake --build build
ctest --test-dir build
```

Planetary ephemerides are installed by `ioccultcalc_setup.py`. It fetches
`de440s.bsp` (31 MB, covering 1849–2150), which is ample for occultation work.
`de441.bsp` spans far longer but is 3 GB and will stall machines with limited
memory.

---

## Configuration

Configuration files are **YAML** or **JSON**; both parse into the same internal
representation, and settings are addressed by dot-path.

```yaml
integrator:
  type: RKF78
  tolerance: 1.0e-11
physics:
  relativity: true
  asteroids:
    enabled: true
```

> The earlier OrbFit-style braced format (`integrator { type = RKF78 }`) is no
> longer supported: such files are rejected rather than silently ignored.

The full key reference is in the
[user manual](astdyn/docs/IOCCULTCALC_MANUAL.md#appendix-a--complete-key-reference).

---

## Status

| area | state |
|------|-------|
| N-body dynamics, planetary perturbations | working |
| asteroid perturbations (AST17 / BC405) | working |
| stellar corrections (proper motion, parallax, aberration, deflection) | working |
| occultation discovery and refinement | working |
| orbit fitting with own covariance | working |
| binary systems from SPK or mutual orbit | working |
| output: Occult4 XML, JSON, SVG, KML | working, under review |
| data installer | working; downloadable packages in preparation |

The stellar catalogue is a processed extract of Gaia DR3 and is not available
from a public source; distributable packages limited by magnitude are being
prepared.

---

## License

MIT — see [LICENSE](LICENSE).
