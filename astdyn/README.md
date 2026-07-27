# AstDyn — library source

This directory holds the library, the `ioccultcalc` tool and their tests.
For an overview of the project, installation and usage, see the
[main README](../README.md).

## Layout

```
astdyn/
├── include/astdyn/     public headers
│   ├── core/           frames, constants, tagged types
│   ├── math/           vectors, rotations
│   ├── time/           epochs and time scales
│   ├── physics/        states, units
│   ├── propagation/    integrators, force model, propagators
│   ├── ephemeris/      SPK reader, planetary ephemerides
│   ├── observations/   observation records, observatory database
│   ├── astrometry/     reduction, occultation logic, Chebyshev ephemerides
│   ├── orbit_determination/  residuals, STM, differential correction
│   └── catalog/        Gaia DR3 access, Chebyshev fitting
├── src/                implementations, mirroring the header layout
├── tools/              ioccultcalc, orchestrator, setup, maintenance scripts
├── tests/              GoogleTest suite (177 tests)
├── examples/           standalone programs, mostly diagnostic
├── docs/               manuals and design notes
└── data/               manifest of external data sources
```

## Building

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

Options: `ASTDYN_BUILD_TESTS`, `ASTDYN_BUILD_EXAMPLES`, `ASTDYN_BUILD_DOCS`,
`ASTDYN_USE_OPENMP` (all ON except docs).

Some tests need external data and skip with an explanatory message when it is
absent. To run them, install the data with `tools/ioccultcalc_setup.py` and set:

```bash
export ASTDYN_EPHEMERIS_PATH=~/.ioccultcalc/ephemerides/de440s.bsp
export ASTDYN_TEST_DATA=$PWD/tests/data
```

A skipped test is not a passing test: the external-validation group is where the
library is checked against AstDyS and JPL, and it is the group that has found
every real defect so far.

## Design notes

**Frames and units are types.** A position in GCRF and one in ECLIPJ2000 cannot
be added; neither can a distance in AU and one in metres. Conversions are
explicit and checked at compile time. This is the choice that most shapes the
code, and it exists because frame confusion is the classic source of silent
errors in this field.

**One implementation per operation.** Where the same computation existed twice —
Earth rotation, integrator tolerance, astrometric reduction — one of the two
copies always turned out to be wrong. New code should reuse rather than
reimplement.

**Failures are loud.** A missing catalogue, an unsupported integrator, an
observation that cannot be reduced: each says so. Falling back silently to a
degraded behaviour has repeatedly cost more than it saved.

## Documentation

| document | content |
|----------|---------|
| [ioccultcalc User Manual](docs/IOCCULTCALC_MANUAL.md) | complete user guide |
| [API Reference](docs/API_REFERENCE.md) | library interfaces |
| [Integrators](docs/integratori.md) | measured comparison and diagnosis |
| [Scientific Reference Manual](docs/manual/) | LaTeX, in progress |
| [Configurable parsers](docs/CONFIGURABLE_PARSERS.md) | observation format definitions |

Design notes, roadmaps and validation reports are also in `docs/`.

## Conventions

- **C++23**; classes `PascalCase`, functions `snake_case`, constants
  `UPPER_SNAKE_CASE`, namespaces lowercase.
- New features: header in `include/astdyn/<module>/`, implementation in
  `src/<module>/`, test in `tests/`. Sources are globbed, so no CMake edit is
  needed for a new file.
- Tests that depend on external data must **skip** with a message, not fail.

## License

MIT — see [LICENSE](../LICENSE).
