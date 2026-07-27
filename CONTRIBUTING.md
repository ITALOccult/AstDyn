# Contributing

Contributions are welcome, whether they are bug reports, questions,
documentation fixes or code.

## Reporting a problem

Open an issue at
<https://github.com/ITALOccult/AstDyn/issues>. What helps most:

- what you ran — the command line or configuration file;
- what you expected and what happened instead;
- for a prediction that looks wrong, the asteroid number, the date, and if
  possible a reference prediction to compare against.

That last point matters more than it might seem. Every real defect found in this
software so far has surfaced from comparison with an external reference — AstDyS
residuals, JPL Horizons positions, published predictions — and none from the
test suite alone. A report that says "this event comes out 40 km from where
Occult4 puts it" is far more useful than one that says the result looks odd.

## Asking a question

Questions about how to use the software are welcome as issues too; if the
answer is not in the [user manual](astdyn/docs/IOCCULTCALC_MANUAL.md), that is
worth knowing.

## Contributing code

1. Fork the repository and branch from `main`.
2. Build and run the tests: `cmake -S astdyn -B build && cmake --build build &&
   ctest --test-dir build`.
3. Add a test for what you changed.
4. Open a pull request describing what the change does and why.

### Conventions

- **C++23**. Classes `PascalCase`, functions `snake_case`, constants
  `UPPER_SNAKE_CASE`, namespaces lowercase.
- Headers go in `astdyn/include/astdyn/<module>/`, implementations in
  `astdyn/src/<module>/`, tests in `astdyn/tests/`. Sources are globbed, so no
  CMake edit is needed for a new file.
- **Frames and units are types.** A position in GCRF and one in ECLIPJ2000
  cannot be added, and neither can a distance in AU and one in metres.
  Conversions are explicit. Please do not work around this with raw doubles.
- **One implementation per operation.** Where the same computation has existed
  twice in this codebase — Earth rotation, integrator tolerance, astrometric
  reduction — one of the two copies always turned out to be wrong. Reuse rather
  than reimplement.
- **Failures are loud.** A missing catalogue, an unsupported integrator, an
  observation that cannot be reduced: each says so. Silently falling back to a
  degraded behaviour has repeatedly cost more than it saved.

### Tests

Some tests need external data — ephemerides, the observatory catalogue, test
fixtures — and **skip with an explanatory message** when it is absent. That is
deliberate: a test that cannot run should say so rather than fail. Install the
data with `astdyn/tools/ioccultcalc_setup.py` and set `ASTDYN_EPHEMERIS_PATH`
and `ASTDYN_TEST_DATA` to run them.

A skipped test is not a passing test. The external-validation group is where
this software is checked against something other than itself, and it is worth
running before submitting a change to the physics.

## Scope

This project predicts stellar occultations and determines asteroid orbits. It is
not a general-purpose N-body integrator; for that, REBOUND and similar packages
are better suited and better tested.
