# AstDyn Library - Scientific Manual

## Overview

This directory contains the source files for the **AstDyn Scientific Manual**, a comprehensive documentation of the celestial mechanics and orbit determination library.

## Structure

```
manual/
├── en/                 # English version
│   ├── main.tex       # Main LaTeX document
│   ├── 00_*.tex       # Front matter (title, preface)
│   ├── 01-07_*.tex    # Part I: Theoretical Foundations
│   ├── 08-11_*.tex    # Part II: Numerical Methods
│   ├── 12-15_*.tex    # Part III: Orbit Determination
│   ├── 16-20_*.tex    # Part IV: Library Implementation
│   ├── 21-23_*.tex    # Part V: Validation
│   ├── 24-25_*.tex    # References and Appendices
│   └── Makefile       # Build system
│
├── it/                 # Italian version
│   ├── main_it.tex    # Main LaTeX document (Italian)
│   ├── 00_*.tex       # Materia preliminare
│   ├── 01-25_*.tex    # Capitoli (struttura parallela)
│   └── Makefile       # Build system
│
└── figures/            # Shared figures and diagrams
    └── (TikZ and image files)
```

## Manual Contents

### Part I: Theoretical Foundations
1. **Introduction** - Overview of celestial mechanics and AstDyn
2. **Time Systems** - UTC, TAI, TT, TDB conversions
3. **Coordinate Systems** - Equatorial, ecliptic, and transformations
4. **Reference Frames** - J2000.0, ICRS, precession
5. **Orbital Elements** - Keplerian, Cartesian, equinoctial
6. **Two-Body Problem** - Analytical solutions and Kepler's equation
7. **Perturbations** - N-body forces, relativistic effects

### Part II: Numerical Methods
8. **Numerical Integration** - Runge-Kutta, Adams methods
9. **Orbit Propagation** - Forward/backward integration
10. **State Transition Matrix** - Variational equations
11. **Ephemeris** - Planetary positions (VSOP87, DE440)

### Part III: Orbit Determination
12. **Observations** - Astrometric measurements, uncertainties
13. **Initial Orbit Determination** - Gauss's method
14. **Differential Correction** - Least-squares fitting
15. **Residuals** - O-C analysis, statistical tests

### Part IV: Library Implementation
16. **Architecture** - Design patterns, module organization
17. **Core Modules** - Time, coordinates, propagation
18. **Parsers** - OrbFit (.eq1, .rwo), MPC formats
19. **API Reference** - Complete class/function documentation
20. **Examples** - Practical usage patterns

### Part V: Validation and Applications
21. **Validation** - Comparison with OrbFit, JPL Horizons
22. **Case Studies** - Asteroid 203 Pompeja analysis
23. **Performance** - Benchmarks, accuracy assessment

## Building the Manual

### Requirements

- **LaTeX distribution**: TeX Live 2020+ or MikTeX
- **Packages**: classico, tikz, pgfplots, listings, hyperref
- **Make**: Build automation

### Compilation

**English version:**
```bash
cd en/
make
# Output: ../../build/manual/main.pdf
```

**Italian version:**
```bash
cd it/
make
# Output: ../../build/manual/main_it.pdf
```

**Quick test (single pass):**
```bash
cd en/
make quick
```

**Clean build files:**
```bash
make clean      # Remove .aux, .log, etc.
make distclean  # Remove everything including PDFs
```

## Typography

- **Font**: URW Classico (Optima-like), 12pt
- **Line spacing**: 1.3× base
- **Page size**: A4 (210 × 297 mm)
- **Margins**: Top/bottom 3cm, left 3.5cm, right 2.5cm
- **Code**: Monospaced, syntax highlighted
- **Math**: Computer Modern (LaTeX default)

## Content Status

| Chapter | English | Italian | Status |
|---------|---------|---------|--------|
| 00 (Title/Preface) | ✅ Complete | ✅ Complete | Ready |
| 01 (Introduction) | ✅ Complete | 🔄 Stub | In progress |
| 02 (Time Systems) | ✅ Complete | 🔄 Stub | In progress |
| 03 (Coordinates) | ✅ Complete | 🔄 Stub | In progress |
| 04-25 | 🔄 Stub | 🔄 Stub | Planned |

**Legend:**
- ✅ Complete: Fully written with figures/examples
- 🔄 Stub: Placeholder, needs content
- ⚠️ Draft: Initial content, needs review

## Contributing

To add or expand content:

1. Edit the relevant `.tex` file in `en/` or `it/`
2. Add figures to `figures/` directory
3. Test compilation: `make quick`
4. Commit changes with descriptive message

### Writing Guidelines

- **Mathematical rigor**: Provide derivations, not just formulas
- **Code examples**: Include working C++ snippets
- **Figures**: Use TikZ for diagrams when possible
- **References**: Cite sources for algorithms and formulas
- **Clarity**: Write for advanced undergraduates/graduate students

## Figures and Diagrams

### Creating TikZ Figures

Example orbit diagram:
```latex
\begin{tikzpicture}
    \draw[thick] (0,0) ellipse (3 and 2);
    \fill (0.8,0) circle (0.1) node[below] {Sun};
    \fill (3,0) circle (0.08) node[right] {Perihelion};
\end{tikzpicture}
```

### External Images

Place PNG/PDF images in `figures/`:
```latex
\includegraphics[width=0.8\textwidth]{figures/pompeja_residuals.pdf}
```

## License

The manual is distributed under the same license as the AstDyn library.

## Contact

- **Author**: Michele Bigi
- **Repository**: https://github.com/manvalan/ITALOccultLibrary
- **Issues**: Use GitHub issue tracker for errors/suggestions

---

**Last Updated**: November 26, 2025
