# AstDyn C++ Implementation Plan & Status

## 🚀 Overview
AstDyn C++ is a modern, high-precision orbital mechanics library. This plan tracks the conversion from the legacy Fortran codebase to a high-performance C++23 system with full parallelization support.

## 🏁 Completed Phases

### Phase 1: Infrastructure & Build System
- [x] CMake 3.15+ integration
- [x] Modern C++23 standard
- [x] Unit testing framework (GTest)
- [x] Continuous Integration setup

### Phase 2: Core Math & Utilities
- [x] Eigen3 integration for linear algebra
- [x] Type-safe physics quantities (Distance, Velocity, Time)
- [x] Frame transformation engine (GCRF, ECLIPJ2000, ITRF)

### Phase 3: Ephemeris Engine
- [x] Native C++ SPK (BSP) reader (Stateless & Thread-safe)
- [x] JPL DE441 integration
- [x] Asteroid perturbation set (AST17)
- [x] VSOP87 analytical fallback

### Phase 4: Observations & Catalogs
- [x] Gaia DR3 local/online integration
- [x] Support for MPC and ADES observation formats
- [x] Stellar correction engine (PM, Parallax, Aberration)

### Phase 5: Orbital Dynamics
- [x] Keplerian, Cartesian, and Equinoctial state conversions
- [x] Initial Orbit Determination (Gauss, Gooding, OrbFit)

### Phase 6: Propagation Core
- [x] Multi-algorithm integrator suite (RKF78, Gauss, SABA, RADAU)
- [x] Encke-type (AAS) high-precision propagation
- [x] General Relativity (PPN) and J2 harmonics

### Phase 7: Orbit Determination
- [x] Batch Least Squares (Iterative Differential Correction)
- [x] Covariance propagation and mapping
- [x] Outlier rejection (Sigma-clipping)

- [x] Bessel Fundamental Plane occultation engine
- [x] Optimized Close Approach detector (Numerical TCA refinement)
- [x] High-precision MOID calculator (2-stage grid refinement)
- [x] Validated against Occult4 benchmarks (Vesta 2026, Nireus 2026)

## 🛠️ Upcoming Milestones

### Phase 9: Advanced Fitter (In Progress)
- [ ] Global batch fitting for multiple objects
- [ ] Non-gravitational force modeling (Yarkovsky)
- [ ] Systematic error (bias) estimation for catalogs

### Phase 10: Scalability & Cloud
- [ ] MPI support for distributed cluster computing
- [ ] AWS/Azure serverless deployment for Gaia mass-scans
- [ ] GPU acceleration for N-body force evaluations

### Phase 11: Release & ecosystem
- [ ] Python bindings (PyAstDyn) refinement
- [ ] Integrated documentation (Doxygen + MKDocs)
- [ ] Beta release (v1.0.0)

## 📈 Validation Status
| Benchmark | Target | Status | Error (vs Ref) |
| :--- | :--- | :--- | :--- |
| N-Body Stability | Jupiter | ✅ PASS | < 1e-12 AU |
| Occultation PA | Vesta 2026 | ✅ PASS | < 0.05° |
| Occultation Time | Vesta 2026 | ✅ PASS | < 0.001s |
| IOD Accuracy | Apophis | ✅ PASS | < 1e-6 (Residuals) |

### Refactoring Status (17/Mar/26)

| Target | Rule | Status | Action Taken |
| :--- | :--- | :--- | :--- |
| `OccultationLogic` | R1, R4 | ✅ COMPLETED | Split monolyth into 3 sub-functions. |
| `OccultationLogic` | R6 | ✅ COMPLETED | Replaced raw doubles with `Angle`. |
| `OccultationMapper` | R1, R4 | ✅ COMPLETED | Extracted projection and marker logic. |
| `OccultationMapper` | R5 | ✅ COMPLETED | Named constants for ERA coefficients. |
| `ioccultcalc.cpp` | Boy Scout | ✅ COMPLETED | Refactored main flow and logic. |
| `CloseApproach` | CFIYH | ✅ COMPLETED | Added MOID refinement and golden search for TCA. |
| `LeastSquaresOD` | R2, R6 | 🕒 PENDING | Convert residuals to `Angle` and matrix SI. |

---
*Last Updated: 2026-03-17*
