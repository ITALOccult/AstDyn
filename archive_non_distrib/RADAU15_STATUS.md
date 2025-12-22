# Radau15 Integrator - Status and Limitations

## ⚠️ WARNING: NOT OPTIMIZED FOR PRODUCTION

**Current Status:** EXPERIMENTAL - Use with caution

## Performance Issues

### Benchmark Results (Kepler 2-body, 10 days)

| Integrator | Time | Speedup vs Radau15 |
|:-----------|-----:|-------------------:|
| **RKF78** | 0.004 ms | **1000×** faster |
| **Gauss** | 0.5 ms | **100×** faster |
| **Radau15** | ~500 ms | baseline |

**Radau15 is 100-1000× SLOWER than alternatives!**

## Why So Slow?

1. **Newton Solver:** Each step requires 4-7 Newton iterations
2. **8 Stages:** Each iteration evaluates 8 implicit equations
3. **LU Decomposition:** 8 matrix inversions per step
4. **No Prediction:** Poor initial guess for Newton solver
5. **Tight Tolerance:** Convergence criterion too strict

## When to Use (Rare Cases)

✅ **Use Radau15 ONLY if:**
- Problem is **genuinely stiff** (e.g., chemical kinetics, not orbits)
- Absolute precision > 1e-13 required
- You've tried RKF78 and it fails
- Computational cost is not a concern

❌ **DO NOT use for:**
- ❌ Asteroid orbit propagation → Use **RKF78**
- ❌ Orbit determination → Use **RKF78 + STM**
- ❌ Long-term evolution → Use **Gauss**
- ❌ Real-time applications → Use **RK4** or **RKF78**
- ❌ Occultation prediction → Use **RKF78**

## Recommended Alternatives

### For Your Use Cases

| Application | Recommended | Why |
|:------------|:------------|:----|
| **Orbit Determination** | RKF78 + STM | Fast, precise, tested |
| **Ephemeris Generation** | RKF78 | Adaptive, efficient |
| **Long-term (> 100 days)** | Gauss | Symplectic, stable |
| **Testing/Debug** | RK4 | Simple, predictable |

## Future Optimization Plan

To make Radau15 competitive, we need:

### Phase 1: Newton Solver (Priority: HIGH)
- [ ] Implement step size prediction
- [ ] Use previous step as initial guess
- [ ] Reduce iterations to 2-3
- [ ] Optimize LU decomposition

### Phase 2: Adaptive Order (Priority: MEDIUM)
- [ ] Implement variable-order scheme (Radau5/9/13)
- [ ] Switch order based on problem stiffness
- [ ] Reduce to 4 stages for non-stiff regions

### Phase 3: Benchmarking (Priority: LOW)
- [ ] Compare with DOPRI8
- [ ] Compare with RADAU5 (Hairer's implementation)
- [ ] Validate on stiff test problems

**Estimated time:** 2-3 weeks of optimization work

## Current Implementation

### What Works
- ✅ Correct Radau IIA coefficients
- ✅ A-stable (good for stiff problems)
- ✅ 15th order accuracy (when it converges)
- ✅ Adaptive step size

### What Doesn't Work Well
- ❌ Newton solver too slow
- ❌ No step prediction
- ❌ Poor initial guess
- ❌ Overkill for non-stiff problems

## Example: When Radau15 Might Be Needed

```cpp
// Hypothetical stiff problem (NOT typical orbit propagation)
// Example: Orbit with very strong drag near perihelion

auto derivative = [](double t, const Eigen::VectorXd& y) {
    // Kepler + very strong atmospheric drag
    // This creates stiffness: fast drag vs slow orbit
    // ...
};

// RKF78 would take tiny steps → very slow
// Radau15 can handle stiffness → potentially faster

RadauIntegrator radau(0.1, 1e-13);
auto result = radau.integrate(derivative, y0, t0, tf);
```

**But even then, test RKF78 first!**

## Summary

**Bottom line:** 
- Radau15 is **not ready** for production
- Use **RKF78** for everything
- Use **Gauss** for long-term
- Radau15 needs **major optimization** before it's useful

**For orbit determination:** RKF78 + STMPropagator is the way to go!

---

**Last Updated:** 2025-12-09  
**Status:** Experimental / Not Recommended
