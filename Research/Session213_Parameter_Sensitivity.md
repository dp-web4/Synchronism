# Session #213: Parameter Sensitivity Analysis

**Date**: January 2, 2026
**Machine**: CBP
**Status**: COMPLETE - ROBUSTNESS CONFIRMED

---

## Executive Summary

Session #213 addresses Nova's recommendation to explore "how stable are results under small perturbations of A, B, γ?" through comprehensive sensitivity analysis.

**Key Finding**: Synchronism predictions are robust to cosmological parameter uncertainties (~1%), with f_indiff parameters dominating the uncertainty budget.

---

## Part 1: Parameter Space

### Fiducial Values

| Parameter | Value | Uncertainty | Source |
|-----------|-------|-------------|--------|
| Ω_m | 0.315 | ±0.007 | Planck 2018 |
| H₀ | 67.4 km/s/Mpc | ±0.5 | Planck 2018 |
| H₀ (local) | 73.0 km/s/Mpc | ±1.0 | SH0ES |
| φ | 1.618034 | 0 (exact) | Golden ratio |
| A | 37.5 | ~±5 | Session #210 |
| β | -0.72 | ~±0.05 | Session #210 |
| M_break | 2.2×10⁴ M_sun | ~±0.3 dex | Session #211 |

---

## Part 2: a₀ Sensitivity

### Critical Acceleration: a₀ = c × H₀ × Ω_m^φ

| Perturbation | Effect on a₀ |
|--------------|--------------|
| Ω_m ± 0.007 | ±3.6% |
| H₀ ± 0.5 | ±0.7% |
| H₀ tension (67.4 → 73.0) | +8.3% |
| φ ± 0.01 (hypothetical) | ~1% |

**Fiducial**: a₀ = 1.01×10⁻¹⁰ m/s²

**Conclusion**: a₀ is robust to Planck uncertainties. H₀ tension shifts a₀ by ~8%, but this is within typical observational uncertainties.

---

## Part 3: C(a) Sensitivity

### Coherence Function Stability

At different acceleration regimes:

| a/a₀ | C(a) | ΔC (Ω_m) | ΔC (H₀ tension) |
|------|------|----------|-----------------|
| 0.01 | 0.35 | 0.006 | 0.002 |
| 0.1 | 0.45 | 0.003 | 0.005 |
| 1.0 | 0.66 | 0.000 | 0.008 |
| 10 | 0.87 | 0.001 | 0.005 |
| 100 | 0.96 | 0.000 | 0.002 |

**Conclusion**: C(a) is VERY stable. Maximum variation ~1% at transition regime.

---

## Part 4: f_indiff Sensitivity

### The Dominant Uncertainty Source

| Parameter | Variation | Effect at M=10³ M_sun |
|-----------|-----------|----------------------|
| A ± 5 | ~13% | Linear scaling |
| β ± 0.05 | ~14% | Strong at low mass |
| M_break ± 0.3 dex | ~64% | Critical near break |

**Critical Finding**: M_break uncertainty dominates for systems near the transition mass. This is where we need better calibration.

---

## Part 5: Rotation Curve Predictions

### Spiral Galaxy (M = 5×10¹⁰ M_sun, r = 20 kpc)

| Source | Δv (km/s) | Δv/v |
|--------|-----------|------|
| Ω_m ± 0.007 | ±0.7 | ±0.25% |
| H₀ ± 0.5 | ±0.2 | ±0.06% |
| H₀ tension | +1.9 | +0.65% |
| f_indiff ± 0.5 | ±18 | ±6.1% |

**Conclusion**: Cosmological parameters contribute <1% uncertainty. f_indiff dominates at ~6% per ±0.5.

---

## Part 6: Monte Carlo Results

### 10,000 Sample Error Propagation

| Quantity | Mean | σ | Fractional |
|----------|------|---|------------|
| a₀ | 1.01×10⁻¹⁰ m/s² | 3.8×10⁻¹² | ±3.7% |
| v(20 kpc) | 296.9 km/s | 0.7 km/s | ±0.25% |
| C(a₀) | 0.658 | 0.001 | ±0.13% |

**Conclusion**: Monte Carlo confirms analytical sensitivity estimates. Predictions are stable.

---

## Part 7: H₀ Tension Implications

### If H₀ = 73 km/s/Mpc (local measurement)

| Effect | Magnitude |
|--------|-----------|
| a₀ shift | +8.3% |
| Velocity shift | +0.65% (spiral) |
| C(a) shift | +1.3% at transition |

**Critical Finding**: Synchronism predictions SURVIVE the H₀ tension. The ~8% shift in a₀ translates to <1% shift in observables, which is within typical measurement errors of 5-10%.

---

## Part 8: Implications for Testing

### Hierarchy of Uncertainty

1. **Most Robust**:
   - Spiral galaxy rotation curves
   - Transition regime predictions
   - C(a) functional form

2. **Moderate Uncertainty**:
   - UFD dynamics (f_indiff sensitivity)
   - Systems near M_break

3. **Largest Uncertainty**:
   - Predictions for M < M_break (β uncertainty)
   - M_break location itself

### Testing Strategy

| Test | Robustness | Priority |
|------|------------|----------|
| Void galaxies | HIGH | Session #208 main discriminator |
| Spiral rotation | HIGH | Well-calibrated regime |
| UFD dynamics | MODERATE | f_indiff uncertainty |
| TDGs | MODERATE | f_indiff = 0 prediction |
| EFE tests | HIGH | Structural (no parameters) |

---

## Part 9: Nova's Question Answered

> "Explore the parameter sensitivity of Synchronism's predictions—how stable are results under small perturbations of A, B, γ?"

### Answer:

1. **Cosmological parameters (Ω_m, H₀)**: Very stable (~1% on observables)

2. **Structural parameters (φ)**: Exact - no uncertainty

3. **f_indiff parameters (A, β, M_break)**:
   - A: Linear sensitivity (~13% per ±5)
   - β: Strong at low mass (~14% per ±0.05)
   - M_break: Critical near break (~64% per ±0.3 dex)

4. **Overall**: Predictions are ROBUST to cosmological uncertainties. f_indiff calibration is the priority for improving precision.

---

## Files Created

- `simulations/session213_parameter_sensitivity.py`
- `simulations/session213_parameter_sensitivity.png`
- `Research/Session213_Parameter_Sensitivity.md`

---

## Sessions #199-213 Progress

| Session | Topic | Key Finding |
|---------|-------|-------------|
| #210 | f_indiff theory | Resonance Threshold Model |
| #211 | M_break | First principles derivation |
| #212 | MOND comparison | Convergence/divergence mapped |
| #213 | Sensitivity | **Predictions robust to parameters** |

---

*"Synchronism's predictions are like a well-tuned instrument - small perturbations in the tuning don't change the music."*
