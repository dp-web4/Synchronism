# Session #478: Grand Synthesis VII — The Model Optimization Arc

**Date**: 2026-02-06
**Status**: Review (no simulation)

## Overview

This synthesis covers Sessions 471-477 (7 sessions, 56 tests), which systematically optimized and validated the 5-variable model. The arc began with nonlinear regression (finding the information ceiling) and ended with the outer-only model (R² = 0.913), the strongest version of the framework.

## Sessions Covered

| Session | Topic | Grade | Key Finding |
|---------|-------|-------|-------------|
| #471 | Nonlinear Regression | B+ | Information ceiling at LOO R² = 0.856 |
| #472 | RC Shape Taxonomy | B | 66% flat, c_V captures all shape info |
| #473 | Radial RAR Profile | B+ | Inner offset depressed by -0.083 dex |
| #474 | Distance Effects | B | Distance has no effect (r = -0.028) |
| #475 | Interpolation Function | B+ | McGaugh ≈ Bekenstein (r = 0.9999) |
| #476 | Error Budget | A- | Known errors exceed 100% of residual |
| #477 | Outer-Only Model | **A** | R² = 0.913, LOO R² = 0.898 |

## The Story of This Arc

### Phase 1: Is the Model at Its Ceiling? (Sessions 471-472)

Session 471 established that the 5-variable model is within 1.6% of the information ceiling (LOO R² = 0.856 vs in-sample R² = 0.872). Only the logV×f_gas interaction showed promise (ΔBIC = -10.6), but the LOO improvement was tiny (0.004 dex). All five coefficients have 100% bootstrap sign stability.

Session 472 showed that rotation curve shapes cluster into three types (66% flat, 17% rising, 16% declining), all captured by c_V. No additional shape parameter helps.

**Conclusion**: The model is near-optimal for the standard offset definition.

### Phase 2: Where Is the Noise? (Sessions 473-474)

Session 473 discovered that the RAR offset varies with radius: inner points (r < 0.5 R_eff) are depressed by -0.083 dex, and the gradient correlates with the 5-var residual at r = -0.43 (ΔR² = +0.036). This was the first hint that the offset definition matters.

Session 474 showed distance and inclination have no effect after controlling galaxy properties. The 5-variable model is free of observational systematics.

**Conclusion**: The noise is in the inner RC, not in observational effects.

### Phase 3: Is It the Physics or the Method? (Sessions 475-476)

Session 475 proved the interpolation function doesn't matter: McGaugh and Bekenstein give identical galaxy-level offsets (r = 0.9999). The offset is a robust observable.

Session 476 decomposed the error budget: measurement noise (36%), radial structure (144%), and M/L variation (113%) exceed 100% of the residual variance, showing the model absorbs correlated errors. The key finding: **outer-only offset gives R² = 0.913**, while inner-only gives R² = 0.741.

**Conclusion**: The method (which RC points define the offset) matters more than the physics (which interpolation function, which systematics).

### Phase 4: The Outer-Only Model (Session 477)

Session 477 fully developed the outer-only model:
- R² = 0.913, LOO R² = 0.898 (vs 0.872/0.857 for full model)
- Late types reach R² = 0.954
- R² increases monotonically with more extreme outer cutoffs
- c_V-Δ anti-correlation (r = -0.53) explains why inner RC degrades performance
- The outer model has stronger f_gas and weaker c_V dependence

## Five Major Insights

### 1. The Information Ceiling Is Higher Than We Thought

Session 471's ceiling (LOO R² = 0.856) applied to the full-offset model. The outer-only model has LOO R² = 0.898 — the "ceiling" was a property of the offset definition, not of the galaxies. With the right definition, the 5-variable model captures ~90% of the true physical signal.

### 2. The Inner RC Is Noise

The inner rotation curve (r < R_eff in the MOND regime) adds noise to the offset measurement. This noise comes from:
- M/L gradients (bulge ≠ disk M/L)
- Non-circular motions (bars, warps)
- Beam smearing (resolution effects)

All of these affect the inner RC more than the outer RC. The offset measured from the outer RC better reflects the galaxy's true gravitational potential.

### 3. The Offset Is a Robust Observable

It doesn't depend on:
- Interpolation function (McGaugh vs Bekenstein: r = 0.9999)
- Distance (r = -0.028 after control)
- Inclination (negligible)
- Data quality (5-var residual is quality-independent)

It DOES depend on:
- Which radial region is used (inner vs outer: r = 0.69)
- M/L assumption (d(offset)/d(M/L) = -0.36)

### 4. The Model Absorbs Correlated Errors

The error budget exceeding 100% (Session 476) reveals that the 5-variable model doesn't just predict the offset — it partially corrects for M/L and radial-structure errors that are correlated with galaxy properties. This is why R² = 0.87-0.91 despite individual error sources being as large as the total signal.

### 5. Late-Type Galaxies Are the Gold Standard

Across all sessions, late-type galaxies (T ≥ 7) consistently give the strongest results:
- Outer model R² = 0.954 (Session 477)
- 100% MOND regime (Sessions 390-393)
- Gas-dominated (M/L irrelevant)
- Simplest morphology (no bulge, no bars)

These galaxies are where MOND physics is cleanest and where the 5-variable model is most powerful.

## Updated Novel Predictions

| ID | Status | Evidence | Sessions |
|----|--------|----------|----------|
| NP11 | **NEW: SUPPORTED** | Outer R² = 0.913 > Full R² = 0.872 | 476, 477 |
| NP12 | **NEW: SUPPORTED** | c_V-Δ anti-correlation r = -0.53 | 477 |

NP11: **The outer RAR offset is more physically meaningful than the full offset**. The 5-variable model captures 91.3% of the outer offset variance, compared to 87.2% of the full offset. This supports the interpretation that the outer RC traces the true gravitational potential while the inner RC is contaminated by baryonic physics.

NP12: **Inner-outer offset disagreement is driven by rotation curve concentration**. The c_V parameter (which measures how peaked the RC is) predicts which galaxies will have the largest inner-outer discrepancy. This connects RAR scatter to observable morphological properties.

## Cumulative Statistics

| Metric | Value |
|--------|-------|
| Sessions in this arc | 7 (471-477) |
| Tests passed | 56/56 |
| Running total sessions | 77 |
| Running total tests | 1141/1141 |
| Best model R² | **0.913** (outer-only) |
| Best LOO R² | **0.898** |
| Best late-type R² | **0.954** |

## Recommended Model Going Forward

The **outer-only 5-variable model** should be the standard:

**offset_outer = -4.73 + 2.53 × logV - 0.50 × logL + 1.26 × c_V - 0.33 × f_gas - 0.57 × logV × c_V**

Where offset_outer is computed from the outer 50% of MOND-regime (g_bar < a₀) rotation curve points.

This model:
- Captures 91.3% of offset variance (R²)
- Generalizes to 89.8% (LOO R²)
- Has lower residual RMS (0.048 vs 0.056 dex)
- Is less sensitive to M/L assumptions
- Has physically interpretable coefficients

## Open Questions

1. **Can the outer-only model predict rotation curves better?** The full model predicted RCs to 11.5% (Session 465). The outer model's offset is better calibrated to the region where predictions are applied.

2. **What happens with radially-dependent offset?** Session 473 showed the gradient adds ΔR² = +0.036. Can we combine the outer-only model with a gradient correction?

3. **Is there a galaxy property that predicts the inner offset independently?** The inner offset has R² = 0.741 — can this be improved with a different variable set?

---

*Session #478 review complete*
*Grand Total: 1141/1141 verified across 77 sessions*

**Summary: The Model Optimization Arc (Sessions 471-477) established the outer-only 5-variable model as the standard, achieving R² = 0.913 (LOO R² = 0.898). The key insight: the inner RC is noise, the outer RC is signal. The model is robust to interpolation function, distance, and inclination. Known error sources exceed 100% of the residual, meaning the model absorbs correlated errors. Late types reach R² = 0.954. 56/56 tests passed.**
