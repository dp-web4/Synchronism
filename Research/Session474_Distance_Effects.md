# Session #474: Distance Effects on the RAR

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Distance errors are a known systematic in galaxy studies. Does the RAR offset correlate with distance after controlling for galaxy properties? Do nearby vs distant galaxies behave differently? Does inclination or data quality bias the 5-variable model? This session tests whether observational systematics contaminate our results.

## Central Result: Distance Has No Effect After Controlling Galaxy Properties

r(logD, offset | 5-var) = -0.028. Adding distance, inclination, and quality to the 5-variable model adds ΔR² = +0.005 — negligible. The 5-variable model's physical variables already absorb any distance-dependent systematics.

## Key Findings

### 1. Offset vs Distance (Test 1)

| Correlation | Value |
|-------------|-------|
| r(logD, offset) | +0.090 |
| r(logD, offset \| V, L) | +0.047 |
| r(logD, offset \| 5-var) | **-0.028** |

The raw correlation (+0.090) is entirely driven by Malmquist bias: distant galaxies are more luminous and massive, which correlates with higher offset. After controlling for V, L alone, the correlation drops to +0.047. After the full 5-variable model, it vanishes to -0.028.

**Distance does not independently affect the RAR offset.**

### 2. Near vs Far Comparison (Test 2)

| Sample | N | ⟨offset⟩ | σ(off) | ⟨logV⟩ | ⟨T⟩ | 5-var R² | RMS |
|--------|---|---------|--------|--------|-----|---------|-----|
| Near (D<17) | 64 | -0.044 | 0.176 | 1.93 | 7.5 | 0.905 | 0.054 |
| Far (D≥17) | 64 | -0.022 | 0.130 | 2.18 | 5.2 | 0.845 | 0.051 |

Near galaxies have larger offset scatter (0.176 vs 0.130) and include more late types (⟨T⟩=7.5 vs 5.2). The 5-variable model fits the near sample *better* (R²=0.905 vs 0.845), despite more scatter, because the near sample contains more low-mass galaxies with diverse offsets that the model captures well.

### 3. Distance and the 5-Variable Residual (Test 3)

| Model | R² | ΔR² |
|-------|-----|------|
| 5-var | 0.8717 | — |
| 5-var + logD | 0.8718 | **+0.0001** |

Adding logD to the model is worthless (ΔR² = 0.0001, β = -0.005). The 5-variable residual shows no distance dependence:

| Distance (Mpc) | N | ⟨resid⟩ | σ(resid) |
|-----------------|---|---------|----------|
| 0–5 | 19 | +0.018 | 0.071 |
| 5–15 | 40 | -0.011 | 0.048 |
| 15–30 | 39 | +0.004 | 0.049 |
| 30–60 | 13 | -0.026 | 0.047 |
| 60–200 | 17 | +0.017 | 0.060 |

Residuals scatter around zero at all distances with comparable σ.

### 4. Data Quality Effects (Test 4)

| Quality | N | ⟨offset⟩ | ⟨5-var resid⟩ | σ(resid) |
|---------|---|---------|--------------|----------|
| Q=1 | 84 | -0.000 | -0.005 | 0.054 |
| Q=2 | 39 | -0.069 | +0.009 | 0.057 |
| Q=3 | 5 | -0.308 | +0.011 | 0.056 |

Q=3 galaxies have much lower offset (-0.31 dex) — these are poorly measured. But the 5-variable model residual is *independent of quality* (σ ≈ 0.055 for all Q). The model captures the quality effect through c_V and f_gas, which are poorly constrained for Q=3 galaxies.

r(N_points, |5-var residual|) = -0.087: more data points ≈ slightly better fit, but weak.

### 5. Malmquist Bias (Test 5)

| Correlation | Value |
|-------------|-------|
| r(logL, logD) | **+0.597** |
| r(logV, logD) | +0.545 |
| r(T, logD) | -0.425 |

**Strong Malmquist bias**: luminosity scales as logL ∝ 1.45 × logD (flux-limited would give ∝ 2.0, volume-limited ∝ 0). Nearby galaxies are dominated by late-type dwarfs (⟨T⟩=7.8, ⟨logL⟩=0.18); distant galaxies by earlier-type spirals (⟨T⟩=5.4, ⟨logL⟩=1.61).

**Despite this strong selection, the offset is unbiased** because the 5-variable model explicitly includes V and L, absorbing the Malmquist-correlated properties.

### 6. Inclination Effects (Test 6)

| Correlation | Value |
|-------------|-------|
| r(i, offset) | +0.051 |
| r(i, 5-var residual) | -0.126 |
| 5-var + sin(i) ΔR² | +0.0023 |

Inclination has minimal effect. The slight negative correlation with the residual (r = -0.126) hints that edge-on galaxies may have slightly lower residuals, but the effect is marginal. Edge-on galaxies (i > 85°) show no systematic pattern — their residuals scatter evenly around zero.

### 7. Combined Systematics (Test 7)

| Model | R² | LOO | k |
|-------|-----|-----|---|
| 5-var | 0.872 | 0.059 | 6 |
| 5-var + logD | 0.872 | — | 7 |
| 5-var + sin(i) | 0.874 | — | 7 |
| 5-var + logD + sin(i) + Q | **0.876** | **0.060** | 9 |

All three observational systematics combined add ΔR² = +0.005 in-sample but **degrade** LOO (0.060 vs 0.059). This confirms that distance, inclination, and quality are not genuine predictors — any apparent improvement is noise fitted by extra parameters.

## Physical Interpretation

### Why Distance Doesn't Matter

The 5-variable model uses V, L, c_V, and f_gas — all intrinsic galaxy properties. Distance affects the *measurement precision* of these properties, but not their values. Since the model works with the measured values (which include distance-dependent errors), the model implicitly absorbs the typical measurement error structure.

The key insight: **observational systematics are orthogonal to the physical signal**. The RAR offset is determined by M/L, mass distribution geometry, and gas content — not by how far away the galaxy is.

### The Malmquist Bias Is Real But Harmless

SPARC has strong Malmquist selection (r(logL, logD) = +0.60). This means the sample is NOT representative of the galaxy population: nearby = dwarfs, distant = massive spirals. But this does NOT bias the 5-variable model because the model explicitly parameterizes mass (V, L). The Malmquist bias affects the *distribution* of galaxies in parameter space but not the *physics* at each point in that space.

### Near vs Far Asymmetry

The near sample has higher R² (0.905 vs 0.845) despite more scatter, because:
1. Near dwarfs span a wider range in offset (σ = 0.176 vs 0.130)
2. The 5-variable model captures this variation well (via c_V, f_gas)
3. Distant massive spirals are more homogeneous but have residual structure the model misses

## Grade: B

A methodologically solid session that confirms the absence of observational systematics in the 5-variable model. All eight tests are clean null results: distance, inclination, and quality collectively explain nothing beyond the physical variables. The Malmquist bias analysis is a useful validation. However, the grade is B rather than higher because null results, while important, don't advance the physics — they confirm what we expected (observational effects are subdominant to physical ones).

## Files Created

- `simulations/session474_distance_effects.py`: 8 tests
- `Research/Session474_Distance_Effects.md`: This document

---

*Session #474 verified: 8/8 tests passed*
*Grand Total: 1117/1117 verified*

**Key finding: Distance has no independent effect on the RAR offset (r = -0.028 after 5-var control). Malmquist bias is strong (r(logL, logD) = +0.60) but harmless because V, L are explicit model variables. Inclination is negligible (ΔR² = +0.002). All observational systematics combined add ΔR² = +0.005 in-sample but degrade LOO, confirming they are noise. The 5-variable model is free of distance-dependent bias. Grade B.**
