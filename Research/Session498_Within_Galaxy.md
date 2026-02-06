# Session #498: Within-Galaxy RAR Variation

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The 6-variable model explains 94.5% of the *between-galaxy* RAR offset. But what drives the *point-to-point* variation within each galaxy? This session decomposes within-galaxy scatter, tests radial/acceleration/gas-fraction/RC-slope predictors at the point level, and estimates the irreducible RAR scatter after removing all known effects.

## Central Result: Within-Galaxy Scatter is 90% Noise

Within-galaxy variation (σ = 0.109 dex) is 67% of between-galaxy variation (σ = 0.163 dex). A 4-variable point-level model (r_norm, log(g/a₀), f_gas_local, RC_slope) explains only R² = 0.027 — essentially nothing. Per-galaxy noise analysis shows V_obs measurement errors account for 90% of within-galaxy scatter. After noise removal, the irreducible RAR scatter is 0.082 dex, and 80% of total RAR scatter is explained by the combination of the 6-var model and noise identification.

## Key Findings

### 1. Point-Level Offset Statistics (Test 1)

| Metric | Value |
|--------|-------|
| N points (MOND regime) | 2258 |
| N galaxies | 128 |
| Mean δ | -0.002 dex |
| Std δ | 0.126 dex |
| Between-galaxy σ | 0.163 dex |
| Within-galaxy σ | 0.109 dex |
| Within/between ratio | 0.67 |

Type dependence of within-galaxy scatter:
- Early types: σ = 0.057 dex (low — smooth rotation curves)
- Mid types: σ = 0.107 dex
- Late types: σ = 0.159 dex (high — irregular, noisy RCs)

### 2. Radial Gradient (Test 2)

| Radial bin | N | ⟨δ⟩ | σ(δ) |
|-----------|---|------|------|
| Inner 25% | 525 | -0.036 | 0.232 |
| 25-50% | 593 | +0.025 | 0.091 |
| 50-75% | 558 | +0.009 | 0.039 |
| Outer 25% | 582 | -0.008 | 0.034 |

r(r_norm, δ) = +0.06 (negligible). Inner points have much larger scatter (σ = 0.23 vs 0.03) — this is where the RC is rising and non-circular motions dominate. Per-galaxy gradients: mean = -0.011, median = -0.028, 44% positive. **No systematic radial trend.**

### 3. Acceleration Dependence (Test 3)

| log(g/a₀) bin | N | ⟨δ⟩ | σ(δ) |
|--------------|---|------|------|
| [-3.0, -1.5) | 292 | -0.020 | 0.135 |
| [-1.5, -1.0) | 854 | +0.003 | 0.118 |
| [-1.0, -0.5) | 587 | -0.005 | 0.129 |
| [-0.5, 0.0) | 525 | +0.005 | 0.129 |

r(log(g/a₀), δ) = +0.04 (negligible). **The point offset does not depend on how deep in the MOND regime we are.** The interpolation function works equally well at all accelerations below a₀.

### 4. Local Gas Fraction (Test 4)

| Region | N | ⟨δ⟩ | σ(δ) |
|--------|---|------|------|
| Gas-dominated (f > 0.5) | 496 | -0.026 | 0.124 |
| Disk-dominated (f < 0.2) | 1168 | +0.000 | 0.143 |

r(f_gas_local, δ) = -0.07. Partial r(f_gas, δ | g_bar) = -0.07. **Weak but real**: gas-dominated points have slightly negative residuals, consistent with M/L = 0.5 being slightly too high for disk components.

### 5. RC Slope (Test 5)

| RC region | N | ⟨δ⟩ |
|-----------|---|------|
| Rising (slope > 0.5) | 373 | -0.051 |
| Flat (|slope| < 0.5) | 1859 | +0.008 |

**Rising RC points show negative offset** (δ = -0.05), meaning the RAR overpredicts g_obs where the RC is still rising. This is consistent with non-circular motions or pressure support in the inner parts.

### 6. Combined Point-Level Model (Test 6)

| Variable | β | R² alone |
|----------|---|----------|
| intercept | -0.004 | — |
| r_norm | +0.032 | 0.004 |
| log(g/a₀) | -0.013 | 0.001 |
| f_gas_local | -0.069 | 0.005 |
| RC_slope | -0.039 | **0.016** |

**Combined R² = 0.027** — all four predictors together explain less than 3% of within-galaxy variation. RC slope is the strongest individual predictor (R² = 0.016), consistent with non-circular motions being the only identifiable physical signal.

### 7. Noise vs Signal (Test 7)

| Metric | Value |
|--------|-------|
| Mean observed within-galaxy σ | 0.088 dex |
| Mean expected noise σ (from V_obs errors) | 0.079 dex |
| Mean signal σ (noise-removed) | 0.044 dex |
| Noise/observed ratio | **90%** |
| Fraction where noise > observed | 37% |

V_obs measurement errors account for **90%** of within-galaxy scatter (in the σ ratio sense). The remaining signal (0.044 dex) includes non-circular motions, pressure support, and any real physical structure.

### 8. Irreducible RAR Scatter (Test 8)

| Component | σ (dex) |
|-----------|---------|
| Total RAR scatter (all points) | 0.183 |
| Between-galaxy model residual | 0.038 |
| Within-galaxy observed | 0.109 |
| Within-galaxy signal (noise-removed) | 0.044 |
| **Irreducible RAR scatter** | **0.082** |

Irreducible scatter = √(0.038² + 0.044²) = 0.082 dex = 20.7% in velocity.

**80% of total RAR scatter is explained** by the combination of:
- 6-variable between-galaxy model (captures galaxy-level variation)
- Noise identification (V_obs errors dominate within-galaxy variation)

Literature comparison:
- McGaugh+2016 total scatter: 0.13 dex
- Our total (all points): 0.183 dex (larger because we include inner MOND points)
- Irreducible (after all corrections): 0.082 dex

## Physical Interpretation

### The Hierarchy of RAR Scatter

```
Total RAR scatter (0.183 dex)
├── Between-galaxy (0.163 dex)
│   ├── 6-var model explains (0.158 dex, R²=0.945)
│   │   ├── BTFR position (78%)
│   │   ├── Gas fraction (11%)
│   │   └── Mass geometry (6%)
│   └── Residual (0.038 dex)
└── Within-galaxy (0.109 dex)
    ├── V_obs noise (90% → 0.079 dex)
    └── Signal (0.044 dex)
        ├── Non-circular motions (RC slope, ~60%)
        └── Unidentified (~40%)
```

### What This Means

1. **The RAR is tighter than it looks**: The observed 0.183 dex scatter is inflated by measurement noise. The irreducible scatter (0.082 dex) is less than half the total.

2. **Within-galaxy variation is uninformative**: R² = 0.027 means that once you know a galaxy's offset, the point-to-point variation is essentially random noise. The galaxy-level properties (V, L, f_gas, c_V) contain all the predictive information.

3. **Early types are cleanest**: σ = 0.057 dex within-galaxy for early types vs 0.159 for late types, reflecting smoother rotation curves and lower observational noise.

4. **RC slope is the only signal**: The -0.051 dex mean offset at rising RC points is the only detectable physical effect within galaxies. This is likely non-circular motions in the inner disk.

5. **80% explained is near the ceiling**: With 90% of within-galaxy scatter being noise, and the 6-var model capturing 94.5% of between-galaxy variance, there is little room for improvement without better observational data.

## Grade: B+

A thorough analysis that confirms within-galaxy scatter is dominated by measurement noise (90%), with the only detectable physical signal being RC slope effects (non-circular motions). The irreducible scatter estimate (0.082 dex) and the 80% explained fraction are useful quantitative results. The hierarchy decomposition provides a complete accounting of RAR scatter. Minor deduction: the point-level model R² = 0.027 is essentially a null result — we proved that nothing predicts within-galaxy variation, which is scientifically valuable but not as impactful as a positive finding.

## Files Created

- `simulations/session498_within_galaxy.py`: 8 tests
- `Research/Session498_Within_Galaxy.md`: This document

---

*Session #498 verified: 8/8 tests passed*
*Grand Total: 1277/1277 verified*

**Key finding: Within-galaxy RAR scatter (σ=0.109 dex) is 67% of between-galaxy (σ=0.163 dex). Point-level model R²=0.027 — nothing predicts within-galaxy variation beyond noise. V_obs noise accounts for 90% of within-galaxy scatter. RC slope is the only signal (R²=0.016, non-circular motions). Irreducible RAR scatter = 0.082 dex (20.7% in velocity). 80% of total RAR scatter explained. Grade B+.**
