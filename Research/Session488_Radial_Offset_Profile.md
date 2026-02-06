# Session #488: The Radial Offset Profile — Within-Galaxy RAR Structure

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The 6-variable model predicts galaxy-level outer offsets at R² = 0.945. But within each galaxy, the offset varies with radius. This session characterizes the radial structure and tests whether it contains independent information.

## Central Result: The 6-Var Model R² Increases Monotonically with Radius (0.67 → 0.94)

The 6-variable model's prediction quality improves systematically from inner to outer radii: R² = 0.67 at 0-25% of R_max, rising to R² = 0.94 at 75-100% of R_max. The radial gradient is dominated by c_V (r = -0.53) and adds nothing to the outer model.

## Key Findings

### 1. Radial Profiles by Type (Test 1)

**Early types**: Offset is slightly positive at all radii (0.02-0.08 dex), peaks in the innermost bin (0.077 at 0-2 R_eff) — bulge contamination.

**Mid types**: Offset is near zero inner (0.02) and slightly negative outer (-0.04) — a mild decreasing trend.

**Late types**: Offset is negative inner (-0.08) and positive outer (+0.15 at 8-10 R_eff) — a STRONG increasing trend. The outermost late-type points deviate furthest from the RAR.

Mean gradient: -0.009 dex per R_max (essentially flat overall). 44.5% have positive gradient.

### 2. Inner-Outer Difference (Test 2)

| Metric | Value |
|--------|-------|
| ⟨outer - inner⟩ | -0.008 dex |
| σ(outer - inner) | 0.133 dex |
| Median |Δ| | 0.074 dex |
| **r(Δ, c_V)** | **-0.526** |

**c_V is the dominant predictor of inner-outer disagreement** (r = -0.53). Concentrated galaxies (high c_V, typically early types with bulges) have inner offset > outer offset. This is the M/L gradient effect: the inner stellar disk has higher M/L than assumed, inflating the inner offset.

### 3. Gradient Predictors (Test 3)

| Property | r(gradient, X) |
|----------|---------------|
| **c_V** | **-0.530** |
| logL | -0.178 |
| f_gas | +0.094 |
| logV | -0.068 |
| T | +0.032 |

**c_V explains 28% of gradient variance.** The 6-variable model predicts the gradient at R² = 0.53 (LOO R² = 0.46). The gradient does NOT predict the 6-var model residual (r = +0.12).

### 4. Model Performance by Radial Bin (Test 4)

| Radial bin | ⟨offset⟩ | σ | 6-var R² |
|-----------|----------|---|----------|
| 0.00-0.25 R_max | -0.043 | 0.235 | **0.666** |
| 0.25-0.50 R_max | -0.005 | 0.164 | **0.758** |
| 0.50-0.75 R_max | -0.028 | 0.153 | **0.888** |
| 0.75-1.00 R_max | -0.049 | 0.163 | **0.936** |

**The 6-variable model's predictive power increases monotonically from inner to outer radii.** The outermost 25% of the MOND RC is predicted at R² = 0.936 — nearly as well as the outer-half offset (R² = 0.945). The inner 25% is at R² = 0.67.

This confirms the outer-only approach: the signal-to-noise improves outward because M/L gradients, beam smearing, and non-circular motions all diminish.

### 5. Within-Galaxy Scatter (Test 5)

| Type | ⟨within σ⟩ |
|------|-----------|
| Early | 0.036 dex |
| Mid | 0.075 dex |
| Late | **0.114 dex** |
| All | 0.087 dex |

Within-galaxy scatter is 53% of between-galaxy scatter. Late types have 3× the within-scatter of early types. This is because:
1. Late types have chaotic inner kinematics (non-circular motions, bars, warps)
2. Their MOND-regime RC spans a wider radial range (deeper into MOND)
3. Gas distribution irregularities create point-to-point scatter

### 6. Inner-Outer Cross-Prediction (Test 6)

| Target | 6-var R² | LOO R² |
|--------|----------|--------|
| **Outer offset** | **0.945** | **0.938** |
| Full offset | 0.882 | 0.868 |
| Inner offset | 0.741 | 0.707 |

r(inner, outer) = 0.692 — moderate correlation. The 6-var model predicts the outer offset 33% better than the inner offset (LOO: 0.938 vs 0.707). Adding the outer offset as a 7th predictor barely helps predict the inner (ΔR² = +0.011).

### 7. Gradient-Residual Connection (Test 7)

The gradient adds only ΔLOO = +0.001 to the 6-var model. Galaxies with positive and negative gradients have indistinguishable residuals (⟨resid⟩ ≈ 0 for both groups).

**The radial gradient is NOT an independent source of information for the outer-offset model.** The 6-var model already captures everything the gradient would tell us, primarily through c_V.

## Physical Interpretation

### Why R² Increases Outward

The monotonic R² increase (0.67 → 0.94 from inner to outer) has four physical causes:

1. **M/L gradients**: Inner disks have variable M/L (star formation history, dust, metallicity). At our fixed M/L = 0.5, the inner offset is contaminated by systematic M/L errors.

2. **Beam smearing**: HI observations have finite angular resolution (~10"). The inner RC is most affected, artificially flattening rotation velocities.

3. **Non-circular motions**: Bars, spiral arms, and warps create significant velocity perturbations in the inner disk (5-20 km/s). These diminish outward.

4. **MOND depth**: At outer radii, galaxies are deeper into the MOND regime (g/g† is smaller), making the RAR prediction more sensitive to the physics and less to baryonic details.

### The c_V-Gradient Anti-Correlation

r(gradient, c_V) = -0.53 means:
- **High c_V** (concentrated RC, typically early types): offset DECREASES outward
- **Low c_V** (flat RC, typically late types): offset is approximately flat or increasing

This is consistent with the M/L gradient effect: concentrated galaxies have bulge-dominated centers where our fixed M/L underestimates the stellar mass, creating an artificially positive inner offset. The outer RC, dominated by the flat disk + gas, has a more accurate offset.

### The 53% Within-Galaxy Scatter

The within-galaxy scatter (0.087 dex mean) is roughly half the between-galaxy scatter (0.163 dex). This means:
- **~47% of total variance** is between galaxies (the "galaxy-level offset" we model)
- **~53% of total variance** is within galaxies (radial structure + noise)

Our 6-var model explains 94.5% of the between-galaxy variance. The within-galaxy variance is dominated by measurement noise and M/L structure that varies with radius.

## Grade: B+

A thorough exploration that beautifully confirms the outer-only approach. The monotonic R² increase (0.67 → 0.94) is the cleanest demonstration yet of why the outer RC is superior. The c_V-gradient correlation (r = -0.53) provides physical interpretation. The within/between scatter decomposition (53%/47%) characterizes the fundamental structure of the data. Slightly lower grade because the session is primarily confirmatory — the radial gradient contains no independent information beyond what's in the 6-var model.

## Files Created

- `simulations/session488_radial_offset_profile.py`: 8 tests
- `Research/Session488_Radial_Offset_Profile.md`: This document

---

*Session #488 verified: 8/8 tests passed*
*Grand Total: 1213/1213 verified*

**Key finding: 6-var R² increases monotonically from 0.67 (inner 25%) to 0.94 (outer 25%). c_V dominates the radial gradient (r = -0.53). Within-galaxy scatter = 0.087 dex (53% of between-galaxy). Inner offset poorly predicted (LOO = 0.71 vs outer 0.94). Gradient adds nothing to 6-var model (ΔLOO ≈ +0.001). Grade B+.**
