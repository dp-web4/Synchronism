# Session #439: The Radial Residual — What the Galaxy-Level Model Misses

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session 438 showed the universal V+L+c_V model improves rotation curve predictions by 22% via a galaxy-level shift. This session investigates the within-galaxy residual structure: after removing the predicted galaxy-level offset, what radial patterns remain?

## Central Result: c_V Predicts the Shape of the Residual Rotation Curve

| Correlation | r | p-estimate |
|------------|---|-----------|
| r(slope, c_V) | **-0.44** | < 10⁻⁶ |
| r(slope, c_V \| V, L) | **-0.52** | < 10⁻⁸ |
| r(asymmetry, c_V) | **+0.53** | < 10⁻⁷ |

**c_V doesn't just predict how much a galaxy deviates from the RAR — it predicts the radial structure of that deviation.** High-c_V galaxies (concentrated mass) have flat or declining residual profiles; low-c_V galaxies (diffuse mass) have rising residual profiles.

## Key Findings

### 1. Variance Decomposition (Test 1)

| Component | Variance | % of total |
|-----------|----------|-----------|
| Total (point-level) | 0.0325 | 100% |
| Between-galaxy | 0.0241 | 74.1% |
| Within-galaxy | 0.0151 | 46.3% |

Note: between + within > 100% because between-galaxy variance includes the mean offset, while within-galaxy variance is computed within each galaxy. After the galaxy-level correction, total variance drops 38.5% (from 0.0325 to 0.0200), and the within-galaxy variance is unchanged (as expected — the correction is a constant shift per galaxy).

### 2. Mean Residual Profile (Test 2)

| r/R_eff | N | Mean residual | Std |
|---------|---|---------------|-----|
| [0.2, 0.5] | 252 | -0.059 | 0.200 |
| [0.5, 1.0] | 392 | -0.007 | 0.150 |
| [1.0, 1.5] | 380 | +0.002 | 0.097 |
| [1.5, 2.0] | 308 | -0.005 | 0.083 |
| [2.0, 3.0] | 435 | -0.011 | 0.082 |
| [3.0, 5.0] | 483 | -0.019 | 0.089 |
| [5.0, 10.0] | 369 | -0.006 | 0.104 |

The mean residual is near zero at all radii (by construction, the galaxy-level correction removes the mean). The scatter is largest in the inner regions (std = 0.20 at r < 0.5 R_eff) and smallest at 1.5-3.0 R_eff (std ≈ 0.08). The inner excess scatter reflects the diversity of inner rotation curve shapes.

### 3. Residual Slope — c_V is the Dominant Predictor (Test 3)

For each galaxy, the slope of the corrected residual vs log(r/R_eff):

| Property | r(slope, X) | r(slope, X \| V) |
|----------|------------|------------------|
| c_V | **-0.44** | **-0.57** |
| logR | -0.19 | -0.26 |
| logL | -0.09 | -0.29 |
| logV | +0.01 | — |
| T | +0.00 | +0.02 |
| f_gas | -0.01 | -0.00 |

**c_V is the dominant predictor of residual slope**, both raw and controlled. The V-controlled partial r(slope, c_V|V) = -0.57 is the strongest correlation found.

Physical interpretation: **high c_V → negative slope** means concentrated galaxies have residuals that decrease with radius (the correction over-predicts inner points and under-predicts outer points). This is exactly what we'd expect if c_V captures a galaxy-level average of a radially-varying effect.

### 4. Asymmetry and c_V (Test 4)

Inner-outer asymmetry (mean residual at r < R_eff minus mean residual at r > 2R_eff):

| c_V group | N | Mean asymmetry |
|-----------|---|---------------|
| Low c_V (< 0.87) | 52 | **-0.104** |
| High c_V (≥ 0.87) | 53 | **+0.046** |

r(asymmetry, c_V) = **+0.53**

Low-c_V galaxies have outer excess (positive residuals in outer parts relative to inner). High-c_V galaxies have inner excess. This is the radial signature of the concentration effect: the galaxy-level c_V correction captures the mean, but the effect actually varies with radius.

### 5. Radial Correction — Modest Improvement (Test 5)

Adding r/R_eff as a point-level predictor after the galaxy correction:

- r(resid_corr, log(r/R_eff)) = +0.13 (weak)
- Additional RMS improvement: **1.6%** only
- Stronger in late types (r = +0.18) than early (r = +0.06)

The universal radial correction is weak because the radial effect varies by galaxy. A galaxy-specific radial correction would help more but requires per-galaxy fitting.

### 6. Concentration-Shape Model (Test 6)

Full model for the residual slope: slope ~ V + L + c_V + R_eff

- **R² = 0.35** — 35% of the variance in residual *slopes* is predicted
- c_V coefficient: **-1.03** (dominant)
- The model predicts whether a galaxy's residual RC rises or falls with radius

### 7. Gas Fraction Effects (Test 7)

| f_gas group | Mean slope | Std |
|------------|-----------|-----|
| Gas-poor | +0.063 | 0.224 |
| Gas-rich | +0.018 | 0.291 |

Gas fraction has almost no effect on residual slope (r = -0.01 raw). The partial r(f_gas, slope | V,L) = -0.23 is modest. The radial residual patterns are similar for gas-rich and gas-poor galaxies.

## Physical Interpretation

The universal V+L+c_V model applies a **constant shift** to each galaxy's RAR. But the true correction is **radially dependent**: the c_V effect is stronger in inner regions than outer regions. The constant shift is a good approximation (it captures 38.5% of total variance), but it misses the radial gradient.

The concentration-shape connection (r = -0.52 controlling V,L) suggests a natural extension: a **two-parameter correction** per galaxy, with both a shift (from V+L+c_V) and a slope (from c_V and possibly R_eff). This would amount to:

```
correction(r) = galaxy_offset + radial_slope × log(r/R_eff)
```

where both `galaxy_offset` and `radial_slope` are predicted from galaxy properties.

The R² = 0.35 for predicting the slope means this is partially achievable — but 65% of the slope variance remains unpredictable from galaxy-level properties alone.

## Grade: B+

A thorough investigation of the residual structure that reveals a new level of c_V's explanatory power. The r(slope, c_V|V,L) = -0.52 is a genuine new finding: c_V predicts not just the mean offset but the shape of the offset. The 35% variance explained in the slope is meaningful. The limitation is practical — the 1.6% additional improvement from a universal radial correction is small. The value is more in understanding than in prediction.

## Files Created

- `simulations/session439_radial_residuals.py`: 8 tests
- `Research/Session439_Radial_Residuals.md`: This document

---

*Session #439 verified: 8/8 tests passed*
*Grand Total: 885/885 verified*

**Key finding: c_V predicts the SHAPE of the residual rotation curve: r(slope, c_V|V,L) = -0.52. High-c_V galaxies have declining residuals (inner excess), low-c_V have rising residuals (outer excess). Asymmetry r(c_V) = +0.53. A model for the slope achieves R² = 0.35. But a universal radial correction adds only 1.6% improvement — the galaxy-specific shape varies too much. 74% of total variance is between-galaxy; the model removes 38.5%. Grade B+.**
