# Session #500: Rotation Curve Shape Classification

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The 6-variable model uses c_V (velocity concentration) as a single-parameter summary of rotation curve shape. But rotation curves come in distinct morphological types: rising, flat, declining, and humpy. Does the full RC shape carry information about the RAR offset beyond what c_V captures?

## Central Result: c_V Captures All Shape Information

RC shape classification yields 4 classes (Rising=39, Flat=74, Declining=10, Humpy=5). Multiple shape parameters (asymmetry, peak ratio, outer slope, baryonic peak position) correlate with the offset, but after controlling for the 6-var model predictors, only 3 show significant partial correlations: asymmetry (r = +0.24), r_bar_max_frac (r = -0.26), and bar_peak_ratio (r = +0.19). However, adding the best shape parameter (asymmetry) to the 6-var model improves LOO R² by only +0.003. Adding the top 2 shapes yields ΔLOO = +0.005. **c_V already encodes the essential shape information.**

## Key Findings

### 1. RC Shape Parameters (Test 1)

| Parameter | Mean | Std | r(offset) | r(c_V) |
|-----------|------|-----|-----------|--------|
| peak_ratio | 1.06 | 0.07 | -0.06 | **+0.62** |
| asymmetry | 1.15 | 0.89 | +0.21 | -0.08 |
| outer_slope_norm | 0.03 | 0.20 | +0.06 | **-0.54** |
| r_max_frac | 0.65 | 0.34 | +0.10 | **-0.68** |
| outer_roughness | 0.01 | 0.01 | -0.03 | -0.16 |
| dm_gradient | 0.51 | 0.74 | +0.19 | +0.03 |
| bar_peak_ratio | 1.40 | 0.49 | +0.17 | **+0.53** |
| r_bar_max_frac | 0.51 | 0.34 | -0.21 | **-0.47** |

Most shape parameters correlate strongly with c_V (|r| > 0.5 for peak_ratio, outer_slope, r_max_frac, bar_peak_ratio). **c_V is a good summary of RC shape.** The exceptions — asymmetry and dm_gradient — carry independent information (r(c_V) < 0.1).

### 2. Shape Classification (Test 2)

| Class | N | % | ⟨offset⟩ | ⟨c_V⟩ | ⟨Type⟩ |
|-------|---|---|----------|-------|--------|
| Rising | 39 | 31% | -0.012 | 0.69 | 7.3 |
| Flat | 74 | 58% | -0.049 | 0.86 | 6.5 |
| Declining | 10 | 8% | -0.050 | 1.15 | 3.4 |
| Humpy | 5 | 4% | -0.047 | 1.11 | 2.8 |

Most galaxies have flat RCs (58%). Declining/humpy RCs are rare (12%) and concentrated in early types. Rising RCs are late-type, low-concentration galaxies.

6-var model residuals by class: DECLINING shows the only notable mean residual (+0.022 dex, i.e., the model underpredicts declining-RC galaxies). But N=10 is too small for strong conclusions.

### 3. Shape → 6-var Residuals (Test 3)

Only 2 shape parameters significantly correlate with 6-var model residuals:
- **asymmetry**: r = +0.23 (p < 0.05) — galaxies with steeper inner rise have positive residuals
- **r_bar_max_frac**: r = -0.19 (p < 0.05) — galaxies where baryonic peak is at large radius have negative residuals

### 4. Shape Beyond c_V (Test 4)

Partial correlations after controlling for all 6 model predictors:

| Parameter | Partial r | Significant? |
|-----------|-----------|-------------|
| **asymmetry** | **+0.24** | **YES** |
| **r_bar_max_frac** | **-0.26** | **YES** |
| **bar_peak_ratio** | **+0.19** | **YES** |
| peak_ratio | +0.13 | no |
| outer_slope_norm | -0.13 | no |
| dm_gradient | +0.09 | no |

Three shape parameters carry information beyond the 6-var model. Asymmetry (how much steeper the inner rise is versus the outer profile) and baryonic RC peak position are the most informative shape features not captured by c_V.

### 5. Inner vs Outer RC Slope (Test 5)

| Correlation | Value |
|------------|-------|
| r(inner_slope, offset) | +0.22 |
| r(outer_slope, offset) | +0.06 |
| r(inner_slope, c_V) | -0.24 |
| r(outer_slope, c_V) | **-0.68** |
| Partial r(inner_slope, offset \| c_V) | +0.20 |
| Partial r(outer_slope, offset \| c_V) | +0.01 |

**c_V is essentially the outer slope** (r = -0.68). The inner slope carries independent information (partial r = +0.20 after controlling for c_V). Inner and outer slopes are weakly correlated (r = +0.20), confirming they are distinct features.

### 6. Baryonic vs Observed RC Shape (Test 6)

| Comparison | r |
|-----------|---|
| Baryonic vs observed peak ratio | **+0.74** |
| Baryonic vs observed peak position | **+0.73** |

Baryonic and observed RC shapes are tightly correlated (r ≈ 0.74). The shape discrepancy (observed - baryonic peak ratio) averages -0.34 — observed RCs are flatter than baryonic (dark matter/MOND flattens the RC). This discrepancy weakly anticorrelates with offset (r = -0.20) and 6-var residual (r = -0.14).

### 7. Shape by Morphological Type (Test 7)

| Type | Best shape predictor of residual | r |
|------|--------------------------------|---|
| Early (T<4) | dm_gradient | -0.35 |
| Mid (4≤T<7) | asymmetry | +0.19 |
| Late (T≥7) | bar_peak_ratio | +0.34 |

Different types have different shape-residual connections:
- **Early types**: Dark matter gradient matters — galaxies with steeper DM dominance gradients have negative residuals
- **Late types**: Baryonic peak ratio matters — more concentrated baryonic distributions correlate with positive residuals
- **Late types have NO declining or humpy RCs** — only Rising (25) and Flat (35)

### 8. Augmented Model (Test 8)

| Addition | R² | ΔR² | LOO R² | ΔLOO |
|----------|-----|------|--------|------|
| 6-var baseline | 0.945 | — | 0.938 | — |
| + asymmetry | 0.948 | +0.003 | **0.941** | +0.003 |
| + r_bar_max_frac | 0.949 | +0.004 | 0.940 | +0.002 |
| + bar_peak_ratio | 0.947 | +0.002 | 0.939 | +0.001 |
| + asymmetry + r_bar_max_frac | 0.951 | +0.006 | **0.943** | **+0.005** |

Adding asymmetry alone: ΔLOO = +0.003. Adding the top 2 shape parameters: ΔLOO = +0.005. This is marginal improvement — the 6-var model with c_V already captures the essential shape information.

**The 8-variable model (6-var + asymmetry + r_bar_max_frac)** achieves R² = 0.951, LOO R² = 0.943, but at the cost of 2 extra parameters. Given N=128 and the noise floor (28% of variance from Session #491), this is not a meaningful improvement.

## Physical Interpretation

### Why c_V Works

c_V = V(R_eff) / V_flat captures the essence of RC shape because:
1. **R_eff marks the transition**: Inside R_eff, the RC is rising (baryons dominate). Outside, it flattens (MOND/DM dominates).
2. **The ratio compresses shape**: A high c_V means the RC reaches V_flat early (compact baryonic distribution). A low c_V means a slowly rising RC (extended baryons).
3. **Outer slope ≈ c_V**: r(outer_slope, c_V) = -0.68 confirms c_V encodes the outer RC behavior.

### What c_V Misses

The only independent shape information is:
1. **Asymmetry** (steepness of inner rise vs outer profile): Partial r = +0.24. This captures how abruptly the RC rises — a feature not determined by c_V alone.
2. **Baryonic peak position** (r_bar_max_frac): Partial r = -0.26. Where the baryonic RC peaks determines the MOND transition radius.

Both are physically meaningful but contribute only ΔLOO = 0.005 — at the noise floor.

### Shape Class Implications

The dominance of Flat (58%) and Rising (31%) RCs reflects the SPARC sample composition (field galaxies, mostly spirals). Declining and humpy RCs (12%) are exclusively early types with concentrated baryonic distributions — they have already reached the MOND regime at smaller radii.

## Grade: B+

A thorough shape analysis confirming that c_V is an effective single-parameter summary of RC shape. The identification of asymmetry and baryonic peak position as independent (but marginal) predictors is interesting. The type-specific shape analysis reveals different physics in early vs late types. However, the practical conclusion — c_V captures all shape information for the 6-var model — is a null result in terms of model improvement. The +0.005 ΔLOO from the best augmented model is at the noise floor. Fitting Session #500 as a comprehensive but ultimately confirmatory result.

## Files Created

- `simulations/session500_rc_shape.py`: 8 tests
- `Research/Session500_RC_Shape.md`: This document

---

*Session #500 verified: 8/8 tests passed*
*Grand Total: 1293/1293 verified*

**Key finding: c_V captures essentially all RC shape information relevant to the RAR offset. Shape classification: Rising=39, Flat=74, Declining=10, Humpy=5. Three shape parameters (asymmetry r=+0.24, r_bar_max_frac r=-0.26, bar_peak_ratio r=+0.19) carry independent partial correlations beyond the 6-var model. But the best augmented model improves LOO R² by only +0.005 (from 0.938 to 0.943). Inner slope carries information beyond c_V (which is essentially the outer slope, r=-0.68). Late types have only Rising/Flat RCs. Grade B+.**
