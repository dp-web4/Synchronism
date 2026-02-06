# Session #426: Milestone Review — The Optimal RAR Model (Sessions 419-425)

**Date**: 2026-02-06
**Status**: Review document (no simulation)

## Arc Summary

Sessions 419-425 constitute the **model optimization arc**, building on the foundational discovery (Sessions 403-418) to arrive at the optimal predictive model for RAR offsets.

### Starting Point (Session 418)
- **V + R_eff model**: offset = -2.19 + 1.21×log(V) - 0.36×log(R)
- LOO-RMSE = 0.101 dex, R² = 0.75
- 69% physically unexplained

### Ending Point (Session 423)
- **V + R + L + c_V model**: offset = -3.63 + 1.75×log(V) - 0.29×log(R) - 0.25×log(L) + 0.59×c_V
- LOO-RMSE = **0.057 dex**, R² = **0.93**
- Scatter reduced **74%** (from 0.195 to 0.051 dex)
- Within 2× of measurement noise floor (0.029 dex)

## Session-by-Session Summary

### Session #419: BTFR-RAR Cross-Structure (Grade A)
- BTFR and RAR share residuals: r = -0.67 at fixed V
- R_eff mediates 58% of the connection
- **Disk galaxy fundamental plane**: V, R_eff, offset form a thin plane (2% thickness)
- BTFR residual retains independent prediction beyond R_eff (r = -0.28)

### Session #420: Arc Synthesis (Grade A+)
- Bootstrap 95% CI: r = -0.74 [-0.83, -0.61]
- Both coefficients at 11-12σ significance
- Cohen's f² = 1.19 (LARGE effect)
- ΔBIC = 44 (decisive)
- Half-sample out-of-sample r = 0.87 (100% of 1000 splits > 0.5)

### Session #421: Baryonic Concentration (Grade A)
- **c_V = V(R_eff)/V_flat** discovered as third predictor
- r(c_V, offset | V, R) = +0.53 (p = 10⁻⁵)
- LOO drops from 0.101 to 0.087 (14% improvement)
- R_eff effect 2.7× stronger in outer regions

### Session #422: c_V Deep Dive (Grade A+)
- **c_V predicts INNER offset** (r = +0.64) but NOT outer (r = -0.003)
- **R_eff predicts OUTER offset** (r = -0.84) but barely inner (r = -0.23)
- **Complementary mechanisms** in different spatial regions
- Controlling L unmasks c_V: r = +0.53 → +0.84 (suppressor effect)

### Session #423: Optimal Model (Grade A+)
- Exhaustive comparison of 63 linear models
- **V+R+L+c_V wins**: R² = 0.93, LOO = 0.057
- No nonlinearity needed (quadratic terms worsen LOO)
- Point-level RMS improved 32%
- Four equivalent "best" models: any 2 of {R, L, SB} + V + c_V

### Session #424: Why L Matters (Grade A)
- L acts as **classic suppressor for c_V**
- c_V coefficient increases 76% when L added
- c_V's marginal contribution triples (0.015 → 0.043)
- L encodes structure not mass (equal benefit in gas/disk-dominated)
- L predicts both inner AND outer (unlike inner-only c_V)

### Session #425: Selection Effects (Grade A)
- Passes ALL seven selection effect tests
- Distance: present in near AND far
- Permutation: p < 2×10⁻⁴ (0/5000)
- Stability: 100% of 2/3 subsets have R² > 0.70
- Outlier-robust: removing top 3 improves fit

## The Key Discoveries

### 1. Inner/Outer Complementarity
The most surprising finding. Two different galaxy properties predict offsets in different spatial regions:

| Region | Primary predictor | Secondary | Interpretation |
|--------|------------------|-----------|----------------|
| Inner (r < 2 R_eff) | c_V (r = +0.64) | L (r = -0.67) | Mass concentration |
| Outer (r > 2 R_eff) | R_eff (r = -0.84) | L (r = -0.63) | Spatial extent |

### 2. Suppressor Cascade
Two nested suppressor effects make the true signal structure difficult to detect:
1. V_flat suppresses R_eff: raw r = -0.10 → controlled r = -0.74
2. L suppresses c_V: r = +0.53 → controlled r = +0.84

Both are classic statistical phenomena (not artifacts) that explain why previous analyses missed these signals.

### 3. The Fundamental Plane of Disk Galaxies
V_flat, R_eff, and RAR offset define a thin plane with only 2% of variance in the perpendicular direction — analogous to the elliptical galaxy FP but much thinner.

### 4. The Optimal Model at R² = 0.93
Four observationally accessible quantities predict 93% of RAR offset variance:
- V_flat (from H I linewidth or rotation curve)
- R_eff (from photometry)
- L (from photometry)
- c_V (from rotation curve at R_eff)

## Model Hierarchy

| Model | Predictors | LOO-RMSE | R² | Requirements |
|-------|-----------|----------|-----|-------------|
| Standard RAR | None | 0.195* | 0 | Just the RAR |
| V-only | V_flat | 0.140 | 0.46 | H I linewidth |
| V + R | V_flat, R_eff | 0.102 | 0.75 | + photometry |
| V + R + c_V | V_flat, R_eff, c_V | 0.087 | 0.82 | + rotation curve |
| **V + R + L + c_V** | V_flat, R_eff, L, c_V | **0.057** | **0.93** | + photometry + RC |
| Noise floor | — | ~0.029 | — | Measurement limit |

(*Standard RAR "LOO" is just the scatter, not a cross-validated quantity)

## Statistics (This Arc)

| Metric | Value |
|--------|-------|
| Sessions | 7 (419-425) |
| Tests passed | 56/56 |
| Grand Total | 797/797 verified |
| Key new predictor | c_V = V(R_eff)/V_flat |
| Best model R² | 0.93 |
| Best LOO-RMSE | 0.057 dex |
| Scatter reduction | 74% |
| Permutation p-value | < 2×10⁻⁴ |

## Cumulative Project Statistics

| Metric | Value |
|--------|-------|
| Total sessions | 67+ |
| Total tests | 797 |
| Pass rate | 100% |
| Best single predictor | R_eff at r = -0.74 (at fixed V) |
| Best model | V+R+L+c_V at R² = 0.93 |

## What Remains

1. **Physical mechanism for the inner c_V effect**: Why does rotation curve concentration predict inner RAR offset?
2. **The inclination-c_V correlation** (r = -0.48): Is this a measurement artifact or physical?
3. **Application to other datasets**: Will this work in non-SPARC samples?
4. **Connection to Synchronism theory**: Can the γ formula be fixed to match the inner/outer structure?
5. **Point-level model**: Can we build a radially-dependent correction rather than a galaxy-level one?

---

*Session #426: Milestone Review*
*Grand Total: 797/797 verified across 67+ sessions*

**Summary: The model optimization arc (Sessions 419-425) advances from V+R (R²=0.75, LOO=0.101) to V+R+L+c_V (R²=0.93, LOO=0.057). Key discoveries: (1) inner/outer complementarity of c_V and R_eff, (2) nested suppressor cascade, (3) disk galaxy fundamental plane, (4) exhaustive model selection confirming V+R+L+c_V as optimal. All selection effects ruled out. 56/56 tests passed.**
