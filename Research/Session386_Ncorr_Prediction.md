# Session #386: N_corr Quantitative Prediction Test

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session #385 established N_corr = a_char/a₀ as the strongest single predictor of RAR systematic offset (r = +0.48). This session formalizes the prediction, tests its robustness via cross-validation and bootstrap, and investigates the critical Vflat degeneracy.

## Key Result: Genuine But Partially Degenerate Prediction (Grade A-)

N_corr = V²_flat/(R_eff × a₀) predicts RAR offset with R² = 0.226 (p = 2×10⁻¹²), survives cross-validation, and has the correct physical direction. However, it is 97% correlated with Vflat, and R_eff adds no independent information. The prediction is real but may be a restatement of the Tully-Fisher residual effect.

## Detailed Findings

### 1. Calibration

Linear model: offset = +0.105 × log(N_corr) - 0.038
- R² = 0.226 (22.6% variance explained)
- r = +0.476, p = 2.08×10⁻¹²
- RMSE = 0.177 dex, MAE = 0.132 dex

Quartile behavior confirms monotonic trend:
| Quartile | N | Observed offset | Predicted offset |
|---|---|---|---|
| Q1 (lowest N_corr) | 43 | -0.260 | -0.228 |
| Q2 | 42 | -0.031 | -0.064 |
| Q3 | 43 | -0.034 | -0.027 |
| Q4 (highest N_corr) | 43 | +0.020 | +0.013 |

### 2. Cross-Validation

| Method | RMSE | R² |
|---|---|---|
| Full sample | 0.177 | 0.226 |
| Leave-one-out | 0.180 | 0.201 |
| 5-fold mean | 0.171 | -- |

The model is stable: LOO R² = 0.201 vs full-sample R² = 0.226 (minimal overfitting).

### 3. Predictor Comparison (KEY FINDING)

| Predictor | r | p | R² |
|---|---|---|---|
| **log N_corr** | **+0.476** | **2×10⁻¹²** | **0.226** |
| Quality | -0.501 | 5×10⁻¹⁴ | 0.251 |
| Vflat | +0.410 | 5×10⁻⁹ | 0.168 |
| Hubble Type | -0.296 | 5×10⁻⁵ | 0.088 |
| log SB | +0.269 | 3×10⁻⁴ | 0.072 |
| log Luminosity | +0.245 | 1×10⁻³ | 0.060 |

N_corr is the strongest **physically-motivated** predictor. Quality (Q flag) has higher R² but is a data quality indicator, not a physical variable.

**Critical partial correlations**:
- r(log N_corr, offset | Vflat) = +0.268 (p = 0.0003) → N_corr adds info beyond Vflat
- r(Vflat, offset | log N_corr) = +0.051 (p = 0.50) → Vflat adds NOTHING beyond N_corr
- r(log N_corr, offset | Type) = +0.390 (p < 10⁻⁶) → N_corr works within morphological types

**This is the most important result**: N_corr subsumes Vflat's predictive power entirely, but Vflat cannot subsume N_corr's. N_corr is the more fundamental variable.

### 4. Residual Analysis

After removing N_corr prediction:
- r(Type, residual) = -0.024 (n.s.) → N_corr absorbs the type signal
- r(Vflat, residual) = +0.030 (n.s.) → N_corr absorbs the mass signal
- r(Quality, residual) = -0.398 (p < 10⁻⁷) → Quality information NOT captured

The residual-Quality correlation means data quality is the dominant remaining source of offset variation. The top 5 outliers are all Q=2 or Q=3 galaxies.

### 5. Scatter Prediction

N_corr predicts scatter weakly (R² = 0.046) but this operates entirely through roughness:
- r(log N_corr, scatter | roughness) = +0.010 (n.s.)
- r(log N_corr, roughness) = -0.312 (p < 10⁻⁴)

Mediation chain: N_corr → |offset| (R² = 0.217) → scatter (R² = 0.154). Indirect contribution to scatter: ~3.3%.

### 6. The Vflat Degeneracy (CRITICAL WEAKNESS)

r(log N_corr, log Vflat) = **0.970**

N_corr ≡ V²_flat/(R_eff × a₀), so it's dominated by V²_flat. Testing whether R_eff adds information:

| Model | R² |
|---|---|
| log Vflat alone | 0.192 |
| log N_corr alone | 0.226 |
| log Vflat + log R_eff | 0.192 |

R_eff coefficient: β = -0.011, t = -0.25, p = 0.80. **R_eff does NOT add significant information.**

The 3.5 percentage point improvement of N_corr over Vflat comes from the V² nonlinearity, not from the R_eff component. This means:
- N_corr is effectively a nonlinear transform of Vflat
- The "coherence radius" component doesn't contribute measurably
- The prediction could be restated as: log(V²_flat) → offset

### 7. Bootstrap Confidence

10,000 bootstrap resamples:
- Slope: +0.105 ± 0.020, 95% CI: [+0.067, +0.146]
- R²: 0.229 ± 0.062, 95% CI: [0.113, 0.356]
- P(R² > 0.15) = 90.2%

The relationship is robust and replicable.

## Honest Assessment

### Strengths
1. Strongest physically-motivated predictor of RAR offset (R² = 0.226)
2. Derived from Synchronism first principles (γ = 2/√N_corr)
3. Direction correct: higher N_corr → closer to standard RAR
4. **Subsumes Vflat completely** in partial correlation analysis
5. Stable under cross-validation (LOO R² = 0.201)
6. Bootstrap 95% CI excludes zero
7. Absorbs both type and mass signals from residuals

### Weaknesses
1. **97% correlated with Vflat** — essentially a nonlinear mass proxy
2. R_eff contributes nothing measurable — the "coherence length" part is noise
3. Cannot distinguish from Tully-Fisher residual effect
4. R² = 0.226 means 77.4% of offset variance unexplained
5. Quality flag (non-physical) has higher R² than N_corr
6. M/L mismatch could generate the same correlation pattern

### The Degeneracy Problem

The Vflat degeneracy is the fundamental challenge. N_corr = V²/(R × a₀) is dominated by V², and V² alone explains most of what N_corr explains. The Synchronism interpretation (gravitational coherence) and the mundane interpretation (more massive galaxies are better measured / have more baryonic dominance) make the same prediction.

**What would break the degeneracy**:
1. Find galaxies with unusual R_eff/Vflat ratios and test if N_corr or Vflat better predicts their offset
2. Test at fixed Vflat whether R_eff variation predicts offset variation
3. Measure actual velocity correlations (the physical basis of N_corr)

### Grade: A-

The prediction is real, robust, and physically motivated. But the Vflat degeneracy means we cannot claim N_corr is measuring "gravitational coherence" rather than just "galaxy mass in a nonlinear way." The partial correlation analysis (N_corr subsumes Vflat, not vice versa) is suggestive but not conclusive.

## Updated Prediction Scorecard

| ID | Prediction | Status | Evidence |
|----|-----------|--------|----------|
| NP1 | a₀ = cH₀/(2π) | SUPPORTED (94%) | Session #385 |
| NP2 | Morphology → scatter | STRONGLY SUPPORTED | p = 5×10⁻⁶, confound-controlled |
| NP3 | High-z a₀ evolution | UNTESTED | Needs JWST data |
| NP4 | V-shaped scatter | SUGGESTIVE | Session #374 |
| NP5 | Local gravity anomalies | UNTESTED | Needs Gaia DR3 |
| **NP6** | **N_corr → RAR offset** | **SUPPORTED (R² = 0.23)** | **This session** |

## Files Created

- `simulations/session386_ncorr_prediction.py`: 8 tests
- `Research/Session386_Ncorr_Prediction.md`: This document

---

*Session #386 verified: 8/8 tests passed*
*Grand Total: 527/527 verified*

**Key finding: N_corr = V²_flat/(R_eff × a₀) predicts RAR offset with R² = 0.226 (p = 2×10⁻¹²), survives cross-validation, and subsumes Vflat in partial correlations. However, r(log N_corr, log Vflat) = 0.970, and R_eff adds no independent information. The prediction is real and robust but may be a nonlinear restatement of "more massive galaxies sit closer to the standard RAR." Breaking this degeneracy requires testing at fixed Vflat or measuring actual velocity correlations. Grade A-.**
