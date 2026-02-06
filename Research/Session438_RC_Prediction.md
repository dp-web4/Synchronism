# Session #438: Predicting Individual Rotation Curves

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The universal model (V+L+c_V, R²=0.75) predicts galaxy-level RAR offset. This session tests the strongest possible application: can it improve predictions of individual rotation curves V(r)?

The standard RAR predicts V_pred(r) = √(r × g_RAR(g_bar(r))). Our correction modifies this to V_pred_corr(r) = √(r × g_RAR × 10^correction), where correction is the model-predicted offset for that galaxy.

## Central Result: 22% Improvement in Rotation Curve Prediction

| Metric | Standard RAR | Corrected RAR | Improvement |
|--------|-------------|---------------|-------------|
| RMS(log V) | 0.0880 dex | 0.0688 dex | 21.9% |
| RMS(V) | 22.9 km/s | 18.4 km/s | 19.7% |
| RMS(ΔV/V) | 25.5% | 20.0% | 21.6% |
| RMS(log g) | 0.1817 | 0.1425 | 21.5% |
| MOND-only RMS(log V) | 0.0881 dex | 0.0638 dex | 27.5% |

The corrected RAR reduces rotation curve scatter by ~22% across 2850 data points from 128 galaxies. The improvement is larger in the MOND regime (27.5%).

## Key Findings

### 1. Galaxy-by-Galaxy (Test 3)

87/128 galaxies improved (68%), 41 worsened (32%). Median improvement: +11.8%.

**Top 5 most improved:**
| Galaxy | Type | N | RMS_std | RMS_corr | Improvement |
|--------|------|---|---------|----------|-------------|
| F561-1 | 9 | 6 | 0.230 | 0.044 | +80.8% |
| NGC4559 | 6 | 32 | 0.067 | 0.014 | +79.9% |
| NGC5055 | 4 | 28 | 0.060 | 0.012 | +79.5% |
| NGC2841 | 3 | 50 | 0.087 | 0.022 | +75.1% |
| UGC07125 | 9 | 13 | 0.177 | 0.045 | +74.6% |

**Top 5 most worsened:**
| Galaxy | Type | N | RMS_std | RMS_corr | Worsening |
|--------|------|---|---------|----------|-----------|
| NGC3521 | 4 | 41 | 0.016 | 0.046 | -195% |
| NGC0300 | 7 | 25 | 0.032 | 0.076 | -136% |
| NGC6503 | 6 | 31 | 0.019 | 0.039 | -108% |
| NGC3741 | 10 | 21 | 0.042 | 0.077 | -83% |
| NGC0247 | 7 | 26 | 0.028 | 0.051 | -81% |

Notable: the galaxies that worsen had very low standard RMS (< 0.05 dex) — they were already well-predicted by the standard RAR, and the correction over-adjusts them. The galaxies that improve the most had large standard deviations (> 0.06 dex).

### 2. Residual Rotation Curves (Test 4)

Mean residual log(V_obs/V_pred) by radial bin:

| r/R_eff | N | Standard | Corrected |
|---------|---|----------|-----------|
| [0.5, 1.0] | 392 | -0.008 | -0.003 |
| [1.0, 1.5] | 380 | -0.002 | +0.001 |
| [1.5, 2.0] | 308 | -0.007 | -0.002 |
| [2.0, 3.0] | 435 | -0.014 | -0.005 |
| [3.0, 5.0] | 483 | -0.017 | -0.009 |
| [5.0, 10.0] | 369 | +0.010 | -0.003 |

The correction flattens the residual rotation curve, reducing the systematic under-prediction at r = 2-5 R_eff and the over-prediction at r > 5 R_eff. The residuals are closer to zero at all radii.

### 3. LOO Rotation Curve Prediction (Test 5)

Leave-one-out: train the V+L+c_V model on 127 galaxies, predict the 128th galaxy's rotation curve.

- 87/128 improved (68%) — identical to in-sample
- Median improvement: +11.4% (vs +11.8% in-sample)
- Mean improvement: +9.9%

The near-identical LOO performance confirms this is genuine predictive power, not overfitting.

### 4. Inner vs Outer (Test 7)

| Region | N | Standard | Corrected | Improvement |
|--------|---|----------|-----------|-------------|
| Inner (r < 2 R_eff) | 1441 | 0.217 | 0.176 | 18.8% |
| Outer (r > 2 R_eff) | 1409 | 0.136 | 0.097 | 29.1% |

The improvement is larger in outer regions (29% vs 19%). This is consistent with the correction capturing a galaxy-level offset that matters more where the scatter is dominated by between-galaxy variation (outer, flat part of RC) rather than within-galaxy variation (inner, rising part).

## Physical Interpretation

The universal model correction is a **galaxy-level shift** applied uniformly to all points within a galaxy. It captures:

1. **M/L correction** (L component): Each galaxy's assumed M/L=0.5 is wrong. The V+L combination estimates each galaxy's true M/L, shifting its entire RC up or down.

2. **Geometry correction** (c_V component): The algebraic RAR assumes spherical symmetry. c_V captures how concentrated the mass distribution is, correcting for the geometry error.

The 22% improvement represents the between-galaxy variance that the standard RAR ignores. The remaining ~78% of the scatter is within-galaxy (radial profile structure, measurement noise) and cannot be captured by a galaxy-level correction.

The galaxies that worsen are those where the model's prediction disagrees with the galaxy's actual offset. These are typically galaxies already well-matched by the standard RAR (low RMS) — the model over-corrects them.

## Grade: B+

A solid practical demonstration that the universal model improves rotation curve predictions. The 22% improvement, replicated in LOO, is genuine and useful. The inner/outer decomposition and residual RC analysis add insight. The limitation is that this is fundamentally a galaxy-level shift, not a radial correction — the residual RC still has structure that a constant offset cannot capture. A future session could explore radial-dependent corrections.

## Files Created

- `simulations/session438_rc_prediction.py`: 8 tests
- `Research/Session438_RC_Prediction.md`: This document

---

*Session #438 verified: 8/8 tests passed*
*Grand Total: 877/877 verified*

**Key finding: The universal V+L+c_V model improves individual rotation curve predictions by 22% (log V: 0.088→0.069 dex), with 87/128 galaxies improving. LOO confirms: 87/128 improved, median +11.4%. Improvement larger in outer regions (29%) than inner (19%). Galaxies that worsen had very small standard errors (<0.05 dex). The correction is a galaxy-level shift capturing M/L and geometry — the residual RC still has radial structure. Grade B+.**
