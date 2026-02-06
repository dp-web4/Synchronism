# Session #440: Two-Parameter Correction — Shift + Slope

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session 439 showed c_V predicts the radial slope of residuals (r = -0.52). This session tests whether a two-parameter per-galaxy correction (shift + radial slope) improves rotation curve predictions beyond the one-parameter shift.

## Central Result: The Slope Correction Hurts — Noise Exceeds Signal

| Correction | RMS(log V) | Improvement | LOO Improvement |
|-----------|-----------|-------------|----------------|
| Standard RAR | 0.0908 | — | — |
| One-param (shift) | 0.0713 | **21.5%** | **26.4%** |
| Two-param (shift+slope) | 0.0757 | 16.7% | 23.1% |
| Three-param (+curvature) | 0.0752 | 17.2% | — |

**The two-parameter correction is worse than the one-parameter correction.** Adding the predicted radial slope *increases* the RMS by 6.2%. The slope prediction (R²=0.35) is not accurate enough — the prediction noise outweighs the signal.

## Key Findings

### 1. Slope Model (Test 1)

| Model | Target | R² | LOO-RMSE |
|-------|--------|-----|---------|
| V+L+c_V | Offset | 0.754 | 0.080 |
| V+L+c_V+R | Slope | 0.351 | 0.220 |
| c_V only | Slope | 0.195 | — |

The offset model is strong (R²=0.75) but the slope model is mediocre (R²=0.35). The LOO-RMSE for the slope (0.220) is large relative to the slopes themselves (std = 0.261). The model predicts only the grossest features of the radial gradient.

### 2. Galaxy-by-Galaxy (Test 3)

- Two-param better: 68/128 (53%) — barely above chance
- LOO two-param better: 66/128 (52%) — essentially coin-flip
- Median marginal improvement: +0.9% (near zero)
- Mean marginal improvement: -4.2% (negative due to outlier damage)

The correction helps roughly half the galaxies and hurts the other half, with the hurt galaxies losing more than the helped galaxies gain.

### 3. BIC Comparison (Test 7)

| Model | k | BIC |
|-------|---|-----|
| 1-param | 4 | **-12174** |
| 2-param | 9 | -11793 |
| 3-param | 14 | -11786 |

BIC strongly favors the one-parameter model. The additional slope and curvature parameters do not justify their complexity.

### 4. Which Galaxies Benefit? (Test 6)

The strongest predictor of benefit is |slope| (r = +0.22): galaxies with large actual radial gradients benefit modestly from the slope correction. But even among large-|slope| galaxies, only 62% benefit (vs 44% for small-|slope|). The prediction is too noisy for reliable galaxy-specific improvements.

### 5. Curvature (Test 7)

The curvature model (quadratic in log(r/R_eff)) achieves only R² = 0.11 — curvature is essentially unpredictable from galaxy properties. Adding it does nothing (3-param ≈ 2-param in performance).

## Physical Interpretation

This is an important **negative result** that clarifies the limits of galaxy-property-based corrections:

1. **The galaxy-level shift captures most of the predictable variance.** The V+L+c_V model explains 75% of galaxy-level offsets. The remaining 25% is noise/unpredictable.

2. **The radial structure is partly systematic (c_V predicts it) but too noisy to use.** With R²=0.35, the slope prediction has a prediction error of σ=0.22, which is comparable to the signal itself (std=0.26).

3. **The lesson: 35% R² is not enough for point-level correction.** When applying a radial correction, the noise from the wrong slope actively degrades predictions. The galaxy-level shift is robust because it's strongly predicted (R²=0.75); the radial slope is not.

4. **The optimal correction is the simplest one**: a single number per galaxy, predicted from V_flat, L, and c_V.

## Grade: B

A valuable negative result that definitively answers whether the radial structure can be exploited. The answer is no — the slope prediction is too noisy to help. The BIC analysis is clean, the LOO confirms the negative result, and the physical interpretation is clear. This closes the "can we do better than a galaxy-level shift?" question.

## Files Created

- `simulations/session440_two_param_correction.py`: 8 tests
- `Research/Session440_Two_Param_Correction.md`: This document

---

*Session #440 verified: 8/8 tests passed*
*Grand Total: 893/893 verified*

**Key finding: The two-parameter correction (shift + radial slope) is WORSE than the one-parameter shift. Marginal: -6.2% (LOO: -4.5%). Slope R²=0.35 too noisy for point-level use. BIC strongly favors 1-param (-12174 vs -11793). Only 53% of galaxies benefit. The one-parameter V+L+c_V shift is optimal. Grade B.**
