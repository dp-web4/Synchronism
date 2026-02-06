# Session #505: MOND Boost Model — γ as Primary Predictor

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Following Session #504's finding that γ has partial r = +0.57 with the MOND boost controlling V and L, this session builds a γ-based MOND boost model, compares it with the empirical 6-var approach, and establishes the exact mathematical connection between boost and offset.

## Central Result: Offset = Boost - MOND Correction (r = 0.998)

The boost (log(g_obs/g_bar)) and offset (log(g_obs/g_rar)) are connected through the MOND interpolation function with r = 0.998 and boost coefficient = 0.986 (theory: 1.0). However, γ fails as a boost predictor in practice (R² = 0.03 raw, though partial r = +0.79 after controlling all other variables). The V, L model dominates the boost (R² = 0.64), and the full 6-var model achieves 0.75 (LOO 0.71). The γ-based route to the offset fails because γ's boost predictions are too noisy.

## Key Findings

### 1. MOND Boost Statistics (Test 1)

| Metric | Value |
|--------|-------|
| Mean boost | 0.588 dex |
| Std boost | 0.221 dex |
| Deep MOND theory mean | 0.552 dex |
| r(theory, observed) | +0.67 |
| RMS(theory - observed) | 0.173 |
| r(boost, offset) | **+0.70** |

The deep MOND prediction (boost = 0.5 × log(a₀/g_bar)) captures the broad trend but has large residuals (RMS = 0.17) because many galaxies are not in deep MOND.

### 2. γ as Boost Predictor (Test 2)

| Control variables | Partial r(γ, boost) |
|------------------|---------------------|
| None | +0.16 |
| logV | +0.09 |
| logV, logL | **+0.57** |
| logV, logL, f_gas, c_V | **+0.79** |

**The partial correlation structure is dramatic**: γ's raw correlation with boost is only +0.16, but controlling for confounders reveals a strong signal. At fixed V, L, f_gas, and c_V, γ explains 63% (r² = 0.79²) of the remaining boost variation. This confirms γ carries genuine MOND regime information, but it's masked by the dominant V and L dependence.

### 3. γ-Based Boost Model (Test 3)

| Model | R² | RMS |
|-------|-----|-----|
| γ alone | 0.03 | 0.218 |
| γ + f_gas + c_V | 0.28 | 0.188 |
| log(g_bar/a₀) trivial | 0.45 | 0.164 |
| Deep MOND theory | — | 0.173 |

γ alone is nearly useless for boost prediction (R² = 0.03). Even with corrections (f_gas, c_V), it reaches only R² = 0.28. The trivial predictor (log(g_bar/a₀)) does better at 0.45.

### 4. Model Comparison (Test 4)

| Model | R² | LOO R² |
|-------|-----|--------|
| γ alone | 0.03 | -0.03 |
| γ + f_gas + c_V | 0.28 | 0.21 |
| V, L | **0.64** | **0.62** |
| V, L + γ | 0.76 | — |
| 6-var | **0.75** | **0.71** |

V and L dominate the boost prediction. Adding γ to V, L improves R² by +0.12 — significant but smaller than the +0.19 for the offset (Session #503). The 6-var model captures 75% of boost variance (LOO 71%).

### 5. LOO Validation (Test 5)

| Model | LOO R² |
|-------|--------|
| γ alone | -0.03 (fails) |
| γ + f_gas + c_V | 0.21 |
| V, L | 0.62 |
| 6-var | **0.71** |
| log(g/a₀) trivial | 0.43 |

γ alone actually has negative LOO R² — it overfits the 2-parameter model on boost data. The 6-var model is the clear winner at LOO R² = 0.71.

### 6. Boost vs Offset Coefficients (Test 6)

| Variable | Boost β | Offset β | Ratio |
|----------|---------|----------|-------|
| logV | +1.38 | +1.90 | 0.73 |
| logL | -0.62 | -0.55 | 1.14 |
| c_V | **-1.25** | -0.22 | 5.72 |
| f_gas | -0.10 | -0.45 | 0.22 |
| logV×c_V | **+0.67** | +0.15 | 4.58 |
| logL×f_gas | +0.37 | +0.18 | 2.06 |

**c_V is 5.7× more important for boost than offset.** The c_V/logV×c_V terms dominate the boost model because mass concentration directly determines the g_obs/g_bar ratio — concentrated baryons create larger accelerations at the center, affecting the overall boost amplitude.

### 7. Residual Analysis (Test 7)

The 6-var boost model residual still correlates with log(γ) at r = +0.38 — γ carries information beyond the 6-var predictors for the boost. This is consistent with γ encoding the MOND regime depth, which the 6-var model doesn't directly capture.

r(boost 6var-residual, offset residual) = +0.60 — the residual structures are correlated, confirming that the same underlying physics drives both.

### 8. Boost-Offset Connection (Test 8)

The exact relationship: **offset = boost - log(ν(g_bar/a₀))**

| Connection | Value |
|-----------|-------|
| offset = a + b×boost + c×log(g/a₀) | R² = **0.995** |
| Boost coefficient | 0.986 (theory: 1.0) |
| Theoretical connection | r = **0.998** |
| RMS of theoretical prediction | **0.014 dex** |

**The offset IS the boost minus the MOND interpolation function correction.** This 0.014 dex RMS is nearly exact — the tiny residual comes from using the mean g_bar rather than the point-by-point transformation.

However, the γ-based route (predict boost from γ → convert to offset) gives R² = -0.32 because γ's boost predictions are too noisy for the nonlinear MOND transformation to preserve accuracy.

## Physical Interpretation

### Why the Boost Is Harder to Predict

The boost (g_obs/g_bar) varies from 0 to 1.1 dex — it depends on how deep in MOND the measurement point is (which varies across galaxies). The offset (g_obs/g_rar) is a DEVIATION from the mean RAR, centered near zero. The MOND interpolation function absorbs much of the boost variation, making the offset more predictable.

This is why the 6-var model gets R² = 0.94 for offset but only 0.75 for boost — the interpolation function provides an additional ~20% of explained variance for free.

### The c_V Amplification

c_V is 5.7× more important for boost than offset because:
- Boost depends on the actual g_obs/g_bar ratio, which is set by mass concentration
- Offset is a deviation from the RAR prediction, which already accounts for much of the concentration effect
- The MOND interpolation function partially cancels the c_V effect when converting boost → offset

### γ's Role: Real but Indirect

γ has genuine predictive power for the MOND boost (partial r = +0.79 after controlling all confounders), but this power cannot be extracted in practice because:
1. The confounders (V, L, f_gas, c_V) explain most of the variance first
2. γ's independent signal is too small (R² = 0.03 raw) to compete
3. The nonlinear boost → offset transformation amplifies errors

## Grade: B+

A thorough investigation establishing the exact boost-offset connection (r = 0.998) and revealing γ's dramatic partial correlations (+0.79) versus raw correlations (+0.16). The coefficient comparison (c_V 5.7× more important for boost) provides physical insight. The finding that γ fails in practice despite strong partial correlations is an important lesson about the limits of partial correlation interpretation. Minor deductions for the negative LOO R² of γ alone and the failed γ-based offset route.

## Files Created

- `simulations/session505_mond_boost_model.py`: 8 tests
- `Research/Session505_MOND_Boost_Model.md`: This document

---

*Session #505 verified: 8/8 tests passed*
*Grand Total: 1325/1325 verified*

**Key finding: offset = boost - log(ν(g/a₀)) with r=0.998, coefficient=0.986. γ has partial r=+0.79 with boost controlling all confounders, but raw R²=0.03. V,L dominate boost (R²=0.64, LOO 0.62). 6-var boost model R²=0.75, LOO=0.71. c_V is 5.7× more important for boost than offset. γ-based offset prediction fails (R²=-0.32). The boost is harder to predict than offset because the MOND interpolation function provides ~20% of explained variance for free. Grade B+.**
