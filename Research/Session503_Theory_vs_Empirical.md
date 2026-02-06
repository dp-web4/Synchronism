# Session #503: Synchronism Theory vs Empirical 6-Variable Model

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The Synchronism framework predicts γ = 2/√N_corr where N_corr = V²/(R×a₀). The empirical 6-variable model achieves R² = 0.945. This session tests whether the theoretical γ prediction matches the empirical offset, whether γ can replace the 6-var model, and what the discrepancy reveals about the theory.

## Central Result: γ Explains Only 28% — The Theory Needs Refinement

The theoretical γ = 2/√N_corr correlates with the offset at r = -0.57 (note: NEGATIVE, opposite to the naive prediction). γ alone achieves R²(LOO) = 0.28, compared to the 6-var model's 0.94. Adding γ to the 6-var model improves LOO by only 0.001. The theoretical prediction does not match the empirical model quantitatively, though it captures the correct qualitative dependence on V and R.

## Key Findings

### 1. γ_theory Statistics (Test 1)

| Metric | Value |
|--------|-------|
| N_corr range | [0.03, 1.01] |
| N_corr median | 0.2 |
| γ(R_max) mean | 4.4 |
| γ(R_max) range | [2.0, 12.5] |
| log₁₀(γ) mean | +0.625 |
| Offset mean | -0.038 |
| Offset range | [-0.77, +0.30] |

γ ranges from 2 to 12, all > 1. log₁₀(γ) ranges from +0.30 to +1.10. The offset ranges from -0.77 to +0.30. **The scales don't match**: log(γ) is always positive while the offset is mostly negative.

### 2. γ vs Offset (Test 2)

| Metric | Value |
|--------|-------|
| r(log γ, offset) | **-0.57** |
| r(log γ_reff, offset) | -0.60 |
| Regression slope | **-0.66** |
| Theory prediction | +1.0 |
| R² (γ alone) | 0.32 |

**The correlation is NEGATIVE**: larger γ → more negative offset. This is the opposite of the naive expectation that offset = log(γ). The regression gives offset ≈ +0.38 - 0.66 × log(γ), with a negative slope instead of the predicted +1.0.

### 3. 6-Var Model Adds Massively Beyond γ (Test 3)

| Model | R² | LOO R² | RMS |
|-------|-----|--------|-----|
| γ alone | 0.32 | 0.28 | 0.134 |
| 6-var (no γ) | **0.945** | **0.938** | 0.038 |
| γ + 6-var | 0.947 | 0.939 | 0.037 |

Adding the 6-var predictors to γ improves LOO by +0.66. Adding γ to the 6-var model improves LOO by only +0.001. **γ is almost entirely redundant with the 6-var model.**

The γ-model residual correlates with c_V (r = -0.31), f_gas (r = +0.33), and logL (r = -0.26) — these are exactly the variables the 6-var model uses to correct γ's deficiencies.

### 4. γ Cannot Replace the 6-Var Model (Test 4)

Even with corrections: γ + f_gas + logL×f_gas + c_V achieves LOO R² = 0.43 — less than half the 6-var model's 0.94. The γ-based model captures only the broad trend (large dwarf galaxies have more negative offsets), not the precise galaxy-by-galaxy variation.

log(γ) correlates with all 6-var predictors:
- r(log γ, logV) = -0.66
- r(log γ, logL) = -0.54
- r(log γ, f_gas) = +0.65
- r(log γ, c_V) = -0.30

### 5. N_corr Decomposition (Test 5)

log(N_corr) = -3.82 + 1.95×logV - 0.44×logL - 0.11×c_V - 1.10×f_gas (R² = 0.78)

The logV coefficient (1.95) is close to the theoretical 2.0 — confirming N_corr ∝ V². The logL coefficient (-0.44) reflects the correlation of R with L (larger galaxies → larger R → lower N_corr). The f_gas coefficient is large (-1.10) because gas-rich galaxies have larger R_max at fixed V.

### 6. Theoretical Coefficient Prediction (Test 6)

The theory predicts offset coefficients via log(γ) decomposition:

| Coefficient | Theory predicts | 6-var observed |
|------------|----------------|----------------|
| β(logV) | **-0.73** | **+1.90** |
| β(logL) | **+0.09** | **-0.55** |

**The signs are wrong.** Theory predicts negative β(logV) and positive β(logL), while the 6-var model has the opposite signs. This is because:

1. The offset is NOT equal to log(γ). It's the deviation of g_obs from the mean RAR.
2. At fixed g_bar (not fixed V), the MOND prediction already absorbs the V dependence.
3. The 6-var model's β(logV) = +1.90 reflects the BTFR (M ∝ V⁴), while γ's dependence on V is through the different quantity N_corr = V²/R.

### 7. γ by Galaxy Type (Test 7)

| Type | ⟨log γ⟩ | ⟨offset⟩ |
|------|---------|----------|
| Early (T<4) | +0.556 | +0.001 |
| Mid (4≤T<7) | +0.556 | -0.026 |
| Late (T≥7) | **+0.702** | **-0.061** |

Late types have the largest γ (smallest N_corr) and the most negative offsets — consistent with the negative correlation.

### 8. Optimal γ-Based Model (Test 8)

| Model | R² | LOO R² | # params |
|-------|-----|--------|----------|
| γ alone | 0.32 | 0.28 | 2 |
| γ + f_gas | 0.45 | 0.41 | 3 |
| γ + f_gas + logL×f_gas | 0.48 | 0.43 | 4 |
| γ + f_gas + logL×f_gas + c_V | 0.49 | 0.43 | 5 |
| **γ + all 6-var** | **0.95** | **0.94** | 8 |
| **6-var (no γ)** | **0.94** | **0.94** | 7 |

The γ-based model saturates at LOO R² ≈ 0.43 with corrections. Only when the full 6-var predictors are added does it reach 0.94 — but at that point γ contributes nothing incremental.

## Physical Interpretation

### Why γ = 2/√N_corr Doesn't Work as a Direct Predictor

1. **Different quantities**: The RAR offset measures log(g_obs/g_rar), which is a DEVIATION from the mean relation at each g_bar. γ measures an absolute scale (2/√N_corr). These live on different scales and have different dependencies on V and R.

2. **The V dependence is inverted**: γ ∝ 1/V (through N_corr ∝ V²), so large galaxies have small γ. But the offset depends on where a galaxy sits on the BTFR, where large V → positive offset (more mass discrepancy at given luminosity). These are opposite trends.

3. **R_max introduces noise**: N_corr uses R_max (the outermost observed radius), which depends on survey depth, not just galaxy physics. This adds scatter that the 6-var model avoids by using c_V instead.

4. **The theory needs a transformation**: γ might predict something related to the offset but through a non-trivial mapping, not a simple log(γ) = offset relationship.

### What γ Does Capture

Despite the poor quantitative match, γ captures 32% of offset variance — more than any single 6-var predictor except logV (which captures ~65% alone). γ encodes a genuine physical quantity (the number of correlated MOND domains), and its negative correlation with offset is consistent with: more MOND domains → smoother acceleration field → smaller deviations from the mean RAR.

### Implications for the Synchronism Framework

The γ = 2/√N_corr prediction, as formulated, does not directly predict the RAR offset. The framework would need to:

1. **Define what γ predicts**: It may predict something other than the RAR offset (e.g., the scatter at each point, or the deviation from the deep-MOND limit specifically).
2. **Account for the BTFR**: The offset is dominated by the galaxy's position on the BTFR (78%, Session #496), which γ doesn't encode.
3. **Use a different radius**: R_eff or a dynamically motivated radius might work better than R_max (r(log γ_reff, offset) = -0.60 vs r(log γ_rmax, offset) = -0.57).

## Grade: B

An important negative result that honestly confronts the gap between theory and data. The γ = 2/√N_corr prediction captures only 28% of offset variance and has the wrong sign structure. However, the analysis is thorough: the N_corr decomposition (Test 5) confirms the V² dependence, the coefficient comparison (Test 6) precisely identifies where the theory diverges, and the model-building exercise (Test 8) quantifies γ's contribution. The session clarifies what the Synchronism framework needs to address to make quantitative predictions.

## Files Created

- `simulations/session503_theory_vs_empirical.py`: 8 tests
- `Research/Session503_Theory_vs_Empirical.md`: This document

---

*Session #503 verified: 8/8 tests passed*
*Grand Total: 1309/1309 verified*

**Key finding: γ = 2/√N_corr correlates with offset at r = -0.57 (NEGATIVE) and explains only 28% of variance vs 6-var model's 94%. Slope is -0.66 (theory: +1.0). Coefficient signs are inverted (theory: β(logV) = -0.73, observed: +1.90). γ + corrections saturate at LOO R² = 0.43. Adding γ to 6-var model improves LOO by only 0.001. The Synchronism prediction needs refinement to connect N_corr to the RAR offset. Grade B.**
