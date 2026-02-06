# Session #480: N_corr Revisited — The γ = 2/√N_corr Prediction

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The Synchronism framework predicts γ = 2/√N_corr where N_corr = V²/(R×a₀). Previous tests (Sessions 385-389) found R² = 0.23. This session revisits N_corr with the outer-only offset and decomposition into physical components.

## Central Result: N_corr Predicts Offset at r = 0.55, But the Sign Is Wrong for γ = 2/√N

N_corr correlates positively with the offset (r = +0.55), meaning higher N_corr → higher offset. The γ = 2/√N prediction gives a *negative* slope (higher N → smaller γ → lower offset). This sign discrepancy means the Synchronism γ formula does not directly predict the observed offset direction. However, N_corr does predict the *magnitude* of the offset (r(logN, |offset|) = -0.48), and the late-type correlation reaches r = +0.85 for the outer offset.

## Key Findings

### 1. N_corr vs Offsets (Test 1)

| N_corr version | r(N, full) | r(N, outer) | R²(full) | R²(outer) |
|----------------|-----------|-------------|---------|----------|
| R_eff-based | +0.546 | **+0.604** | 0.298 | **0.365** |
| R_max-based | +0.531 | +0.567 | 0.282 | 0.321 |

N_corr predicts the outer offset better than the full offset (R² = 0.37 vs 0.30), consistent with Session 477's finding that the outer offset is cleaner. The R_eff-based version is slightly better than R_max-based.

### 2. The γ Prediction — Sign Problem (Test 2)

| Metric | Predicted | Observed |
|--------|-----------|----------|
| Slope (offset vs logN) | -0.217 | **+0.235** |
| r(γ prediction, offset) | expected positive | **-0.571** |

**The observed slope is positive, not negative.** This means:
- High N_corr galaxies (massive, compact, high acceleration) have *positive* offsets (g_obs > g_RAR)
- Low N_corr galaxies (dwarf, extended, low acceleration) have *negative* offsets (g_obs < g_RAR)

The γ = 2/√N prediction says the opposite: low N should have larger deviations *above* the RAR. The sign flip means the offset is driven by M/L effects (massive galaxies have higher true M/L → positive offset) rather than by statistical fluctuations as the Synchronism framework proposes.

### 3. N_corr in the 5-Variable Model (Test 3)

| Model | k | R² | LOO |
|-------|---|-----|-----|
| logN alone | 2 | 0.298 | 0.133 |
| logN + f_gas | 3 | 0.309 | 0.133 |
| 5-var | 6 | 0.872 | 0.059 |
| 5-var + logN | 7 | 0.872 | 0.060 |

**N_corr adds nothing to the 5-variable model** (ΔR² = 0.0001). The 5 variables (V, L, c_V, f_gas, V×c_V) completely absorb the N_corr information. This is because N_corr = V²/(R×a₀), and R is captured by the combination of L and c_V.

### 4. N_corr Decomposition (Test 4)

| Component | r(X, full) | r(X, outer) | Partial r(\|logV) |
|-----------|-----------|-------------|-------------------|
| logV | +0.422 | +0.389 | — |
| logR_eff | +0.030 | -0.085 | **-0.382** |
| logN_corr | +0.546 | +0.604 | — |

**R_eff adds significant information beyond V** (partial r = -0.38 for full, -0.51 for outer). At fixed V, smaller galaxies have lower offset — consistent with the size effect (Sessions 390-393). The raw r(logR, offset) ≈ 0 because R and V are correlated; the partial correlation reveals the true size effect.

### 5. N_corr by Galaxy Type (Test 5)

| Type | N | r(N, offset) | R² |
|------|---|-------------|-----|
| S0-Sb | 22 | +0.504 | 0.254 |
| Sbc-Sd | 46 | -0.062 | 0.004 |
| **Sdm-Im** | **60** | **+0.793** | **0.629** |

**Late types show a dramatic correlation: r = +0.79 (R² = 0.63)**, while Sbc-Sd types show essentially zero correlation. The late-type outer offset correlation reaches **r = +0.85**. This is the strongest single-predictor result in the entire research program for late types.

Gas-dominated late types (f_gas > 0.5): r = +0.70 — the signal persists even in the M/L-independent regime.

### 6. N_corr vs 5-Var Residual (Test 6)

r(logN, 5-var residual) = +0.013 for full model, +0.064 for outer model. **N_corr is completely captured by the 5-variable model.** Adding it to the outer model gives ΔR² = +0.001.

### 7. Physical Interpretation (Test 7)

N_corr = g_cent(R_eff)/a₀ = the centripetal acceleration at the effective radius in units of a₀. Physically:
- N_corr > 1 (57% of sample): Newtonian regime at R_eff
- N_corr < 1 (43% of sample): Deep MOND at R_eff

**r(logN, |offset|) = -0.48**: High N_corr → small |offset|, meaning galaxies deep in MOND (low N_corr) have larger absolute offsets. This is consistent with the γ = 2/√N prediction for *magnitude* but not *direction*.

## Physical Interpretation

### The Sign Problem Explained

The γ = 2/√N prediction assumes the offset is a random fluctuation around zero, with amplitude ∝ 1/√N. But the observed offset has systematic direction: it correlates with M/L-related properties (V, L, f_gas). The positive slope means:

1. High-V, compact galaxies (high N_corr) → high M/L → offset > 0
2. Low-V, extended galaxies (low N_corr) → low M/L → offset < 0

This is an M/L effect, not a statistical fluctuation effect. The Synchronism γ prediction works for offset *magnitude* (r = -0.48 for |offset|) but not offset *sign*.

### Late Types: The Best Test

In late types (T ≥ 7), which are gas-dominated and have negligible M/L uncertainty, N_corr predicts the offset at r = 0.79. This is the cleanest test because M/L is irrelevant. The signal here may genuinely reflect the N_corr physics rather than M/L contamination. The R² = 0.63 for late types from a single parameter (N_corr) is remarkable.

### Why N_corr Is Redundant With 5-Var

N_corr = V²/(R × a₀) ∝ V² / R. In the 5-variable model:
- V is directly included (logV)
- R ∝ √(L/SB) is captured by logL and c_V (which relates V(R_eff) to V_flat)
- So logN ≈ 2logV - log(√(L/SB)) is a linear combination of the model variables

## Grade: B+

A valuable revisit that sharpens the N_corr findings. The sign discrepancy with γ = 2/√N is an important negative result that constrains the Synchronism framework. The late-type r = 0.85 (outer offset) is a strong positive finding. The decomposition showing R_eff adds beyond V (partial r = -0.38) confirms the size effect. Slightly lower grade because the sign problem is a significant challenge for the theoretical framework.

## Files Created

- `simulations/session480_ncorr_revisited.py`: 8 tests
- `Research/Session480_Ncorr_Revisited.md`: This document

---

*Session #480 verified: 8/8 tests passed*
*Grand Total: 1157/1157 verified*

**Key finding: N_corr = V²/(R×a₀) correlates with RAR offset at r = +0.55, but the sign is OPPOSITE to γ = 2/√N prediction. Late types reach r = +0.85 for outer offset — the strongest single-predictor result. R_eff adds beyond V at partial r = -0.38 (confirming size effect). N_corr adds nothing to 5-var model (ΔR² = 0.000). The prediction works for |offset| magnitude (r = -0.48) but not direction. Grade B+.**
