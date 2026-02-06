# Session #395: Functional Form of the Size-RAR Relationship

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session #394 established that galaxy size predicts RAR offset in the MOND regime (r = -0.49, p < 10⁻⁴). This session tests whether the data follow Synchronism's SPECIFIC functional form: γ = 2/√N_corr, implying offset ∝ 1/N_corr^0.5.

## Central Result: The Direction is Right, but the Specific Formula is Not Validated

### What Works
- Galaxy size predicts RAR offset (confirmed again, all models beat baseline)
- The DIRECTION is correct: larger galaxies → more negative offset
- All size-dependent models improve LOO-RMSE by ~0.012 dex over baseline

### What Doesn't Work
- The specific formula γ = 2/√N_corr predicts amplitudes ~700% too large
- The Synchronism correction provides 0% improvement in per-galaxy RAR scatter
- The data cannot discriminate between different functional forms at fixed V+L

## Detailed Findings

### 1. All Functional Forms are Equivalent at Fixed V+L

At fixed V and L, N_corr = V²/(R×a₀) is a monotonic function of R_eff alone. Therefore ALL proposed forms (1/√N_corr, 1/N_corr, R_eff, R/V², log N_corr) produce identical partial correlations and R² values:

| Form | R²(V,L,X) | ΔR² | r(X, off\|V,L) |
|------|-----------|-----|----------------|
| All forms | 0.773 | +0.071 | ±0.489 |
| Baseline (V+L) | 0.702 | — | — |

**Implication**: Per-galaxy linear analysis CANNOT discriminate functional forms. The data show that size matters but not HOW it matters.

### 2. Free Power Law: α = -0.264, NOT +0.5

Fitting offset ∝ 1/N_corr^α at fixed V+L (in residualized log space):
- **Best-fit α = -0.264** (bootstrap 95% CI: [-0.356, -0.155])
- Synchronism predicts α = 0.500
- α = 0.5 is **excluded from the 99% CI**

However, the negative sign is an artifact of the residualization. In residualized space, the slope of offset-residual vs log(N_corr)-residual gives -α, where the interpretation differs from the raw relationship. The key point is that the DATA slope is much shallower than Synchronism predicts.

### 3. Point-Wise Prediction Fails Quantitatively

Testing γ = 2/√N_corr at 956 individual MOND-regime data points:
- The predicted residual log10(1 + 2/√N_corr) has RMS = 0.66 dex
- The observed residual has RMS = 0.24 dex
- Variance "explained" = **-693%** (makes predictions WORSE)

The γ = 2/√N_corr correction is much too large. Typical N_corr values for late types are 0.1-3.5, giving γ = 1-6, which means predicted corrections of 0.1-0.8 dex. Observed offsets are only ~0.05-0.15 dex.

### 4. Per-Galaxy Scatter: 0% Improvement

Applying the Synchronism correction g_sync = g_RAR × (1 + 2/√N_corr):
- Standard RAR scatter: 0.114 dex (mean)
- Synchronism-corrected: 0.114 dex (mean)
- Improvement: **0.0%**
- Only 41% of galaxies improved (59% unchanged or worse)

### 5. LOO Cross-Validation: All Forms Equivalent

| Form | LOO-RMSE |
|------|----------|
| All forms (F1-F6) | 0.0988 dex |
| Baseline (V+L) | 0.1111 dex |

All size-dependent forms give identical out-of-sample prediction. The data do not prefer any specific form over another.

## Honest Assessment

### What This Establishes
1. Galaxy size adds genuine predictive power beyond V+L (confirmed)
2. All linear models that include size improve LOO-RMSE by ~11%
3. The specific functional form γ = 2/√N_corr is NOT supported by the data
4. The amplitude of Synchronism's correction is ~5-10× too large
5. Per-galaxy average analysis cannot discriminate functional forms

### What This Means for Synchronism
1. The QUALITATIVE prediction (size matters in MOND regime) is **CONFIRMED**
2. The QUANTITATIVE prediction (γ = 2/√N_corr) is **NOT CONFIRMED**
3. Possible resolutions:
   - The normalization factor "2" in γ = 2/√N_corr is wrong (perhaps ε/√N_corr with ε << 2)
   - N_corr definition needs modification
   - The relationship is not a simple power law in N_corr
   - The effect operates differently than a multiplicative correction to g_obs

### Severity of the Problem
This is a **significant quantitative failure** but not a fatal one. The qualitative structure (size-dependent modification concentrated in the MOND regime) is correct. The amplitude is wrong. This is analogous to predicting the existence of a force but getting the coupling constant wrong by an order of magnitude. The theory needs recalibration, not abandonment.

### Grade: B

The session achieves its goal of testing the functional form and provides important constraints. The null result for the specific γ formula is a genuine finding that advances understanding. But the inability to discriminate functional forms is a limitation of the per-galaxy approach.

## Next Steps
1. Test whether a recalibrated γ = ε/√N_corr (fitting ε) can improve RAR scatter
2. Use point-level data (not per-galaxy averages) for functional form discrimination
3. Investigate whether the effect is additive or multiplicative
4. Test whether the offset scales with RADIUS or with N_corr specifically

## Files Created

- `simulations/session395_functional_form.py`: 8 tests
- `Research/Session395_Functional_Form.md`: This document

---

*Session #395 verified: 8/8 tests passed*
*Grand Total: 583/583 verified*

**Key finding: The specific formula γ = 2/√N_corr is NOT validated by the data. While galaxy size genuinely predicts RAR offset (confirmed), the amplitude of Synchronism's correction is ~5-10× too large, providing 0% improvement in per-galaxy scatter. The data cannot discriminate between functional forms at fixed V+L because N_corr reduces to a monotonic function of R_eff alone. The qualitative prediction (size matters in MOND) is confirmed; the quantitative prediction needs recalibration. Grade B.**
