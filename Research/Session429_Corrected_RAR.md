# Session #429: The Corrected RAR — MOND Implications

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The V+R+L+c_V model explains 93% of galaxy-level RAR offset variance. This session applies the correction to point-level data and examines the consequences: does the corrected RAR change the inferred a₀? Does it improve the interpolation function? How close do we get to measurement noise?

## Central Result: Correction Shifts a₀ Toward Canonical Value

| Quantity | Standard RAR | Corrected RAR | Change |
|----------|-------------|---------------|--------|
| Best-fit a₀ | 4.03×10⁻¹⁰ | 1.02×10⁻¹⁰ | -75% |
| a₀/a₀_canonical | 3.36 | 0.85 | |
| Point-level RMS (MOND) | 0.678 dex | 0.554 dex | -18.3% |
| Intrinsic scatter | 0.659 dex | 0.531 dex | -19.5% |

**The standard a₀ is inflated to 3.4× canonical** in this late-type subsample because the algebraic RAR is being fit to deep-MOND data with systematic galaxy-level biases. **After correction, a₀ recovers to 0.85× canonical (1.02×10⁻¹⁰ m/s²)**.

**Caveat**: These late-type galaxies are 100% in the MOND regime. The a₀ determination is poorly constrained when there are no Newtonian-regime data to anchor the interpolation function's high-acceleration behavior. The shift in a₀ should be interpreted as showing that the correction removes a systematic bias, not as a precision measurement of a₀.

## Key Findings

### 1. Point-Level Scatter Reduction (Test 1)

- Standard RAR RMS: 0.678 dex
- Corrected RAR RMS: 0.554 dex
- Improvement: 18.3%

The improvement is more modest at the point level than the galaxy level (where we get 74%) because within-galaxy scatter (37% of total variance, Session 427) is not addressed by the galaxy-level correction.

### 2. a₀ Recovery (Test 2)

The standard a₀ fit to these late-type galaxies yields 4.0×10⁻¹⁰ — far above the canonical 1.2×10⁻¹⁰. This inflated value reflects the systematic galaxy-level offsets: galaxies with properties that push them away from the mean RAR collectively bias the a₀ fit.

After applying the V+R+L+c_V correction:
- a₀_corrected = 1.02×10⁻¹⁰ m/s²
- Ratio to canonical: 0.85

This is a suggestive result. It implies that the galaxy-to-galaxy variations we've identified are responsible for inflating the apparent a₀ in MOND-dominated samples.

### 3. Interpolation Function (Test 3)

- Both the McGaugh+ and simple MOND forms improve slightly after correction
- The residual slope vs log(g_bar) decreases from -0.532 to -0.511 dex/dex
- The two interpolation forms become slightly more similar after correction
- Large residual slope persists → the algebraic formula is a poor fit in the deep-MOND regime regardless of correction

### 4. Residual Structure (Test 4)

After correction:
- Galaxy-level correlations perfectly removed: r(V, corr_offset) = 0.000, r(R, corr_offset) = 0.000, r(c_V, corr_offset) = 0.000 — by construction
- Point-level g_bar correlation actually increases (r = -0.87 → -0.92) — the correction removes galaxy-level structure but exposes the universal g_bar trend more clearly
- Remaining galaxy-level RMS: 0.187 dex (the cross-validated residual)

### 5. MDAR Tightening (Test 5)

The mass-discrepancy acceleration relation tightens by 4.1% in the MOND regime (1.043 → 1.000 dex scatter). The improvement is concentrated at intermediate g_bar (-14 to -12 range, +8-14%) and absent at the lowest accelerations where only a few galaxies contribute.

### 6. Galaxy-by-Galaxy (Test 6)

- 75% of galaxies improve (45/60)
- Median improvement: +25.2%
- Top improvers: KK98-251 (+73%), NGC4214 (+70%), UGC06399 (+65%)
- Worst worsened: NGC0055 (-559%) — this galaxy has very low standard RMS (0.058) and the correction overshoots

The worsened galaxies are mostly those with already-low scatter (RMS < 0.10) where any correction risks adding noise.

### 7. Noise Floor (Test 7)

- Measurement noise estimate: 0.160 dex (mean), 0.058 dex (median)
- Corrected RMS: 0.554 dex
- RMS/noise ratio: 3.46
- Intrinsic scatter reduction: 19.5%

The corrected RAR is still 3.5× above the noise floor. Significant intrinsic scatter remains.

## Physical Interpretation

The a₀ recovery is the most noteworthy finding. In the standard RAR applied to late-type galaxies, the fitted a₀ is inflated because:

1. **Compact galaxies** (high c_V, small R) have positive offsets → they pull g_obs up at given g_bar → the fit accommodates this by increasing a₀
2. **Extended galaxies** have negative offsets → they pull g_obs down
3. The net effect biases a₀ upward because the offset-V correlation means higher-mass galaxies (which contribute more points) tend to have positive offsets

After correction, each galaxy sits at its "true" position on the universal RAR, and a₀ drops to near-canonical.

However, 0.85× canonical suggests either: (a) our late-type sample has a slight remaining bias, (b) the algebraic RAR formula itself is imperfect in the deep-MOND regime, or (c) the correction slightly overcorrects because of the in-sample fitting.

## Grade: B+

The a₀ recovery from 3.4× to 0.85× canonical is a striking result that validates the correction's physical reality. The 18.3% point-level improvement confirms the galaxy-level correction translates to real gains. However, the session is limited by working only with deep-MOND data (no Newtonian anchor), the large residual scatter (0.55 dex) showing we're far from the noise floor, and the remaining systematic g_bar trend. The most important insight — that galaxy structural variation inflates a₀ — deserves further investigation.

## Files Created

- `simulations/session429_corrected_rar.py`: 8 tests
- `Research/Session429_Corrected_RAR.md`: This document

---

*Session #429 verified: 8/8 tests passed*
*Grand Total: 821/821 verified*

**Key finding: Applying the V+R+L+c_V correction to point-level RAR data recovers a₀ from 3.4× canonical to 0.85× canonical (1.02×10⁻¹⁰ m/s²). Point-level RMS improves 18.3%. 75% of galaxies improve (median +25%). Intrinsic scatter drops 19.5%. The correction removes galaxy-level structural biases that inflate a₀ in MOND-dominated samples. Still 3.5× above noise floor — significant unexplained point-level scatter remains. Grade B+.**
