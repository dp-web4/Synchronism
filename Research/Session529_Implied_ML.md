# Session #529: The Implied M/L — Extracting Mass-to-Light from the Model

**Date**: 2026-02-07
**Status**: 8/8 verified

## Overview

Session #526 found logL×f_gas implies M/L ∝ L^0.36. Session #528 showed the BTFR IS MOND when f_gas is included. The 6-var model predicts the per-galaxy RAR offset, which measures the M/L correction. This session extracts the implied M/L for each galaxy and tests whether it matches stellar population synthesis expectations.

## Central Result: Implied M/L Is Physically Sensible But Nearly Flat with Luminosity

The median implied M/L is 0.44 (3.6μm), close to the assumed 0.5. The M/L correlates with morphological type (r=-0.21, later types have lower M/L) and surface brightness (r=+0.21, HSB have higher M/L), both in the expected direction from stellar populations. However, the M/L-luminosity slope is only 0.027 — far flatter than the 0.36 predicted from the logL×f_gas coefficient. The "0.36" was predominantly gas-luminosity covariance, not a true M/L-luminosity relation.

## Key Findings

### 1. Implied M/L from the Offset (Test 1)

In the deep MOND limit: offset ≈ 0.5 × log(M/L_true / M/L_assumed). Therefore:

M/L_implied = 0.5 × 10^(2 × offset)

| Statistic | Simple | Gas-corrected |
|-----------|--------|---------------|
| Mean | 0.523 | 0.631 |
| Median | 0.439 | 0.468 |
| Std | 0.336 | 0.529 |
| Range | [0.014, 2.03] | [0.027, 2.80] |

The median M/L of 0.44 is close to the assumed 0.5 and within the SPS expectation of 0.5 ± 0.1 for 3.6μm disk galaxies. The gas-corrected values are slightly higher (median 0.47) because the correction accounts for gas contributing to g_bar without an M/L uncertainty.

### 2. M/L vs Luminosity: Flatter Than Predicted (Test 2)

| Metric | Value |
|--------|-------|
| r(logL, log M/L) | +0.089 (p=0.32) — NOT significant |
| Slope log(M/L) vs logL | 0.027 ± 0.027 |
| Session #526 prediction | 0.36 |

**The M/L-luminosity slope is 0.027, not 0.36.** The discrepancy of 0.33 is explained by the fact that the logL×f_gas coefficient (β=0.181) captures the joint effect of M/L AND gas fraction covarying with luminosity. When we extract the pure M/L (removing the gas effect from the offset), the luminosity dependence nearly vanishes.

This means:
- The logL×f_gas interaction is mostly about **gas-luminosity covariance** (gas-rich galaxies are low-L)
- True M/L varies little with luminosity (consistent with a "universal" M/L at 3.6μm)
- The ~0.36 slope is ~0.03 from M/L variation and ~0.33 from f_gas-L covariance

### 3. M/L vs Galaxy Type (Test 3)

| Type | N | Median M/L | Mean M/L |
|------|---|-----------|----------|
| Sa-Sb (T=1-3) | 21 | 0.415 | 0.537 |
| Sbc-Sc (T=4-5) | 30 | 0.433 | 0.488 |
| Scd-Sd (T=6-7) | 29 | 0.491 | 0.573 |
| Sm-Im/Irr (T≥8) | 47 | 0.383 | 0.496 |

r(T, log M/L) = -0.214 (p=0.015): **later types have lower M/L**, consistent with younger stellar populations. The trend is weak but in the expected direction.

Interestingly, Scd-Sd galaxies have the highest median M/L (0.49), not early types (0.42). This may reflect distance or inclination systematics — early-type spirals in SPARC tend to be at larger distances with more uncertain parameters.

### 4. M/L vs Surface Brightness (Test 4)

| Correlation | r | p |
|------------|---|---|
| r(log SB_eff, log M/L) | +0.211 | 0.017 |
| r(log SB_disk, log M/L) | +0.244 | 0.006 |
| r_partial(log SB, log M/L \| logL) | +0.237 | 0.007 |

**Higher surface brightness galaxies have higher implied M/L**, even after controlling for luminosity. This is consistent with HSB galaxies having older, redder stellar populations with higher M/L. This is the strongest M/L correlation after controlling for luminosity.

### 5. M/L vs Gas Fraction (Test 5)

| Correlation | r | p |
|------------|---|---|
| r(f_gas, log M/L) | -0.096 | 0.28 |
| r(f_gas, log M/L_corrected) | +0.049 | 0.61 |
| r_partial(f_gas, log M/L \| logL) | -0.043 | 0.63 |

**M/L is NOT significantly correlated with gas fraction** after proper extraction. This is surprising — gas-rich galaxies should have younger stars and lower M/L. The lack of correlation means either: (a) the extraction already absorbs the f_gas-M/L link through the offset, or (b) the true M/L-f_gas correlation is weak at 3.6μm where stellar populations have relatively uniform M/L.

### 6. SPS Comparison (Test 6)

- 90.6% of galaxies have M/L in [0.2, 1.5] — physically reasonable
- 11 galaxies have M/L < 0.2 (including known problematic cases: F561-1, PGC51017, UGC04305, UGC06628)
- 1 galaxy has M/L > 1.5 (UGC07399)
- Gas-poor median M/L = 0.45, gas-rich median M/L = 0.41 — small difference

The low-M/L outliers are mostly dwarf irregulars or known problematic galaxies, consistent with distance/inclination errors rather than truly anomalous M/L values.

### 7. The f_gas-M/L Degeneracy (Test 7)

Variance decomposition of the 6-var model prediction:

| Component | Variance | % (raw) |
|-----------|----------|---------|
| BTFR (V, L) | 0.0496 | 198% |
| Gas (f_gas, L×f_gas) | 0.0214 | 86% |
| Structure (c_V, V×c_V) | 0.0019 | 8% |
| BTFR-Gas covariance | -0.0241 | — |
| Total predicted | 0.0250 | 100% |

The BTFR and gas components have enormous variance individually (198% and 86%) but large negative covariance (-0.024), because luminosity correlates with both the BTFR and gas fraction. The net effect: BTFR predicts mass, gas corrects the luminosity-mass mapping.

The BTFR residual has σ=0.105 dex (65% of total offset σ). This residual correlates strongly with f_gas (r=-0.73) and with implied M/L (r=+0.65), confirming that what the model corrects is a mix of gas effects and M/L variation.

### 8. Synthesis (Test 8)

The 6-var model predicts M/L to 0.038 dex (factor 1.09 accuracy). Key findings:

1. **Median M/L = 0.44** — close to the assumed 0.5, consistent with SPS at 3.6μm
2. **M/L-luminosity slope = 0.027** — much flatter than the 0.36 from logL×f_gas coefficient
3. **The 0.36 is gas-luminosity covariance**, not true M/L variation
4. **r(type, M/L) = -0.21** — later types have lower M/L (correct direction)
5. **r(SB, M/L) = +0.21** — HSB have higher M/L (correct direction)
6. **r(f_gas, M/L) = -0.10** — weak, not significant
7. **90.6% of implied M/L values are physically reasonable** (0.2-1.5)

The offset IS the M/L correction, and the 6-var model IS a M/L predictor. But the "M/L ∝ L^0.36" from Session #526 was misleading — it's predominantly gas-luminosity covariance, not a stellar population gradient. The true M/L at 3.6μm is nearly universal (~0.44), with weak dependence on type and surface brightness.

## Grade: A-

A valuable session that tests and partially refutes a prediction from Session #526. The finding that the M/L-luminosity slope is 0.027 (not 0.36) corrects an overinterpretation of the logL×f_gas coefficient. The SPS comparison and type/SB correlations are in the right direction, providing physical validation. The variance decomposition clarifies the BTFR-gas degeneracy. Minor deductions: the extraction assumes deep MOND (not always valid), and the gas-corrected M/L formula is approximate. Could have tested whether using the implied M/L as input improves the model iteratively.

## Files Created

- `simulations/session529_implied_ml.py`: 8 tests
- `Research/Session529_Implied_ML.md`: This document

---

*Session #529 verified: 8/8 tests passed*
*Grand Total: 1469/1469 verified*

**Key finding: Median implied M/L = 0.44 (close to assumed 0.5, consistent with SPS). M/L-luminosity slope = 0.027 (NOT 0.36 from Session #526 — the difference is gas-luminosity covariance). r(type, M/L) = -0.21 (later→lower M/L, correct). r(SB, M/L) = +0.21 (HSB→higher M/L, correct). 90.6% of M/L physically reasonable. Model predicts M/L to 0.038 dex (factor 1.09). True M/L at 3.6μm is nearly universal. Grade A-.**
