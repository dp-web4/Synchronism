# Session #418: Modified Interpolating Function with R_eff Dependence

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Building on Session 417's finding that the R_eff effect amplifies outward, this session tests whether modifying the RAR's characteristic acceleration g† to depend on galaxy structure can capture the effect.

## Central Result: Additive Correction Beats Modified g†

| Approach | Point-level RMS improvement | Galaxy-level RMS improvement |
|----------|---------------------------|----------------------------|
| Modified g†_eff | **14.4%** | 25.8% |
| Additive correction (V+R) | **28.6%** | 50.3% (LOO) |

The additive correction (constant offset per galaxy in log space) works 2× better than modifying the interpolating function's transition acceleration.

## The Modified g†

**g†_eff = g† × (V_flat/75 km/s)^{+1.10} × (R_eff/2.27 kpc)^{-0.34}**

- c1 (V coefficient) = +1.10 (compare empirical: +1.21)
- c2 (R coefficient) = -0.34 (compare empirical: -0.36)
- These match the galaxy-level model closely

g†_eff range: [0.26, 1.84] × g† — a factor of 7 across galaxies

Interpretation:
- Compact fast galaxies: g†_eff up to 1.84 × g† → stronger MOND effect
- Extended slow galaxies: g†_eff down to 0.26 × g† → weaker MOND effect

## Detailed Findings

### 1. R_eff Absorption (Test 2)

The modified g† reduces but doesn't eliminate the R_eff correlation:
- Standard: r(R_eff, offset | V) = **-0.74**
- Modified: r(R_eff, offset | V) = **-0.55**

Only ~26% of the R_eff predictive power is absorbed.

### 2. Residual Structure (Test 5)

After modified g†, the remaining offsets still correlate with:
- V_flat: r = +0.58 (strong — the V coefficient in the interpolation doesn't fully capture the V dependence)
- SB: r = +0.29 (moderate)
- R_eff: r = -0.07 (mostly absorbed for this specific quantity)

### 3. Outward Amplification (Test 6)

Compact-extended difference by radius:

| r/R_eff | Standard | Modified g† |
|---------|----------|------------|
| [0, 2) | +0.076 | +0.017 |
| [2, 5) | +0.185 | +0.126 |
| [5, ∞) | +0.199 | +0.165 |

The inner part is well corrected (0.076 → 0.017) but the outer amplification persists (0.199 → 0.165). The g†_eff modification is too simple to capture the radius-dependent structure.

### 4. Cross-Validation (Test 4)

LOO-CV comparison:
- Standard RAR: 0.203 dex
- Additive V+R model: **0.101 dex** (50.3% improvement)

This confirms the galaxy-level additive correction is robust and not overfitted.

## Why Additive > Multiplicative

Modifying g† changes the entire shape of the interpolating function. But the data shows:
1. The offset is mostly **uniform** within galaxies (Session 417)
2. A constant shift per galaxy (additive in log space) captures this better than reshaping the interpolation
3. The modified g† over-corrects the inner regions and under-corrects the outer regions

This suggests the true correction is NOT simply a modified acceleration scale, but rather a galaxy-dependent offset that may reflect something beyond the standard interpolating function framework.

## Grade: B+

The modified g† approach is a natural idea but is outperformed by the simpler additive correction. The g†_eff coefficients closely matching the empirical model is reassuring but the 14% vs 29% improvement gap shows the limitation of this approach. Good negative result.

## Files Created

- `simulations/session418_modified_interpolation.py`: 8 tests
- `Research/Session418_Modified_Interpolation.md`: This document

---

*Session #418 verified: 8/8 tests passed*
*Grand Total: 741/741 verified*

**Key finding: Modified g†_eff = g† × (V/75)^{1.1} × (R/2.27)^{-0.34} improves RAR by 14.4% at point level but is outperformed 2× by simple additive correction (28.6%). The g†_eff coefficients match the galaxy-level model. R_eff correlation reduces from -0.74 to -0.55 (only 26% absorbed). Inner regions well corrected but outer amplification persists. Additive > multiplicative. Grade B+.**
