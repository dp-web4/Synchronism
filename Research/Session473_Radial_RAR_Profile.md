# Session #473: The Radial RAR Profile — Offset as a Function of Radius

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The 5-variable model uses a single galaxy-level offset. But does the RAR offset vary systematically with radius within a galaxy? If so, this would reveal M/L gradients, interpolation function errors, or physics beyond the algebraic RAR.

## Central Result: The Offset Gradient Correlates With the 5-Var Residual at r = -0.43

The radial gradient of the RAR offset — measuring how the offset changes from inner to outer regions — correlates with the 5-variable model residual at r = -0.43. Adding the gradient as a 6th variable gives ΔR² = +0.036, the strongest improvement found since f_gas. However, this is a derived quantity requiring the full rotation curve, not a simple galaxy property.

## Key Findings

### 1. Radial Offset Profile (Test 1)

| r/R_eff | N | ⟨offset⟩ | σ | ⟨log g_bar⟩ |
|---------|---|---------|---|-------------|
| 0-0.5 | 361 | **-0.083** | 0.305 | -10.0 |
| 0.5-1 | 392 | -0.016 | 0.202 | -10.2 |
| 1-1.5 | 380 | -0.005 | 0.156 | -10.3 |
| 1.5-2 | 308 | -0.015 | 0.146 | -10.5 |
| 2-3 | 435 | -0.029 | 0.138 | -10.7 |
| 3-5 | 483 | -0.035 | 0.143 | -10.9 |
| >5 | 491 | **+0.018** | 0.118 | -11.1 |

**The innermost points (r < 0.5 R_eff) have a systematically negative offset (-0.083 dex)**, meaning the RAR over-predicts g_obs in galaxy centers. The outermost points (r > 5 R_eff) have a slightly positive offset (+0.018 dex). The scatter decreases monotonically from 0.305 dex (inner) to 0.118 dex (outer).

**Physical interpretation**: The inner negative offset reflects M/L being too low (M/L_disk = 0.5 underestimates the central stellar mass), or beam smearing making V_obs too low in the center. The scatter decrease with radius reflects the outer rotation curve being smoother and better measured.

### 2. Inner vs Outer Offset (Test 2)

| Region | ⟨offset⟩ | σ |
|--------|---------|---|
| Inner (r < R_eff) | -0.062 | 0.226 |
| Outer (r > 2 R_eff) | -0.034 | 0.145 |
| Δ(outer - inner) | +0.028 | 0.198 |

Outer > Inner in only 54% of galaxies — barely significant. The inner and outer offsets are moderately correlated (r = 0.50), meaning they share a common galaxy-level component but also have independent radial variation.

### 3. Gradient vs Properties (Test 3)

| Property | r(gradient, X) |
|----------|---------------|
| c_V | **-0.357** |
| offset | -0.088 |
| T | -0.072 |
| logV | +0.066 |
| logL | -0.044 |
| f_gas | +0.026 |

The gradient correlates most strongly with c_V (r = -0.36): galaxies with concentrated rotation curves (high c_V) have *negative* gradients (offset decreases outward). This makes sense: concentrated galaxies have high M/L centers (bulge) where the RAR offset is high, and low M/L outskirts where the offset is lower.

### 4. Gradient by Hubble Type (Test 4)

| Type | N | ⟨gradient⟩ |
|------|---|-----------|
| S0-Sa | 12 | +0.006 |
| Sab-Sb | 26 | +0.008 |
| Sbc-Sc | 30 | +0.002 |
| Scd-Sm | 19 | +0.014 |
| Im-BCD | 41 | -0.009 |

The gradient is nearly zero for all types (±0.01 dex/R_eff). Im-BCD galaxies have a weak negative gradient, possibly because their inner regions are affected by stochastic star formation.

### 5. Gradient Predicts 5-Var Residual (Test 5)

| Metric | Value |
|--------|-------|
| r(gradient, 5-var residual) | **-0.431** |
| 5-var R² | 0.872 |
| 5-var + gradient R² | **0.908** |
| ΔR² | **+0.036** |

**This is a significant finding.** The gradient adds 3.6% to R², the largest single-variable improvement since f_gas. The negative correlation means: galaxies with *steeper* (more negative) gradients have *higher* 5-var residuals. This suggests the gradient captures M/L radial structure that the uniform offset misses.

However, the gradient is not a simple galaxy property — it requires computing the radial RAR profile. It's also the outer offset that better predicts the galaxy offset (r = +0.90) than the inner offset (r = +0.76).

### 6. Offset by Acceleration Regime (Test 6)

The radial gradient differs by acceleration regime:
- **Deep MOND**: Δoffset = -0.006 ± 0.012 (essentially zero)
- **Transition**: Δoffset = +0.001 ± 0.014 (zero)
- **Newtonian**: Δoffset = +0.036 ± 0.034 (marginal positive)

The gradient is only detectable in the Newtonian regime, where M/L gradients matter most (baryonic mass dominates). In deep MOND, the gradient is washed out because the acceleration is too low for baryonic structure to matter.

### 7. The M/L Gradient (Test 7)

If M/L varies with radius: inner M/L > outer M/L (due to bulge, older stars in center), then the offset should be more negative in the inner regions (where M/L is under-estimated by our fixed value) and more positive in outer regions → positive gradient.

The observed mean gradient (+0.002 dex/R_eff) is slightly positive but essentially zero. Strong gradients exist in individual galaxies (33% positive, 30% negative, 37% weak), but they nearly cancel in the mean.

## Physical Interpretation

The radial RAR profile reveals three effects:

1. **Inner offset depression**: The first 0.5 R_eff has systematically low offset (-0.08 dex), reflecting M/L underestimation or beam smearing in galaxy centers.

2. **Scatter gradient**: Inner points have 3× more scatter than outer points (0.30 vs 0.12 dex), because inner rotation curves are affected by non-circular motions, beam smearing, and uncertain M/L.

3. **Gradient-residual correlation**: The r = -0.43 correlation between gradient and 5-var residual shows that radial M/L structure contains information not captured by the galaxy-level model. This could potentially be exploited with radius-dependent offset corrections.

## Grade: B+

A valuable session revealing that the RAR offset has systematic radial structure. The ΔR² = +0.036 from the gradient is the strongest single-variable improvement found since f_gas, though it requires per-galaxy radial analysis. The inner offset depression (-0.08 dex) and the scatter gradient (0.30 → 0.12 dex) provide insight into where M/L and observational errors are concentrated. The gradient-c_V anti-correlation (r = -0.36) connects radial structure to galaxy concentration.

## Files Created

- `simulations/session473_radial_rar_profile.py`: 8 tests
- `Research/Session473_Radial_RAR_Profile.md`: This document

---

*Session #473 verified: 8/8 tests passed*
*Grand Total: 1109/1109 verified*

**Key finding: The RAR offset varies with radius: inner (r < 0.5 R_eff) is -0.083 dex below zero. The gradient correlates with the 5-var residual at r = -0.43, adding ΔR² = +0.036 (the strongest improvement since f_gas). Inner scatter (0.30 dex) is 3× higher than outer (0.12 dex). Gradient anti-correlates with c_V (r = -0.36). The outer offset predicts galaxy offset better (r = 0.90) than inner (r = 0.76). Grade B+.**
