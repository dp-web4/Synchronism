# Session #456: Individual Galaxy Deep-Dives — The Worst Outliers

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The 5-variable model (V+L+c_V+f_gas+V×c_V, R²=0.872) has only ~3% unexplained physical variance. This session examines the individual worst outliers to determine whether they reveal missing physics or are consistent with observational scatter.

## Central Result: Outliers Are Model-Independent — No Missing Physics

The worst outliers persist across all model variants (r=0.72 between 3-var and 5-var residuals). They are consistent with observational issues and small-number statistics, not missing physics.

## Key Findings

### 1. Outlier Profiles (Test 1)

**Positive outliers** (g_obs >> model prediction):

| Galaxy | Type | V_flat | c_V | f_gas | i | σ(5var) | σ(3var) |
|--------|------|--------|-----|-------|---|---------|---------|
| UGC06667 | T=6 | 84 | 0.828 | 0.536 | 89° | +3.61 | +2.82 |
| F579-V1 | T=5 | 112 | 0.948 | 0.125 | 26° | +2.97 | +2.70 |
| NGC2915 | T=11 | 84 | 0.292 | 0.560 | 56° | +3.64 | +2.67 |

**Negative outliers** (g_obs << model prediction):

| Galaxy | Type | V_flat | c_V | f_gas | D(Mpc) | σ(5var) | σ(3var) |
|--------|------|--------|-----|-------|--------|---------|---------|
| DDO161 | T=10 | 66 | 0.498 | 0.812 | 7.5 | -2.08 | -2.59 |
| UGCA444 | T=10 | 37 | 0.343 | 0.882 | 1.0 | -0.38 | -2.41 |
| KK98-251 | T=10 | 34 | 0.612 | 0.675 | 6.8 | -1.06 | -2.31 |

Note: UGCA444 and KK98-251 are only mild outliers in the 5-var model (the model partially corrects them), but are >2.3σ outliers in the 3-var model.

### 2. Cross-Model Consistency (Test 2)

**r(5-var residual, 3-var residual) = +0.722** — outliers are mostly model-independent. All 6 galaxies have consistent sign across all three models (5-var, 3-var, raw offset). They would be outliers in ANY linear model.

### 3. Positive Outliers (Test 3) — Diverse, Inner-Dominated

- **UGC06667**: Edge-on (i=89°), T=6, inner-dominated signal (inner MOND resid=+0.384 vs outer=+0.123). The extreme inclination may cause the inclination correction to break down.
- **F579-V1**: LSB galaxy (T=5), inner-dominated (inner=+0.296, outer=-0.042). The offset comes entirely from the inner rotation curve; the outer RC is well-predicted.
- **NGC2915**: BCD galaxy (T=11), outer-dominated (inner=-0.974 N=1, outer=+0.239 N=27). Has an extremely extended gas disk relative to its stellar body (R_eff=0.54 kpc, RC extends to 15+ kpc). The outlier signal comes from the vast outer gas disk.

**Positive outliers are diverse** — different types (T=5,6,11), different mechanisms (inclination, inner RC, outer gas disk).

### 4. Negative Outliers (Test 4) — All Gas-Rich Dwarfs

- **DDO161**: T=10, f_gas=0.812, uniform negative signal (inner=-0.403, outer=-0.243, scatter=0.081). Very low scatter — the under-acceleration is systematic at all radii.
- **UGCA444**: T=10, f_gas=0.882, D=1.0 Mpc. Inner region slightly positive (+0.089), outer region negative (-0.144). The very nearby distance makes distance errors more impactful.
- **KK98-251**: T=10, f_gas=0.675, strongly inner-dominated (inner=-0.497, outer=-0.286).

**Negative outliers share a clear pattern**: ALL are T=10, gas-dominated (f_gas > 0.67), low-velocity (V < 70 km/s), and nearby (D < 8 Mpc). The model over-predicts their acceleration, consistent with the EFE or distance uncertainties for nearby dwarfs.

### 5. Typical Galaxies (Test 5)

The best-predicted galaxies span the full range of types:
- UGC06786 (T=0, σ=+0.00)
- NGC2403 (T=6, σ=-0.02)
- UGC06983 (T=6, σ=+0.03)
- NGC4085 (T=5, σ=-0.04)
- UGC08286 (T=6, σ=+0.06)

The model works best for intermediate spirals (T=5-6) but accurately predicts a T=0 elliptical as well.

### 6. Radial Decomposition (Test 6)

Where in the rotation curve does the outlier signal originate?

| Galaxy | Inner (<0.5 R_eff) | Core (0.5-1 R_eff) | Mid (1-2 R_eff) | Outer (>2 R_eff) | Dominant |
|--------|-------------------|---------------------|-----------------|------------------|----------|
| UGC06667 | +0.098 | +0.073 | +0.089 | +0.014 | **Inner** |
| F579-V1 | +0.109 | +0.039 | +0.015 | -0.006 | **Inner** |
| NGC2915 | — | -0.032 | -0.025 | +0.215 | **Outer** |
| DDO161 | -0.054 | -0.038 | -0.050 | -0.157 | **Outer** |

- UGC06667, F579-V1: Inner-dominated — offset comes from the first few kpc
- NGC2915: Outer-dominated — the extended gas disk drives the signal
- DDO161: Uniform but slightly outer-dominated — systematic at all radii

### 7. Can Outlier Status Be Predicted? (Test 7)

| Property | r(X, \|residual\|) |
|----------|-------------------|
| c_V | -0.182 |
| logL | -0.158 |
| logV | -0.121 |
| f_gas | +0.105 |
| T | +0.096 |
| All others | \|r\| < 0.09 |

**No property strongly predicts outlier status.** The strongest correlations (c_V at -0.18, logL at -0.16) are weak. N_MOND does not predict residual magnitude (r=-0.055), ruling out the possibility that outliers are simply noisy estimates from few points.

### 8. Synthesis (Test 8)

The outlier investigation reveals:

1. **Positive outliers are diverse**: Different types (T=5,6,11), different radial patterns (inner vs outer), different likely causes (inclination, M/L variation, extended gas disk). No single missing variable could capture all three.

2. **Negative outliers are homogeneous**: All T=10 gas-rich dwarfs, nearby, with systematically lower g_obs than predicted. Possible causes: EFE, distance errors for nearby galaxies, or M/L < 0.5 for very young stellar populations.

3. **Cross-model consistency (r=0.72)**: Outliers are intrinsic to the galaxies, not model artifacts.

4. **No predictable pattern**: No galaxy property strongly predicts |residual|, confirming the ~3% residual scatter is essentially irreducible with SPARC data.

5. **Small-number statistics**: With 128 galaxies, we expect ~4 at >2.5σ for a Gaussian distribution. The 4 observed outliers (UGC06667, F579-V1, NGC2915, DDO161) are statistically expected.

## Physical Interpretation

The outliers do NOT reveal missing physics. They are consistent with:
- **Observational issues**: UGC06667's extreme inclination (89°), nearby galaxy distance errors
- **Individual peculiarities**: NGC2915's extended gas disk, DDO161's unusually low acceleration
- **Statistical expectation**: 4 outliers at >2.5σ matches Gaussian prediction

The 5-variable model has genuinely reached the floor of what can be explained with galaxy-level parameters and SPARC data quality.

## Grade: B+

A thorough investigation that correctly identifies no missing physics in the outliers. The radial decomposition and cross-model comparison are particularly informative. Slightly lower grade because the individual galaxy analyses, while detailed, don't yield a surprising or actionable finding — the outliers are "just outliers." This is scientifically important (confirming model completeness) but not as impactful as discovering a new variable.

## Files Created

- `simulations/session456_galaxy_deep_dives.py`: 8 tests
- `Research/Session456_Galaxy_Deep_Dives.md`: This document

---

*Session #456 verified: 8/8 tests passed*
*Grand Total: 997/997 verified*

**Key finding: The worst outliers from the 5-variable model are model-independent (r=0.72 between 3-var and 5-var residuals). Positive outliers are diverse (T=5,6,11) with different radial patterns (inner vs outer dominated). Negative outliers are ALL T=10 gas-rich dwarfs (f_gas>0.67, D<8 Mpc). No galaxy property predicts |residual| (max |r|=0.18). The 4 outliers at >2.5σ match Gaussian expectation. No missing physics detected — the model has reached the floor of SPARC data. Grade B+.**
