# Session #482: Residual Forensics — What the Model Misses

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The outer-only 5-variable model has R² = 0.911 and RMS = 0.048 dex. This session forensically examines the residuals to understand what creates them, which galaxies are persistent outliers, and whether hidden patterns remain.

## Central Result: 41% of Galaxies Within Measurement Noise; 7 Persistent Outliers

41% of galaxies (52/128) have outer-model residuals smaller than the estimated measurement noise (0.025 dex). Seven galaxies are persistent outliers in both full and outer models. The irreducible scatter is ~0.041 dex (~10% in velocity), consistent with residual M/L variation.

## Key Findings

### 1. The Outlier Hall of Fame (Test 1)

| Galaxy | Full resid | Outer resid | T | V | Notes |
|--------|-----------|------------|---|---|-------|
| NGC2915 | +0.202 | +0.137 | 11 | 84 | Extended gas disk |
| UGC06667 | +0.201 | +0.140 | 6 | 84 | High c_V, gas-rich |
| PGC51017 | -0.091 | -0.110 | 11 | 19 | Q=3, poorly measured |
| DDO161 | -0.116 | -0.075 | 10 | 66 | Gas-rich irregular |
| UGC01281 | -0.134 | +0.074 | 8 | 55 | Sign FLIPS between models |
| UGC07603 | +0.081 | +0.079 | 7 | 62 | Stable outlier |
| NGC0891 | -0.070 | -0.078 | 3 | 216 | Edge-on, massive |

**UGC01281 flips sign** between full and outer models (+0.074 outer vs -0.134 full), meaning its inner and outer offsets disagree strongly. NGC2915 is the classic extended-gas-disk outlier that no galaxy-level model can fix.

### 2. Residual Autocorrelation (Test 2)

| k-NN | Full r | Outer r |
|------|--------|---------|
| k=1 | +0.035 | **+0.463** |
| k=3 | +0.016 | +0.417 |
| k=5 | +0.064 | +0.405 |
| k=10 | +0.050 | +0.370 |

**The outer model shows significant residual autocorrelation (r = +0.46)**: similar galaxies tend to have similar residuals. The full model does NOT show this (r ≈ 0). This means the outer offset contains structure not captured by the 5 variables that is correlated with galaxy properties. However, this may also reflect the outer offset being noisier (fewer points) in correlated ways.

### 3. Residual Predictors (Test 3)

| Property | r(X, \|resid_outer\|) | r(X, resid_outer) |
|----------|---------------------|-------------------|
| RC roughness | +0.19 | +0.15 |
| RC gradient | -0.07 | **-0.20** |
| N_mond | -0.01 | **+0.29** |
| N_points | -0.03 | +0.22 |

**N_mond correlates with signed residual (r = +0.29)**: galaxies with more MOND points have more positive residuals. This could mean the offset measurement becomes biased with more points (averaging over a larger radial range) or that galaxies with extended MOND regions have systematically different offsets.

**RC gradient correlates at r = -0.20**: galaxies with steeper outward decline in the RAR residual have more negative model residuals.

### 4. Discrepant Pairs (Test 4)

UGC06667 appears in 7 of the top 10 most discrepant close pairs — it is the most anomalous galaxy in parameter space. It has similar properties to many Sbc-Sd galaxies but an unusually high offset (+0.14).

The top pair (NGC2366 vs UGC00731) has near-identical properties but opposite residuals (+0.044 vs -0.144). NGC2366 is at 3.3 Mpc (nearby) while UGC00731 is at 12.5 Mpc — distance may play a subtle role through resolution effects.

### 5. Outlier Removal (Test 5)

| Sample | R² | RMS | N |
|--------|-----|-----|---|
| All | 0.911 | 0.048 | 128 |
| Remove top 5 | 0.927 | 0.040 | 123 |
| Remove top 10 | 0.939 | 0.037 | 118 |
| Q=1 only | 0.855 | 0.044 | 84 |
| i < 80° | 0.918 | 0.048 | 104 |

Removing top 10 outliers pushes R² to 0.939. Q=1-only galaxies actually show LOWER R² (0.855) — this is because Q=1 galaxies are more homogeneous (less offset variance to explain).

### 6. Residual Stability (Test 6)

r(resid_full, resid_outer) = 0.57 — moderate correlation. The same galaxies tend to be outliers, but the outer model reshuffles the rankings significantly. 5 galaxies improve dramatically (F579-V1: 0.165 → 0.049), while 5 degrade (UGCA444: 0.021 → 0.146).

### 7. Irreducible Scatter (Test 7)

| Component | Value |
|-----------|-------|
| Total RMS | 0.048 dex |
| Measurement noise | ~0.025 dex |
| **Irreducible** | **~0.041 dex** |
| In velocity | ~10% |

41% of galaxies have residuals within measurement noise (< 0.025 dex). The irreducible ~10% velocity scatter sets the fundamental limit for galaxy-level RAR prediction.

## Physical Interpretation

### The 10% Velocity Floor

The irreducible scatter (~0.041 dex, ~10% in V) has three possible origins:
1. **M/L variation**: Even at fixed V, L, c_V, f_gas, galaxies have different stellar populations (age, metallicity) that create ~5% V scatter
2. **Mass geometry**: The disk thickness, warp, and asymmetry create ~5% V scatter that c_V doesn't capture
3. **Gas distribution**: The gas mass profile varies between galaxies with the same integrated f_gas

These sources add in quadrature to ~10%.

### The Autocorrelation Signal

The NN autocorrelation (r = 0.46 in outer model) suggests a "6th variable" exists that is correlated with the 5 model variables. Candidate: the **radial scale length** or **concentration parameter** that was partially identified in Sessions 390-393. The SB partial correlation (Session 481, r = -0.49) also points to this.

## Grade: B+

A thorough forensic analysis that establishes the 10% velocity floor, identifies 7 persistent outliers, and reveals NN autocorrelation in the outer model. The discrepant-pair analysis fingering UGC06667 is informative. The outlier removal test (R² → 0.939 without 10 galaxies) shows the model is excellent for 90% of the sample. Slightly lower grade because no new variable was discovered, though the autocorrelation finding is suggestive.

## Files Created

- `simulations/session482_residual_forensics.py`: 8 tests
- `Research/Session482_Residual_Forensics.md`: This document

---

*Session #482 verified: 8/8 tests passed*
*Grand Total: 1173/1173 verified*

**Key finding: 41% of galaxies within measurement noise; 7 persistent outliers (NGC2915, UGC06667 dominant). NN autocorrelation r = +0.46 in outer model (structure remains). Irreducible scatter ~0.041 dex (~10% velocity). Removing 10 outliers: R² → 0.939. RC gradient and N_mond weakly predict signed residual. No single hidden variable explains residuals. Grade B+.**
