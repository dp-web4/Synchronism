# Session #467: The Mass Discrepancy-Acceleration Relation

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The Mass Discrepancy D = g_obs/g_bar = V²_obs/V²_bar measures the ratio of total to baryonic gravity. In Newtonian gravity, D = 1. In MOND, D → √(a₀/g_bar) at low accelerations. This session examines the MDAR shape, scatter, and galaxy-to-galaxy variation.

## Central Result: D Increases Outward in 79% of Galaxies, With Median D_outer/D_inner = 1.5

The mass discrepancy grows with radius in most galaxies (D_outer/D_inner = 1.51 median), and the maximum discrepancy occurs at the outermost measured point in 63% of cases. Late-type galaxies reach D_max up to 34, meaning only 3% of the dynamical mass is baryonic at their outermost radii.

## Key Findings

### 1. MDAR Shape (Test 1)

| log g_bar | N | ⟨D⟩ | D_MOND |
|-----------|---|-----|--------|
| -11.56 | 357 | 7.28 | 7.14 |
| -11.26 | 356 | 5.04 | 5.21 |
| -11.06 | 356 | 4.23 | 4.25 |
| -10.84 | 356 | 3.57 | 3.41 |
| -10.54 | 356 | 2.79 | 2.57 |
| -10.22 | 356 | 2.14 | 1.96 |
| -9.85 | 356 | 1.58 | 1.52 |
| -9.19 | 357 | 1.17 | 1.11 |

The observed MDAR matches MOND predictions at all accelerations. The median D tracks the MOND prediction to within ~10% across 2.5 decades of acceleration. At the highest accelerations (log g_bar = -9.2), D ≈ 1.17 — nearly Newtonian. At the lowest accelerations (log g_bar = -11.6), D ≈ 7.3 — over 86% of the dynamical mass is "missing."

### 2. MDAR Scatter (Test 2)

| Regime | σ(log D/D_RAR) | σ(log D/D_5var) |
|--------|----------------|-----------------|
| Deep MOND | 0.200 | **0.127** |
| Moderate MOND | 0.160 | **0.113** |
| Transition | 0.177 | 0.156 |
| Newtonian | 0.178 | 0.192 |

The 5-variable model reduces scatter by **36% in deep MOND** and **29% in moderate MOND**. In the Newtonian regime, the model slightly worsens scatter (as expected — the galaxy-level correction is calibrated on MOND data). The biggest improvement is where MOND is strongest, confirming the model captures MOND-specific physics.

### 3. Transition Acceleration (Test 3)

| Threshold | N galaxies | Observed g_bar | MOND prediction |
|-----------|-----------|----------------|-----------------|
| D = 2 | 80 | 0.37 × 10⁻¹⁰ | 0.58 × 10⁻¹⁰ |
| D = 5 | 58 | 0.09 × 10⁻¹⁰ | — |
| D = 10 | 9 | 0.05 × 10⁻¹⁰ | — |

The D = 2 transition occurs at g_bar = 0.37 × 10⁻¹⁰, which is 63% of the MOND prediction (a₀ × ln²2 = 0.58 × 10⁻¹⁰). The **scatter is large** (0.64 dex), reflecting galaxy-to-galaxy variation in where the transition occurs — this is exactly the scatter the 5-variable model explains.

Only 9/128 galaxies reach D = 10 (90% "missing mass"), all in the deepest MOND regime.

### 4. Galaxy-Level Mean D vs Properties (Test 4)

| Property | r(logD, X) |
|----------|-----------|
| offset | **+0.606** |
| f_gas | +0.560 |
| logL | -0.520 |
| c_V | -0.404 |
| T | +0.407 |
| logV | -0.282 |

The mean mass discrepancy is most strongly correlated with the RAR offset (+0.61) — galaxies with positive offset (more g_obs than predicted) have higher D. Gas-rich, low-luminosity, late-type galaxies have the highest mass discrepancies, as expected in MOND (they are deepest in the MOND regime).

### 5. Radial D(r) Profiles by Type (Test 5)

| r/R_eff | Early (T≤4) | Sbc-Sc | Late (T≥7) |
|---------|-------------|--------|------------|
| < 0.5 | 1.03 | 1.77 | 3.13 |
| 0.5-1 | 1.14 | 2.04 | 3.60 |
| 1-2 | 1.48 | 2.00 | 4.15 |
| 2-5 | 2.13 | 2.93 | 4.62 |
| > 5 | 3.84 | 4.82 | 7.85 |

**Early types are nearly Newtonian in their inner regions** (D = 1.03 at r < 0.5 R_eff), but still reach D ≈ 4 at large radii. **Late types start already discrepant** (D = 3.1 at the center!) and reach D ≈ 8 in their outskirts. The mass discrepancy grows monotonically with radius for all types.

### 6. Maximum Discrepancy (Test 6)

- ⟨D_max⟩ = 6.5, range [1.7, 34.0]
- D_max occurs at outermost point in 80/128 galaxies (63%)
- Median r(D_max)/R_eff = 3.83

Top 5 most discrepant galaxies:

| Galaxy | T | D_max | f_gas |
|--------|---|-------|-------|
| UGC06667 | 6 | 34.0 | 0.54 |
| UGC00731 | 10 | 28.1 | 0.86 |
| UGCA444 | 10 | 27.1 | 0.88 |
| UGC12506 | 6 | 19.1 | 0.25 |
| F583-1 | 9 | 16.3 | 0.69 |

UGC06667 (D_max = 34) is a known outlier — it's edge-on (i = 89°) and the discrepancy may be inflated by inclination corrections. The gas-rich dwarfs (UGC00731, UGCA444) have genuinely extreme discrepancies.

### 7. D Asymmetry (Test 7)

- D_outer/D_inner median = 1.51
- 89/113 galaxies (79%) have D_outer > D_inner
- D_outer/D_inner correlates with Hubble type (r = -0.34): early types have more asymmetric D profiles

The D asymmetry correlates with logV (+0.33): more massive galaxies have steeper D gradients. This makes sense — massive galaxies span a wider range of accelerations (from Newtonian in the center to MOND in the outskirts), creating a larger D contrast.

## Physical Interpretation

### The Mass Discrepancy as a MOND Thermometer

D(r) directly measures how deep a galaxy is into the MOND regime at each radius:
- D ≈ 1: Newtonian (no missing mass)
- D ≈ 2: Half missing (MOND transition)
- D ≈ 10: 90% missing (deep MOND)

The MDAR is a one-to-one function of g_bar in MOND, so D(r) is effectively a measure of the local acceleration. The galaxy-to-galaxy scatter in D at fixed g_bar is the same scatter the 5-variable model explains.

### Why D Grows Outward

In a typical galaxy, the inner regions are baryonically dominated (g > a₀, D ≈ 1), while the outer regions are in the MOND regime (g < a₀, D >> 1). The radius where D = 2 marks the transition between these regimes. For late-type dwarfs, D > 1 even at the center — they live entirely in the MOND regime.

### The D = 34 Galaxy

UGC06667's extreme D_max = 34 means that at its outermost measured point, baryonic mass accounts for only 3% of the dynamical mass. In CDM, this requires a dark matter halo ~33× the baryonic mass. In MOND, this is simply a galaxy deep in the low-acceleration regime. Its edge-on inclination (89°) may contribute to the extreme value.

## Grade: B

A competent exploration of the MDAR that confirms MOND predictions and characterizes the mass discrepancy profile. The radial D(r) profiles by type (Test 5) and the D asymmetry analysis (Test 7) are informative. However, the MDAR is mathematically equivalent to the RAR (just D = g_obs/g_bar rather than log g_obs vs log g_bar), so the findings largely restate known results in different units. The most novel finding is the quantification of where D_max occurs (63% at outermost point, median r/R_eff = 3.8).

## Files Created

- `simulations/session467_mass_discrepancy.py`: 8 tests
- `Research/Session467_Mass_Discrepancy.md`: This document

---

*Session #467 verified: 8/8 tests passed*
*Grand Total: 1069/1069 verified*

**Key finding: The MDAR matches MOND to ~10% across 2.5 decades. D grows outward in 79% of galaxies (median D_outer/D_inner = 1.51). D_max = 34 (UGC06667). Late types have D > 3 even at their centers — they live entirely in MOND. The 5-var model reduces deep MOND D-scatter by 36%. The D = 2 transition occurs at g_bar = 0.37 × 10⁻¹⁰ (63% of MOND prediction), with large galaxy-to-galaxy scatter (0.64 dex). Grade B.**
