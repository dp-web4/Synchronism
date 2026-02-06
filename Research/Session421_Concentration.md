# Session #421: Baryonic Concentration and Mass Profile Shape

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The 69% unexplained fraction of the R_eff effect (Session 415) demands explanation. The outward amplification (Session 417) suggests the baryonic mass DISTRIBUTION matters, not just total mass or size. This session tests whether concentration and profile shape metrics carry additional information beyond V_flat and R_eff.

## Central Result: c_V Carries Strong Information Beyond V + R_eff

| Metric | r(metric, offset \| V, R) | p | ΔR² over V+R | LOO-RMSE |
|--------|--------------------------|---|---------------|----------|
| **c_V = V(R_eff)/V_flat** | **+0.53** | **1×10⁻⁵** | **+0.070** | **0.087** |
| RC slope (outer) | -0.39 | 2×10⁻³ | +0.038 | 0.089 |
| g_bar gradient | -0.27 | 0.04 | +0.018 | 0.092 |
| g_bar range (MOND) | -0.22 | 0.08 | +0.012 | 0.094 |
| c_disk (r_half/R_eff) | -0.10 | 0.44 | +0.003 | 0.095 |

**c_V is a new third predictor** that reduces LOO-RMSE from 0.101 to 0.087 — a 14% further improvement beyond V + R_eff.

## Detailed Findings

### 1. Concentration Metrics (Test 1)

Five concentration metrics were computed from rotation curve data:
- **c_V** = V(R_eff)/V_flat: how quickly the RC reaches its asymptote (mean 0.71)
- **g_bar gradient**: log(g_bar_inner) - log(g_bar_outer) (mean +0.30)
- **RC slope**: dlog(V)/dlog(r) in outer region (mean +0.24, most are still rising)
- **c_disk**: r_half_disk/R_eff (mean 0.44)
- **g_bar range**: dynamic range of g_bar in MOND regime (mean 0.64 dex)

### 2. c_V: The Strongest New Predictor (Tests 2, 5)

At fixed V_flat and R_eff, **r(c_V, offset) = +0.53 (p = 10⁻⁵)**. This means:
- **High c_V** (RC reaches V_flat quickly → concentrated mass) → **more positive offset** (more acceleration than standard RAR)
- **Low c_V** (RC still rising at R_eff → diffuse mass) → **more negative offset**

Adding c_V to V + R_eff:
- ΔR² = +0.070 (7% additional variance)
- F-test: p = 2×10⁻⁵
- LOO-RMSE: 0.101 → **0.087** (14% improvement)

Crucially, c_V does NOT mediate R_eff — it's independent information. Controlling c_V actually STRENGTHENS the R_eff effect (mediation = -12%), meaning R_eff and c_V capture different aspects of galaxy structure.

### 3. RC Slope: Significant Third Predictor (Test 4)

The outer rotation curve slope carries information beyond V + R:
- r(RC slope, offset | V, R) = **-0.39** (p = 0.002)
- Galaxies with steeper (more positive) outer slopes have more negative offsets

RC slope is weakly correlated with R_eff (r = -0.14 at fixed V), so it captures mostly independent information.

### 4. Outward Amplification Confirmed and Quantified (Test 7)

The R_eff effect in inner vs outer regions (N = 55):

| Region | r(R_eff, offset \| V) | R coefficient | RMS |
|--------|----------------------|---------------|-----|
| Inner (r < 2 R_eff) | -0.23 (p = 0.08) | -0.165 | 0.183 |
| Outer (r > 2 R_eff) | **-0.84** (p = 10⁻¹⁵) | **-0.444** | 0.077 |
| Overall | -0.74 (p = 10⁻¹¹) | -0.365 | 0.096 |

The R coefficient is **2.7× stronger** in the outer region. The inner signal is marginally significant at best. The V+R model achieves its best performance in the outer region (RMS = 0.077).

The offset gradient (outer - inner) is predicted by R_eff at fixed V: r = -0.39 (p = 0.003).

### 5. g_bar Gradient and Dynamic Range (Tests 3, 6)

- g_bar gradient weakly predicts offset beyond V+R: r = -0.27
- g_bar gradient does NOT explain the offset gradient once R_eff is controlled (r drops from 0.28 to 0.08)
- g_bar dynamic range mediates only 7% of the R_eff effect

These quantities are collinear with R_eff and carry little independent information.

## Physical Interpretation

The c_V result reveals a new dimension of the R_eff-dependent RAR. At fixed V_flat and R_eff, galaxies with more centrally concentrated mass (high c_V) show more observed acceleration than predicted. This makes physical sense:

1. **High c_V** → mass is centrally concentrated → at outer radii, g_bar drops rapidly → the RAR formula, which is calibrated on average profiles, overestimates how much the acceleration should drop → apparent positive offset
2. **Low c_V** → mass is diffuse → g_bar profile is flatter → less extreme outer behavior

This is consistent with the outward amplification: the inner region (where most galaxies have similar g_bar regardless of c_V) shows a weak R_eff effect, while the outer region (where profile shape matters most) shows a strong effect.

The three-parameter model (V, R_eff, c_V) achieves LOO-RMSE = 0.087, reducing scatter by 57% compared to the standard RAR.

## Grade: A

A genuine new finding. c_V at r = +0.53 beyond V+R is the strongest third predictor we've found. The 14% LOO improvement is real and cross-validated. The physical interpretation is clear and connects naturally to the outward amplification. The inner/outer decomposition (2.7× amplification of R coefficient) is quantitatively sharp. Slightly below A+ because c_V may partly reflect the same physics as R_eff (both are structural measures).

## Files Created

- `simulations/session421_concentration.py`: 8 tests
- `Research/Session421_Concentration.md`: This document

---

*Session #421 verified: 8/8 tests passed*
*Grand Total: 765/765 verified*

**Key finding: c_V = V(R_eff)/V_flat is a strong third predictor beyond V + R_eff: r = +0.53 (p = 10⁻⁵), reducing LOO-RMSE from 0.101 to 0.087 (14% improvement). RC slope also significant (r = -0.39). The R_eff effect is 2.7× stronger in outer regions (R coeff: -0.165 inner vs -0.444 outer). Inner signal is marginally significant; the effect is dominated by outer radii. c_V and R_eff are independent — controlling c_V strengthens the R_eff effect. Grade A.**
