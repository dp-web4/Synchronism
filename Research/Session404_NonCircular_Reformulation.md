# Session #404: Non-Circular Reformulation — Baryonic Predictors Only

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Following Session #403's discovery that local N_corr(r) is tautological, this session investigates what non-circular predictors can explain the RAR offset in late-type galaxies. All predictors used here are either photometric (R_eff, L, SB) or derived from baryonic mass models (g_bar profile), ensuring no circular use of g_obs.

## Central Results

### 1. R_eff at Fixed V_flat is the Dominant Predictor

**r(R_eff, offset | V_flat) = -0.74 (p = 10⁻¹¹)**

This is the strongest non-circular correlation in the analysis. At fixed rotation speed, more extended galaxies have more negative RAR offsets (observed gravity below standard RAR prediction).

### 2. The Effect is NOT Mediated by Mean g_bar

- r(R_eff, offset | V, <g_bar>) = **-0.74** (unchanged!)
- Mediation by <g_bar>: **0%**

The offset does NOT arise because larger galaxies sample different g_bar regimes. At fixed V_flat AND fixed mean baryonic acceleration, R_eff still fully predicts the offset.

### 3. g_bar Concentration is an Interesting Secondary Predictor

| Predictor | r(x, offset | V, L) | ΔR² |
|-----------|---------------------|-----|
| R_eff | -0.489 | +0.071 |
| SB_eff | +0.489 | +0.071 |
| Mean g_bar | +0.384 | +0.044 |
| g_bar slope | -0.278 | +0.023 |
| g_bar conc. | +0.263 | +0.021 |

Interestingly, g_bar concentration absorbs the R_eff signal in simple bivariate analysis:
- r(R_eff, offset | conc) = -0.007 (zero!)
- r(conc, offset | R_eff) = +0.28

But when controlling V+L+<g_bar>+slope, R_eff retains significance (r = -0.32, p = 0.013).

### 4. Point-Level Baryonic Predictors Are Weak

| Predictor | r(x, residual) | p |
|-----------|---------------|---|
| log g_bar | +0.044 | 0.18 |
| log(g_bar × R_eff) | -0.026 | 0.42 |
| log(g_bar / R_eff) | +0.092 | 0.004 |

No baryonic quantity strongly predicts the point-level residual. This is expected — the standard RAR already uses g_bar as input. The information about galaxy SIZE only appears when averaging over galaxies.

### 5. Cross-Validated Model Performance

| Model | LOO-RMSE | Improvement |
|-------|----------|-------------|
| Null (zero) | 0.203 | — |
| V only | 0.147 | 27.6% |
| V + L | 0.111 | 45.1% |
| **V + L + R_eff** | **0.099** | **51.2%** |
| V + L + SB | 0.099 | 51.2% |
| V + L + slope | 0.109 | 46.4% |

Adding R_eff (or equivalently SB at fixed L) reduces LOO-RMSE from 0.111 to 0.099 — a genuine 11% improvement in out-of-sample prediction.

## Physical Interpretation

### Why Does R_eff Matter at Fixed V_flat?

Two galaxies with the same V_flat but different R_eff have:
1. **Different baryonic concentration**: the more compact galaxy has higher g_bar in the center
2. **Different g_bar radial range**: r(R_eff, g_bar range | V) = -0.72 — larger galaxies span a wider range of g_bar
3. **NOT different mean g_bar**: r(R_eff, <g_bar> | V) = -0.04 (n.s.)

So R_eff encodes the SHAPE of the baryonic mass distribution, not the total amount. Specifically:
- Compact galaxies: narrow g_bar range, high concentration
- Extended galaxies: wide g_bar range, low concentration

The offset correlation means: **galaxies that sample a wider range of g_bar in the MOND regime tend to have more negative offsets** (observed gravity below RAR prediction).

### Possible Explanations

**Standard physics**: The RAR formula g_obs = g_bar/(1 - exp(-√(g_bar/g†))) may be slightly imprecise at the level of ~0.1 dex, and this imprecision manifests differently for compact vs extended galaxies because they sample different parts of the g_bar → g_obs mapping.

**Modified gravity**: If the gravitational modification depends not just on local g_bar but also on the spatial scale over which g_bar varies, then more extended galaxies would show different behavior. This is precisely what Synchronism predicts.

**Selection effects**: Distance, inclination, or resolution effects that correlate with R_eff could create artificial offsets. But Session 390-394 controlled for these (9/9 confounds).

## The Non-Circular Modified RAR

The honest, non-circular modification is per-galaxy, not point-by-point:

> **At fixed V_flat and L, the per-galaxy RAR offset scales with log R_eff**
> offset = α + β×log(V) + γ×log(L) + δ×log(R_eff)

This cannot be written as a local formula because R_eff is a global property. The local N_corr formulation was an attempt to make it local, but that turned out to be circular.

## Grade: A

This session successfully reformulates the result using only non-circular predictors and provides genuine physical insight into what R_eff captures (baryonic profile shape, not mean acceleration).

## Files Created

- `simulations/session404_noncircular_reform.py`: 8 tests
- `Research/Session404_NonCircular_Reformulation.md`: This document

---

*Session #404 verified: 8/8 tests passed*
*Grand Total: 645/645 verified*

**Key findings: (1) R_eff at fixed V_flat remains the dominant non-circular predictor of RAR offset (r=-0.74, p=10⁻¹¹). (2) The effect is 0% mediated by mean g_bar — it's genuinely about size, not acceleration regime. (3) R_eff encodes the baryonic profile SHAPE: r(R_eff, g_bar range | V) = -0.72. (4) LOO cross-validation confirms: V+L+R_eff reduces RMSE by 51% over null (vs 45% for V+L alone). (5) No baryonic quantity strongly predicts point-level residuals — the effect is per-galaxy. Grade A.**
