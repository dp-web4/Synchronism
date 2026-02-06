# Session #465: Predicted Rotation Curves From the 5-Variable Model

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

The 5-variable model predicts each galaxy's RAR offset from (V, L, c_V, f_gas). This session asks: can we use this to predict full rotation curves? By applying the galaxy-level offset correction to the standard RAR, we generate predicted rotation curves without fitting any per-galaxy parameters.

## Central Result: Median 12% Accuracy Without Per-Galaxy Fitting

The 5-variable model predicts individual rotation curves to a median fractional accuracy of 11.5%. 80% of galaxies are predicted to better than 20%. The correlation r(V_obs, V_pred) = 0.978 improves over the raw RAR (r = 0.963).

## Key Findings

### 1. Best-Predicted Rotation Curves (Test 1)

| Galaxy | Type | V_flat | RC RMS (RAR) | RC RMS (5-var) | Improvement |
|--------|------|--------|--------------|----------------|-------------|
| UGC06786 | T=0 | 219 | 31.2 km/s | 15.1 km/s | +52% |
| NGC2403 | T=6 | 131 | 6.4 km/s | 4.4 km/s | +31% |
| UGC06983 | T=6 | 109 | 9.2 km/s | 8.5 km/s | +8% |
| NGC4085 | T=5 | 132 | 22.1 km/s | 14.7 km/s | +34% |
| UGC08286 | T=6 | 82 | 11.2 km/s | 5.6 km/s | +50% |

The best-predicted galaxies span the full mass range. NGC2403 (a benchmark galaxy in SPARC) is predicted to 4.4 km/s RMS — 3% of V_flat.

### 2. Worst-Predicted Rotation Curves (Test 2)

| Galaxy | Type | V_flat | Offset resid | RC RMS (5-var) |
|--------|------|--------|-------------|----------------|
| NGC2915 | T=11 | 84 | +0.202 | 23.0 km/s |
| UGC06667 | T=6 | 84 | +0.201 | 12.6 km/s |
| F579-V1 | T=5 | 112 | +0.165 | 18.8 km/s |

The worst outliers are the same galaxies identified in Session 456 (deep-dive outliers). NGC2915's extended gas disk remains unpredictable from global properties. The 5-var model actually **worsens** NGC2915 slightly (RAR: 20.6 → 5-var: 23.0 km/s) because the offset correction is in the wrong direction for this atypical galaxy.

### 3. Full-Sample Statistics (Test 3)

| Metric | RAR-only | 5-variable |
|--------|----------|------------|
| Point-level RMS | 22.9 km/s | 18.1 km/s |
| Point-level frac RMS | 28.6% | 22.1% |
| Median galaxy frac RMS | 15.2% | 8.7% |
| Galaxies improved | — | 97/128 (76%) |

The 5-variable model improves 76% of galaxies and reduces the median fractional RMS from 15.2% to 8.7%.

### 4. Radial Performance (Test 4)

| r/R_eff | N | RMS (RAR) | RMS (5-var) | Improvement |
|---------|---|-----------|-------------|-------------|
| < 0.5 | 361 | 27.9 | 32.4 | **-16%** |
| 0.5-1 | 392 | 26.1 | 22.9 | +12% |
| 1-2 | 688 | 23.0 | 15.3 | +34% |
| 2-5 | 918 | 19.8 | 11.9 | +40% |
| > 5 | 491 | 21.0 | 10.9 | +48% |

**Critical finding**: The 5-var model **degrades** the innermost regions (r < 0.5 R_eff) by 16%. This is because the galaxy-level offset is applied uniformly at all radii, but the inner RC has different physics (baryonic dominance, beam smearing, non-circular motions). The improvement increases monotonically with radius, reaching +48% at r > 5 R_eff.

**Implication**: A radially-varying correction would further improve predictions, but would require per-galaxy fitting rather than global properties alone.

### 5. Fractional Velocity Residual Profile (Test 5)

| r/R_eff | ⟨ΔV/V⟩ RAR | ⟨ΔV/V⟩ 5-var |
|---------|-----------|-------------|
| < 0.5 | -12.1% | -12.1% |
| 0.5-1 | -4.8% | -1.5% |
| 1-2 | -2.8% | +0.1% |
| 2-5 | -5.3% | -0.8% |
| > 5 | +1.1% | +2.0% |

The RAR systematically underestimates velocities at all inner radii (⟨ΔV/V⟩ < 0). The 5-var correction largely removes this bias at r > 0.5 R_eff, bringing the mean residual to near zero. The inner region (r < 0.5 R_eff) is unchanged — the model cannot correct inner RC shapes.

### 6. Velocity Correlation (Test 6)

| Metric | RAR | 5-variable |
|--------|-----|------------|
| r(V_obs, V_pred) | 0.9625 | 0.9776 |

By velocity range:
| Range | N | RMS (RAR) | RMS (5-var) | Frac RMS |
|-------|---|-----------|-------------|----------|
| V < 50 | 395 | 15.3 | 11.0 | 50.6% |
| 50 ≤ V < 150 | 1129 | 16.4 | 12.4 | 14.7% |
| V ≥ 150 | 1326 | 28.7 | 23.2 | 10.3% |

Low-velocity points have the highest fractional residuals (50.6%) even after correction — these are the rising portions of dwarf galaxy rotation curves, where beam smearing and non-circular motions dominate.

### 7. Predicted RC Quality (Test 7)

| Quality | Frac RMS | Count | Fraction |
|---------|----------|-------|----------|
| Excellent | < 10% | 49 | 38% |
| Good | < 20% | 103 | 80% |
| Fair | < 30% | 118 | 92% |
| Poor | ≥ 30% | 10 | 8% |

**80% of rotation curves are predicted to better than 20% without any per-galaxy fitting.** This is remarkable — given only four global numbers (V_flat, L, c_V, f_gas), the model generates a reasonable rotation curve.

## Physical Interpretation

### Why Does a Galaxy-Level Correction Work?

The 5-variable offset acts as a vertical shift of the entire RAR for each galaxy. Since the RAR maps g_bar → g_obs, shifting the RAR uniformly shifts log(g_obs) at all radii. This translates to a multiplicative correction to g_obs, which becomes an additive correction to V²(r). At the flat part of the rotation curve (where the offset is measured), this correction is exact. At inner radii, the correction is approximate because:

1. The inner RC is dominated by baryonic physics (disk/bulge mass distribution)
2. The MOND modification varies with radius (deep MOND in outer, Newtonian in inner)
3. The galaxy-level offset averages over radial structure that varies point-to-point

### The -16% Inner Degradation

The 5-var model worsens the inner RC prediction because it applies a uniform offset correction that was calibrated on the MOND regime (g < g†). For the innermost points (r < 0.5 R_eff), many galaxies are in the Newtonian regime (g > g†) where the offset has a different meaning. Applying the MOND-calibrated correction to Newtonian-regime points systematically over-corrects.

### The 48% Outer Improvement

At large radii (r > 5 R_eff), galaxies are deep in the MOND regime, the RAR mapping is most uncertain, and the galaxy-level correction is most effective. This is where the 5-variable model's knowledge of stellar populations (M/L via V+L), mass geometry (c_V), and gas content (f_gas) provides the most information.

## Grade: B+

A solid session that demonstrates the 5-variable model's predictive power extends from galaxy-level offsets to full rotation curves. The 80% of galaxies predicted to <20% is a strong result. The discovery that inner regions are degraded (-16% at r < 0.5 R_eff) while outer regions are dramatically improved (+48% at r > 5 R_eff) provides physical insight into where the galaxy-level correction works and fails. Slightly lower grade because the approach is straightforward (apply offset uniformly) and the inner-region degradation limits practical applicability.

## Files Created

- `simulations/session465_predicted_rotation_curves.py`: 8 tests
- `Research/Session465_Predicted_Rotation_Curves.md`: This document

---

*Session #465 verified: 8/8 tests passed*
*Grand Total: 1053/1053 verified*

**Key finding: The 5-variable model predicts rotation curves to 11.5% median accuracy (80% of galaxies < 20%) from global properties alone, without per-galaxy fitting. The correction improves outer regions by +48% (r > 5 R_eff) but degrades inner regions by -16% (r < 0.5 R_eff) because the galaxy-level offset doesn't capture radial structure. r(V_obs, V_pred) = 0.978, up from 0.963 for the raw RAR. The worst outliers are the same galaxies from Session 456. Grade B+.**
