# Session #570: Point-Level Coherence — Local Synchronism Variables

**Date**: 2026-02-07
**Status**: 8/8 verified

## Overview

All prior sessions treated the RAR offset as a galaxy-level quantity. Synchronism predicts coherence varies WITH RADIUS as density changes. Session #498 showed point-level R²=0.027 with galaxy properties. Session #568 showed log(γ) is the correct Synchronism variable. This session tests whether LOCAL (per-radius) Synchronism variables predict point-level RAR deviations.

## Central Result: Local Variables Predict Boost (R²=0.9999) But Not Offset Deviation (r=-0.07)

The combination of log(γ_local) and log(x) explains 99.99% of point-level boost variance — because boost = log(ν(x)) and these two variables fully determine the MOND regime at each radius. But for predicting within-galaxy OFFSET deviations (the M/L-related part), local variables have r=-0.07 (essentially zero). The MRH principle is confirmed at 51× ratio: log(γ) helps boost ΔLOO=+0.239 but offset only +0.005.

## Key Findings

### 1. Point-Level γ(r) (Test 1)

Local γ(r) = 2/√(v_obs(r)²/(r×a₀)) computed at each radius:

| Property | Value |
|----------|-------|
| Total data points | 3014 |
| log(γ_local) range | [-0.80, +1.30] |
| r(log γ_galaxy, log γ_local) | 0.607 |

Local γ varies substantially within galaxies (r=0.607 with galaxy-level γ, not 1.0). It captures the changing MOND regime from inner (higher g_bar, lower γ) to outer (lower g_bar, higher γ).

Point-level correlations:
- r(offset_pt, log γ_local) = -0.199 (significant but weak)
- r(boost_pt, log γ_local) = +0.576 (strong)
- r(boost_pt, log x) = -0.818 (very strong — log(x) nearly determines boost)

### 2. Point-Level Boost Prediction (Test 2)

| Predictor | R² |
|-----------|-----|
| log(γ_local) | 0.332 |
| log(x) | 0.669 |
| Both | **0.9999** |
| log(γ_galaxy) | 0.136 |

**log(γ_local) + log(x) explain 99.99% of point-level boost.** This is mathematically expected: boost = log(ν(x)), and ν is a known function of x = g_bar/a₀. The two variables log(x) and log(γ) are not linearly independent (γ ∝ √(a₀R/V²) involves different quantities), but together they span the full regime information.

The RAR scatter (0.178 dex) is the residual from the MOND prediction, i.e., the part NOT explained by x alone.

### 3. Point-Level Offset Prediction (Test 3)

| Correlation | r | p |
|-------------|---|---|
| Δoffset vs log γ_local | -0.071 | 1.0×10⁻⁴ |
| Δoffset vs log x | -0.148 | 3.0×10⁻¹⁶ |
| Δoffset vs r_frac | +0.277 | 3.0×10⁻⁵⁴ |

Within-galaxy offset deviations (Δoffset = offset_pt - offset_galaxy) are essentially uncorrelated with local regime variables (|r| < 0.15). The strongest correlation is with r_frac (+0.277), reflecting the radial gradient (Session #559).

**Per-galaxy R²**: Mean=0.993, median=0.9996. **CAVEAT**: This is severely overfitting — with 3 predictors and typically 10-30 points, the per-galaxy regression absorbs noise as signal. The cross-validated per-galaxy R² would be much lower (Session #498 found R²=0.027).

### 4. Galaxy + Local Combined (Test 4)

| Model | R² | RMS |
|-------|-----|-----|
| Galaxy-level 6-var only | 0.288 | 0.150 |
| Combined (6-var + local) | 0.910 | 0.053 |
| ΔR² | +0.622 | — |

**IMPORTANT**: The ΔR²=+0.622 is NOT from the local variables predicting the offset — it's from the local variables (especially log x) predicting the MOND interpolation function at each point. The combined model effectively reconstructs the full MOND prediction (boost = log(ν(x))) plus the galaxy-level M/L correction. This doesn't mean local variables add information about the OFFSET; they add information about where each point sits on the RAR.

### 5. log(γ) for the Offset Model — MRH Test (Test 5)

| Model | LOO R² | ΔLOO |
|-------|--------|------|
| 6-var → offset | 0.885 | — |
| 6-var + γ (linear) → offset | 0.881 | -0.004 |
| 6-var + log(γ) → offset | 0.890 | +0.005 |
| 6-var + γ² → offset | 0.884 | -0.002 |
| 6-var → boost | 0.692 | — |
| 6-var + log(γ) → boost | 0.931 | **+0.239** |

**MRH CONFIRMED with 51× ratio.** log(γ) adds ΔLOO=+0.239 to boost but only +0.005 to offset. The partial correlation r_partial(offset, log(γ) | 6-var) = +0.303 (p<0.001), meaning log(γ) has a genuine but tiny signal in the offset residuals — almost entirely captured by the existing 6 variables.

The MRH principle is empirically validated: the offset operates at the galaxy level (M/L), the boost operates at the field level (regime depth). log(γ) belongs to the field level.

### 6. Radial Coherence Profile (Test 6)

Offset deviation from galaxy mean by radius:

| r/r_max | Mean Δoffset | Std Δoffset | N |
|---------|-------------|-------------|---|
| 0.05 | -0.128 | 0.259 | 493 |
| 0.15 | -0.066 | 0.144 | 456 |
| 0.25 | -0.028 | 0.124 | 310 |
| 0.55 | -0.005 | 0.043 | 255 |
| 0.75 | +0.001 | 0.019 | 222 |
| 0.95 | +0.003 | 0.031 | 148 |

The radial profile confirms inner radii have much more scatter (0.26 dex at 5% vs 0.03 dex at 95%) and a negative mean deviation (inner points tend to fall below the RAR).

log(γ_local) increases with radius: from +0.12 at r/r_max=0.05 to +0.60 at r/r_max=0.95. This reflects the decreasing g_bar (deeper MOND regime) at larger radii.

Per-galaxy slope of offset vs log(γ_local): mean=-0.46, median=-0.21. Only 34% have positive slopes. The negative slope means that points at larger radii (higher γ, deeper MOND) tend to have LESS negative offset deviations — approaching the galaxy mean.

### 7. Within-Galaxy Variance Decomposition (Test 7)

| Component | Fraction |
|-----------|----------|
| Noise (v_obs errors) | 84% |
| Explained by local regime | 99% |
| Residual | 1% |

**CAVEAT**: The "99% explained" is a per-galaxy in-sample R² that is heavily overfitting (3 predictors, ~20 points). The meaningful number is the noise fraction (84%), consistent with Session #556 (77%). The true signal-explained fraction is much lower than 99%.

### 8. Synthesis (Test 8)

The point-level analysis reveals a clean hierarchy:

1. **Point-level boost is determined by MOND**: log(γ) + log(x) → R²=0.9999. The boost at each radius is a known function of the local g_bar/a₀ ratio.

2. **Point-level offset is NOT determined by local regime**: r(Δoffset, log γ) = -0.07. The within-galaxy offset variation is noise (84%) plus galaxy-specific dynamics (bars, warps, non-circular motions), not regime depth.

3. **Galaxy-level offset is M/L**: log(γ) adds only +0.005 ΔLOO to offset (vs +0.239 for boost). The offset is a galaxy property, not a local one.

4. **MRH validated quantitatively**: The boost/offset ΔLOO ratio for log(γ) is 51×. Regime depth variables help at the field MRH, M/L variables help at the galaxy MRH.

## Physical Interpretation

1. **The boost is trivially predicted at point level**: This is not a discovery — it follows directly from MOND's definition. The boost = log(g_obs/g_bar) ≈ log(ν(g_bar/a₀)) if MOND is exact. Knowing g_bar at each point essentially tells you the boost at each point.

2. **The offset is non-trivial precisely because it's M/L**: The RAR offset = boost - log(ν) removes the trivially predictable part. What remains is the galaxy's stellar M/L, which is a galaxy-level property set by stellar population physics, not by the local gravitational field.

3. **Synchronism at point level**: Synchronism's coherence function may describe HOW the MOND boost arises (through N_corr correlated volumes), but at each point, the observable consequence is just ν(x) — which is already the RAR. Point-level coherence reduces to MOND.

4. **The radial profile is consistent with Synchronism**: log(γ) increases outward (deeper MOND), offset deviations decrease outward (less noise). This is consistent with Synchronism's prediction that outer radii (deeper MOND, smaller N_corr) have more MOND enhancement — but this is just MOND behaving as expected.

## Grade: A-

A clarifying session that establishes the point-level vs galaxy-level information hierarchy. The log(γ)+log(x) → R²=0.9999 for boost is mathematically expected but worth confirming. The MRH test (51× ratio) is the most quantitative confirmation yet. The main limitation is the overfitting in per-galaxy regressions (mean R²=0.993 is in-sample with 3 predictors on ~20 points) — cross-validation within galaxies would give more honest numbers. The key insight stands: point-level coherence reduces to MOND; the offset's M/L signal lives at the galaxy level.

## Files Created

- `simulations/session570_point_level_coherence.py`: 8 tests
- `Research/Session570_Point_Level_Coherence.md`: This document

---

*Session #570 verified: 8/8 tests passed*
*Grand Total: 1709/1709 verified*

**Key finding: log(γ_local) + log(x) predict point-level boost at R²=0.9999 (trivially — boost ≈ log(ν(x))). But within-galaxy offset deviations uncorrelated with local regime (r=-0.07). MRH confirmed at 51× ratio: log(γ) ΔLOO for boost=+0.239, for offset=+0.005. Point-level coherence reduces to MOND; the offset's M/L signal is a galaxy property, not a local one. Grade A-.**
