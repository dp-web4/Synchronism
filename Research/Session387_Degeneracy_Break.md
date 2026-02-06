# Session #387: Breaking the N_corr-Vflat Degeneracy

**Date**: 2026-02-06
**Status**: 8/8 verified

## Overview

Session #386 found that N_corr = V²_flat/(R_eff × a₀) is 97% correlated with Vflat (on the full 171-galaxy sample), suggesting it might be just a fancy mass proxy. This session directly attacks the degeneracy by testing whether R_eff contributes independent information at fixed Vflat.

## Key Result: R_eff IS Significant — Degeneracy Partially Broken (Grade B+)

At fixed Vflat, R_eff significantly predicts RAR offset: r(R_eff, offset | Vflat) = -0.306 (p = 0.0002). The sign is correct (more compact → higher offset → closer to standard RAR). The effect is 3x stronger in the MOND regime than the Newtonian regime, consistent with coherence interpretation.

## Detailed Findings

### 1. R_eff Residual Test (THE KEY RESULT)

| Model | R² |
|---|---|
| Vflat alone | 0.161 |
| Vflat + R_eff | 0.240 |
| **ΔR² from R_eff** | **+0.078** |

- r(R_eff, offset | Vflat) = **-0.306** (p = 0.0002)
- r(N_corr, offset | Vflat) = +0.306 (p = 0.0002)

**R_eff adds 7.8 percentage points of explanatory power beyond Vflat.** This is not a trivial contribution. The sign is correct: larger R_eff (lower N_corr at same mass) → more negative offset (further from standard RAR).

### 2. Compact vs Extended at Matched Vflat

After residualizing R_eff on Vflat (controlling for mass):
- **Compact galaxies**: mean offset = -0.003
- **Extended galaxies**: mean offset = -0.070
- **Difference**: +0.067 dex (p = 0.010, permutation test)

Within Vflat tertiles:
| Vflat bin | Compact - Extended | N_comp | N_ext |
|---|---|---|---|
| Low (18-80 km/s) | **+0.209** | 23 | 22 |
| Mid (80-113 km/s) | +0.002 | 20 | 25 |
| High (113+ km/s) | -0.003 | 24 | 21 |

**The effect is concentrated in low-mass galaxies** — exactly where MOND effects dominate and coherence should matter most.

### 3. N_corr Outlier Behavior

Galaxies with unusually high or low N_corr for their Vflat:

| Group | N | Mean offset |
|---|---|---|
| High N_corr outliers (>1σ) | 24 | +0.025 |
| Normal | 94 | -0.036 |
| Low N_corr outliers (<-1σ) | 17 | -0.128 |
| **High - Low** | | **+0.153 dex** |

Theory predicts high-N_corr outliers should sit higher. ✓ Correct.

Notable outliers:
- **NGC2915** (logV=1.92, high N_corr): offset = +0.158 (ABOVE RAR)
- **NGC1705** (logV=1.86, high N_corr): offset = +0.247 (well ABOVE RAR)
- **F561-1** (logV=1.70, low N_corr): offset = -0.459 (well BELOW RAR)
- **UGC06628** (logV=1.62, low N_corr): offset = -0.478 (well BELOW RAR)

### 4. Binned Vflat Analysis

Within Vflat bins, r(R_eff, offset):

| Vflat bin | N | r(R_eff, offset) | p |
|---|---|---|---|
| 18-80 km/s | 34 | **-0.361** | 0.029 |
| 80-113 km/s | 33 | -0.321 | 0.059 |
| 113-183 km/s | 34 | -0.055 | 0.758 |
| 183-340 km/s | 34 | +0.203 | 0.241 |

Weighted average: -0.132. The R_eff effect is **strongest in the lowest Vflat bin** and disappears at high Vflat. This is the expected pattern: coherence effects dominate in the low-acceleration (low-mass) regime.

### 5. Surface Density (Σ) vs N_corr

| Predictor | r(X, offset | Vflat) | p |
|---|---|---|
| log N_corr (∝ V²/R) | +0.306 | 0.0002 |
| log Σ_eff (∝ V²/R²) | -0.240 | 0.0043 |

N_corr (1/R scaling) outperforms Σ (1/R² scaling). The data prefer the linear R⁻¹ dependence over the quadratic R⁻² dependence. This is consistent with N_corr = V²/(R × a₀) being the correct combination.

### 6. Acceleration-Regime Dependence (STRONGEST EVIDENCE)

| Regime | r(R_eff, offset | Vflat) | p |
|---|---|---|
| **MOND (g < g†)** | **-0.366** | **< 10⁻⁴** |
| Newtonian (g ≥ g†) | -0.116 | 0.473 |

The R_eff effect is **3.2x stronger in the MOND regime** and only significant there. This is the prediction of coherence theory: gravitational coherence effects should dominate where g_bar < g†, not in the baryonic-dominated regime.

Also: N_corr is more predictive in MOND (r = 0.539) than Newtonian (r = 0.179), while Vflat shows the opposite pattern (r = 0.428 MOND vs r = 0.326 Newtonian).

### 7. Optimal V^α × R^β

Best fit: V^3.5 × R^-1.5 (R² = 0.240)

| Model | α | β | R² |
|---|---|---|---|
| V^3.5 R^-1.5 (optimal) | 3.5 | -1.5 | 0.240 |
| **N_corr (V² R⁻¹)** | **2** | **-1** | **0.237** |
| Vflat only | 1-4 | 0 | 0.161 |
| Σ (V² R⁻²) | 2 | -2 | 0.110 |

N_corr (α=2, β=-1) achieves R² = 0.237, essentially matching the optimal 0.240. The N_corr combination is near-optimal and physically interpretable.

## Honest Assessment

### Strengths
1. **R_eff IS significant at fixed Vflat** (p = 0.0002) — the degeneracy is partially broken
2. **Correct direction**: compact → higher on RAR → more coherence
3. **MOND-regime dominated**: 3x stronger at g < g† — exactly where coherence should matter
4. **Low-mass dominated**: Strongest in Vflat < 80 km/s bin
5. **N_corr near-optimal**: V²/R is within 0.3 percentage points of the best V^α R^β
6. **Outliers behave as predicted**: 0.15 dex difference between high and low N_corr outliers

### Weaknesses
1. **Sample shrunk to 135 galaxies** (from 171 in Session #386) — different selection
2. **Effect concentrated in low-mass galaxies** — where data quality is worst
3. **Quality still dominates residuals** — r(Q, residual) = -0.40
4. **M/L sensitivity untested** — the R_eff effect could be an M/L artifact for low-mass galaxies
5. **No independent R_eff measurements** — derived from SB and luminosity, not direct
6. **V^3.5 R^-1.5 optimal, not V² R⁻¹** — the exact N_corr formula isn't the mathematical optimum

### What Changed from Session #386?

Session #386 found r(log N_corr, log Vflat) = 0.970 with 171 galaxies. This session finds 0.752 with 135 galaxies. The difference is galaxy selection — this session's prepare_dataset has slightly different filtering. The lower correlation means more dynamic range in R_eff, making the degeneracy test more powerful.

### Grade: B+

The R_eff signal is real, in the correct direction, and acceleration-regime dependent. This is genuinely suggestive of coherence rather than pure mass. But the effect is concentrated in low-mass galaxies where systematics are largest, and M/L sensitivity is untested.

## Implications for Synchronism

1. **N_corr is NOT just a Vflat proxy** — R_eff contributes 7.8 pp of R²
2. **The coherence length (R_eff) matters** — but only in the MOND regime
3. **The low-mass regime is key** — coherence effects dominate for Vflat < 80 km/s
4. **Next test**: M/L sensitivity of the R_eff effect. If it survives M/L variation, it's robust.

## Files Created

- `simulations/session387_degeneracy_break.py`: 8 tests
- `Research/Session387_Degeneracy_Break.md`: This document

---

*Session #387 verified: 8/8 tests passed*
*Grand Total: 535/535 verified*

**Key finding: The N_corr-Vflat degeneracy is PARTIALLY BROKEN. R_eff significantly predicts RAR offset at fixed Vflat (r = -0.306, p = 0.0002, correct sign). The effect is 3x stronger in the MOND regime than the Newtonian regime, concentrated in low-mass galaxies (Vflat < 80 km/s), and N_corr (V²/R) is near-optimal among V^α R^β combinations. This is genuine evidence that N_corr measures something beyond galaxy mass — consistent with gravitational coherence. However, the effect is concentrated where data quality is worst, M/L sensitivity is untested, and Quality still dominates residuals. Grade B+.**
