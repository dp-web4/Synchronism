# Session #462: Grand Synthesis V — The a₀ Arc and the End of Synchronism's Last Prediction

**Date**: 2026-02-06
**Status**: Review (no simulation)

## Scope

This document synthesizes Sessions 458-461 (the "a₀ arc") and provides the final assessment of Synchronism's only surviving prediction: a₀ = cH₀/(2π).

## The a₀ Arc in One Paragraph

Session 458 found that SPARC's point-level best-fit a₀ = 1.043 × 10⁻¹⁰, remarkably close to cH₀/(2π) = 1.042 × 10⁻¹⁰ with Planck H₀. Session 459 showed this depends on the assumed M/L: the degeneracy line M/L ≈ -0.315 × a₀ + 0.86 means (M/L=0.50, a₀=1.04) and (M/L=0.48, a₀=1.20) give equal fits. Session 460 revealed that a₀ varies by regime: 0.99 in deep MOND vs 1.27 in moderate MOND (>4σ). Session 461 resolved this: the interpolation function exponent α = 0.458 (not 0.500) is preferred by ΔBIC = -51, and with free α the best-fit a₀ shifts to 1.28. **The a₀ = cH₀/(2π) agreement was an artifact of fixing α = 0.5.**

## The Four Sessions

### Session 458: SPARC Prefers a₀ ≈ cH₀/(2π) — Or Does It?

- Point-level best-fit: a₀ = 1.043 × 10⁻¹⁰ (0.0σ from Planck cH₀/(2π), 1.8σ from MOND)
- Bootstrap: 1.046 ± 0.087
- **But**: Late types prefer 0.89, early types prefer 1.19 (30% subsample variation)
- **Grade: A-**

### Session 459: The M/L-a₀ Degeneracy

- M/L and a₀ are strongly degenerate: M/L ≈ -0.315 × a₀ + 0.861
- Profile likelihood minimum: a₀ = 0.89 at M/L = 0.575
- The RMS surface is flat: only 0.8% variation across [0.7, 1.5] × 10⁻¹⁰
- At a₀ = 1.2 (MOND): optimal M/L = 0.475 (plausible)
- At a₀ = 1.04 (Planck): optimal M/L = 0.525 (plausible)
- **Cannot distinguish theories with M/L free**
- **Grade: B+**

### Session 460: a₀ Varies With Regime

- Deep MOND: a₀ = 0.994 ± 0.027
- Moderate MOND: a₀ = 1.270 ± 0.050
- These differ by >4σ — the McGaugh function is imperfect
- Galaxy-by-galaxy a₀ spans 0.40-2.50 with enormous scatter
- Low-mass galaxies (V < 60) prefer a₀ ≈ 0.50
- "Zero correlation" a₀ = 0.894
- **Grade: A-**

### Session 461: The Interpolation Function Fix

- Generalized exp(-(g_bar/a₀)^α) with α = 0.458 is preferred (ΔBIC = -51)
- With free α: best-fit a₀ = 1.276 (not 1.04)
- The α-a₀ degeneracy: lower α ↔ higher a₀
- Regime dependence partially resolved (gap shrinks from 0.29 to 0.09)
- Gas-rich galaxies still prefer low a₀ regardless of α
- **Grade: A-**

## The Death of a₀ = cH₀/(2π)

The a₀ = cH₀/(2π) prediction was the last surviving element of Synchronism. It is now shown to be an artifact:

### The Artifact Chain

1. **Start**: MOND's standard interpolation uses α = 0.5 (i.e., exp(-√x))
2. **Assumption**: M/L_disk = 0.5 (conventional choice)
3. **Consequence**: At these values, the RAR's point-level scatter is minimized at a₀ ≈ 1.04
4. **Coincidence**: This happens to equal cH₀/(2π) with Planck H₀

But:
- **Free α**: When α is a free parameter, α = 0.46 is preferred and a₀ shifts to 1.28
- **Free M/L**: When M/L is free, a₀ shifts depending on M/L choice
- **Regime dependence**: Different acceleration ranges prefer different a₀
- **Subsample dependence**: Late types and early types prefer different a₀

**The "agreement" is fragile — it depends on three simultaneous assumptions being fixed at conventional values.**

### Could It Still Be Real?

The a₀ = cH₀/(2π) prediction could survive if:
1. Theoretical derivation shows α = 0.5 EXACTLY (not approximately)
2. SPS models firmly establish M/L_disk = 0.50 ± 0.02
3. The gas-rich subsample low-a₀ preference is explained by observational effects

None of these conditions are currently met. The prediction is **not falsified** (it's within the broad uncertainty) but it is **not supported** (the agreement vanishes under free-parameter analysis).

## The Complete Status of Synchronism

| Prediction | Status | Session(s) |
|------------|--------|-----------|
| γ = 2/√N_corr | **FALSIFIED** (wrong sign) | 430 |
| N_corr as key variable | **FALSIFIED** (r=0.01 after M/L) | 444 |
| Geometric offset as Synchronism physics | **FALSIFIED** (it's MOND phantom DM) | 447 |
| a₀ = cH₀/(2π) | **NOT SUPPORTED** (artifact of fixed α, M/L) | 458-461 |

**Synchronism has no surviving testable predictions with SPARC data.**

## What the a₀ Arc Achieved

Despite falsifying the last Synchronism prediction, the a₀ arc produced valuable results:

### 1. The M/L-a₀ Degeneracy (Session 459)
A clean demonstration that rotation curve data alone cannot uniquely determine both M/L and a₀. The degeneracy line M/L ≈ -0.315 × a₀ + 0.86 is a fundamental limitation.

### 2. Interpolation Function Imperfection (Sessions 460-461)
The first clear demonstration that the McGaugh interpolation function is imperfect at the ~2% level, with α = 0.458 preferred over 0.500 by ΔBIC = -51. This is a testable prediction for MOND theory.

### 3. Regime-Dependent a₀ (Session 460)
A₀ varies by regime: 0.99 (deep MOND), 1.27 (moderate MOND), 0.68 (Newtonian). This systematic pattern constrains any interpolation function and rules out some theoretical forms.

### 4. Gas-Rich Galaxies Are Special (Sessions 458-461)
Gas-dominated galaxies consistently prefer lower a₀ across all analyses. This likely reflects their deep-MOND nature and residual M/L sensitivity, but deserves further study with better data.

## Research Program Statistics

| Metric | Value |
|--------|-------|
| Total sessions | 62 (403-462, including 5 synthesis reviews) |
| Verified tests | **1029** |
| Python simulations | 57 |
| Synthesis documents | 5 |
| Sample | ~130 galaxies, ~3000 data points |
| Original predictions | 0/4 surviving |
| Novel methodological contributions | ~10 |

## The Best Interpolation Function

```
g_obs = g_bar / (1 - exp(-(g_bar/a₀)^α))

Best fit: a₀ = 1.276 × 10⁻¹⁰ m/s², α = 0.458
ΔBIC = -51 vs standard McGaugh (α = 0.5)
RMS = 0.1768 dex
```

This is the best-fitting single-function description of the SPARC RAR.

## What Comes Next

The SPARC RAR is now comprehensively characterized:
- **5-variable galaxy-level model** (R² = 0.872) decomposes the scatter
- **Generalized interpolation function** (α = 0.46) better fits the point-level data
- **All Synchronism predictions are falsified or unsupported**

Further progress requires:
1. **New data**: MeerKAT, SKA, and future surveys with larger samples, better distances, and independent M/L estimates
2. **Numerical MOND**: Full Poisson-equation solutions for SPARC galaxies to verify the phantom DM interpretation
3. **Redshift evolution**: Testing whether a₀ (or α) varies with cosmic time

The Synchronism research program is complete.

---

*Session #462: Grand Synthesis V — The a₀ Arc and the End of Synchronism's Last Prediction*
*Grand Total: 1029/1029 verified across 62 sessions*

**The a₀ = cH₀/(2π) agreement is an artifact of fixing α=0.5 and M/L=0.5. With free interpolation exponent, ΔBIC=-51 favors α=0.458 and a₀ shifts to 1.28. With free M/L, the a₀ constraint spans [0.7,1.5] depending on M/L choice. Synchronism has no surviving testable predictions with SPARC data. The best interpolation function is exp(-(g_bar/a₀)^0.458) with a₀=1.276×10⁻¹⁰. The research program is complete: 62 sessions, 1029 verified tests, comprehensive RAR decomposition.**
