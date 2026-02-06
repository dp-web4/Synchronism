# Session #487: Grand Synthesis VIII — The Model Improvement Arc

**Date**: 2026-02-06
**Status**: Review (synthesizes Sessions #482-486)

## Arc Overview

Sessions #482-486 form the **Model Improvement Arc**, discovering and validating the logL×f_gas interaction that pushes the 6-variable outer-only model to R² = 0.945, LOO R² = 0.938. This is the largest single advance in the research program since the outer-only model itself (Session #477).

## The Arc's Narrative

| Session | Title | Key Finding | Grade |
|---------|-------|-------------|-------|
| #482 | Residual Forensics | 41% within noise; NN autocorrelation r=+0.46; 7 persistent outliers | B+ |
| **#483** | **The Sixth Variable** | **logL×f_gas pushes LOO R² from 0.896 to 0.938** | **A** |
| **#484** | **6-Var Validation** | **F=73.6; autocorrelation ELIMINATED; bootstrap 100%** | **A** |
| #485 | Type-Specific Models | Cross-prediction fails (R²=0.61-0.67); late-type 3-var LOO=0.957 | A- |
| #486 | M/L Sensitivity | logL×f_gas stable (β≈0.18, |t|>6.7) at all M/L values | B+ |

**Cumulative tests**: 40/40 verified (Sessions #482-486)
**Grand total**: 1205/1205 verified

## The Central Discovery

### The logL×f_gas Interaction

The gas fraction's effect on the RAR offset is luminosity-dependent:
- **Low luminosity** (dwarfs): effective f_gas coefficient ≈ -0.63 (strong)
- **High luminosity** (spirals): effective f_gas coefficient ≈ -0.09 (weak)

Gas fraction has **7× stronger effect** in dwarfs than in giants.

### Evidence Strength

| Test | Result | Source |
|------|--------|--------|
| F-statistic | 73.6 (p ≈ 10⁻¹⁴) | #484 |
| Bootstrap P(ΔR² > 0) | 100% | #484 |
| Bootstrap 95% CI on ΔR² | [0.013, 0.058] | #484 |
| Jackknife min ΔR² | 0.028 (all positive) | #484 |
| NN autocorrelation eliminated | r: +0.46 → +0.005 | #484 |
| M/L invariant | |t| > 6.7 at all M/L | #486 |
| Persistent after outlier removal | ΔLOO = +0.027 (top 5 removed) | #484 |

## The 6-Variable Model

```
offset_outer = -3.379
  + 1.897 × logV
  - 0.548 × logL
  - 0.218 × c_V        [non-significant, t = -0.89]
  - 0.451 × f_gas
  + 0.147 × logV×c_V   [non-significant, t = +1.24]
  + 0.181 × logL×f_gas  [t = +8.58, F = 73.6]
```

| Metric | 5-var (old) | 6-var (new) | Improvement |
|--------|-------------|-------------|-------------|
| R² | 0.911 | **0.945** | +0.034 |
| LOO R² | 0.896 | **0.938** | +0.042 |
| RMS | 0.048 dex | **0.038 dex** | -21% |
| % within noise (0.025 dex) | 41% | **50%** | +9% |
| % within 0.05 dex | 74% | **84%** | +10% |
| NN autocorrelation (k=5) | +0.41 | **+0.005** | eliminated |

## Architecture of the Model

### The 4 Significant Terms

Though the model has 6 predictors (+ intercept), only 4 are significant:

1. **logV** (t = +15.4): Velocity scale — the dominant predictor
2. **logL** (t = -36.2): Luminosity — at fixed V, more luminous galaxies have lower offsets
3. **f_gas** (t = -14.6): Gas fraction — more gas → lower offset
4. **logL×f_gas** (t = +8.6): The interaction — gas matters more in dwarfs

The c_V and logV×c_V terms (t = -0.89, +1.24) are retained for backward compatibility but contribute negligibly in the full-sample model.

### The Late-Type Minimal Model

For late types (T ≥ 7), a 4-variable model (logV, logL, f_gas, logL×f_gas) achieves LOO R² = 0.954 — c_V is irrelevant (t = -0.14). The even simpler 3-variable model (logV + logL + f_gas²) gives LOO R² = 0.957.

### Type-Dependent Structure

| Type | Dominant physics | Best minimal model | LOO R² |
|------|-----------------|-------------------|--------|
| Early (T<4) | Bulge concentration | logV + logL + c_V (3 vars) | 0.877 |
| Mid (4≤T<7) | Gas/disk balance | logV + logL + f_gas (3 vars) | 0.846 |
| Late (T≥7) | Gas dynamics | logV + logL + f_gas² (3 vars) | **0.957** |

Cross-prediction fails across the type boundary: Late→Early R² = 0.61, Early→Late R² = 0.67.

## Novel Prediction Updates

| ID | Previous | Current | Change |
|----|----------|---------|--------|
| NP11 | — | **NEW**: logL×f_gas interaction | R² → 0.945 |
| NP12 | — | **NEW**: Late 3-var LOO = 0.957 | Simplest strong model |
| NP1 | SUPPORTED | Unchanged | a₀ = cH₀/(2π) |
| NP6 | SUPPORTED | Enhanced | N_corr + logL×f_gas |
| NP7 | SUPPORTED | Unchanged | R_eff MOND-dominated |
| NP8 | SUPPORTED | Confirmed at all M/L | M/L-independent |

## Remaining Questions

### 1. What Explains the Last 5.5%?
The model explains 94.5% of offset variance. The remaining 5.5% consists of:
- ~2.4% measurement noise (0.025 dex)
- ~3.1% intrinsic variation (M/L scatter, mass geometry)

The noise ceiling R² = 0.976. Gap remaining: 3.1%.

### 2. Why M/L_disk ≈ 0.8-0.9 Is Optimal
Higher M/L_disk reduces the offset scatter and improves R² up to M/L ≈ 0.8-0.9. This suggests the standard M/L = 0.5 at 3.6 μm slightly underestimates stellar masses, consistent with recent population synthesis models (Chabrier IMF gives M/L ≈ 0.6-0.8).

### 3. The 7 Persistent Outliers
NGC2915, UGC06667, PGC51017, DDO161, UGC01281, UGC07603, NGC0891 persist across all models and M/L values. These likely represent genuine physical edge cases (extended gas disks, poorly resolved kinematics, non-equilibrium dynamics).

## The Research Trajectory

### Cumulative Progress

| Phase | Sessions | Key advance | R² |
|-------|----------|-------------|-----|
| Initial discovery | 385-389 | N_corr prediction | 0.23 |
| 5-variable model | 395-470 | Full regression | 0.87 |
| Outer-only model | 476-477 | Outer radial cut | 0.91 |
| **logL×f_gas** | **482-486** | **Interaction term** | **0.945** |

### Remaining Opportunities

1. **Radial profile analysis**: The 6-var model predicts galaxy-level offsets but not the radial profile of the offset within each galaxy
2. **N_corr sign problem**: The Synchronism prediction γ = 2/√N has the wrong sign (r = +0.55 vs predicted negative)
3. **External validation**: Test on non-SPARC datasets (e.g., THINGS, LITTLE THINGS)
4. **Theoretical derivation**: Can the logL×f_gas interaction be derived from the Synchronism framework?

## Summary Statistics

| Metric | Value |
|--------|-------|
| Sessions in arc | 5 (#482-486) |
| Tests verified | 40/40 |
| Grand total | 1205/1205 |
| Key grades | A, A, A-, B+, B+ |
| Model R² | 0.945 |
| Model LOO R² | 0.938 |
| Model RMS | 0.038 dex |
| Improvement over 5-var | +0.042 LOO R² |

---

*Session #487: Grand Synthesis VIII*
*Grand Total: 1205/1205 verified across 87 sessions*

**The Model Improvement Arc (#482-486) discovers and validates the logL×f_gas interaction as the dominant missing term. R² jumps from 0.911 to 0.945, LOO R² from 0.896 to 0.938. The NN autocorrelation is eliminated (r: +0.46 → +0.005). The term is M/L-invariant (|t| > 6.7 at all M/L). Late types reach LOO R² = 0.957 with just 3 variables. Cross-prediction fails across type boundaries. The 6-variable model is the new standard.**
