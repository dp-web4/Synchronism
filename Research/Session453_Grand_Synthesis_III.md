# Session #453: Grand Synthesis III — The Complete RAR Decomposition

**Date**: 2026-02-06
**Status**: Review (no simulation)

## Scope

This document synthesizes Sessions 449-452 (the "extension arc") and provides the definitive summary of the entire research program (Sessions 403-452, 52 sessions, 973 verified tests).

## The Complete Result in One Paragraph

The Radial Acceleration Relation's galaxy-to-galaxy scatter is now decomposed into five sources with a single equation: **offset = -5.51 + 2.77×logV - 0.49×logL + 2.29×c_V - 0.18×f_gas - 0.92×logV×c_V** (R² = 0.872, LOO RMS = 0.059 dex). The five components are: stellar M/L variation (62%, via BTFR residual), MOND phantom dark matter (13%, via velocity concentration), gas fraction correction (6%, via M/L modulation), mass-dependent geometry (6%, via V×c_V interaction), and irreducible residual (~13%, of which ~10% is measurement noise and ~3% is genuine physics). The geometry correction vanishes at V ≈ 305 km/s (typical L* galaxy mass), matching MOND's prediction that phantom DM is absent in the Newtonian regime. The original Synchronism prediction γ = 2/√N_corr is falsified.

## The Final Model

### Equation
```
offset = -5.51 + 2.77×logV - 0.49×logL + 2.29×c_V - 0.18×f_gas - 0.92×logV×c_V
```

### Equivalent Form
```
offset = β₀ + β_BTFR × (logV - 0.25×logL) + β_L × logL
         + β_geom × c_V × log(305/V_flat) + β_gas × f_gas
```

The geometry correction is **c_V × log(V₀/V)** — velocity concentration scaled by distance to the crossover velocity.

### Performance

| Metric | V+L+c_V (3-var) | V+L+c_V+f_gas+V×c_V (5-var) | Improvement |
|--------|-----------------|-------------------------------|-------------|
| R² (galaxy-level) | 0.754 | **0.872** | +15.7% |
| LOO RMS (galaxy) | 0.080 | **0.059** | -26.3% |
| BIC | -636.8 | **-710.6** | -73.8 |
| Point-level RMS | 0.143 dex | **0.136 dex** | -4.9% |
| RC prediction RMS | 0.072 dex | **0.069 dex** | -4.2% |
| Deep MOND RMS | 0.204 dex | **0.126 dex** | -38.1% |
| Galaxies improved (LOO) | — | 89/128 (70%) | — |

### Variance Budget

| Source | % of total | Nature | Key variable | Discovery |
|--------|-----------|--------|-------------|-----------|
| V (mass scale) | 17.8% | Dynamical | V_flat | Session 436 |
| L (M/L correction) | 44.4% | Stellar population | BTFR residual | Session 436 |
| c_V (geometry) | 13.1% | MOND phantom DM | V(R_eff)/V_flat | Session 436 |
| f_gas (gas correction) | 6.0% | M/L modulation | Gas mass fraction | Session 450 |
| V×c_V (mass-dep. geometry) | 5.8% | MOND regime dependence | logV × c_V | Session 451 |
| Velocity noise | ~4% | Measurement | e_v | Session 442 |
| Distance noise | ~6% | Measurement | Distance errors | Session 442 |
| True scatter | **~3%** | Unknown | — | Session 452 |

## The Extension Arc (Sessions 449-452): Four Key Findings

### 1. f_gas Is the Fourth Variable (Session 449-450)

**Discovery**: The V+L+c_V residuals correlate strongly with gas fraction: r(f_gas, resid | V,L,c_V) = -0.494.

**Physical mechanism**: The V+L correction calibrates stellar M/L, but gas (M/L = 1 exactly) doesn't need calibration. Gas-rich galaxies are over-corrected.

**Impact**: R² improves from 0.754 to 0.814 (ΔR² = 0.060). BIC improves by -31.

### 2. V×c_V Interaction Reveals Mass-Dependent Geometry (Session 451)

**Discovery**: The c_V effect depends on galaxy mass. Effective c_V coefficient:
- V = 30 km/s (deep MOND): +0.93
- V = 120 km/s (mild MOND): +0.37
- V = 300 km/s (Newtonian): +0.01

**Physical mechanism**: Phantom dark matter only exists where MOND modifications are active. Low-mass galaxies are deep in the MOND regime → maximum phantom DM → maximum c_V effect.

**Impact**: R² improves from 0.814 to 0.872 (ΔR² = 0.058). Equivalent to c_V × log(305/V).

### 3. Only ~3% True Physics Remains (Session 452)

**Discovery**: After subtracting ~10% measurement noise from the 12.8% unexplained, only ~3% of variance reflects genuine unmodeled physics.

**Implication**: The 5-variable model captures ~97% of the physical signal. Any "new physics" contribution to RAR scatter is constrained to ≤3%.

### 4. Synchronism's Status (Session 452)

**Falsified**: γ = 2/√N_corr, N_corr as key variable, geometric component as Synchronism physics.

**Surviving**: a₀ ≈ cH₀/(2π) (13% agreement with H₀ = 67.4; better with H₀ ≈ 74).

**Outlook**: The ~3% residual is too small to distinguish competing theories. Future tests require redshift evolution or larger datasets.

## Physical Story

The five model variables tell a coherent story:

1. **V** (17.8%): The galaxy's mass scale determines its position on the BTFR and its overall acceleration regime.

2. **L** (44.4%): At fixed V, luminosity reveals the stellar mass-to-light ratio. Brighter galaxies at fixed V have lower M/L, meaning the standard M/L = 0.5 assumption over-estimates their baryonic mass.

3. **c_V** (13.1%): Velocity concentration measures how non-spherical the effective mass distribution is. The algebraic RAR assumes spherical symmetry; concentrated galaxies (high c_V) have steeper inner potentials → phantom dark matter → positive offset.

4. **f_gas** (6.0%): Gas fraction modulates the M/L calibration. Gas has M/L = 1 exactly, so gas-dominated galaxies need less stellar M/L correction.

5. **V×c_V** (5.8%): The phantom DM effect depends on the acceleration regime. At low V (deep MOND), the modified Poisson equation maximally deviates from the algebraic RAR. At high V (Newtonian), the algebraic RAR is exact and no geometry correction is needed.

**All five components are physically motivated**: BTFR calibration (V+L) + MOND phantom DM (c_V, V×c_V) + gas correction (f_gas). There is no ad hoc variable.

## Connections to MOND Theory

| Empirical finding | MOND prediction | Match? |
|-------------------|----------------|--------|
| c_V → positive offset | Phantom DM is positive | YES |
| Effect ~20% (high vs low c_V) | N-body: ~10-20% | YES |
| Inner-dominated (25× at r<R_eff) | Inner disk = non-spherical | YES |
| Vanishes at V≈305 km/s | No phantom DM in Newtonian regime | YES |
| c_V² ≈ enclosed mass fraction (r=0.94) | Geometry from mass distribution | YES |
| f_gas → negative coefficient | Gas doesn't need M/L calibration | YES |

The 5-variable model is essentially an **empirical implementation of full MOND with M/L calibration**. A numerical solution of the MOND modified Poisson equation with variable M/L should reproduce these results.

## Research Program Statistics

| Metric | Value |
|--------|-------|
| Total sessions | 52 (403-452, plus reviews) |
| Total verified tests | 973 |
| Python simulations | 52 |
| Sample size | 128 galaxies, 2850 points |
| Galaxy-level R² | 0.872 |
| LOO RMS | 0.059 dex |
| Physical variance explained | ~97% |
| Original prediction (γ = 2/√N_corr) | Falsified |

## What This Research Achieved

### For Galaxy Physics
1. **Quantified RAR scatter decomposition** at unprecedented detail (5 components)
2. **Connected empirical scatter to MOND phantom DM** (sign, magnitude, radial profile)
3. **Discovered mass-dependent geometry**: c_V effect vanishes at V≈305 km/s
4. **Established gas fraction as M/L modulator**: explains 6% of variance
5. **Set the noise floor**: ~3% genuine physics remains unexplained
6. **Provided 8 testable constraints** for any modified gravity theory

### For Synchronism Theory
1. **Falsified γ = 2/√N_corr** (wrong sign, irrelevant after M/L)
2. **Falsified N_corr as key variable** (r=0.01 after V+L)
3. **Identified geometry component as MOND**, not Synchronism
4. **Constrained new physics to ≤3%** of RAR variance
5. **Preserved a₀ = cH₀/(2π)** as the surviving prediction

### Methodological Contributions
1. **Suppressor variable analysis**: V+L reveals structure hidden in raw correlations
2. **Galaxy-level vs point-level**: constant galaxy shift outperforms variable corrections
3. **LOO as essential validation**: in-sample R² can be misleading
4. **Permutation testing**: necessary for significance with small samples
5. **Sequential model building**: adding one variable at a time reveals contributions

## The Definitive Model

```
RAR offset = -5.51 + 2.77 log V - 0.49 log L + 2.29 c_V - 0.18 f_gas - 0.92 log V × c_V

Where:
  V = V_flat (km/s)
  L = luminosity (10⁹ L_sun)
  c_V = V(R_eff) / V_flat
  f_gas = V²_gas / (V²_gas + V²_disk) at flat region

R² = 0.872, LOO RMS = 0.059 dex, BIC = -710.6
```

---

*Session #453: Grand Synthesis III — The Complete RAR Decomposition*
*Grand Total: 973/973 verified across 52 sessions*

**The RAR scatter is now fully decomposed: M/L calibration (62%) + MOND phantom DM (19%) + gas correction (6%) + noise (10%) + true scatter (3%). The 5-variable model achieves R²=0.872 and captures ~97% of physical variance. The original Synchronism prediction is falsified, but a₀ = cH₀/(2π) survives. The research program is complete for the SPARC dataset.**
