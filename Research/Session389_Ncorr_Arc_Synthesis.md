# Session #389: N_corr Prediction Arc Synthesis

**Date**: 2026-02-06
**Status**: Synthesis (no tests)

## Arc Overview: Sessions #385-388

This arc derived and tested the N_corr = V²_flat/(R_eff × a₀) prediction from Synchronism's γ = 2/√N_corr framework. It is the most quantitatively successful arc to date.

## Arc Timeline

| Session | Topic | Grade | Key Finding |
|---|---|---|---|
| #385 | g† First Principles | B+ | a₀ = cH₀/(2π) = 94% of MOND; N_corr strongest offset predictor (r = +0.48) |
| #386 | N_corr Prediction Test | A- | R² = 0.226, survives CV; but 97% correlated with Vflat |
| #387 | Degeneracy Breaking | B+ | R_eff significant at fixed Vflat (r = -0.31, p = 0.0002); MOND-regime dominated |
| #388 | M/L Robustness | A- | Signal survives all M/L tests; STRONGER in gas-dominated galaxies |

## The N_corr Story

### What We Derived
From γ = 2/√N_corr and the Hubble frequency connection:
- N_corr = V²_flat/(R_eff × a₀) measures gravitational coherence
- Higher N_corr → lower γ → less departure from standard gravity → less RAR offset
- The transition a₀ = cH₀/(2π) connects the coherence scale to the Hubble volume

### What We Found

```
N_corr = V²_flat / (R_eff × a₀)
       ↓
RAR offset = 0.105 × log(N_corr) - 0.038
       ↓
R² = 0.226, p = 2×10⁻¹², survives cross-validation
```

The prediction:
1. **Direction**: Correct (higher N_corr → higher on RAR)
2. **Magnitude**: R² = 0.23, explaining 23% of offset variance
3. **Specificity**: Subsumes Vflat entirely in partial correlations
4. **R_eff component**: Significant at fixed Vflat (p = 0.0002)
5. **M/L robustness**: Survives all 7 M/L tests
6. **Gas-dominated**: STRONGER where M/L is irrelevant (r = -0.59)
7. **MOND-regime**: 3x stronger at g < g† than g > g†
8. **SB suppressor**: Controlling SB unmasks r = -0.68

### The Evidence Pyramid

```
Level 5 (Strongest):
  ✓ Gas-dominated galaxies show strongest R_eff effect
  ✓ M/L variation never eliminates the signal

Level 4:
  ✓ SB control reveals suppressor (r doubles)
  ✓ Monte Carlo M/L 95% CI excludes zero
  ✓ Effect is MOND-regime dominated

Level 3:
  ✓ R_eff significant at fixed Vflat
  ✓ N_corr near-optimal among V^α R^β combinations
  ✓ Cross-validation stable

Level 2:
  ✓ Correct prediction direction
  ✓ Reasonable R² (0.23)

Level 1 (Supporting):
  ✓ Three independent a₀ derivation routes agree
  ✓ Measured γ_eff at transition = 1.24 (theory: 1.23)
```

## What This Arc Establishes

1. **N_corr is a genuine physical predictor** — not just a mass proxy
2. **The R_eff (size) component carries independent information**
3. **The effect is strongest where modified gravity matters** (low acceleration, gas-dominated)
4. **M/L systematics cannot explain the signal**
5. **a₀ = cH₀/(2π) provides a physical connection** to cosmology

## What This Arc Does NOT Establish

1. **Causation**: N_corr could be a correlate, not the cause
2. **Uniqueness**: Other size-dependent quantities might work equally well
3. **External validation**: All analysis uses one dataset (SPARC)
4. **Mechanism**: The word "coherence" is a label, not an explanation
5. **QFT grounding**: No field-theoretic derivation of the coherence mechanism

## Updated Prediction Scorecard

| ID | Prediction | Status | Evidence |
|----|-----------|--------|----------|
| NP1 | a₀ = cH₀/(2π) | SUPPORTED (94%) | Session #385 |
| NP2 | Morphology → scatter | STRONGLY SUPPORTED | p = 5×10⁻⁶ |
| NP3 | High-z a₀ evolution | UNTESTED | Needs JWST |
| NP4 | V-shaped scatter | SUGGESTIVE | Session #374 |
| NP5 | Local gravity anomalies | UNTESTED | Needs Gaia DR3 |
| **NP6** | **N_corr → RAR offset** | **SUPPORTED** | **R² = 0.23, M/L-robust** |
| **NP7** | **R_eff effect MOND-dominated** | **SUPPORTED** | **3x stronger at g < g†** |
| **NP8** | **R_eff effect M/L-independent** | **SUPPORTED** | **7/7 robustness tests** |

## Overall Synchronism Status After 389 Sessions

### What's Working
- **Empirical predictions**: NP2, NP6, NP7, NP8 all supported with SPARC data
- **Dimensional analysis**: a₀ = cH₀/(2π) correct to 6%
- **Quantitative framework**: N_corr = V²/(R×a₀) → offset works
- **Three-component scatter model**: Roughness 51%, offset 11%, unexplained 38%

### What's Not Working
- **First principles**: No derivation of WHY γ = 2/√N_corr
- **The 2π factor**: Assumed, not derived
- **N_corr < 1 for late types**: Theoretically problematic
- **No QFT grounding**: "Coherence" remains undefined at the fundamental level

### What's Next
The natural progression from this arc:

1. **Alternative R_eff measures**: Use HI radius or optical radius instead of SB-derived R_eff
2. **External field effect**: Compare N_corr predictions with Chae+ (2020) EFE estimates
3. **The SB suppressor**: Understand WHY controlling for SB doubles the R_eff signal
4. **High-z prediction**: a₀(z) = cH(z)/(2π) should be testable with JWST rotation curves
5. **Theoretical grounding**: What physical mechanism produces the coherence function?

## Arc Statistics

- Sessions: 4 (#385-388) + synthesis (#389)
- Tests: 32/32 verified
- Grand Total: 543/543 verified
- Arc grade: **A-** (strongest arc yet — genuine quantitative prediction with M/L robustness)

## Files Created

- `Research/Session389_Ncorr_Arc_Synthesis.md`: This document

---

*Grand Total: 543/543 verified (32 from this arc)*

**Arc summary: N_corr = V²_flat/(R_eff × a₀) is Synchronism's strongest quantitative prediction. It explains 23% of RAR offset variance, the R_eff component is significant at fixed Vflat (p = 0.0002), the effect is 3x stronger in the MOND regime, and survives every M/L robustness test including being STRONGER in gas-dominated galaxies. This arc establishes that N_corr measures genuine physics beyond galaxy mass — consistent with the gravitational coherence interpretation. The main gaps remain: no first-principles derivation, no external validation, and no field-theoretic grounding.**
